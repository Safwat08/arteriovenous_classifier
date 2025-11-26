import scanpy as sc
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from sklearn.preprocessing import LabelEncoder
import wandb


# inputs:

project_name = "classify_arteriovenous_ECs"

model_name = "SubtypeClassifier_ver1"

fold = 1

batch_size = 64

num_epochs = 10
lr = 1e-3
weight_decay = 1e-4


class SubtypeClassifier(nn.Module):
    def __init__(
        self,
        tissue_input_n: int, # len(adata_norm.obs["Tissue"].unique())
        class_input_n: int, # len(adata.obs["vascular_subclass"].unique)
        tissue_output_dim: int = 8,
        feature_input_dim: int = 50, # adata.varm['PCs'].shape[1]
        encoder_hidden_dim: int = 64,
        encoder_output_dim: int = 32,
    ):
        super().__init__()

        # Get Tissue Embedding
        self.tissue_embedder = nn.Embedding(
            num_embeddings = tissue_input_n,
            embedding_dim = tissue_output_dim
        )

        # Concat tissue embedding and input feature dim
        encoder_input_dim = tissue_output_dim + feature_input_dim

        # Cell-tissue encoder:
        self.cell_tissue_encoder = nn.Sequential(
            nn.Linear(encoder_input_dim, encoder_hidden_dim),
            nn.ReLU(),
            nn.BatchNorm1d(encoder_hidden_dim), # Batch Normalization 1D
            nn.Linear(encoder_hidden_dim, encoder_output_dim)
        )

        # Class embeddings:
        self.class_embedder = nn.Parameter(
            torch.randn(class_input_n, encoder_output_dim) * 0.1 # Same as output of cell_tissue_encoder
        )

    def forward(self, x, tissue_idx):

        t_emb = self.tissue_embedder(tissue_idx)

        h = torch.cat([x,t_emb], dim=1)

        z = self.cell_tissue_encoder(h)

        z = F.normalize(z, dim=1)

        c_emb = F.normalize(self.class_embedder, dim=1)

        logits = z @ c_emb.T

        return logits, z


class ScanpyECDataset(Dataset):
    def __init__(
        self,
        adata,
        tissue_encoder: LabelEncoder = None,
        class_encoder: LabelEncoder = None,
        tissue_var_name: str = "Tissue",
        class_var_name: str = "vascular_subclass",
    ):
        # Use PCs
        X = adata.obsm["X_pca"]
        self.tissue_encoder = tissue_encoder
        self.class_encoder = class_encoder

        tissue_ids = self.tissue_encoder.transform(adata.obs[tissue_var_name])
        class_ids = self.class_encoder.transform(adata.obs[class_var_name])

        self.X = torch.tensor(X, dtype=torch.float32)
        self.tissue_ids = torch.tensor(tissue_ids, dtype=torch.long)
        self.class_ids  = torch.tensor(class_ids, dtype=torch.long)

        # Ensure lenght of samples/cells is equal
        assert self.X.shape[0] == self.tissue_ids.shape[0] == self.class_ids.shape[0]

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        # Single sample
        x = self.X[idx]              # (n_features,)
        tissue_id = self.tissue_ids[idx]  # scalar long
        class_id  = self.class_ids[idx]   # scalar long

        return x, tissue_id, class_id

# Load adata
adata_val = sc.read_h5ad('data/splits/adata_test_fold1.h5ad')
adata_train = sc.read_h5ad('data/splits/adata_train_fold1.h5ad')

tissue_encoder = LabelEncoder()
class_encoder  = LabelEncoder()

tissue_encoder.fit(adata_train.obs["Tissue"])
class_encoder.fit(adata_train.obs["vascular_subclass"])

# Create Dataset
train_dataset = ScanpyECDataset(adata_train,
                                tissue_encoder = tissue_encoder,
                                class_encoder = class_encoder)

val_dataset = ScanpyECDataset(adata_val,
                              tissue_encoder = tissue_encoder,
                              class_encoder = class_encoder)

train_loader = DataLoader(
    train_dataset,
    batch_size=batch_size,
    shuffle=True,      # shuffle for training
    drop_last=False,
)

val_loader = DataLoader(
    val_dataset,
    batch_size=len(adata_val),
    shuffle=False,     # no need to shuffle for validation
    drop_last=False,
)

tissue_len =  len(adata_train.obs["Tissue"].unique())
class_len = len(adata_train.obs["vascular_subclass"].unique())
feature_len = adata_train.varm['PCs'].shape[1]

model = SubtypeClassifier(
    tissue_input_n = tissue_len,
    class_input_n = class_len,
    tissue_output_dim = 8,
    feature_input_dim = feature_len,
    encoder_hidden_dim = 64,
    encoder_output_dim = 32
)

optimizer = optim.Adam(
    model.parameters(),
    lr=lr,
    weight_decay=weight_decay)

criterion = nn.CrossEntropyLoss()

# wandb init
run = wandb.init(
    project=project_name,
    name=f"{model_name}_fold{fold+1}",
    group=model_name,
    config={
        "fold": fold,
        "model": model_name,
        "lr": lr,
        "weight_decay": weight_decay,
        "batch_size": batch_size,
        "epochs": num_epochs,
        "loss": "CrossEntropyLoss",
    },
    reinit=True,
)

for epoch in range(num_epochs):
    model.train()
    running_train_loss = 0.0
    n_train_batches = 0
    for x_batch, tissue_batch, class_batch in train_loader:

        # print(f"Epoch {epoch+1}/{num_epochs}, Batch {n_train_batches+1}/{len(train_loader)}")

        logits, z = model(x_batch, tissue_batch)

        loss = criterion(logits, class_batch)

        optimizer.zero_grad()

        loss.backward()

        optimizer.step()

        # Add to running loss and batch
        running_train_loss += loss.item()

        n_train_batches += 1

    # Calculate average loss across batch
    avg_train_loss = running_train_loss / max(n_train_batches, 1)
    print(f"avg_train_loss = {avg_train_loss}")
    # wandb summary

    # Validation loop
    model.eval()
    running_val_loss = 0.0
    n_val_batches = 0

    # Disable gradient computation for validation
    with torch.no_grad():
        for x_batch, tissue_batch, class_batch in train_loader:

            # print(f"Epoch {epoch+1}/{num_epochs}, Batch {n_train_batches+1}/{len(train_loader)}")

            logits, z = model(x_batch, tissue_batch)

            loss = criterion(logits, class_batch)

            # Add to running loss and batch
            running_val_loss += loss.item()

            n_val_batches += 1

            probs = F.softmax(logits, dim=1)  # (batch, n_classes)
            # predicted class IDs:
            pred_ids = probs.argmax(dim=1)

    # Calculate average loss across batch
    avg_val_loss = running_val_loss / max(n_val_batches, 1)
    print(f"avg_val_loss = {avg_val_loss}")

# Implement metrics
# Macro F1
macro_f1 = f1_score(y_test, y_pred, average="macro")

# Balanced Accuracy
bal_acc = balanced_accuracy_score(y_test, y_pred)

# Classificaiton Report
report = classification_report(y_test, y_pred, output_dict=True)
df_report = pd.DataFrame(report).transpose()

# Confusion Matrix
conf_mat = confusion_matrix(y_test, y_pred)

# AUROC
roc_auc = roc_auc_score(y_test, y_prob, multi_class="ovr")

plt.figure(figsize = (6,5))
sns.heatmap(conf_mat,
            annot = True,
            cmap = "Blues",
            fmt=".2f",
           xticklabels = label_enc.classes_,
           yticklabels = label_enc.classes_)
plt.xlabel("Predicted labels")
plt.ylabel("True labels")
plt.title("Confusion Matrix")
plt.tight_layout()
plt.show()

print(macrof1,"\n",
     bal_acc, "\n",
     df_report, "\n",
     roc_auc)
