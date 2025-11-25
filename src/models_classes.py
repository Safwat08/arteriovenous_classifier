import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import LabelEncoder

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
        tissue_var_name: str = "Tissue",
        class_var_name: str = "vascular_subtype",
    ):
        # Use PCs 
        X = adata.obsm["X_pca"]
        tissue_labels = adata.obs[tissue_var_name]
        class_labels = adata.obs[class_var_name]

        self.X = torch.tensor(X, dtype=torch.float32)

        tissue_encoder = LabelEncoder()
        class_encoder  = LabelEncoder()

        tissue_ids = tissue_encoder.fit_transform(tissue_labels)   # numpy array of ints
        class_ids  = class_encoder.fit_transform(class_labels)

        self.tissue_ids = torch.tensor(tissue_ids, dtype=torch.long)
        self.class_ids  = torch.tensor(class_ids, dtype=torch.long) 

        # Ensure lenght of samples/cells is equal
        assert self.X.shape[0] == self.tissue_ids.shape[0] == self.class_ids.shape[0]

        def __len__(self):
            return self.X.shape[0]
            
        def __getidem__(self, idx):
            # Single sample
            x = self.X[idx]              # (n_features,)
            tissue_id = self.tissue_ids[idx]  # scalar long
            class_id  = self.class_ids[idx]   # scalar long

            return x, tissue_id, class_id


