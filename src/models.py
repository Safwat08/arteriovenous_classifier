import torch
import torch.nn as nn
import torch.nn.functional as F


class SubtypeClassifier(nn.Module):
    def __init__(
        self,
        tissue_input_n, # len(adata.obs["Tissue"].unique)
        tissue_output_dim: int = 8,
        feature_input_dim: int = 50, # PCA analysis
        encoder_hidden_dim: int = 64,
        encoder_output_dim: int = 32
        class_input_n, # len(adata.obs["vascular_subclass"].unique)
    )
    super().__init__()

    # Get Tissue Embedding
    self.tissue_embedder = nn.Embedding(
        num_embeddings = tissue_input_n,
        embeddings_dim = tissue_output_dim
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

    def forward(self, x, tissue_input_n):

        t_emb = self.tissue_embedder(tissue_input_n)

        h = torch.cat([x,t_emb], dim=1)

        z = self.encoder(h)

        z = F.normalize(z, dim=1)

        c_emb = self.class_embedder

        c_emb = F.normalize(c_emb, dim=1)

        logits = z @ c_emb.T

        return logits, z





