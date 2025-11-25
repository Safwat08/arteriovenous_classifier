# normalization_log1p.py
# Script takes in a adata_qc.h5ad file, conducts log normalization and then saves a adata_norm.h5ad

import scanpy as sc
from pathlib import Path
import argparse
import numpy as np

def main(adata_qc_path: str, adata_norm_path: str):

    # Ensure output path is okay
    Path(adata_qc_path).parent.mkdir(parents=True, exist_ok=True)
    Path(adata_norm_path).parent.mkdir(parents=True, exist_ok=True)

    print("loading adata_qc from", adata_qc_path)
    
    # Read data
    adata = sc.read_h5ad(adata_qc_path)

    adata = adata[adata.obs["Endothelial cell"] == "Yes"].copy()

    labels = adata.obs["Cluster"].str.lower()

    adata.obs["vascular_class"] = (
        np.select(
            [
                labels.str.contains("artery") | labels.str.contains("arterial") | labels.str.contains("arteriole"),
                labels.str.contains("vein") | labels.str.contains("venous") | labels.str.contains("venule"),
                labels.str.contains("proliferating"),
                labels.str.contains("capillary") & ~labels.str.contains("arterial") & ~labels.str.contains("arteriole") & ~labels.str.contains("venous") & ~labels.str.contains("venule"),
                labels.str.contains("lymphatic") | labels.str.contains("lypmhatic"),
                labels.str.contains("angiogenic"),
                labels.str.contains("choroidplexus") | labels.str.contains("choroid plexus") | labels.str.contains("glomeruli") | labels.str.contains("colon_...7", regex=True),
                labels.str.contains("interferon"),
            ],
            [
                "Arterial",
                "Venous",
                "Proliferating",
                "Capillary",
                "Lymphatic",
                "Angiogenic",
                "Capillary",
                "Immunogenic",
            ],
            default="Unassigned",
        )
    )

    adata.obs["vascular_subclass"] = (
        np.select(
            [
                labels.str.contains("large") & labels.str.contains("artery"),
                labels.str.contains("artery") & ~labels.str.contains("large"),
                labels.str.contains("arterial"),
                labels.str.contains("arteriole"),
                labels.str.contains("vein"),
                labels.str.contains("venous"),
                labels.str.contains("venule"),
                labels.str.contains("proliferating"),
                labels.str.contains("capillary") & ~labels.str.contains("arterial") & ~labels.str.contains("arteriole") & ~labels.str.contains("venous") & ~labels.str.contains("venule"),
                labels.str.contains("lymphatic") | labels.str.contains("lypmhatic"),
                labels.str.contains("angiogenic"),
                labels.str.contains("choroidplexus") | labels.str.contains("choroid plexus") | labels.str.contains("glomeruli") | labels.str.contains("colon_...7", regex=True),
                labels.str.contains("interferon"),
            ],
            [
                "Large Artery",
                "Artery",
                "Arterial",
                "Arteriole",
                "Vein",
                "Venous",
                "Venule",
                "Proliferating",
                "Capillary",
                "Lymphatic",
                "Angiogenic",
                "Capillary",
                "Immunogenic",
            ],
            default="Unassigned",
        )
    )

    # Log normalize
     # Normalizing to median total counts
    sc.pp.normalize_total(adata)

    # Logarithmize the data
    sc.pp.log1p(adata)

    # Get highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # Conduct PCA
    sc.tl.pca(adata)

    # Get nearest neighbor graph
    sc.pp.neighbors(adata)

    # Get UMAP
    sc.tl.umap(adata)

    # Conduct clustering 
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    # QC Plot 2: Scatter of ngenes by count & total counts, colored by outlier
    sc.pl.umap(
        adata,
        color=["Tissue", "vascular_class", "vascular_subclass", "leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
        save = "umap1_metrics"
    )

    print("saving adata_norm to", adata_norm_path)
   
    # Save file
    adata.write(adata_norm_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata_qc_path", required=True)
    ap.add_argument("--adata_norm_path", required=True)
    args = ap.parse_args()
    main(args.adata_qc_path, args.adata_norm_path)