# qc_data
# This script loads the adata_raw.h5ad file, conducts QC, generates QC figures and saves the QC'd AnnData

import scanpy as sc 
import numpy as np 
from pathlib import Path
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
from upsetplot import from_memberships, UpSet
import matplotlib.pyplot as plt
import os

def is_outlier(adata, 
               metric: str, 
               nmads: int):
    """
    Identify outliers in the AnnData object based on the specified metric and 
    the number of median absolute deviations (MADs). 
    Anna Schaar. Theiss Lab

    Parameters:
    -----------
    adata : AnnData
        The AnnData object containing the observation data.
    metric : str
        The name of the column in `adata.obs` on which to detect outliers.
    nmads : int
        The number of median absolute deviations (MADs) to use as the threshold 
        for determining outliers.

    Returns:
    --------
    outlier : pandas.Series
        A boolean Series where `True` indicates that the corresponding 
        observation is an outlier.
    """

    
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def qc_adata(adata, 
         species: str = "mouse"):
    """
    Load an AnnData file and preprocesss using various QC metrics.
    
    Parameters:
    -----------
    adata 
        The adata file to be processeed. Assumes that rawcounts are the default 

    species : str
        The species of the data, default is "mouse" with options for "human"

    Returns:
    --------
    adata : adata file
        The adata file with the qc metrics
    """

    # Make a new layer 'rawcounts' to keep the raw data
    adata.layers['rawcounts'] = adata.X.copy()

    # Handle species-specific mitochondrial gene detection
    if species.lower() == "human":
        mito_prefix = "MT-"
        ribo_prefix = ("RPS", "RPL")
        hemo_prefix = "^HB[^(P)]"
    elif species.lower() == "mouse":
        mito_prefix = "mt-"
        ribo_prefix = ("Rps", "Rpl")
        hemo_prefix = "^Hb[^(P)]"
    else:
        raise ValueError("Species not recognized. Use 'human' or 'mouse'.")
    
    # Insert a column 'gene_names' with the values from the index of adata.var
    adata.var['gene_names'] = adata.var.index       
    
    # Identify mitochondrial, ribosomal and hemoglobin genes based on species-specific prefix
    adata.var["mt"] = adata.var["gene_names"].str.startswith(mito_prefix)
    adata.var["ribo"] = adata.var["gene_names"].str.startswith(ribo_prefix)
    adata.var["hb"] = adata.var["gene_names"].str.contains(hemo_prefix)
    
    # Calculate qc_metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ribo", "hb"], 
        inplace=True, 
        log1p=True, 
        percent_top=[20]
    )
    
    # Filter by 5 MADs depending on counts distributions
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    
    # Filter based on mitochondrial, ribosomal and hemoglobin gene expression
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 5)
    adata.obs["ribo_outlier"] = is_outlier(adata, "pct_counts_ribo", 5)
    adata.obs["hb_outlier"] = is_outlier(adata, "pct_counts_hb", 5)

    return adata

def main(adata_raw_path: str,
         adata_qc_path: str,
         figures_path: str):
    
    # Ensure path of QC'd AnnData exists
    Path(adata_qc_path).parent.mkdir(parents=True, exist_ok=True)

    # Ensure path of QC figures exists
    os.makedirs(figures_path, exist_ok=True)

    adata = sc.read_h5ad(adata_raw_path)

    adata = qc_adata(adata)

    # QC Plot 1: Overlap of all outliers
    
    # build memberships
    memberships = []
    for cell in adata.obs.index:
        mem = []
        if cell in outliers: mem.append("Counts")
        if cell in mt_outliers: mem.append("MT")
        if cell in ribo_outliers: mem.append("Ribo")
        if cell in hb_outliers: mem.append("Hb")
        memberships.append(mem)

    # convert to UpSet data
    data = from_memberships(memberships)
    UpSet(data).plot()

    plt.savefig(figures_path + "qc1_upset.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Change boolean to categorical for easier visualization
    adata.obs["outlier"] = pd.Categorical(adata.obs["outlier"])
    adata.obs["mt_outlier"] = pd.Categorical(adata.obs["mt_outlier"])
    adata.obs["ribo_outlier"] = pd.Categorical(adata.obs["ribo_outlier"])
    adata.obs["hb_outlier"] = pd.Categorical(adata.obs["hb_outlier"])
    
    # QC Plot 2: Scatter of ngenes by count & total counts, colored by outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="outlier",
        save= figures_path + "qc2_scatter_counts.png"
    )

    # QC Plot 3: Scatter of ngenes by count & total counts, colored by mt_outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="mt_outlier",
        save=figures_path + "qc3_scatter_mt.png"
    )

    # QC Plot 4: Scatter of ngenes by count & total counts, colored by mt_outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="ribo_outlier",
        save=figures_path + "qc4_scatter_ribo.png"
    )

    # QC Plot 5: Scatter of ngenes by count & total counts, colored by mt_outlier
    sc.pl.scatter(
        adata,
        x="n_genes_by_counts", 
        y="total_counts",     
        color="hb_outlier",
        save=figures_path + "qc5_scatter_hb.png"
    )

    # QC Plot 6: Violin plot of n genes by count, total counts and pct_counts_mt, colored by outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'outlier',
        save = figures_path + "qc6_violin_counts.png"
    )

    # QC Plot 7: Violin plot of n genes by count, total counts and pct_counts_mt, colored by mt_outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'mt_outlier',
        save = figures_path + "qc7_violin_mt.png"
    )

    # QC Plot 8: Violin plot of n genes by count, total counts and pct_counts_mt, colored by mt_outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'ribo_outlier',
        save = figures_path + "qc8_violin_ribo.png"
    )

    # QC Plot 9: Violin plot of n genes by count, total counts and pct_counts_mt, colored by mt_outlier
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, groupby= 'hb_outlier',
        save = figures_path + "qc9_violin_hb.png"
    )

    # Change back to bool
    adata.obs["outlier"] = adata.obs["outlier"].astype(bool)
    adata.obs["mt_outlier"] = adata.obs["mt_outlier"].astype(bool)
    adata.obs["ribo_outlier"] = adata.obs["ribo_outlier"].astype(bool)
    adata.obs["hb_outlier"] = adata.obs["hb_outlier"].astype(bool)

    # Filter based on outlier
    adata_filtered = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier) & (~adata.obs.ribo_outlier) & (~adata.obs.hb_outlier)].copy()

    adata_filtered.write(adata_qc_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata_raw_path", required=True)
    ap.add_argument("--adata_qc_path", required=True)
    ap.add_argument("--figures_path", required=True)
    args = ap.parse_args()
    main(args.adata_raw_path, args.adata_qc_path, args.figures_path)