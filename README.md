# Single Cell RNA-seq arteriovenous classifier

A pipeline to assess different machine learning (ML) architectures in classifying arteriovenous ECs using the "Murine EC Atlas" from Kalucka et al. (Kalucka, J. et al. Cell. 2020)

# Purpose
To establish a model that can accurately identify arterial, venous and capillary population in single cell RNA-sequencing data.

# Methodology

## Workflow
Docker and snakemake used for workflow management.

## Data acquisition
Raw data counts and metadata was downloaded from the Peter Carmeliet lab shinyapp https://endotheliomics.shinyapps.io/ec_atlas/ . Files downloaded include Data.csv and Metadata.csv and placed under data/raw

## Load data
Data.csv was loaded using the polars package to expedite speed. Data and metadata combined into an AnnData object for downstream analysis. Saved under data/processed as adata_raw.h5ad. See load_data.py

## Quality Control
adata_raw.h5ad was used to establish several QC metrics including but not limited to percent expression from mitochondiral, ribosomal and hemoglobin genes. Outliers identified as samples beyond 5 Median Absolute Deviations (MADs) and subsetted. QC plots can be found in figures/. See qc_data.py

## Scanpy Pipeline
QC'd adata was then processed through the generic Scanpy pipeline. Briefly, cells not annotated as "Endothelial cell" were excluded. "Cluster" annotations were remapped to "vascular_subtype" vairable based on vascular hieararchy. Counts were scaled to the median count depth followed by log1p transformation. Top 2000, highly variable genes were identified and used for downstream principal component analysis (PCA). PCA was used for constructured nearest neighbor graph, which was then embedded into two dimensions using UMAP for visualization. Cells were clustered using Leiden clustering. UMAP plots can be found in figures/. See scanpy_pipeline.py

## Train/Test splits
Data was subsetted further to only include features that are "Highly Variable". Data was then split into train/test datasets using stratified 5fold splitting using sklearn's "StratifiedKFold" and splits were saved into data/splits.

## Evaluation Metrics
Macro F1 score, balanced_accuracy_score, per-class precision, recall, F1, AUROC


# Results
## Substantial class imbalance








