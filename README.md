# Single Cell RNA-seq arteriovenous classifier

A pipeline to assess different machine learning (ML) architectures in classifying arteriovenous ECs using the "Murine EC Atlas" from Kalucka et al. (Kalucka, J. et al. Cell. 2020)

Purpose
To establish a model that can accurately identify arterial, venous and capillary population in single cell RNA-sequencing data.

Methodology

Workflow
Docker and snakemake used for workflow management.

Data acquisition
Raw data counts and metadata was downloaded from the Peter Carmeliet lab shinyapp https://endotheliomics.shinyapps.io/ec_atlas/ . Files downloaded include Data.csv and Metadata.csv and placed under data/raw

Load data
Data.csv was loaded using the polars package to expedite speed. Data and metadata combined into an AnnData object for downstream analysis. Saved under data/processed as adata_raw.h5ad. See load_data.py

Quality Control
adata_raw.h5ad was used to establish several QC metrics including but not limited to percent expression from mitochondiral, ribosomal and hemoglobin genes








