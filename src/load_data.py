#load_data
# This script loads the raw data and metadata using polars package, saves the raw data as parquet, then creates an anndata object and saves it in the output folder in data

import polars as pl
import pandas as pd 
import anndata as ad
import scipy.sparse as sp
from pathlib import Path
import argparse

def main(rawdata_path: str, 
         metadata_path: str, 
         adata_path: str):
        
    """
    Load a csv file using polars and convert it into a an AnnData object
    
    Parameters:
    -----------
    rawdata_path : str
        The file path to the data .csv file to be loaded.
    
    metadata_path : str
        The file path to the metadata .csv file to be loaded.

    adata_path : str
        The file path to where the AnnData will be saved
    """

    # Ensure path of AnnData exists
    Path(adata_path).parent.mkdir(parents=True, exist_ok=True)

    # Read raw data using polars (faster for larger files) using lazy/streaming method
    rawdata = pl.scan_csv(rawdata_path).collect(streaming=True)

    # Cell names from index column
    gene_names = rawdata[:, 0].to_list()
    cell_names = rawdata.columns[1:] 

    print(gene_names[0:5])
    print(cell_names[0:5])

    print("loading raw data from", rawdata_path)
    print(rawdata.head())

    # Read metadata with index 0, to read cell names as index
    metadata = pd.read_csv(metadata_path, index_col=0)

    # Reindex incase of mismatch cell_names
    metadata = metadata.reindex(cell_names)

    print("loading metadata from", metadata_path)
    print(metadata.head())

    # Extract numeric values as NumPy
    X = rawdata[:, 1:].to_numpy()

    # Convert to sparse and transpose because should be observations X variables
    X_sparse = sp.csr_matrix(X.T)

    # Build AnnData
    adata = ad.AnnData(X_sparse)

    # Add cell and gene names
    adata.obs_names = cell_names
    adata.var_names = gene_names

    # Add all metadata
    adata.obs = metadata

    # Save file
    adata.write(adata_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rawdata_path", required=True)
    ap.add_argument("--metadata_path", required=True)
    ap.add_argument("--adata_path", required=True)
    args = ap.parse_args()
    main(args.rawdata_path, args.metadata_path, args.adata_path)