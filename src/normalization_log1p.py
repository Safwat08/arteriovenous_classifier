# normalization_log1p.py
# Script takes in a adata_qc.h5ad file, conducts log normalization and then saves a adata_norm.h5ad

import scanpy as sc
from pathlib import Path
import argparse

def main(adata_qc_path: str, adata_norm_path: str):

    # Ensure output path is okay
    Path(adata_qc_path).parent.mkdir(parents=True, exist_ok=True)
    Path(adata_norm_path).parent.mkdir(parents=True, exist_ok=True)

    print("loading adata_qc from", adata_qc_path)
    
    # Read data
    adata = sc.read_h5ad(adata_qc_path)

    # Log normalize
    sc.pp.log1p(adata)  

    print("saving adata_norm to", adata_norm_path)
   
    # Save file
    adata.write(adata_norm_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata_qc_path", required=True)
    ap.add_argument("--adata_norm_path", required=True)
    args = ap.parse_args()
    main(args.adata_qc_path, args.adata_norm_path)