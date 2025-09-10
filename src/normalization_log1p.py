# normalization_log1p.py
# Script takes in a adata_raw.h5ad file, conducts log normalization and then saves a adata_norm.h5ad

import scanpy as sc

def main(rawdata_path, adata_path):

    # Ensure output path is okay
    Path(adata_path).parent.mkdir(parents=True, exist_ok=True)

    print("loading raw data from", rawdata_path)

    

    # Save file
    adata.write(adata_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rawdata_path", required=True)
    ap.add_argument("--metadata_path", required=True)
    ap.add_argument("--adata_path", required=True)
    args = ap.parse_args()
    main(args.rawdata_path, args.metadata_path, args.adata_path)