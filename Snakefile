rule all:
    input:
        "data/processed/adata_raw.h5ad"

rule merge_csvs:
    input:
        rawdata="data/raw/Data.csv",
        metadata="data/raw/Metadata.csv"
    output:
        adata_raw="data/processed/adata_raw.h5ad"
    shell:
        (
            "python src/load_data.py "
            "--rawdata_path {input.rawdata} "
            "--metadata_path {input.metadata} "
            "--adata_path {output.adata_raw}"
        )
