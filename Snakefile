configfile: "config/config.yaml"

raw_dir  = config["paths"]["raw"]
out_dir = config["paths"]["processed"]
figs_dir = config["paths"]["figs"]

rule all:
    input:
        f"{out_dir}/adata_qc.h5ad"

rule load_data:
    input:
        rawdata=f"{raw_dir}/Data.csv",
        metadata=f"{raw_dir}/Metadata.csv"
    output:
        adata_raw=f"{out_dir}/adata_raw.h5ad"
    shell:
        """
        python src/load_data.py \
            --rawdata_path {input.rawdata} \
            --metadata_path {input.metadata} \
            --adata_path {output.adata_raw}
        """

rule qc_data:
    input:
        adata_raw = f"{out_dir}/adata_raw.h5ad"
    output:
        adata_qc  = f"{out_dir}/adata_qc.h5ad"
    params:
        figs      = figs_dir,
        species   = config.get("species", "mouse")
    shell:
        """
        python src/qc_data.py \
            --adata_raw_path {input.adata_raw} \
            --adata_qc_path {output.adata_qc} \
            --figures_path {params.figs} \
            --species {params.species}
        """
