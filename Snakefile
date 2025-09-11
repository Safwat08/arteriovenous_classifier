configfile: "config/config.yaml"

raw_dir  = config["paths"]["raw"]
out_dir = config["paths"]["processed"]
figs_dir = config["paths"]["figs"]

rule all:
    input:
        "data/processed/adata_raw.h5ad"

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
