configfile: "config/config.yaml"

raw_dir  = config["paths"]["raw"]
out_dir = config["paths"]["processed"]
figs_dir = config["paths"]["figs"]
results_dir = config["paths"]["results"]
split_dir = config["paths"]["splits"]

species_type = config["params"]["species"]

rawdata_file = config["files"]["rawdata"]
metadata_file = config["files"]["metadata"]
adata_raw_file = config["files"]["adata_raw"]
adata_qc_file = config["files"]["adata_qc"]
adata_norm_file = config["files"]["adata_norm"]

# Test/Train splits
adata_test_fold1_file = config["files"]["adata_test_fold1"]
adata_test_fold2_file = config["files"]["adata_test_fold2"]
adata_test_fold3_file = config["files"]["adata_test_fold3"]
adata_test_fold4_file = config["files"]["adata_test_fold4"]
adata_test_fold5_file = config["files"]["adata_test_fold5"]
adata_train_fold1_file = config["files"]["adata_train_fold1"]
adata_train_fold2_file = config["files"]["adata_train_fold2"]
adata_train_fold3_file = config["files"]["adata_train_fold3"]
adata_train_fold4_file = config["files"]["adata_train_fold4"]
adata_train_fold5_file = config["files"]["adata_train_fold5"]

rule all:
    input:
        adata_norm = f"{out_dir}/{adata_norm_file}"

rule load_data:
    input:
        rawdata = f"{raw_dir}/{rawdata_file}",
        metadata = f"{raw_dir}/{metadata_file}"
    output:
        adata_raw = f"{out_dir}/{adata_raw_file}"
    shell:
        """
        python src/load_data.py \
            --rawdata_path {input.rawdata} \
            --metadata_path {input.metadata} \
            --adata_path {output.adata_raw}
        """

rule qc_data:
    input:
        adata_raw = f"{out_dir}/{adata_raw_file}"
    output:
        adata_qc  = f"{out_dir}/{adata_qc_file}"
    params:
        figs      = figs_dir,
        results   = results_dir,
        species    = species_type
    shell:
        """
        python src/qc_data.py \
            --adata_raw_path {input.adata_raw} \
            --adata_qc_path {output.adata_qc} \
            --figures_path {params.figs} \
            --results_path {params.results} \
            --species {params.species}
        """

rule scanpy_pipeline:
    input:
        adata_qc = f"{out_dir}/{adata_qc_file}"
    output:
        adata_norm = f"{out_dir}/{adata_norm_file}"
    shell:
        """
        python src/scanpy_pipeline.py \
            --adata_qc_path {input.adata_qc} \
            --adata_norm_path {output.adata_norm}
        """

rule data_split:
    input:
        adata_path = f"{out_dir}/{adata_norm_file}",
        split_path = split_dir
    output:
        adata_test_fold1 = f"{split_dir}/{adata_test_fold1_file}",
        adata_test_fold2 = f"{split_dir}/{adata_test_fold2_file}",
        adata_test_fold3 = f"{split_dir}/{adata_test_fold3_file}",
        adata_test_fold4 = f"{split_dir}/{adata_test_fold4_file}",
        adata_test_fold5 = f"{split_dir}/{adata_test_fold5_file}",
        adata_train_fold1 = f"{split_dir}/{adata_train_fold1_file}",
        adata_train_fold2 = f"{split_dir}/{adata_train_fold2_file}",
        adata_train_fold3 = f"{split_dir}/{adata_train_fold3_file}",
        adata_train_fold4 = f"{split_dir}/{adata_train_fold4_file}",
        adata_train_fold5 = f"{split_dir}/{adata_train_fold5_file}"

    shell:
        """
        python src/data_split.py \
            --adata_path {input.adata_path} \
            --split_path {input.split_path}
        """