# data_split.py
# Splits adata file into training and test splits and outputs class distribution

from sklearn.model_selection import StratifiedKFold
import scanpy as sc
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from pathlib import Path
import argparse


def main(adata_path: str, split_path: str):
    
    adata_path = "data/processed/adata_norm.h5ad"
    split_path = "data/splits"

    # Path exists
    Path(adata_path).parent.mkdir(parents=True, exist_ok=True)
    Path(split_path).parent.mkdir(parents=True, exist_ok=True)

    # Load adata
    adata = sc.read_h5ad(adata_path)

    # subset only highly variable genes
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Get vascular_subtype labels
    labels = adata.obs['vascular_subtype']

    # StratifiedKFold function
    skf = StratifiedKFold(n_splits = 5, shuffle=True, random_state=42)

    # Get an array of zero to input in skf function
    zero_array = np.zeros(len(labels))

    for fold, (train_idx, test_idx) in enumerate(skf.split(zero_array, labels), start=1):

        # subset_adata based on train-test idx
        adata_train = adata[train_idx]
        adata_test = adata[test_idx]

        adata_train.write(f"{split_path}/adata_train_fold{fold}.h5ad")
        adata_test.write(f"{split_path}/adata_test_fold{fold}.h5ad")

        # Get train counts and values for plotting
        train_counts = adata_train.obs['vascular_subtype'].value_counts()
        train_counts_values = train_counts.values
        train_counts_index = train_counts.index

        # Get test counts and values for plotting
        test_counts = adata_test.obs['vascular_subtype'].value_counts()
        test_counts_values = test_counts.values
        test_counts_index = test_counts.index

        # Plot figures
        fig, axes = plt.subplots(1,2, figsize=(10,5))

        # Train Pie Chart
        axes[0].pie(train_counts_values,
            #  labels = train_counts_index,
                autopct = '%1.1f%%',
                startangle=90,
                pctdistance=1.1,
                colors=plt.cm.Set3.colors,
                wedgeprops={'linewidth': 1, 
                            'edgecolor': 'white'})
        axes[0].set_title("Train Dataset")
        axes[0].legend(train_counts_index,
                loc="center left", 
                bbox_to_anchor=(1, 0.5))

        # Test Pie Chart
        axes[1].pie(test_counts_values,
        #   labels = test_counts_index,
            autopct = '%1.1f%%',
            startangle=90,
            pctdistance=1.1,
            colors=plt.cm.Set3.colors,
            wedgeprops={'linewidth': 1, 
                        'edgecolor': 'white'})
        axes[1].set_title("Test Dataset")

        axes[1].legend(test_counts_index,
                    loc="center left", 
                    bbox_to_anchor=(1, 0.5))

        plt.tight_layout()

        Path("figures/").parent.mkdir(parents=True, exist_ok=True)

        # Figure path for each fold
        fig_path = f"figures/piechart_train_test_split_fold{fold}.png"

        plt.savefig(fig_path)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--adata_path", required=True)
    ap.add_argument("--split_path", required=True)
    args = ap.parse_args()
    main(args.adata_path, args.split_path)


