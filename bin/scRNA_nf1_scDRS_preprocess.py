import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend (no GUI)
import scdrs
import copy
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
import os
import sys

def export_processed_data(adata, output_file="processed_data.h5ad"):
    # Save the processed data
    adata.write_h5ad(output_file)
    df_cov = pd.DataFrame(index=adata.obs.index)
    df_cov["const"] = 1
    df_cov["n_genes"] = adata.obs["n_genes"]
    df_cov.to_csv(output_file[:-18]+"_covariates.tsv", sep="\t")
    # Clear intermediate variables
    del df_cov  # Free memory
    print(f"Data preprocessing complete. Processed data saved to {output_file}.")

def sanitize_anndata(adata):
    adata_sanitized = adata.copy()

    # Clean up raw attribute if present
    if adata_sanitized.raw is not None:
        raw_clean = adata_sanitized.raw.to_adata()
        if "_index" in raw_clean.var.columns:
            raw_clean.var = raw_clean.var.drop(columns=["_index"])
        adata_sanitized.raw = raw_clean

    def convert_series_in_dict(d):
        for key in list(d.keys()):
            if isinstance(d[key], pd.Series):
                d[key] = d[key].to_dict()
            elif isinstance(d[key], dict):
                convert_series_in_dict(d[key])
        return d

    # Sanitize .uns by converting pandas.Series to dict
    if isinstance(adata_sanitized.uns, dict):
        adata_sanitized.uns = convert_series_in_dict(adata_sanitized.uns)

    return adata_sanitized

def restore_series_in_dict(d):
    for key in list(d.keys()):
        if isinstance(d[key], dict):
            # Recursively process nested dicts first
            d[key] = restore_series_in_dict(d[key])
            try:
                # Check if it's likely a Series (heuristic: all values are scalar and not dicts/lists)
                if all(isinstance(v, (int, float, str, bool, type(None), np.number)) for v in d[key].values()):
                    d[key] = pd.Series(d[key])
            except Exception:
                pass  # Leave as dict if conversion fails
    return d


######### Execute script #########
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--batch_column', required=True)
    parser.add_argument('--cell_type_column', required=True)
    parser.add_argument('--output_dir', required=True)
    args = parser.parse_args()

    input_file = args.input_file
    batch_column = args.batch_column
    cell_type_column = args.cell_type_column

    output_prefix = os.path.basename(input_file).replace(".h5ad", "")
    adata = scdrs.util.load_h5ad(h5ad_file=input_file, flag_filter_data=True, flag_raw_count=True)
    adata_cov = pd.DataFrame(index=adata.obs.index)
    adata_cov["n_genes"] = adata.obs["n_genes"]
    if batch_column in adata.obs and adata.obs[batch_column].nunique() > 1:
        adata_cov["animal_ID"] = adata.obs[batch_column]

    adata_processed = scdrs.preprocess(
        data=adata,
        cov=adata_cov,
        adj_prop=cell_type_column,
        n_mean_bin=20,
        n_var_bin=20,
        n_chunk=None,
        copy=True
    )
    adata_sanitized = sanitize_anndata(adata_processed)
    output_file = os.path.join(args.output_dir, output_prefix + "_preprocessed.h5ad")
    export_processed_data(adata_sanitized, output_file=output_file)
