#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import scanpy as sc
import pandas as pd
import glob

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", required=True, help="Path to folder containing scDRS result files (*.full_score.gz, *.scdrs_group.*, *.scdrs_gene)")
    parser.add_argument("--scRNA_input", required=True, help="Annotated .h5ad file")
    parser.add_argument("--groupAnalysis_column", required=True, help="Column name in .obs to group by for downstream stats")
    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    if not input_folder.exists():
        raise FileNotFoundError(f"Input folder not found: {input_folder}")

    # Load AnnData
    adata = sc.read_h5ad(args.scRNA_input)

    # Load all .full_score.gz files and assign norm_score to adata.obs
    trait_scores = {}
    score_files = sorted(input_folder.glob("*.full_score.gz"))
    if not score_files:
        raise FileNotFoundError(f"No .full_score.gz files found in {input_folder}")

    for f in score_files:
        trait = f.name.split(".")[0]
        df_score = pd.read_csv(f, sep="\t", index_col=0)
        trait_scores[trait] = df_score
        adata.obs[trait] = df_score.loc[adata.obs_names, "norm_score"]

    # Save the annotated adata
    out_adata_file = (Path(args.scRNA_input).stem + "_final.h5ad")
    adata.write(out_adata_file)

    # Group-level statistics
    dict_df_stats = {}
    group_pattern = f"*.scdrs_group.{args.groupAnalysis_column}"
    group_files = sorted(input_folder.glob(group_pattern))
    for f in group_files:
        trait = f.name.split(".")[0]
        df_stats = pd.read_csv(f, sep="\t", index_col=0)
        dict_df_stats[trait] = df_stats

    if dict_df_stats:
        pd.concat(dict_df_stats).to_csv(f"{input_folder.name}_group_final.txt", sep="\t")

    # Gene-level analysis files
    gene_files = sorted(input_folder.glob("*.scdrs_gene"))
    if gene_files:
        df_gene_all = pd.concat(
            [pd.read_csv(f, sep="\t").assign(trait=f.name.split(".")[0]) for f in gene_files]
        )
        df_gene_all.to_csv(f"{input_folder.name}_gene_final.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()
