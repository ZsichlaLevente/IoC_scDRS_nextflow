import os
import argparse
import pandas as pd
import scanpy as sc
from gprofiler import GProfiler

def map_gene_names(input_file, species):
    # Load the h5ad file into memory with write access
    ad = sc.read_h5ad(input_file).copy()

    # Only map if species is mmulatta
    if species == "mmulatta":
        print("Mapping gene names from Macaca mulatta to Homo sapiens...")

        # Use gene names from ad.var_names
        gene_names = ad.var_names.tolist()

        # Initialize GProfiler instance
        gp = GProfiler(return_dataframe=True)

        # Map genes from Macaca mulatta to Homo sapiens
        conversion_result = gp.orth(
            organism=species,
            query=gene_names,
            target="hsapiens"
        )

        # Create dictionary of mapped names
        matched_names = {
            row['incoming']: row['ortholog_ensg']
            for _, row in conversion_result.iterrows()
            if pd.notna(row['ortholog_ensg'])
        }

        # Replace gene names
        mapped_names = [matched_names.get(g, g) for g in gene_names]
        ad.var_names = mapped_names

        print("Gene names mapped successfully.")
    else:
        print("No gene name mapping required.")

    # Save the modified file
    output_file = input_file.replace(".h5ad", "_mapped.h5ad")
    ad.write(output_file)

    return output_file

def main():
    parser = argparse.ArgumentParser(description="Map gene names from Macaca mulatta to Homo sapiens.")
    parser.add_argument('--input_file', type=str, required=True, help='Input h5ad file')
    parser.add_argument('--species', type=str, required=True, help='Species name (e.g., mmulatta)')
    args = parser.parse_args()

    output_file = map_gene_names(args.input_file, args.species)
    print(f"Processed gene names saved to: {output_file}")

if __name__ == "__main__":
    main()
