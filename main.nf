#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process scRNA_nf1_scDRS_preprocess {
    publishDir "results/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_h5ad), val(batch_column), val(cell_type_column), path(script_preprocess)

    output:
    path "*.h5ad", emit: processed_h5ad
    path "*.tsv",  emit: covariate_tsv

    script:
    """
    python3 $script_preprocess \
        --input_file $input_h5ad \
        --batch_column $batch_column \
        --cell_type_column $cell_type_column \
        --output_dir .
    """
}

workflow {
    script_preprocess = file("bin/scRNA_nf1_scDRS_preprocess.py")

    Channel
        .fromPath("data/scRNA_samples.tsv")
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id,
                file(row.input_file),
                row.batch_column,
                row.cell_type_column,
                script_preprocess
            )
        }
        .set { sample_tuples }

    scRNA_nf1_scDRS_preprocess(sample_tuples)
}

