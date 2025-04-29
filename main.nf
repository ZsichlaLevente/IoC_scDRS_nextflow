#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//// Define processes ////

process scRNA_nf1_preprocess {
    container = 'marcoschinas/scdrs:1.0.2'
    //publishDir "results/scRNA", mode: 'copy'

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

process GWAS_nf1_MAGMA {
    container 'magma:1.10'
    //publishDir "results/magma", mode: 'copy'

    input:
    tuple path(tsv_file), val(N), path(ncbi37_file), path(g1000_file_zip)

    output:
    path("*.genes.out"), emit: gene_out

    script:
    def base_name = tsv_file.simpleName

    """
    export base_name=${base_name}

    # process gene location file
    awk 'BEGIN {OFS="\\t"} {print \$6, \$2, \$3, \$4}' ${ncbi37_file} > NCBI37.3.gene_symbol.loc

    # unzip reference genome file
    unzip ${g1000_file_zip}

    # Annotate SNPs to genes
    magma --annotate \
        --snp-loc ${tsv_file} \
        --gene-loc NCBI37.3.gene_symbol.loc \
        --out \$base_name

    # Extract the first and fourth columns from the tsv_file and save to a new file
    cut -f1,4 ${tsv_file} > ${base_name}.pval.txt

    # Gene-based analysis
    magma \
        --bfile g1000_eur \
        --pval ${base_name}.pval.txt N=${N} \
        --gene-annot ${base_name}.genes.annot \
        --out \$base_name
    """
}

process GWAS_nf2_combineMAGMA {
    container = 'r-base:4.5.0'
    //publishDir "results/magma", mode: 'copy'

    input:
    path(gene_out_files)

    output:
    path("GWAS_allPheno_zscores.tsv"), emit: combined_zscores

    script:
    """
    # Create a temp directory
    mkdir temp_zscores

    # Loop over each .genes.out file in the working directory
    for file in *.genes.out; do
        base_name=\$(basename \$file .genes.out)
        phenotype=\$(echo \$base_name | sed 's/^[^_]*_//')
        [ -z "\$phenotype" ] && phenotype=\$base_name

        awk 'BEGIN {print "GENE\\tZSTAT"} NR>1 {print \$1"\\t"\$8}' \$file > temp_zscores/\${phenotype}.txt
    done

    # Merge all phenotype zscore files by GENE using full join in R
    Rscript -e '
        files <- list.files("temp_zscores", full.names = TRUE)
        if (length(files) == 0) stop("No files found in temp_zscores")

        # Read and merge files in a loop
        combined <- NULL
        for (file in files) {
            data <- read.table(file, header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
            phenotype <- sub("\\\\.txt\$", "", basename(file))  # Corrected escape sequences
            colnames(data)[2] <- phenotype
            if (is.null(combined)) {
                combined <- data
            } else {
                combined <- merge(combined, data, by = "GENE", all = TRUE)
            }
        }

        # Write the result
        write.table(combined, "GWAS_allPheno_zscores.tsv", sep = "\\t", row.names = FALSE, quote = FALSE)
    '
    """
}

process GWAS_nf3_mungeGS {
    container = 'marcoschinas/scdrs:1.0.2'
    publishDir "results/munge_gs", mode: 'copy'

    input:
    path(zscore_file)

    output:
    path "*.gs", emit: gs_file

    script:
    """
    scdrs munge-gs \
        --out-file ${zscore_file}.gs \
        --zscore-file ${zscore_file} \
        --weight zscore \
        --n-max 1000
    """
}

process scDRS_nf1_computeScores {
    container = 'marcoschinas/scdrs:1.0.2'

    input:
    tuple path(scRNA_file), path(covariate_file), path(gs_file)

    output:
    path("${scRNA_file.simpleName}/*"), emit: scDRS_scores_files
    path("${scRNA_file.simpleName}"), emit: scDRS_scores_folder

    script:
    """
    mkdir ${scRNA_file.simpleName}
    scdrs compute-score \
        --h5ad-file ${scRNA_file} \
        --h5ad-species mouse \
        --gs-file ${gs_file} \
        --gs-species human \
        --cov-file ${covariate_file} \
        --out-folder ${scRNA_file.simpleName}
    """
}

process scDRS_nf2_performDownstream {
    container = 'marcoschinas/scdrs:1.0.2'
    publishDir "results/scDRS", mode: 'copy'

    input:
    tuple path(scDRS_scores_folder), path(scRNA_file), val(sample_id), path(input_h5ad), val(batch_column), val(cell_type_column), path(script_preprocess)

    output:
    path("${scRNA_file.simpleName}/*"), emit: scDRS_downstream_folder

    script:
    """
    scdrs perform-downstream \
        --h5ad-file ${scRNA_file} \
        --score-file ${scDRS_scores_folder}/@.full_score.gz \
        --out-folder ${scRNA_file.simpleName} \
        --group-analysis ${cell_type_column} \
        --gene-analysis\
    """
}

//// Execute workflow ////

workflow {
    script_preprocess = file("bin/scRNA_nf1_scDRS_preprocess.py")

    // scRNA preprocessing
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
        .set { scRNA_tuples }

    processed_scRNA = scRNA_nf1_preprocess(scRNA_tuples)

    // GWAS preprocessing
    ncbi37_file = Channel.of(file("data/gene_locations/NCBI37.3.gene.loc"))
    g1000_file_zip = Channel.of(file("data/reference_genomes/g1000_eur.zip"))

    Channel
        .fromPath("data/GWAS_samples.tsv")
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                file(row.input_file),
                row.N as Integer
            )
        }
        .set { GWAS_tuples }

    magma_results = GWAS_tuples
        .combine(ncbi37_file)
        .combine(g1000_file_zip)
        | GWAS_nf1_MAGMA

    magma_results
        .collect()
        .set { all_gene_out_files }

    // scDRS
    gs_file = GWAS_nf2_combineMAGMA(all_gene_out_files) | GWAS_nf3_mungeGS

    scDRS_score_results = processed_scRNA.processed_h5ad
        .combine(processed_scRNA.covariate_tsv)
        .combine(gs_file)
        | scDRS_nf1_computeScores

    scDRS_score_results.scDRS_scores_folder
        .combine(processed_scRNA.processed_h5ad)
        .combine(scRNA_tuples)
        | scDRS_nf2_performDownstream

}
