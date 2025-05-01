#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//// Define processes ////

process scRNA_nf1_preprocess {
    container = 'scdrs-python:1.0.2'
    //publishDir "results/scRNA", mode: 'copy'

    input:
    tuple val(sample_id), path(input_h5ad), val(batch_column), val(cell_type_column), val(species), path(script_preprocess)

    output:
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(species),  path("*.h5ad") , path("*.tsv")

    script:
    """
    python3 $script_preprocess \
        --input_file $input_h5ad \
        --batch_column $batch_column \
        --cell_type_column $cell_type_column \
        --output_dir .
    """
}

process scRNA_nf2_orthomap {
    container = 'scdrs-python:1.0.2'
    publishDir "results/processed_scRNA", mode: 'copy'

    input:
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(species), path(processed_h5ad) , path(covariate_file), path(script_orthomap)

    output:
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(species), path(covariate_file), path("*_mapped.h5ad")

    script:
    """
    python3 ${script_orthomap} --input_file ${processed_h5ad} --species ${species}
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
    # Check if running inside a Docker container
    if [[ -f /.dockerenv ]]; then
        echo "Running inside Docker container"
    else
        echo "Not running inside Docker container"
    fi

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
    container = 'scdrs-python:1.0.2'
    publishDir "results/processed_geneset", mode: 'copy'

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
    container = 'scdrs-python:1.0.2'

    input:
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(h5ad_species), path(covariate_file), path(scRNA_file), path(gs_file)

    output:
    path("${scRNA_file.simpleName}"), emit: scDRS_scores_folder
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(h5ad_species), path(covariate_file), path(scRNA_file), path(gs_file), path("${scRNA_file.simpleName}"), emit: scDRS_scores_passon

    script:
    def effective_species = h5ad_species == 'mmulatta' ? 'hsapiens' : h5ad_species
    """
    mkdir ${scRNA_file.simpleName}
    scdrs compute-score \
        --h5ad-file ${scRNA_file} \
        --h5ad-species ${effective_species} \
        --gs-file ${gs_file} \
        --gs-species human \
        --cov-file ${covariate_file} \
        --out-folder ${scRNA_file.simpleName}
    """
}

process scDRS_nf2_performDownstream {
    container = 'scdrs-python:1.0.2'
    //publishDir "results/scDRS_results", mode: 'copy'

    input:
    tuple val(sample_id), val(batch_column), val(cell_type_column), val(h5ad_species), path(covariate_file), path(scRNA_file), path(gs_file), path(scDRS_scores_folder)

    output:
    path("${scRNA_file.simpleName}"), emit: scDRS_downstream_folder
    output:tuple path("${scRNA_file.simpleName}"), val(cell_type_column), path(scRNA_file), emit: scDRS_results_gather


    script:
    """
    scdrs perform-downstream \
        --h5ad-file ${scRNA_file} \
        --score-file ${scDRS_scores_folder}/@.full_score.gz \
        --out-folder ${scRNA_file.simpleName} \
        --group-analysis ${cell_type_column} \
        --gene-analysis
    """
}


process scDRS_nf3_gatherResults {
    container = 'scdrs-python:1.0.2'
    publishDir "results/scDRS_results_gathered", mode: 'copy'

    input:
    tuple path(scDRS_results_folder), val(cell_type_column), path(scRNA_file), path(script_gather)

    output:
    path("${scRNA_file.simpleName}_final.h5ad"), emit: scRNA_annotated_final
    path("*_gene_final.txt"), emit: geneAnalysis_gathered
    path("*_group_final.txt"), emit: groupAnalysis_gathered

    script:
    """
    python3 ${script_gather} \
        --input_folder ${scDRS_results_folder} \
        --scRNA_input ${scRNA_file} \
        --groupAnalysis_column ${cell_type_column}
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
                row.species,
                script_preprocess
            )
        }
        .set { scRNA_tuples }
    scRNA_preprocessed_mapped = scRNA_nf1_preprocess(scRNA_tuples)
        .combine(Channel.of(file("bin/scRNA_nf2_scDRS_orthomap.py")))
        | scRNA_nf2_orthomap

    // GWAS preprocessing
    if ( params.GWAS_panel ) {
        gs_file = Channel.of(file(params.GWASfile))
    }
    else {
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

        // scDRS with our prepared GWAS
        gs_file = GWAS_nf2_combineMAGMA(all_gene_out_files) | GWAS_nf3_mungeGS
    }

    // perform scDRS analysis
    scDRS_scores = scRNA_preprocessed_mapped
        .combine(gs_file)
        | scDRS_nf1_computeScores

    scDRS_downstream = scDRS_scores.scDRS_scores_passon
        | scDRS_nf2_performDownstream

    // gather and export results
    script_gather = Channel.of(file("bin/scDRS_nf1_scDRS_gather.py"))
    scDRS_downstream.scDRS_results_gather
        .combine(script_gather)
        | scDRS_nf3_gatherResults

}
