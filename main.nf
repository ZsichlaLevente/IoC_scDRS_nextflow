#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//// Define processes ////

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

process GWAS_nf1_scDRS_MAGMA {
    container 'magma:1.10'
    publishDir "results/magma", mode: 'copy'

    input:
    tuple path(tsv_file), val(N), path(ncbi37_file), path(g1000_file_bim), path(g1000_file_bed), path(g1000_file_fam), path(g1000_file_synonims)

    output:
    path("*.genes.out"), emit: gene_out

    script:
    def base_name = tsv_file.simpleName

    """
    export base_name=${base_name}

    # process gene location file
    awk 'BEGIN {OFS="\\t"} {print \$6, \$2, \$3, \$4}' ${ncbi37_file} > NCBI37.3.gene_symbol.loc

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

process GWAS_nf2_scDRS_mungeGS {
    publishDir "results/munge_gs", mode: 'copy'

    input:
    path(genes_out_file)

    output:
    path "*.gs", emit: gs_file

    script:
    """
    # Extract the base name of the genes_out_file
    base_name=\$(basename ${genes_out_file} .genes.out)

    # Ensure the base_name is passed into the shell script
    echo "Base name is: \$base_name"

    # Extract GENE and ZSTAT columns (tab-separated) before munge-gs
    awk 'NR==1 {print \$1"\\t"\$8} NR>1 {print \$1"\\t"\$8}' ${genes_out_file} > \${base_name}.zscore.txt

    # Run munge-gs
    scdrs munge-gs \
        --out-file \${base_name}.gs \
        --zscore-file \${base_name}.zscore.txt \
        --weight zscore \
        --n-max 1000
    """
}


//// Execute workflow ////

workflow {
    script_preprocess = file("bin/scRNA_nf1_scDRS_preprocess.py")

    // Define a channel for scRNA inputs
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

    scRNA_nf1_scDRS_preprocess(scRNA_tuples)

    // Define a channel for GWAS inputs
    ncbi37_file = Channel.of(file("data/gene_locations/NCBI37.3.gene.loc"))
    g1000_file_bim = Channel.of(file("data/reference_genomes/g1000_eur.bim"))
    g1000_file_bed = Channel.of(file("data/reference_genomes/g1000_eur.bed"))
    g1000_file_fam = Channel.of(file("data/reference_genomes/g1000_eur.fam"))
    g1000_file_synonims = Channel.of(file("data/reference_genomes/g1000_eur.synonims"))

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

    GWAS_tuples
        .combine(ncbi37_file)
        .combine(g1000_file_bim)
        .combine(g1000_file_bed)
        .combine(g1000_file_fam)
        .combine(g1000_file_synonims)
        | GWAS_nf1_scDRS_MAGMA
        | GWAS_nf2_scDRS_mungeGS
}
