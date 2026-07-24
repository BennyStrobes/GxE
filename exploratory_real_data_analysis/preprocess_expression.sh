#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB


gtex_v8_gene_reads_file="$1"
gtex_v8_eqtl_expression_matrices_dir="$2"
xcell_ct_proportions_file="$3"
gtex_sample_attributes_file="$4"
gtex_covariate_dir="$5"
processed_expression_dir="$6"
tissue_name="$7"
genotype_psam_sample_file="$8"


python preprocess_expression.py \
    --gtex_v8_gene_reads_file ${gtex_v8_gene_reads_file} \
    --gtex_v8_eqtl_expression_matrices_dir ${gtex_v8_eqtl_expression_matrices_dir} \
    --xcell_ct_proportions_file ${xcell_ct_proportions_file} \
    --gtex_sample_attributes_file ${gtex_sample_attributes_file} \
    --gtex_covariate_dir ${gtex_covariate_dir} \
    --processed_expression_dir ${processed_expression_dir} \
    --tissue_name ${tissue_name} \
    --valid_individuals_file ${genotype_psam_sample_file}

