#!/bin/bash
#SBATCH -t 0-6:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=50GB


expression_bed_file="$1"
plink2_genotype_stem="$2"
E_var_file="${3}"
output_stem="${4}"
PA_H2_code_dir="${5}"
gene_category_file="${6}"
covariate_file="${7}"

source ~/.bashrc
conda activate PA-h2


# Default run with variance component for covariates
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--output-stem $output_stem"_with_covariates"

# Default run
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--output-stem $output_stem 


# Permuted run
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--output-stem $output_stem"_permuted" \
    --permute "True"


