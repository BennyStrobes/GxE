#!/bin/bash
#SBATCH -t 0-20:00                         # Runtime in D-HH:MM format
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
if false; then
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--output-stem $output_stem

# Permuted run
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--output-stem $output_stem"_permuted" \
    --permute "True"
fi

if false; then
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--residualize-genotype \
	--cross-fit-folds "3" \
	--output-stem $output_stem"_geno_resid_cross_val"

# Permuted run
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--residualize-genotype \
	--cross-fit-folds "3" \
	--output-stem $output_stem"_geno_resid_cross_val_permuted" \
    --permute "True"
fi

python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--cross-fit-folds "3" \
	--output-stem $output_stem"_cross_val"

# Permuted run
python ${PA_H2_code_dir}PA_h2.py \
	--expression-bed $expression_bed_file \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--gene-list $gene_category_file \
	--cross-fit-folds "3" \
	--output-stem $output_stem"_cross_val_permuted" \
    --permute "True"


if false; then
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
fi
