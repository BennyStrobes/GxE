#!/bin/bash
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=20GB


expression_bed_file="$1"
plink2_genotype_stem="$2"
E_var_file="${3}"
covariate_file="${4}"
variant_gene_pairs_file="${5}"
output_stem="${6}"
PA_mixture_model_code_dir="${7}"

source ~/.bashrc
conda activate PA-h2

# Stop at the first failing step (otherwise Step 2 runs on missing input files)
set -e


###################
# Step 1: Generate mixture model input
# (if variant_gene_pairs_file is none, pairs are the top variant per gene with a significant cis interaction QTL)
###################
if false; then
python ${PA_mixture_model_code_dir}make_pa_mixture_model_input_from_data.py \
	--expression-bed $expression_bed_file \
	--plink2-per-chrom-stem $plink2_genotype_stem \
	--binary-E-interaction-covariate-file $E_var_file \
	--covariate-file $covariate_file \
	--variant-gene-pairs-file $variant_gene_pairs_file \
	--output-stem $output_stem
fi

###################
# Step 2: Run mixture model inference
###################
# Genotype-residualized S (recommended)
if false; then
python ${PA_mixture_model_code_dir}run_pa_mixture_model_inference.py \
	--mixture-model-input-file ${output_stem}"_pa_mixture_model_geno_resid_S_input.txt" \
	--output-stem ${output_stem}"_geno_resid_mixture_model_results" \
	--n-warmup-iters 100
fi
# Raw S (for comparison)
if false; then
python ${PA_mixture_model_code_dir}run_pa_mixture_model_inference.py \
	--mixture-model-input-file ${output_stem}"_pa_mixture_model_input.txt" \
	--output-stem ${output_stem}"_mixture_model_results" \
	--n-warmup-iters 100
fi

# Standard Gibbs cross-check (genotype-residualized S)
python ${PA_mixture_model_code_dir}run_pa_mixture_model_inference_standard_gibbs.py \
	--mixture-model-input-file ${output_stem}"_pa_mixture_model_geno_resid_S_input.txt" \
	--output-stem ${output_stem}"_standard_gibbs_geno_resid_mixture_model_results"
fi