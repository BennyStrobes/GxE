#!/bin/bash
#SBATCH -t 0-24:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 



simulation_number="${1}"
simulation_type="${2}"
simulated_data_dir="${3}"
plink2_genotype_stem="${4}"
simulation_results_dir="${5}"
pa_h2_code_dir="${6}"

source ~/.bashrc
conda activate PA-h2

# Simulation sample sizes
sample_size_arr=("100" "300")

date
for sample_size in "${sample_size_arr[@]}"; do
    echo $sample_size"_"${simulation_number}


	simulation_name_string="simulation_"${simulation_number}"_sim_type_"${simulation_type}"_default_""N_"${sample_size}
	output_stem=${simulation_results_dir}${simulation_name_string}


	expression_bed_file=${simulated_data_dir}${simulation_name_string}".expression.bed"
	E_var_file=${simulated_data_dir}${simulation_name_string}".interaction_covariates.txt"

	if false; then
	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--output-stem $output_stem

	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--output-stem $output_stem"_perm" \
		--permute True
	fi


	if false; then
	# Run PA-H2 with residualized genotype
	geno_resid_output_stem=${simulation_results_dir}${simulation_name_string}"_geno_resid"
	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--residualize-genotype \
		--output-stem $geno_resid_output_stem

	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--residualize-genotype \
		--output-stem $geno_resid_output_stem"_perm" \
		--permute True
	fi


	# Run PA-H2 without residualized genotype and cross-validation
	cross_val_output_stem=${simulation_results_dir}${simulation_name_string}"_cross_val"
	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--cross-fit-folds "4" \
		--output-stem $cross_val_output_stem

	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--cross-fit-folds "4" \
		--output-stem $cross_val_output_stem"_perm" \
		--permute True

	# Run PA-H2 with residualized genotype and cross-validation
	geno_resid_cross_val_output_stem=${simulation_results_dir}${simulation_name_string}"_geno_resid_cross_val"
	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--residualize-genotype \
		--cross-fit-folds "4" \
		--output-stem $geno_resid_cross_val_output_stem

	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--residualize-genotype \
		--cross-fit-folds "4" \
		--output-stem $geno_resid_cross_val_output_stem"_perm" \
		--permute True


done
date