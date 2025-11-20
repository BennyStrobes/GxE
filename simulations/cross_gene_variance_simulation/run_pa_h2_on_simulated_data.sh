#!/bin/bash
#SBATCH -t 0-9:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=3GB 



simulation_number="${1}"
simulation_type="${2}"
simulated_data_dir="${3}"
plink2_genotype_stem="${4}"
simulation_results_dir="${5}"
pa_h2_code_dir="${6}"

source ~/.bashrc
conda activate PA-h2

# Simulation sample sizes
sample_size_arr=("100" "300" "450")

for sample_size in "${sample_size_arr[@]}"; do
    echo $sample_size"_"${simulation_number}


	simulation_name_string="simulation_"${simulation_number}"_sim_type_"${simulation_type}"_default_""N_"${sample_size}
	output_stem=${simulation_results_dir}${simulation_name_string}


	expression_bed_file=${simulated_data_dir}${simulation_name_string}".expression.bed"
	E_var_file=${simulated_data_dir}${simulation_name_string}".interaction_covariates.txt"

	python ${pa_h2_code_dir}PA_h2.py \
		--expression-bed $expression_bed_file \
		--plink2-per-chrom-stem $plink2_genotype_stem \
		--binary-E-interaction-covariate-file $E_var_file \
		--output-stem $output_stem

done