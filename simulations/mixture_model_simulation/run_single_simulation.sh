#!/bin/bash
#SBATCH -t 0-0:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB


simulation_iter="${1}"
eqtl_sample_size="${2}"
n_tests="${3}"
e_variable_ratio="${4}"
simulated_data_dir="${5}"
simulation_results_dir="${6}"



source ~/.bashrc
conda activate plink_env


# Loop through prop_PA_tests values
prop_PA_tests_arr=( "0.2" "0.4" "0.6" "0.8" )

for prop_PA_tests in "${prop_PA_tests_arr[@]}"; do

simulation_name_stem="simulation_iter_"${simulation_iter}"_eqtl_sample_size_"${eqtl_sample_size}"_n_tests_"${n_tests}"_e_variable_ratio_"${e_variable_ratio}"_prop_PA_tests_"${prop_PA_tests}


###################
# First generate simulated data
###################
# Files to be generated
simulation_oracle_file=${simulated_data_dir}${simulation_name_stem}"_oracle.txt"
pa_mixture_model_input_file=${simulated_data_dir}${simulation_name_stem}"_pa_mixture_model_input.txt"
pa_mixture_model_geno_resid_S_input_file=${simulated_data_dir}${simulation_name_stem}"_pa_mixture_model_geno_resid_S_input.txt"
python simulate_mixture_model_data.py $simulation_iter $eqtl_sample_size $n_tests $e_variable_ratio $prop_PA_tests $simulation_oracle_file $pa_mixture_model_input_file $pa_mixture_model_geno_resid_S_input_file


###################
# Second, run mixture model inference on the simulated data
###################
##
# COmplex inference

output_mixture_model_results_stem=${simulation_results_dir}${simulation_name_stem}"_mixture_model_results"
python run_pa_mixture_model_inference.py \
	--mixture-model-input-file $pa_mixture_model_input_file \
	--output-stem $output_mixture_model_results_stem

output_mixture_model_results_stem=${simulation_results_dir}${simulation_name_stem}"_geno_resid_mixture_model_results"
python run_pa_mixture_model_inference.py \
	--mixture-model-input-file $pa_mixture_model_geno_resid_S_input_file \
	--output-stem $output_mixture_model_results_stem


##
# Standard Gibbs inference
output_mixture_model_results_stem=${simulation_results_dir}${simulation_name_stem}"_standard_gibbs_mixture_model_results"
python run_pa_mixture_model_inference_standard_gibbs.py \
	--mixture-model-input-file $pa_mixture_model_input_file \
	--output-stem $output_mixture_model_results_stem

output_mixture_model_results_stem=${simulation_results_dir}${simulation_name_stem}"_standard_gibbs_geno_resid_mixture_model_results"
python run_pa_mixture_model_inference_standard_gibbs.py \
	--mixture-model-input-file $pa_mixture_model_geno_resid_S_input_file \
	--output-stem $output_mixture_model_results_stem

done