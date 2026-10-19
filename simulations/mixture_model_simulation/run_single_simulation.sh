#!/bin/bash
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB


simulation_iter="${1}"
eqtl_sample_size="${2}"
n_tests="${3}"
e_variable_ratio="${4}"
simulated_data_dir="${5}"
simulation_results_dir="${6}"



source ~/.bashrc
conda activate PA-h2


# Loop through prop_PA_tests values
prop_PA_tests_arr=( "0.2" "0.4" "0.6" "0.8" )

for prop_PA_tests in "${prop_PA_tests_arr[@]}"; do

simulation_name_stem="simulation_iter_"${simulation_iter}"_eqtl_sample_size_"${eqtl_sample_size}"_n_tests_"${n_tests}"_e_variable_ratio_"${e_variable_ratio}"_prop_PA_tests_"${prop_PA_tests}


###################
# First generate simulated data
# (in same format as PA_h2.py input: expression bed, plink2 genotype, binary E file; plus variant-gene pairs file)
###################
simulation_oracle_file=${simulated_data_dir}${simulation_name_stem}"_oracle.txt"
simulated_data_stem=${simulated_data_dir}${simulation_name_stem}
python simulate_mixture_model_data.py $simulation_iter $eqtl_sample_size $n_tests $e_variable_ratio $prop_PA_tests $simulation_oracle_file $simulated_data_stem


###################
# Second, generate mixture model input from simulated data (same script as used on real data)
###################
python make_pa_mixture_model_input_from_data.py \
	--expression-bed ${simulated_data_stem}"_expression.bed.gz" \
	--plink2-per-chrom-stem ${simulated_data_stem}"_genotype_chr" \
	--binary-E-interaction-covariate-file ${simulated_data_stem}"_E_binary.txt" \
	--covariate-file none \
	--variant-gene-pairs-file ${simulated_data_stem}"_variant_gene_pairs.txt" \
	--output-stem $simulated_data_stem
pa_mixture_model_input_file=${simulated_data_stem}"_pa_mixture_model_input.txt"
pa_mixture_model_geno_resid_S_input_file=${simulated_data_stem}"_pa_mixture_model_geno_resid_S_input.txt"

# Delete raw simulated data to save disk space (re-generatable from simulation_iter seed)
rm ${simulated_data_stem}"_expression.bed.gz" ${simulated_data_stem}"_genotype_chr1."{pgen,pvar,psam} ${simulated_data_stem}"_E_binary.txt" ${simulated_data_stem}"_variant_gene_pairs.txt"


###################
# Third, run mixture model inference on the simulated data
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