#!/bin/bash
#SBATCH -t 0-3:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=3GB 



simulation_number="${1}"
sample_size_arr="${2}"
simulation_type="${3}"
simulated_data_dir="${4}"
gene_tss_file="${5}"
protein_coding_gene_file="${6}"
plink2_genotype_stem="${7}"

source ~/.bashrc
conda activate PA-h2


sample_size="450"



simulation_name_string="simulation_"${simulation_number}"_sim_type_"${simulation_type}"_default_""N_"${sample_size}
output_stem=${simulated_data_dir}${simulation_name_string}

python simulate_eqtl_data.py \
  --sample-size "${sample_size}" \
  --output-stem "${output_stem}" \
  --gene-tss-file "${gene_tss_file}" \
  --protein-coding-gene-file "${protein_coding_gene_file}" \
  --simulation-type "${simulation_type}" \
  --plink2-genotype-stem "${plink2_genotype_stem}"