#!/bin/bash
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=3GB 



simulation_number="${1}"
simulation_type="${2}"
simulated_data_dir="${3}"
gene_tss_file="${4}"
protein_coding_gene_file="${5}"
plink2_genotype_stem="${6}"

source ~/.bashrc
conda activate PA-h2



# Simulation sample sizes
sample_size_arr=("100" "300" "450")

for sample_size in "${sample_size_arr[@]}"; do
  echo $sample_size"_"${simulation_number}

  simulation_name_string="simulation_"${simulation_number}"_sim_type_"${simulation_type}"_default_""N_"${sample_size}
  output_stem=${simulated_data_dir}${simulation_name_string}

  python simulate_eqtl_data.py \
    --sample-size "${sample_size}" \
    --output-stem "${output_stem}" \
    --gene-tss-file "${gene_tss_file}" \
    --protein-coding-gene-file "${protein_coding_gene_file}" \
    --simulation-type "${simulation_type}" \
    --plink2-genotype-stem "${plink2_genotype_stem}" \
    --random-seed "${simulation_number}"

done




