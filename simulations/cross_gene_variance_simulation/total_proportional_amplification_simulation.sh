#!/bin/bash
#SBATCH -t 0-2:30                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB 



plink_genotype_stem="${1}"
chrom_num="${2}"
e_variable_proportion="${3}"
proportional_amplification_scaling_factor="${4}"
n_proportional_amplication_causal_snps_per_gene="${5}"
orig_ge_h2="${6}"
sample_size="${7}"
cis_radius="${8}"
simulation_results_dir="${9}"
gene_tss_file="${10}"
protein_coding_gene_file="${11}"

source ~/.bashrc
conda activate plink_env



python total_proportional_amplification_simulation.py $plink_genotype_stem $chrom_num $e_variable_proportion $proportional_amplification_scaling_factor $n_proportional_amplication_causal_snps_per_gene $orig_ge_h2 $sample_size $cis_radius $simulation_results_dir $gene_tss_file $protein_coding_gene_file
