#################
# Input data
#################
# Genotype data
# Hg19 1KG files from o2: /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/
plink_genotype_stem="/lab-share/CHIP-Strober-e2/Public/ben/gxe/simulations/input_data/1000G.EUR.QC.hm3_noMHC."

# Gene tss file
gene_tss_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/gene-ids-and-positions.tsv"

# Protein coding gene file
protein_coding_gene_file="/lab-share/CHIP-Strober-e2/Public/ben/s2e_uncertainty/rsage_testing/input_data/protein_coding_genes.csv"



#######################
# output directories
#######################
output_stem="/lab-share/CHIP-Strober-e2/Public/ben/gxe/simulations/single_gene_variance_simulation/"

simulation_results_dir=${output_stem}"simulation_results/"

visualization_results_dir=${output_stem}"visualization_results/"





#################
# Simulation 1: Total proportional amplification simulation
#################
chrom_num="1"
e_variable_proportion="0.5"
n_proportional_amplication_causal_snps_per_gene="5"
cis_radius="500000"
orig_ge_h2="0.1"

sample_size_arr=("100" "200" "300" "400" "450")
proportional_amplification_scaling_factor_arr=("1.5" "2.0" "4.0" "6.0" "8.0" "10.0")
if false; then
for sample_size in "${sample_size_arr[@]}"
do
for proportional_amplification_scaling_factor in "${proportional_amplification_scaling_factor_arr[@]}"
do
	sbatch total_proportional_amplification_simulation.sh $plink_genotype_stem $chrom_num $e_variable_proportion $proportional_amplification_scaling_factor $n_proportional_amplication_causal_snps_per_gene $orig_ge_h2 $sample_size $cis_radius $simulation_results_dir $gene_tss_file $protein_coding_gene_file
done
done
fi




#################
# Simulation 2: Mixed variance component simulation
#################
chrom_num="1"
e_variable_proportion="0.5"
n_proportional_amplication_causal_snps_per_gene="5"
cis_radius="500000"
orig_ge_h2="0.1"

n_constant_causal_snps_per_gene="5"
constant_ge_h2="0.05"
n_interaction_causal_snps_per_gene="5"
interaction_ge_h2="0.01"


sample_size_arr=("100" "200" "300" "400" "450")
proportional_amplification_scaling_factor_arr=("1.5" "2.0" "4.0" "6.0" "8.0" "10.0")
if false; then
for sample_size in "${sample_size_arr[@]}"
do
for proportional_amplification_scaling_factor in "${proportional_amplification_scaling_factor_arr[@]}"
do
	sbatch mixed_variance_component_simulation.sh $plink_genotype_stem $chrom_num $e_variable_proportion $proportional_amplification_scaling_factor $n_proportional_amplication_causal_snps_per_gene $orig_ge_h2 $sample_size $cis_radius $simulation_results_dir $gene_tss_file $protein_coding_gene_file $n_constant_causal_snps_per_gene $constant_ge_h2 $n_interaction_causal_snps_per_gene $interaction_ge_h2
done
done
fi



#################
# Simulation 3: Interaction effects prior to proportional amplification
#################
chrom_num="1"
e_variable_proportion="0.5"
n_proportional_amplication_causal_snps_per_gene="5"
cis_radius="500000"
orig_ge_h2="0.1"

n_constant_causal_snps_per_gene="5"
constant_ge_h2="0.05"
n_interaction_causal_snps_per_gene="5"
interaction_ge_h2="0.01"


sample_size_arr=("100" "200" "300" "400" "450")
proportional_amplification_scaling_factor_arr=("1.5" "2.0" "4.0" "6.0" "8.0" "10.0")
if false; then
for sample_size in "${sample_size_arr[@]}"
do
for proportional_amplification_scaling_factor in "${proportional_amplification_scaling_factor_arr[@]}"
do
	sbatch pre_scaling_interaction_simulation.sh $plink_genotype_stem $chrom_num $e_variable_proportion $proportional_amplification_scaling_factor $n_proportional_amplication_causal_snps_per_gene $orig_ge_h2 $sample_size $cis_radius $simulation_results_dir $gene_tss_file $protein_coding_gene_file $n_constant_causal_snps_per_gene $constant_ge_h2 $n_interaction_causal_snps_per_gene $interaction_ge_h2
done
done
fi



#################
# Visualize interaction effects
#################
if false; then
Rscript visualize_simulation_results.R $simulation_results_dir $visualization_results_dir
fi






