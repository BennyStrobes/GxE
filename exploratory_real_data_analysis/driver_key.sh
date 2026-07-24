#######################
# Input Data
########################

gtex_v8_tpm_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
gtex_v8_gene_reads_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"

gtex_v8_eqtl_expression_matrices_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_v8_eQTL_expression_matrices/"

xcell_ct_proportions_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz"

gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

gtex_covariate_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_eQTL_covariates/"

gtex_input_genotype_dir="/lab-share/CHIP-Strober-e2/Public/ben/process_gtex_genotype_data/processed_genotype/"

PA_H2_code_dir="/lab-share/CHIP-Strober-e2/Public/ben/gxe/PA-h2/"




########################
# Output Data
########################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gxe/exploratory_real_data_analysis/"

processed_expression_dir="${output_root}processed_expression/"

processed_genotype_dir=${output_root}"processed_genotype/"

pa_h2_results_dir=${output_root}"pa_h2_results/"



########################
# Code
########################


##############################
# Quick reprocessing of genotype data
if false; then
source ~/.bashrc
conda activate plink_env
python quick_reprocessing_of_genotype_data.py $gtex_input_genotype_dir $processed_genotype_dir
fi

echo "CONSIDER FILTERING TO EUROPEAN ANCESTRY INDIVIDUALS"
##############################
# Preprocess expression data (with various normalizations and filtering))
tissue_name="Whole_Blood"
if false; then
sh preprocess_expression.sh \
    ${gtex_v8_gene_reads_file} \
    ${gtex_v8_eqtl_expression_matrices_dir} \
    ${xcell_ct_proportions_file} \
    ${gtex_sample_attributes_file} \
    ${gtex_covariate_dir} \
    ${processed_expression_dir} \
    ${tissue_name} \
    ${processed_genotype_dir}"gtex_v9_eqtl_chr1.psam"
fi


###############################
# Extract per-gene variance stratified by E-variable
tissue_name="Whole_Blood"
cell_type="Neutrophils"
tissue_xcell_ct_proportions_file="${processed_expression_dir}/${tissue_name}.xcell_ct_proportions.txt.gz"
normalization_methods="log_tmm log_tmm.unstandardized inverse_normal_transform"
if false; then
source ~/.bashrc
conda activate plink_env
for normalization_method in ${normalization_methods}; do
    tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
    per_gene_variance_output_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.${cell_type}.per_gene_variance.txt"
    python compute_per_gene_variance_stratified_by_e_variable.py \
        --expression_matrix_file ${tissue_expression_matrix_file} \
        --xcell_ct_proportions_file ${tissue_xcell_ct_proportions_file} \
        --tissue_name ${tissue_name} \
        --cell_type ${cell_type} \
        --per_gene_variance_output_file ${per_gene_variance_output_file}
done
fi



###############################
# Run PA-H2 regression
tissue_name="Whole_Blood"
cell_type="neutrophils"
normalization_method="log_tmm"

if false; then
# Files
tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type
sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir 

tissue_name="Whole_Blood"
cell_type="neutrophils"
normalization_method="inverse_normal_transform"

# Files
tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type
sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir 
fi




if false; then
source ~/.bashrc
conda activate plink_env
# Visualize the per-gene variance stratified by the E-variable, comparing across normalization methods
per_gene_variance_plot_stem="${processed_expression_dir}/${tissue_name}.${cell_type}.per_gene_variance_method_comparison"
Rscript visualize_per_gene_variance.R \
    ${processed_expression_dir} \
    ${tissue_name} \
    ${cell_type} \
    "log_tmm.unstandardized" \
    ${per_gene_variance_plot_stem}
fi