#######################
# Input Data
########################

gtex_v8_tpm_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
gtex_v8_gene_reads_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"

gtex_v8_eqtl_expression_matrices_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/GTEx_Analysis_v8_eQTL_expression_matrices/"

xcell_ct_proportions_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz"

gtex_sample_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

gtex_covariate_dir="/lab-share/CHIP-Strober-e2/Public/GTEx/expression/v8_per_tissue_expression/covariates/GTEx_Analysis_v8_eQTL_covariates/"






########################
# Output Data
########################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gxe/exploratory_real_data_analysis/"

processed_expression_dir="${output_root}/processed_expression/"



########################
# Code
########################

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
    ${tissue_name}
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

source ~/.bashrc
conda activate plink_env
# Visualize the per-gene variance stratified by the E-variable, comparing across normalization methods
per_gene_variance_plot_stem="${processed_expression_dir}/${tissue_name}.${cell_type}.per_gene_variance_method_comparison"
Rscript visualize_per_gene_variance.R \
    ${processed_expression_dir} \
    ${tissue_name} \
    ${cell_type} \
    "$(echo ${normalization_methods} | tr ' ' ',')" \
    ${per_gene_variance_plot_stem}






