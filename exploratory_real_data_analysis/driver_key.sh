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

# Gtex subject attributes (has ancestry)
gtex_subject_attributes_file="/lab-share/CHIP-Strober-e2/Public/GTEx/genotype_dbgap_download/phs000424.v11.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"





########################
# Output Data
########################
output_root="/lab-share/CHIP-Strober-e2/Public/ben/gxe/exploratory_real_data_analysis/"

processed_expression_dir="${output_root}processed_expression/"

downsampling_processed_expression_dir="${output_root}downsampling_processed_expression/"

processed_genotype_dir=${output_root}"processed_genotype/"

gene_categories_dir=${output_root}"gene_categories/"

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
    ${gtex_subject_attributes_file} \
    ${processed_genotype_dir}"gtex_v9_eqtl_chr1.psam"
fi

###############################
# Run exploratory analysis where we down-sample expression by randomly sampling read-counts
tissue_read_count_file=${processed_expression_dir}${tissue_name}".gene_reads.filtered.txt.gz"
if false; then
for downsampling_percentage in "0.2" "0.4" "0.6" "0.8" "1.0"; do
    sh run_downsampling_expression_quantification.sh $tissue_read_count_file $downsampling_percentage $downsampling_processed_expression_dir${tissue_name}
done
fi


###############################
# Extract per-gene variance stratified by E-variable
tissue_name="Whole_Blood"
cell_type="Neutrophils"
tissue_xcell_ct_proportions_file="${processed_expression_dir}/${tissue_name}.xcell_ct_proportions.txt.gz"
normalization_methods="log_tmm log_tmm.unstandardized inverse_normal_transform"
if false; then
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
# Fit the mean-variance relationship separately within each E stratum, and predict each
# gene's per-stratum variance from its per-stratum mean.
# Note: --log_transform_mean is deliberately NOT passed. These expression matrices are
# already on a log2 scale, so their per-gene means can be negative and must not be logged again.
source ~/.bashrc
conda activate plink_env
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm.unstandardized"
per_gene_variance_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.${cell_type}.per_gene_variance.txt"
mean_variance_model_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.${cell_type}.per_gene_variance.mean_variance_model.txt"
mean_variance_coefficient_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.${cell_type}.mean_variance_model_coefficients.txt"
if false; then
python fit_mean_variance_model.py \
    --per_gene_variance_file ${per_gene_variance_file} \
    --mean_variance_model_output_file ${mean_variance_model_file} \
    --coefficient_output_file ${mean_variance_coefficient_file} \
    --polynomial_degree 3

# Visualize the fitted mean-variance relationship, separately for each E stratum
mean_variance_plot_stem="${processed_expression_dir}/${tissue_name}.${normalization_method}.${cell_type}"
Rscript visualize_mean_variance_model.R \
    ${mean_variance_model_file} \
    ${tissue_name} \
    ${cell_type} \
    ${mean_variance_plot_stem}
fi

###############################
# Create gene based groupings
# Built from the mean-variance model file (a superset of the per-gene variance file), so that
# the pred_abs_var_diff categories have the fitted variances available
if false; then
gene_categories_output_stem=${gene_categories_dir}"gene_categories_"${tissue_name}"_"${cell_type}
python make_gene_category_files.py $mean_variance_model_file $gene_categories_output_stem
fi



###############################
# Run PA-H2 regression
## Log-TMM based normalization
if false; then
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
# files
gene_category_file="none"
tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_all_genes"
sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file

## INT based normalization
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="inverse_normal_transform"
# files
gene_category_file="none"
tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_all_genes"
sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file

# Stratified by gene annotations
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
gene_category_statistic="mean"
for gene_category_bin in "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"; do
    # files
    gene_category_file=${gene_categories_dir}"gene_categories_"${tissue_name}"_"${cell_type}"_"$gene_category_statistic"_based_categories_category_"$gene_category_bin".txt"
    tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
    genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
    covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
    E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_"${gene_category_statistic}"_"${gene_category_bin}
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file
done

tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
gene_category_statistic="abs_var_diff"
for gene_category_bin in "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"; do
    # files
    gene_category_file=${gene_categories_dir}"gene_categories_"${tissue_name}"_"${cell_type}"_"$gene_category_statistic"_based_categories_category_"$gene_category_bin".txt"
    tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
    genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
    covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
    E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_"${gene_category_statistic}"_"${gene_category_bin}
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file
done

tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
gene_category_statistic="abs_de_t"
for gene_category_bin in "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"; do
    # files
    gene_category_file=${gene_categories_dir}"gene_categories_"${tissue_name}"_"${cell_type}"_"$gene_category_statistic"_based_categories_category_"$gene_category_bin".txt"
    tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
    genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
    covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
    E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_"${gene_category_statistic}"_"${gene_category_bin}
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file
done

tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
gene_category_statistic="pred_abs_var_diff"
for gene_category_bin in "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"; do
    # files
    gene_category_file=${gene_categories_dir}"gene_categories_"${tissue_name}"_"${cell_type}"_"$gene_category_statistic"_based_categories_category_"$gene_category_bin".txt"
    tissue_expression_matrix_file="${processed_expression_dir}/${tissue_name}.${normalization_method}.txt.gz"
    genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
    covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
    E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_"${normalization_method}"_"$cell_type"_"${gene_category_statistic}"_"${gene_category_bin}
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file
done



###############################
# Run PA-H2 regression on log_tmm downsampled data
## Log-TMM based normalization
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
# files that do not depend on the downsampling percentage
gene_category_file="none"
genotype_stem=${processed_genotype_dir}"gtex_v9_eqtl_chr"
E_var_file=${processed_expression_dir}${tissue_name}".xcell_"${cell_type}"_binary.txt"
covariate_file="${processed_expression_dir}/${tissue_name}.covariates.txt"
for downsampling_percentage in "0.2" "0.4" "0.6" "0.8" "1.0"; do
    # All genes (no lowly-expressed-gene filter)
    tissue_expression_matrix_file=${downsampling_processed_expression_dir}${tissue_name}".downsample_"${downsampling_percentage}"."${normalization_method}".txt.gz"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_downsample_"${downsampling_percentage}"_"${normalization_method}"_"$cell_type"_all_genes"
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file

    # Lowly-expressed genes removed
    tissue_expression_matrix_file=${downsampling_processed_expression_dir}${tissue_name}".downsample_"${downsampling_percentage}".filtered_genes."${normalization_method}".txt.gz"
    pa_h2_output_stem=${pa_h2_results_dir}"pa_h2_results_"${tissue_name}"_downsample_"${downsampling_percentage}"_filtered_genes_"${normalization_method}"_"$cell_type"_all_genes"
    sbatch run_pa_h2.sh $tissue_expression_matrix_file $genotype_stem $E_var_file $pa_h2_output_stem $PA_H2_code_dir $gene_category_file $covariate_file
done
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


###############################
# Visualize PA-H2 results across the gene category bins
tissue_name="Whole_Blood"
cell_type="Neutrophils"
normalization_method="log_tmm"
gene_category_statistics="mean,abs_var_diff,abs_de_t,pred_abs_var_diff"
num_gene_category_bins="10"
downsampling_percentages="0.2,0.4,0.6,0.8,1.0"
comparison_normalization_methods="log_tmm,inverse_normal_transform"
Rscript visualize_pa_h2.R \
    ${pa_h2_results_dir} \
    ${tissue_name} \
    ${normalization_method} \
    ${cell_type} \
    ${gene_category_statistics} \
    ${num_gene_category_bins} \
    ${downsampling_percentages} \
    ${comparison_normalization_methods} \
    ${pa_h2_results_dir}