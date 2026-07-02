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

processed_expression_dir="${output_root}/processed_expression"



########################
# Code
########################

tissue_name="Whole_Blood"
sh preprocess_expression.sh \
    ${gtex_v8_gene_reads_file} \
    ${gtex_v8_eqtl_expression_matrices_dir} \
    ${xcell_ct_proportions_file} \
    ${gtex_sample_attributes_file} \
    ${gtex_covariate_dir} \
    ${processed_expression_dir} \
    ${tissue_name}



