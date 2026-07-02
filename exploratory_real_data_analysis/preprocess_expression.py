from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
import os
import sys
import pdb
import gzip
import re
import rnaseqnorm

def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess GTEx expression data for a single tissue.')
    parser.add_argument('--gtex_v8_gene_reads_file', required=True,
                        help='Path to GTEx v8 gene read counts file (.gct.gz)')
    parser.add_argument('--gtex_v8_eqtl_expression_matrices_dir', required=True,
                        help='Directory containing GTEx v8 eQTL expression matrices')
    parser.add_argument('--xcell_ct_proportions_file', required=True,
                        help='Path to xCell cell-type proportions file')
    parser.add_argument('--gtex_sample_attributes_file', required=True,
                        help='Path to GTEx sample attributes file')
    parser.add_argument('--gtex_covariate_dir', required=True,
                        help='Directory containing GTEx v8 eQTL covariate files')
    parser.add_argument('--processed_expression_dir', required=True,
                        help='Output directory for processed expression data')
    parser.add_argument('--tissue_name', required=True,
                        help='Tissue name (e.g. Whole_Blood)')
    return parser.parse_args()

def load_in_eqtl_genes_and_tissues(filer):
    valid_chroms = {}
    for chrom_num in range(1, 23):
        valid_chroms['chr' + str(chrom_num)] = 1

    f = gzip.open(filer, 'rt')
    ordered_genes = []
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            ordered_indi_ids = np.asarray(data[4:])
            head_count += 1
            continue
        if data[0] not in valid_chroms:
            continue
        gene_id = data[3]
        ordered_genes.append(gene_id)
    f.close()

    gene_dictionary = {}
    for i in range(len(ordered_genes)):
        gene_dictionary[ordered_genes[i]] = i

    indi_dictionary = {}
    for i in range(len(ordered_indi_ids)):
        indi_dictionary[ordered_indi_ids[i]] = i

    return ordered_indi_ids, np.asarray(ordered_genes), gene_dictionary, indi_dictionary

def to_underscore(s: str) -> str:
    s = s.replace(")", "")
    # replace " (" OR any run of spaces and hyphens with "_"
    return re.sub(r"(?:\s*\(|[ -]+)", "_", s)


def create_mapping_from_gtex_sample_id_to_individual_tissue_format(gtex_sample_attributes_file, tissue_type):
    gtex_sample_to_tissue_name = {}
    gtex_tissues = {}
    f = open(gtex_sample_attributes_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        gtex_sample_id = data[0]
        tissue_type = to_underscore(data[6])
        if tissue_type == 'Cells_EBV_transformed_lymphocytes':
            tissue_type = 'Cells_EBV-transformed_lymphocytes'
        elif tissue_type == 'Brain_Spinal_cord_cervical_c_1':
            tissue_type = 'Brain_Spinal_cord_cervical_c-1'

        gtex_sample_to_tissue_name[gtex_sample_id] = tissue_type
        gtex_tissues[tissue_type] = 1
    f.close()
    return gtex_sample_to_tissue_name

def generate_tissue_gene_reads_file(tissue_gene_reads_file, gtex_v8_gene_reads_file, ordered_eqtl_indi_ids, ordered_eqtl_gene_names, eqtl_gene_dictionary, gtex_sample_id_to_individual_tissue_format, eqtl_indi_dictionary, tissue_name):
    f = gzip.open(gtex_v8_gene_reads_file, 'rt')
    t = gzip.open(tissue_gene_reads_file, 'wt')
    head_count = 0
    used_genes = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if line.startswith('#'):
            continue
        if head_count == 0:
            head_count += 1
            continue
        if head_count == 1:
            head_count += 1
            file_sample_ids = np.asarray(data[2:])
            valid_sample_indices = []
            ordered_indi_ids = []
            for i in range(len(file_sample_ids)):
                sample_id = file_sample_ids[i]
                if sample_id not in gtex_sample_id_to_individual_tissue_format:
                    continue
                if gtex_sample_id_to_individual_tissue_format[sample_id] != tissue_name:
                    continue
                individual_id = '-'.join(sample_id.split('-')[:2])
                if individual_id not in eqtl_indi_dictionary:
                    continue
                valid_sample_indices.append(i)
                ordered_indi_ids.append(individual_id)
            ordered_indi_ids = np.asarray(ordered_indi_ids)
            valid_sample_indices = np.asarray(valid_sample_indices)
            if len(ordered_indi_ids) != len(ordered_eqtl_indi_ids):
                print('assumption error: number of individuals in gene reads file does not match number of individuals in eQTL expression matrix')
                pdb.set_trace()
            t.write(data[0] + '\t' + data[1] + '\t' + '\t'.join(ordered_indi_ids) + '\n')
            continue
        ensamble_id = data[0]
        gene_id = data[1]
        if ensamble_id not in eqtl_gene_dictionary:
            continue
        valid_data = np.asarray(data[2:])[valid_sample_indices]
        t.write(ensamble_id + '\t' + gene_id + '\t' + '\t'.join(valid_data) + '\n')
        used_genes += 1
    f.close()
    t.close()
    if used_genes != len(ordered_eqtl_gene_names):
        print('assumption error: number of genes in gene reads file does not match number of genes in eQTL expression matrix')
        pdb.set_trace()
    return

def edger_logcpm(counts_df, prior_count=0.25):
    # edgeR::cpm(..., log=TRUE) with TMM-normalized effective library sizes
    tmm = rnaseqnorm.edgeR_calcNormFactors(counts_df)     # per-sample TMM factors
    lib_size = counts_df.sum(axis=0) * tmm                # effective library sizes, L_j
    pc = prior_count * lib_size / lib_size.mean()         # library-size-scaled prior, pc_j
    lib_size_adj = (lib_size + 2 * pc) * 1e-6             # inflated, per-million, L'_j
    return np.log2((counts_df + pc) / lib_size_adj)       # broadcasts across samples

def standardize_rows(df):
    # Standardize each gene (row) to mean 0 and variance 1 (population std, ddof=0)
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1, ddof=0), axis=0)

def normalize_expression(tissue_gene_reads_file, tissue_log_tmm_file, tissue_int_file, tissue_log_cpm_file):
    # Load filtered read counts into a genes (rows) x samples (cols) DataFrame
    reads_df = pd.read_csv(tissue_gene_reads_file, sep='\t', index_col=0)
    # First column is the gene Description; drop it so only sample columns remain
    reads_df = reads_df.drop(columns=reads_df.columns[0])

    # TMM-normalized CPM (edgeR_cpm computes TMM library-size factors internally)
    tmm_df = rnaseqnorm.edgeR_cpm(reads_df, normalized_lib_sizes=True)

    # Each output is standardized so every gene (row) has mean 0 and variance 1

    # Output 1: edgeR log2 CPM (TMM-normalized effective library sizes)
    log_tmm_df = edger_logcpm(reads_df)
    standardize_rows(log_tmm_df).to_csv(tissue_log_tmm_file, sep='\t', compression='gzip')

    # Output 2: inverse normal transform of TMM-normalized CPM (per gene, across samples)
    int_df = rnaseqnorm.inverse_normal_transform(tmm_df)
    standardize_rows(int_df).to_csv(tissue_int_file, sep='\t', compression='gzip')

    # Output 3: naive log2 of TMM-normalized CPM (log2(cpm + 1))
    log_cpm_df = np.log2(tmm_df + 1)
    standardize_rows(log_cpm_df).to_csv(tissue_log_cpm_file, sep='\t', compression='gzip')
    return

def main():
    args = parse_args()

    gtex_v8_gene_reads_file = args.gtex_v8_gene_reads_file
    gtex_v8_eqtl_expression_matrices_dir = args.gtex_v8_eqtl_expression_matrices_dir
    xcell_ct_proportions_file = args.xcell_ct_proportions_file
    gtex_sample_attributes_file = args.gtex_sample_attributes_file
    gtex_covariate_dir = args.gtex_covariate_dir
    processed_expression_dir = args.processed_expression_dir
    tissue_name = args.tissue_name
    
    # Extract gtex individual IDs for the tissue of interest
    ordered_eqtl_indi_ids, ordered_eqtl_gene_names, eqtl_gene_dictionary, eqtl_indi_dictionary = load_in_eqtl_genes_and_tissues(gtex_v8_eqtl_expression_matrices_dir + tissue_name + '.v8.normalized_expression.bed.gz')

    # Create mapping from sample ID to tissue name
    gtex_sample_id_to_individual_tissue_format = create_mapping_from_gtex_sample_id_to_individual_tissue_format(gtex_sample_attributes_file, tissue_name)
    
    
    # Filter gene reads file to only samples in this tissue, and filter to genes in the eQTL expression matrix
    tissue_gene_reads_file = processed_expression_dir + tissue_name + '.gene_reads.filtered.txt.gz'
    generate_tissue_gene_reads_file(tissue_gene_reads_file, gtex_v8_gene_reads_file, ordered_eqtl_indi_ids, ordered_eqtl_gene_names, eqtl_gene_dictionary, gtex_sample_id_to_individual_tissue_format, eqtl_indi_dictionary, tissue_name)


    # Normalize expression: edgeR logCPM, inverse normal transformed TMM-CPM, and naive log2(cpm+1)
    tissue_log_tmm_file = processed_expression_dir + tissue_name + '.log_tmm.txt.gz'
    tissue_int_file = processed_expression_dir + tissue_name + '.inverse_normal_transform.txt.gz'
    tissue_log_cpm_file = processed_expression_dir + tissue_name + '.log_cpm.txt.gz'
    normalize_expression(tissue_gene_reads_file, tissue_log_tmm_file, tissue_int_file, tissue_log_cpm_file)


if __name__ == '__main__':
    main()
