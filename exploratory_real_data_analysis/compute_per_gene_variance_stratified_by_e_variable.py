from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
import os
import sys
import pdb
import gzip
import scipy.stats as stats


def parse_args():
    parser = argparse.ArgumentParser(description='Compute per-gene expression variance stratified by an E (environmental) variable.')
    parser.add_argument('--expression_matrix_file', required=True,
                        help='Path to the tissue expression matrix (genes x individuals, gzipped)')
    parser.add_argument('--xcell_ct_proportions_file', required=True,
                        help='Path to the tissue xCell cell-type proportions file (cell types x individuals, gzipped)')
    parser.add_argument('--tissue_name', required=True,
                        help='Tissue name (e.g. Whole_Blood)')
    parser.add_argument('--cell_type', required=True,
                        help='Cell type to use as the E variable (e.g. Neutrophils)')
    parser.add_argument('--per_gene_variance_output_file', required=True,
                        help='Output file for per-gene variance stratified by the E variable')
    return parser.parse_args()


def extract_binary_E_variable_from_xcell(xcell_ct_proportions_file, cell_type):
    """
    Extract a binary E variable (cell type) from the xCell cell-type proportions file.
    The E variable is 1 if the cell type proportion is above the median, and 0 otherwise.
    """
    f = gzip.open(xcell_ct_proportions_file, 'rt')
    head_count = 0
    counter = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            indi_ids = np.asarray(data[1:])
            head_count += 1
            continue
        if data[0] != cell_type:
            continue
        if counter > 0:
            print('assumption eroror')
            pdb.set_trace()
        counter += 1
        E_variable = np.asarray(data[1:]).astype(float)
    f.close()

    # Now binarize E-variable: inverse normal transform, then threshold at 0
    n = len(E_variable)
    E_variable_int = stats.norm.ppf((stats.rankdata(E_variable) - 0.5) / n)  # rankdata averages ties
    E_variable_binary = (E_variable_int > 0.0).astype(int)

    return E_variable_binary, indi_ids


def main():
    args = parse_args()

    expression_matrix_file = args.expression_matrix_file
    xcell_ct_proportions_file = args.xcell_ct_proportions_file
    tissue_name = args.tissue_name
    cell_type = args.cell_type
    per_gene_variance_output_file = args.per_gene_variance_output_file

    # First, extract the binary E-variable (cell type) from the xCell cell-type proportions file
    E_variable, E_variable_sample_names = extract_binary_E_variable_from_xcell(xcell_ct_proportions_file, cell_type)

    var_diff = []
    mean_diff = []
    rat = []
    t = open(per_gene_variance_output_file, 'w')
    t.write('gene_name\tvar_E0\tvar_E1\tmean_E0\tmean_E1\n')
    f = gzip.open(expression_matrix_file, 'rt')
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            sample_names = np.asarray(data[4:])
            if np.array_equal(sample_names, E_variable_sample_names) == False:
                print('Sample names in expression matrix and xCell proportions file do not match!')
                pdb.set_trace()
            head_count += 1
            continue
        gene_name = data[3]  # columns are #chr, start, end, gene_id, <samples...>
        expression_values = np.asarray(data[4:]).astype(float)

        # Now compute per-gene variance stratified by the E-variable
        var_E0 = np.var(expression_values[E_variable == 0], ddof=1)
        var_E1 = np.var(expression_values[E_variable == 1], ddof=1)

        # Now compute per-gene variance stratified by the E-variable
        mean_E0 = np.mean(expression_values[E_variable == 0])
        mean_E1 = np.mean(expression_values[E_variable == 1])


        t.write(f'{gene_name}\t{var_E0}\t{var_E1}\t{mean_E0}\t{mean_E1}\n')
        rat.append(var_E1 / var_E0)
        var_diff.append(var_E1 - var_E0)
        mean_diff.append(mean_E1 - mean_E0)

    t.close()
    f.close()

    var_diff = np.asarray(var_diff)
    mean_diff = np.asarray(mean_diff)

    print(np.sort(rat))




if __name__ == '__main__':
    main()
