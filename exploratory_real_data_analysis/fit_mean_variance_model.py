from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
import os
import sys
import pdb


def parse_args():
    parser = argparse.ArgumentParser(description='Fit the mean-variance relationship across genes, separately within each E stratum, and predict each gene per-stratum variance from its per-stratum mean.')
    parser.add_argument('--per_gene_variance_file', required=True,
                        help='Per-gene variance file produced by compute_per_gene_variance_stratified_by_e_variable.py')
    parser.add_argument('--mean_variance_model_output_file', required=True,
                        help='Output per-gene file: the input columns with the fitted variances appended')
    parser.add_argument('--coefficient_output_file', required=True,
                        help='Output file holding the fitted polynomial coefficients for each E stratum')
    parser.add_argument('--polynomial_degree', type=int, default=3,
                        help='Degree of the polynomial in the mean (default 3, ie. cubic)')
    parser.add_argument('--log_transform_mean', action='store_true',
                        help='Take log2 of the per-gene mean before fitting. Leave this off for the log_tmm / '
                             'inverse_normal_transform matrices, whose means are already on a log scale and can be '
                             'negative; only use it if the expression matrix is on a raw (linear) scale.')
    return parser.parse_args()


def fit_mean_variance_model(mean_values, variance_values, polynomial_degree, log_transform_mean):
    """
    Regress log2 per-gene variance on a polynomial in the per-gene mean, across genes within
    a single E stratum. Returns the coefficients (highest power first, as np.polyfit gives
    them), the predicted log2 variance for every gene, and the mask of genes that informed
    the fit.
    """
    # Genes with a non-positive mean have no log2, and genes with a non-positive variance have
    # no log2 variance to regress; both are dropped from the fit rather than silently becoming nan
    with np.errstate(invalid='ignore', divide='ignore'):
        xx = np.log2(mean_values) if log_transform_mean else np.copy(mean_values)
        yy = np.log2(variance_values)
    fittable = np.isfinite(xx) & np.isfinite(yy)

    if np.sum(fittable) <= polynomial_degree:
        print('assumption error: only ' + str(int(np.sum(fittable))) + ' genes available to fit a degree ' + str(polynomial_degree) + ' polynomial')
        pdb.set_trace()

    coefficients = np.polyfit(xx[fittable], yy[fittable], polynomial_degree)

    # Predict for every gene with a usable mean, including genes excluded from the fit because
    # their observed variance was non-positive
    with np.errstate(invalid='ignore'):
        predicted_log2_variance = np.polyval(coefficients, xx)
    predicted_log2_variance[np.isfinite(xx) == False] = np.nan

    return coefficients, predicted_log2_variance, fittable


def variance_explained(observed, predicted, fittable):
    # Fraction of variance in observed log2 variance explained by the fitted curve
    resid = observed[fittable] - predicted[fittable]
    total = observed[fittable] - np.mean(observed[fittable])
    return 1.0 - (np.sum(np.square(resid))/np.sum(np.square(total)))


def main():
    args = parse_args()

    df = pd.read_csv(args.per_gene_variance_file, sep='\t')
    for required_column in ['gene_name', 'var_E0', 'var_E1', 'mean_E0', 'mean_E1']:
        if required_column not in df.columns:
            print('assumption error: column ' + required_column + ' missing from ' + args.per_gene_variance_file)
            pdb.set_trace()

    coefficient_file = open(args.coefficient_output_file, 'w')
    coefficient_file.write('E_stratum\tpolynomial_degree\tlog_transform_mean\tn_genes_fit\tvariance_explained\tcoefficients_highest_power_first\n')

    with np.errstate(invalid='ignore', divide='ignore'):
        observed_log2_variance = {'E0': np.log2(df['var_E0'].values), 'E1': np.log2(df['var_E1'].values)}

    for e_stratum in ['E0', 'E1']:
        coefficients, predicted_log2_variance, fittable = fit_mean_variance_model(
            df['mean_' + e_stratum].values,
            df['var_' + e_stratum].values,
            args.polynomial_degree,
            args.log_transform_mean)

        df['pred_log2_var_' + e_stratum] = predicted_log2_variance
        df['pred_var_' + e_stratum] = np.power(2.0, predicted_log2_variance)

        coefficient_file.write(e_stratum + '\t' + str(args.polynomial_degree) + '\t' + str(args.log_transform_mean) + '\t' +
                               str(int(np.sum(fittable))) + '\t' +
                               str(variance_explained(observed_log2_variance[e_stratum], predicted_log2_variance, fittable)) + '\t' +
                               ','.join(str(coefficient) for coefficient in coefficients) + '\n')

    coefficient_file.close()

    df.to_csv(args.mean_variance_model_output_file, sep='\t', index=False, na_rep='nan')

    print('Fit mean-variance model for ' + str(df.shape[0]) + ' genes; wrote ' + args.mean_variance_model_output_file)


if __name__ == '__main__':
    main()
