import numpy as np
import os
import sys
import pdb
import argparse
import gzip
import pgenlib as pg
from scipy import stats
from mixture_model_input_utils import MIXTURE_MODEL_INPUT_HEADER, compute_test_summary_stats, format_input_line


# Generates PA mixture model input files from the same inputs as PA_h2.py
# (expression bed, plink2 per-chrom genotype, binary E file, covariate file)
# plus a file of variant-gene pairs to test.
# If no variant-gene pairs file is provided, pairs are identified with a cis interaction-QTL scan
# (top variant per gene, for genes with significant interaction QTLs)


def open_expression_bed(filer):
    if filer.endswith('.gz'):
        return gzip.open(filer, 'rt')
    return open(filer)


def load_in_E_variable(E_var_file):
    raw_cov = np.loadtxt(E_var_file, dtype=str, delimiter='\t')
    sample_names = raw_cov[0, 1:]
    interaction_var = raw_cov[1, 1:].astype(float)

    if np.array_equal(np.sort(np.unique(interaction_var)), np.asarray([0.0, 1.0])) == False:
        print('Fatal error: Interaction variable must be binary and encoded as 0 / 1')
        sys.exit(1)

    return interaction_var, sample_names


def load_in_covariates(covariate_file):
    # Returns (covariate matrix (num_samples X num_covariates), sample_names)
    # Returns (None, None) if no covariate file was provided
    if covariate_file == 'none':
        return None, None

    raw_cov = np.loadtxt(covariate_file, dtype=str, delimiter='\t', ndmin=2)
    sample_names = raw_cov[0, 1:]
    CC = np.transpose(raw_cov[1:, 1:].astype(float))

    return CC, sample_names


def load_in_variant_gene_pairs(variant_gene_pairs_file):
    # First two columns are gene id and variant id (header line is skipped)
    # A tensorQTL cis_qtl_top_assoc file (phenotype_id, variant_id, ...) works directly
    pairs = []
    f = open_expression_bed(variant_gene_pairs_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        if head_count == 0:
            head_count = head_count + 1
            continue
        if line == '':
            continue
        data = line.split('\t')
        pairs.append((data[0], data[1]))
    f.close()
    return pairs


def load_in_expression_for_genes(expression_bed, gene_set):
    # Returns (sample_names, dictionary mapping gene_id to (chrom_string, gene_tss, expression vector))
    # If gene_set is None, all genes are loaded
    gene_to_expr = {}
    f = open_expression_bed(expression_bed)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            sample_names = np.asarray(data[4:])
            continue
        gene_id = data[3]
        if gene_set is not None and gene_id not in gene_set:
            continue
        gene_to_expr[gene_id] = (data[0], float(data[1]), np.asarray(data[4:]).astype(float))
    f.close()
    return sample_names, gene_to_expr


def load_in_pvar(pvar_path):
    # Returns (variant positions, variant ids), in pgen row order
    positions = []
    variant_ids = []
    f = open(pvar_path)
    for line in f:
        if line.startswith('##'):
            continue
        data = line.rstrip().split('\t')
        if line.startswith('#'):
            pos_column = data.index('POS')
            id_column = data.index('ID')
            continue
        positions.append(float(data[pos_column]))
        variant_ids.append(data[id_column])
    f.close()
    return np.asarray(positions), np.asarray(variant_ids)


def load_in_psam_sample_names(psam_path):
    f = open(psam_path)
    header = f.readline().rstrip().split('\t')
    iid_column = header.index('IID')
    sample_names = []
    for line in f:
        sample_names.append(line.rstrip().split('\t')[iid_column])
    f.close()
    return np.asarray(sample_names)


def reorder_genotype_samples(geno_samples, E_var_sample_names):
    mapping = {}
    for ii, geno_sample in enumerate(geno_samples):
        mapping[geno_sample] = ii

    reordering_vec = []
    for sample_name in E_var_sample_names:
        if sample_name not in mapping:
            print('Fatal error: Sample in interaction (E) variable that is not found in plink2 IID')
            sys.exit(1)
        reordering_vec.append(mapping[sample_name])
    reordering_vec = np.asarray(reordering_vec)

    return reordering_vec


def residualize_and_standardize(xx, design_mat):
    # Regress out covariates (design_mat includes intercept) and standardize globally (across both E bins)
    if design_mat is not None:
        coef = np.linalg.lstsq(design_mat, xx, rcond=None)[0]
        xx = xx - np.dot(design_mat, coef)
    return (xx - np.mean(xx))/np.std(xx)


def mean_impute_missing_genotypes(G):
    # G: (num_samples X num_variants). Missing genotypes are encoded as -9
    missing = G < 0
    if np.sum(missing) > 0:
        G = np.where(missing, np.nan, G)
        G = np.where(missing, np.nanmean(G, axis=0, keepdims=True), G)
    return G


def benjamini_hochberg(pvalues):
    n_tests = len(pvalues)
    order = np.argsort(pvalues)
    ranked_qvalues = pvalues[order] * n_tests / np.arange(1, n_tests + 1)
    ranked_qvalues = np.minimum(np.minimum.accumulate(ranked_qvalues[::-1])[::-1], 1.0)
    qvalues = np.empty(n_tests)
    qvalues[order] = ranked_qvalues
    return qvalues


def interaction_qtl_test_for_single_gene(yy, GG, EE, Q):
    # Tests y ~ g + E + g*E + covariates for each column (variant) of GG
    # Q: orthonormal basis of covariate design matrix (includes intercept)
    # Returns (interaction betas, interaction beta standard errors, interaction pvalues), one per variant
    dof = len(yy) - Q.shape[1] - 3

    # Regress out covariates (Frisch-Waugh-Lovell)
    yy_r = yy - np.dot(Q, np.dot(np.transpose(Q), yy))
    EE_r = EE - np.dot(Q, np.dot(np.transpose(Q), EE))
    GG_r = GG - np.dot(Q, np.dot(np.transpose(Q), GG))
    GE = GG * EE[:, np.newaxis]
    GE_r = GE - np.dot(Q, np.dot(np.transpose(Q), GE))

    # Per-variant design matrix (num_samples X num_variants X 3): genotype, E, genotype*E
    XX = np.stack((GG_r, np.broadcast_to(EE_r[:, np.newaxis], GG_r.shape), GE_r), axis=2)
    XtX_inv = np.linalg.inv(np.einsum('nmi,nmj->mij', XX, XX))
    coef = np.einsum('mij,mj->mi', XtX_inv, np.einsum('nmi,n->mi', XX, yy_r))
    rss = np.sum(np.square(yy_r[:, np.newaxis] - np.einsum('nmi,mi->nm', XX, coef)), axis=0)

    beta_ge = coef[:, 2]
    beta_ge_se = np.sqrt(rss / dof * XtX_inv[:, 2, 2])
    pvalues = 2.0 * stats.t.sf(np.abs(beta_ge / beta_ge_se), dof)
    return beta_ge, beta_ge_se, pvalues


def run_cis_interaction_qtl_scan(args, EE, E_var_sample_names, design_mat):
    # For every autosomal gene in expression bed, test all cis variants for interaction with E,
    # and keep the top variant per gene. Gene-level pvalue is Bonferroni corrected across the gene's cis variants,
    # and genes are called significant by Benjamini-Hochberg across genes.
    # Returns list of (gene_id, variant_id) for significant genes
    expression_sample_names, gene_to_expr = load_in_expression_for_genes(args.expression_bed, None)
    if np.array_equal(expression_sample_names, E_var_sample_names) == False:
        print('Fatal error: Expression sample names do not match interaction variable sample names (need 1 to 1 match)')
        sys.exit(1)

    if design_mat is None:
        design_mat = np.ones((len(EE), 1))
    Q = np.linalg.qr(design_mat)[0]

    # Keep track of top variant per gene: (gene_id, variant_id, n_cis_variants_tested, beta, beta_se, pvalue)
    top_per_gene = []

    for chrom_num in range(1, 23):
        chrom_string = 'chr' + str(chrom_num)
        chrom_gene_ids = [gene_id for gene_id in gene_to_expr if gene_to_expr[gene_id][0] == chrom_string]
        if len(chrom_gene_ids) == 0:
            continue

        # Load in genotype data for this chromosome
        plink2_genotype_stem = args.plink2_per_chrom_stem + str(chrom_num)
        if os.path.exists(plink2_genotype_stem + '.pgen') == False:
            print('Warning: skipping ' + str(len(chrom_gene_ids)) + ' genes on ' + chrom_string + ' because ' + plink2_genotype_stem + '.pgen does not exist')
            continue
        pgen = pg.PgenReader((plink2_genotype_stem + '.pgen').encode())
        variant_positions, variant_ids = load_in_pvar(plink2_genotype_stem + '.pvar')
        if np.all(np.diff(variant_positions) >= 0) == False:
            print('Fatal error: variants in ' + plink2_genotype_stem + '.pvar are not sorted by position')
            sys.exit(1)
        genotype_sample_reordering = reorder_genotype_samples(load_in_psam_sample_names(plink2_genotype_stem + '.psam'), E_var_sample_names)

        for gene_id in chrom_gene_ids:
            gene_tss = gene_to_expr[gene_id][1]
            yy = gene_to_expr[gene_id][2]
            if np.std(yy) == 0.0:
                continue

            # Extract cis-genotype matrix for this gene (same cis window definition as PA_h2.py)
            cis_start = np.searchsorted(variant_positions, gene_tss - args.cis_radius, side='left')
            cis_end = np.searchsorted(variant_positions, gene_tss + args.cis_radius, side='left')
            if cis_end == cis_start:
                continue
            cis_genotype_mat = np.empty((cis_end - cis_start, pgen.get_raw_sample_ct()), dtype=np.float64)
            pgen.read_dosages_range(int(cis_start), int(cis_end), cis_genotype_mat)
            GG = mean_impute_missing_genotypes(np.transpose(cis_genotype_mat[:, genotype_sample_reordering]))
            cis_variant_ids = variant_ids[cis_start:cis_end]

            # Filter to variants with MAF >= min_maf_per_bin in each E bin
            af_0 = np.mean(GG[EE == 0.0], axis=0) / 2.0
            af_1 = np.mean(GG[EE == 1.0], axis=0) / 2.0
            maf_filter = (np.minimum(af_0, 1.0 - af_0) >= args.min_maf_per_bin) & (np.minimum(af_1, 1.0 - af_1) >= args.min_maf_per_bin)
            if np.sum(maf_filter) == 0:
                continue
            GG = GG[:, maf_filter]
            cis_variant_ids = cis_variant_ids[maf_filter]

            beta_ge, beta_ge_se, pvalues = interaction_qtl_test_for_single_gene(yy, GG, EE, Q)
            top_index = np.argmin(pvalues)
            top_per_gene.append((gene_id, cis_variant_ids[top_index], GG.shape[1], beta_ge[top_index], beta_ge_se[top_index], pvalues[top_index]))

    # Gene-level multiple testing correction
    top_pvalues = np.asarray([gene_tuple[5] for gene_tuple in top_per_gene])
    n_cis_variants_tested = np.asarray([gene_tuple[2] for gene_tuple in top_per_gene])
    gene_bonferroni_pvalues = np.minimum(top_pvalues * n_cis_variants_tested, 1.0)
    gene_qvalues = benjamini_hochberg(gene_bonferroni_pvalues)

    # Print top variant per gene to output (first two columns are a valid variant-gene pairs file)
    pairs = []
    t = open(args.output_stem + '_interaction_qtl_top_per_gene.txt', 'w')
    t.write('gene_id\tvariant_id\tn_cis_variants_tested\tinteraction_beta\tinteraction_beta_se\tinteraction_pvalue\tgene_bonferroni_pvalue\tgene_bh_qvalue\tsignificant\n')
    for gene_index, gene_tuple in enumerate(top_per_gene):
        significant = gene_qvalues[gene_index] < args.interaction_fdr
        if significant:
            pairs.append((gene_tuple[0], gene_tuple[1]))
        t.write(gene_tuple[0] + '\t' + gene_tuple[1] + '\t' + str(gene_tuple[2]) + '\t' + str(gene_tuple[3]) + '\t' + str(gene_tuple[4]) + '\t' + str(gene_tuple[5]) + '\t' + str(gene_bonferroni_pvalues[gene_index]) + '\t' + str(gene_qvalues[gene_index]) + '\t' + str(significant) + '\n')
    t.close()

    print(str(len(pairs)) + ' of ' + str(len(top_per_gene)) + ' tested genes have a significant interaction QTL (BH FDR < ' + str(args.interaction_fdr) + ')')

    return pairs


################################################################################
# main
################################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--expression-bed', default='', type=str,
                        help='Bed file containing gene expression (same format as PA_h2.py)')
    parser.add_argument('--plink2-per-chrom-stem', default='', type=str,
                        help='Path to genotype stem (same format as PA_h2.py)')
    parser.add_argument('--binary-E-interaction-covariate-file', default='', type=str,
                        help='Path to binary E covariate file (same format as PA_h2.py)')
    parser.add_argument('--covariate-file', default='none', type=str,
                        help='File containing additional covariates to regress out of expression and genotype (same format as PA_h2.py). Use none if dont wish to control for additional covariates')
    parser.add_argument('--variant-gene-pairs-file', default='none', type=str,
                        help='File with header whose first two columns are gene id and variant id (one test per line). Use none to identify pairs with a cis interaction-QTL scan (top variant per gene with a significant interaction QTL)')
    parser.add_argument('--output-stem', default='', type=str,
                        help='Path to output file stem')
    # Interaction-QTL scan (only used if --variant-gene-pairs-file is none)
    parser.add_argument('--cis-radius', default=500000, type=int,
                        help='cis window around TSS to consider snps')
    parser.add_argument('--min-maf-per-bin', default=0.05, type=float,
                        help='Minimum minor allele frequency within each E bin for a variant to be tested')
    parser.add_argument('--interaction-fdr', default=0.05, type=float,
                        help='Benjamini-Hochberg FDR threshold (across genes) on gene-level Bonferroni-corrected interaction pvalues')
    args = parser.parse_args()

    ######################
    # Load in data
    ######################
    EE, E_var_sample_names = load_in_E_variable(args.binary_E_interaction_covariate_file)

    CC, covariate_sample_names = load_in_covariates(args.covariate_file)
    design_mat = None
    if CC is not None:
        if np.array_equal(covariate_sample_names, E_var_sample_names) == False:
            print('Fatal error: Covariate sample names do not match interaction variable sample names (need 1 to 1 match)')
            sys.exit(1)
        design_mat = np.hstack((np.ones((len(EE), 1)), CC))

    if args.variant_gene_pairs_file.lower() == 'none':
        pairs = run_cis_interaction_qtl_scan(args, EE, E_var_sample_names, design_mat)
    else:
        pairs = load_in_variant_gene_pairs(args.variant_gene_pairs_file)

    expression_sample_names, gene_to_expr = load_in_expression_for_genes(args.expression_bed, set([pair[0] for pair in pairs]))
    if np.array_equal(expression_sample_names, E_var_sample_names) == False:
        print('Fatal error: Expression sample names do not match interaction variable sample names (need 1 to 1 match)')
        sys.exit(1)

    ######################
    # Compute summary stats for each variant-gene pair
    # Keep track of results in a dictionary mapping pair index to (raw_S_row, geno_resid_S_row)
    ######################
    pair_results = {}

    # Group pairs by chromosome of the gene
    chrom_to_pair_indices = {}
    for pair_index, (gene_id, variant_id) in enumerate(pairs):
        if gene_id not in gene_to_expr:
            print('Warning: skipping ' + gene_id + ':' + variant_id + ' because gene is not in expression bed')
            continue
        chrom_string = gene_to_expr[gene_id][0]
        if chrom_string not in chrom_to_pair_indices:
            chrom_to_pair_indices[chrom_string] = []
        chrom_to_pair_indices[chrom_string].append(pair_index)

    for chrom_string in sorted(chrom_to_pair_indices.keys()):
        chrom_pair_indices = chrom_to_pair_indices[chrom_string]

        # Load in genotype data for this chromosome
        plink2_genotype_stem = args.plink2_per_chrom_stem + chrom_string.replace('chr', '')
        if os.path.exists(plink2_genotype_stem + '.pgen') == False:
            print('Warning: skipping ' + str(len(chrom_pair_indices)) + ' pairs on ' + chrom_string + ' because ' + plink2_genotype_stem + '.pgen does not exist')
            continue
        pgen = pg.PgenReader((plink2_genotype_stem + '.pgen').encode())
        variant_to_index = {}
        for variant_index, variant_id in enumerate(load_in_pvar(plink2_genotype_stem + '.pvar')[1]):
            variant_to_index[variant_id] = variant_index
        genotype_sample_reordering = reorder_genotype_samples(load_in_psam_sample_names(plink2_genotype_stem + '.psam'), E_var_sample_names)

        buf = np.empty(pgen.get_raw_sample_ct(), dtype=np.float64)

        for pair_index in chrom_pair_indices:
            gene_id, variant_id = pairs[pair_index]
            if variant_id not in variant_to_index:
                print('Warning: skipping ' + gene_id + ':' + variant_id + ' because variant is not in ' + plink2_genotype_stem + '.pvar')
                continue

            # Extract genotype for this variant (filtered/reordered to E variable samples)
            pgen.read_dosages(variant_to_index[variant_id], buf)
            genotype = buf[genotype_sample_reordering].copy()
            # Mean impute missing genotypes (encoded as -9)
            missing = genotype < 0
            if np.sum(missing) > 0:
                genotype[missing] = np.mean(genotype[~missing])

            expr = gene_to_expr[gene_id][2]

            if np.std(genotype) == 0.0 or np.std(expr) == 0.0:
                print('Warning: skipping ' + gene_id + ':' + variant_id + ' because genotype or expression has zero variance')
                continue

            # Regress out covariates and standardize expression and genotype globally
            std_genotype = residualize_and_standardize(genotype, design_mat)
            std_expr = residualize_and_standardize(expr, design_mat)

            pair_results[pair_index] = compute_test_summary_stats(std_genotype, std_expr, EE)

    ######################
    # Print to output (tests ordered as in variant-gene pairs file)
    ######################
    t_mixture_model_input = open(args.output_stem + '_pa_mixture_model_input.txt', 'w')
    t_mixture_model_geno_resid_S_input = open(args.output_stem + '_pa_mixture_model_geno_resid_S_input.txt', 'w')
    t_test_key = open(args.output_stem + '_pa_mixture_model_test_key.txt', 'w')
    t_mixture_model_input.write(MIXTURE_MODEL_INPUT_HEADER)
    t_mixture_model_geno_resid_S_input.write(MIXTURE_MODEL_INPUT_HEADER)
    t_test_key.write('test_idx\tgene_id\tvariant_id\n')

    test_idx = 0
    for pair_index, (gene_id, variant_id) in enumerate(pairs):
        if pair_index not in pair_results:
            continue
        raw_S_row, geno_resid_S_row = pair_results[pair_index]
        t_mixture_model_input.write(format_input_line(test_idx, raw_S_row))
        t_mixture_model_geno_resid_S_input.write(format_input_line(test_idx, geno_resid_S_row))
        t_test_key.write(str(test_idx) + '\t' + gene_id + '\t' + variant_id + '\n')
        test_idx = test_idx + 1

    t_mixture_model_input.close()
    t_mixture_model_geno_resid_S_input.close()
    t_test_key.close()

    print(str(test_idx) + ' of ' + str(len(pairs)) + ' variant-gene pairs written to ' + args.output_stem + '_pa_mixture_model_input.txt')

    return



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
