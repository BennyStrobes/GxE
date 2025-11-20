import numpy as np
import os
import sys
import pdb
import argparse
import pgenlib as pg
import pandas as pd
from tqdm import tqdm


def get_gene_tss(protein_coding_gene_file, gene_tss_file):
	# First extract list of protein coding genes
	pc_genes = {}
	f = open(protein_coding_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		gene_id = data[0]
		pc_genes[data[0]] = 1
	f.close()

	# Second get array of tss's of all PC gene on desired chromosome
	tss_arr = []
	f = open(gene_tss_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		gene_id = data[1]
		if len(data) != 15:
			pdb.set_trace()
		if gene_id not in pc_genes:
			continue

		if len(data[7].split('hr')) != 2:
			continue

		line_chrom = data[7].split('hr')[1]
		if line_chrom == 'Y' or line_chrom == 'X':
			continue

		tss = float(data[8])
		tss_arr.append((gene_id, int(line_chrom), tss))
	f.close()
	sorted_list = sorted(tss_arr, key=lambda x: (x[1], x[2]))

	return sorted_list


def load_in_genotype(plink2_genotype_stem, chrom_num):
	pgen_path = plink2_genotype_stem + str(chrom_num) + '.pgen'
	pvar_path = plink2_genotype_stem + str(chrom_num) + '.pvar'
	psam_path = plink2_genotype_stem + str(chrom_num) + '.psam'

	# Load in thing
	pgen = pg.PgenReader(pgen_path.encode())

	n_samples  = pgen.get_raw_sample_ct()
	n_variants = pgen.get_variant_ct()

	pvar = np.loadtxt(pvar_path, dtype=str, delimiter='\t')
	psam = np.loadtxt(psam_path, dtype=str, delimiter='\t')

	if pvar.shape[0] != n_variants:
		print('assumption erroror')
		pdb.set_trace()
	if psam.shape[0] != n_samples:
		print('assumpton error')
		pdb.set_trace()

	return pgen, pvar, psam

def extract_dosages_from_reader(ra, cis_snp_indices, dtype=np.float32):
	"""
	Extract dosage genotypes from an existing pgenlib.PgenReader.

	Parameters
	----------
	ra : pgenlib.PgenReader
		Already opened reader.
	cis_snp_indices : np.ndarray[bool]
		Boolean mask of length K (#variants).
	dtype : numpy dtype
		Output dtype (default float32 for dosage values).

	Returns
	-------
	D : np.ndarray of shape (M, N)
		Dosage matrix for selected variants.
	idx : np.ndarray[int] of length M
		The variant indices selected (to align with .pvar).
	"""
	n_samples  = ra.get_raw_sample_ct()
	n_variants = ra.get_variant_ct()

	# Convert mask → indices
	idx = np.flatnonzero(cis_snp_indices)
	M = idx.size

	# Allocate output matrix
	D = np.empty((M, n_samples), dtype=dtype)

	# Temporary buffer for a single variant
	buf = np.empty(n_samples, dtype=dtype)

	# Read dosages one variant at a time
	for j, v_idx in enumerate(idx):
		ra.read_dosages(int(v_idx), buf)
		D[j] = buf

	return D

def standardize_genotype(G, axis=1, eps=1e-12):
	"""
	Standardize a genotype matrix SNP-wise (rows = SNPs, columns = samples).

	Parameters
	----------
	G : np.ndarray (num_snps, num_samples)
		Genotype or dosage matrix.
	axis : int
		Axis along which to compute mean/variance (default: 1 = per SNP).
	eps : float
		Small constant to avoid division by zero.

	Returns
	-------
	G_std : np.ndarray
		Standardized genotype matrix with zero-variance SNPs removed.
	keep_mask : np.ndarray[bool]
		Boolean mask indicating which SNPs were retained.
	"""

	# Compute mean and std per SNP
	mean = G.mean(axis=axis, keepdims=True)
	std  = G.std(axis=axis, keepdims=True)

	# Identify SNPs with nonzero variance
	keep_mask = (std.squeeze() > eps)

	# Subset to variance>0 SNPs
	G_keep = G[keep_mask]

	# Standardize
	mean_keep = mean[keep_mask]
	std_keep  = std[keep_mask]

	G_std = (G_keep - mean_keep) / std_keep

	return G_std, keep_mask


def simulate_total_proportional_amplification(gene_geno, EE, proportional_amplification_scaling_factor, n_proportional_amplication_causal_snps_per_gene, orig_ge_h2):
	# Remove snps that have no variance across indivudals
	polymorphic_snp_indices = np.var(gene_geno,axis=0) != 0.0
	gene_geno = gene_geno[:, polymorphic_snp_indices]

	# Get sample size and number of snps
	NN,KK = gene_geno.shape

	if np.sum(polymorphic_snp_indices) != len(polymorphic_snp_indices):
		print('assumption erororor')
		pdb.set_trace()

	# Generate original causal effect sizes
	orig_causal_effect_sizes = np.zeros(KK)
	causal_indices = np.random.choice(np.arange(KK), size=n_proportional_amplication_causal_snps_per_gene, replace=False)
	orig_causal_effect_sizes[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(orig_ge_h2/n_proportional_amplication_causal_snps_per_gene), size=n_proportional_amplication_causal_snps_per_gene)


	# Get variance scaling factor
	SS = np.ones(len(EE))
	SS[EE==1.0] = proportional_amplification_scaling_factor
	SS = SS/np.mean(SS)

	PA_geno = gene_geno * np.sqrt(SS[:, np.newaxis])
	mean_effects = np.dot(PA_geno, orig_causal_effect_sizes)

	only_genetic_effects = np.dot(gene_geno, orig_causal_effect_sizes)


	sqrtS = np.sqrt(SS)

	G_const = np.mean(sqrtS) * only_genetic_effects
	G_int = (sqrtS - np.mean(sqrtS)) * only_genetic_effects
	# Note that G_const + G_int == mean_effects


	scaled_resid_var = np.copy(SS)*(1.0 - orig_ge_h2)
	residiuals = np.random.normal(loc=0, scale=np.sqrt(scaled_resid_var))


	YY = mean_effects + residiuals

	return YY, SS, np.var(G_const,ddof=1), np.var(G_int, ddof=1)





parser = argparse.ArgumentParser()
parser.add_argument('--sample-size', default=450, type=int,
					help='QTL sample size to run the simulation on')
parser.add_argument('--output-stem', default='', type=str,
					help='Path to output file stem')
parser.add_argument('--gene-tss-file', default='', type=str,
					help='file containing gene tss')
parser.add_argument('--protein-coding-gene-file', default='', type=str,
					help='file containing list of protein coding genes')
parser.add_argument('--simulation-type', default='total_PA', type=str,
					help='general type of simulation')
parser.add_argument('--plink2-genotype-stem', default='', type=str,
					help='Path to genotype stem')


# Defaults
parser.add_argument('--cis-radius', default=500000, type=int,
					help='cis window around TSS to consider snps')
parser.add_argument('--e-proportion', default=0.3, type=float,
					help='# Proportion of samples where E==1')
parser.add_argument('--orig_ge_h2', default=0.05, type=float,
					help='Heritability in original space')
parser.add_argument('--n-pa-causal-snps', default=5, type=int,
					help='Number of proportinal amplification causal snps per gene')
parser.add_argument('--variance-scaling-factor', default=5.0, type=float,
					help='scaling of variances')
parser.add_argument('--min-snps-per-gene', default=50, type=int,
					help='Minimum number of snps per gene. else we throw out the gene.')
parser.add_argument('--random-seed', default=1, type=int,
					help='Random seed.')
args = parser.parse_args()

raw_pa_scaling_factors = np.asarray([0.0, .05, .1, .2, .5])
pa_scaling_factors = 1.0 + ((args.variance_scaling_factor)*np.copy(raw_pa_scaling_factors))

print('SIMULATING PROPORTIONAL AMPLIFICATION DATA')
print(args)

# Set seed
np.random.seed(args.random_seed)

#######################
# Get list of gene's tss to sample over
# Sorted
#######################
gene_tss_arr = get_gene_tss(args.protein_coding_gene_file, args.gene_tss_file)

#######################
# Load in genotype data for chrom 1
#######################
pgen, pvar, psam = load_in_genotype(args.plink2_genotype_stem, 1)

#######################
# Subset samples to be correct sample size
#######################
sample_indices = np.sort(np.random.choice(np.arange(psam.shape[0]), size=args.sample_size, replace=False))
psam_subset = psam[sample_indices, :]

#######################
# Simulate E variable (fixed across simulation)
#######################
EE = np.random.binomial(n=1, p=args.e_proportion, size=args.sample_size).astype(float)


#######################
# Open output file handles
#######################
# Expression output
expression_output_file = args.output_stem + '.expression.bed'
t_expr = open(expression_output_file,'w')
t_expr.write('#chr\tstart\tend\tgene_id\t' + '\t'.join(psam_subset[:,1]) + '\n')  # NOTE: USE IID NOT FID
# E-variable output
interaction_variable_output_file = args.output_stem + '.interaction_covariates.txt'
t_cov = open(interaction_variable_output_file,'w')
t_cov.write('\t' + '\t'.join(psam_subset[:,1]) + '\n')
t_cov.write('sim_EE\t' + '\t'.join(EE.astype(int).astype(str)) + '\n')
t_cov.close()
# Simulation summary
simulation_summary_output_file = args.output_stem + '.simulation_summary.txt'
t_sim = open(simulation_summary_output_file,'w')
t_sim.write('gene_id\tgene_chrom\tgene_tss\tn_cis_snps\ttotal_genetic_variance\tconstant_genetic_variance\tinteraction_genetic_variance\tpa_interaction_genetic_variance\tprop_amp_0\tprop_amp_1\n')


#######################
# Loop through genes
#######################
aa = []
bb = []
count = 0
cur_chrom = 1
for gene_tuple in tqdm(gene_tss_arr):

	# Extract relevent gene level info
	gene_id = gene_tuple[0]
	gene_chrom = gene_tuple[1]
	gene_tss = gene_tuple[2]

	# Check if we have to update genotype chromosomes
	if cur_chrom != gene_chrom:
		if gene_chrom != cur_chrom + 1:
			print('assumption errorro')
			pdb.set_trace()
		pgen, pvar, psam = load_in_genotype(args.plink2_genotype_stem, cur_chrom + 1)
		cur_chrom = cur_chrom + 1

	# Extract cis snps to the gene
	variant_positions = pvar[:, 1].astype(float)
	# Filter genotype data for this simulation
	cis_snp_indices = (variant_positions >= gene_tss - args.cis_radius) & (variant_positions < gene_tss + args.cis_radius)

	# Skip genes with too few snps
	n_cis_snps = np.sum(cis_snp_indices)
	if n_cis_snps < args.min_snps_per_gene:
		print('skipping gene because too few snps')
		continue

	# Extract cis snp genotype
	cis_genotype_mat = extract_dosages_from_reader(pgen, cis_snp_indices)

	# Filter to sample subset
	cis_genotype_mat = cis_genotype_mat[:, sample_indices]

	# Standardize genotype and filter out snps with 0 variance
	standardized_cis_genotype_mat, no_var_filter_mask = standardize_genotype(cis_genotype_mat)
	# updated number of cis snps
	n_cis_snps = standardized_cis_genotype_mat.shape[0]

	# Randomly pick genes scaling factor
	gene_pa_scaling_factor = np.random.choice(pa_scaling_factors)

	# Simulate Expression data
	if args.simulation_type == 'total_PA':
		YY, SS, constant_genetic_variance, pa_interaction_variance = simulate_total_proportional_amplification(np.transpose(standardized_cis_genotype_mat), EE, gene_pa_scaling_factor, args.n_pa_causal_snps, args.orig_ge_h2)
		interaction_genetic_variance = 0.0
		total_genetic_variance = constant_genetic_variance + pa_interaction_variance + interaction_genetic_variance



	##################
	# Save to output
	##################
	# Expr
	t_expr.write('chr' + str(gene_chrom) + '\t' + str(gene_tss) + '\t' + str(gene_tss+1) + '\t' + gene_id + '\t' + '\t'.join(YY.astype(str)) + '\n')
	# Simulation summary
	SS_0 = SS[np.where(EE==0)[0][0]]
	SS_1 = SS[np.where(EE==1)[0][0]]
	t_sim.write(gene_id + '\t' + 'chr' + str(gene_chrom) + '\t' + str(gene_tss) + '\t' + str(n_cis_snps) + '\t' + str(total_genetic_variance) + '\t' + str(constant_genetic_variance) + '\t' + str(interaction_genetic_variance) + '\t' + str(pa_interaction_variance) + '\t' + str(SS_0) + '\t' + str(SS_1) + '\n')

t_expr.close()
t_sim.close()




