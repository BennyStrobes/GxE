import numpy as np
import os
import sys
import pdb
import gzip
import pgenlib as pg


def write_expression_bed(expression_bed_file, gene_ids, gene_positions, sample_names, expr_mat):
    # Same format as expression bed used by PA_h2.py: chr, start (TSS), end, gene_id, then one column per sample
    t = gzip.open(expression_bed_file, 'wt')
    t.write("#chr\tstart\tend\tgene_id\t" + "\t".join(sample_names) + "\n")
    for gene_idx, gene_id in enumerate(gene_ids):
        t.write("chr1\t" + str(gene_positions[gene_idx]) + "\t" + str(gene_positions[gene_idx] + 1) + "\t" + gene_id + "\t" + "\t".join([str(val) for val in expr_mat[gene_idx]]) + "\n")
    t.close()


def write_plink2_genotype(plink2_genotype_stem, variant_ids, variant_positions, sample_names, genotype_mat):
    # All simulated variants are placed on chromosome 1
    pgen_writer = pg.PgenWriter((plink2_genotype_stem + '1.pgen').encode(), len(sample_names), len(variant_ids), False)
    for variant_idx in range(len(variant_ids)):
        pgen_writer.append_biallelic(genotype_mat[variant_idx].astype(np.int8))
    pgen_writer.close()

    t = open(plink2_genotype_stem + '1.pvar', 'w')
    t.write("#CHROM\tPOS\tID\tREF\tALT\n")
    for variant_idx, variant_id in enumerate(variant_ids):
        t.write("1\t" + str(variant_positions[variant_idx]) + "\t" + variant_id + "\tA\tG\n")
    t.close()

    t = open(plink2_genotype_stem + '1.psam', 'w')
    t.write("#FID\tIID\n")
    for sample_name in sample_names:
        t.write(sample_name + "\t" + sample_name + "\n")
    t.close()


def write_binary_E_file(E_var_file, sample_names, e_variable_binary):
    # Same format as binary E interaction covariate file used by PA_h2.py
    t = open(E_var_file, 'w')
    t.write("covariate\t" + "\t".join(sample_names) + "\n")
    t.write("E_binary\t" + "\t".join([str(val) for val in e_variable_binary]) + "\n")
    t.close()


def write_variant_gene_pairs_file(variant_gene_pairs_file, gene_ids, variant_ids):
    t = open(variant_gene_pairs_file, 'w')
    t.write("gene_id\tvariant_id\n")
    for test_idx in range(len(gene_ids)):
        t.write(gene_ids[test_idx] + "\t" + variant_ids[test_idx] + "\n")
    t.close()










simulation_iter = sys.argv[1]
eqtl_sample_size = int(sys.argv[2])
n_tests = int(sys.argv[3])
e_variable_ratio = float(sys.argv[4])
prop_PA_tests = float(sys.argv[5])
simulation_oracle_file = sys.argv[6]
simulated_data_stem = sys.argv[7]

# Random seed for reproducibility
np.random.seed(int(simulation_iter))

constant_genetic_variance = 0.1 # Variance of the constant genetic effect for non-PA tests
interaction_genetic_variance = 0.05 # Variance of the interaction genetic effect for nonPA tests
constant_e_effect_variance = 0.025
eta_variance = 0.1 # Variance of the eta parameter for PA tests

# First simulate E-variable
e_variable = np.random.normal(0, 1, eqtl_sample_size)
# now binarize E-variable into two groups based on quantiles
e_variable_binary = np.where(e_variable > np.quantile(e_variable, 1 - e_variable_ratio), 1, 0)
# Standardize
e_variable_binary_std = (e_variable_binary - np.mean(e_variable_binary))/np.std(e_variable_binary)


# Get simulated test assignment
test_assignment = np.random.choice([0, 1], size=n_tests, p=[1 - prop_PA_tests, prop_PA_tests]) # 1 indicates PA test, 0 indicates non-PA test

# Open output file handles
t_oracle = open(simulation_oracle_file,'w')
t_oracle.write("test_idx\ttest_type\tS1\tS2\n")

# Keep track of simulated genotype and expression (one row per test)
genotype_mat = []
expr_mat = []



# Now loop through tests and simulate data for each test
for test_idx in range(n_tests):
    # First simulate genotype data for the test
    maf = np.random.uniform(0.05, 0.5) # Randomly select a minor allele frequency between 0.05 and 0.5
    genotype_data_a1 = np.random.choice([0, 1], size=eqtl_sample_size, p=[1 - maf, maf])
    genotype_data_a2 = np.random.choice([0, 1], size=eqtl_sample_size, p=[1 - maf, maf])
    genotype_data = genotype_data_a1 + genotype_data_a2 # Add the two alleles to get genotype data (0, 1, or 2)
    std_genotype_data = (genotype_data - np.mean(genotype_data))/np.std(genotype_data) # Standardize genotype data

    # Now simulate phenotype data based on whether the test is a PA test or not
    if test_assignment[test_idx] == 0: # Non PA test
        # Simulate phenotype data as a linear combination of genotype and E-variable with some noise
        constant_genetic_effect = np.random.normal(0, scale=np.sqrt(constant_genetic_variance))
        interaction_genetic_effect = np.random.normal(0, scale=np.sqrt(interaction_genetic_variance))
        env_effect = np.random.normal(0, scale=np.sqrt(constant_e_effect_variance))

        noise_variance = 1 - (constant_genetic_variance + interaction_genetic_variance + constant_e_effect_variance)
        noise = np.random.normal(0, scale=np.sqrt(noise_variance), size=eqtl_sample_size)
        expr_data = (constant_genetic_effect * std_genotype_data) + (interaction_genetic_effect * std_genotype_data * e_variable_binary_std) + (env_effect * e_variable_binary_std) + noise

        oracle_S1 = 1.0 # no amplification for non-PA tests
        oracle_S2 = 1.0
    elif test_assignment[test_idx] == 1: # PA test
        log_ratio = np.random.normal(0.3, 0.15)   # e.g. amp_mean=0.3, amp_sd=0.15
        if np.random.rand() < 0.5:
            S1 = 1.0
            S2 = np.exp(log_ratio)
        else:
            S1 = np.exp(log_ratio)
            S2 = 1.0
        eta = np.random.normal(0, np.sqrt(eta_variance))          # e.g. eta_variance = 0.1
        per_person_S = np.where(e_variable_binary == 1, S2, S1)   # bin-specific scale
        # Normalize scales so mean(per_person_S^2) == 1; preserves the S2/S1 amplification ratio
        per_person_S = per_person_S / np.sqrt(np.mean(np.square(per_person_S)))
        genetic_effect = eta * per_person_S * std_genotype_data   # effect scales with S
        env_effect = np.random.normal(0, np.sqrt(constant_e_effect_variance))
        noise_variance = 1 - (eta_variance + constant_e_effect_variance)
        noise = np.random.normal(0, np.sqrt(noise_variance), eqtl_sample_size) * per_person_S
        expr_data = genetic_effect + env_effect * e_variable_binary_std + noise

        oracle_S1 = S1
        oracle_S2 = S2

    genotype_mat.append(genotype_data)
    expr_mat.append(expr_data)

    # Write oracle information to file
    t_oracle.write(str(test_idx) + "\t" + str(test_assignment[test_idx]) + "\t" + str(oracle_S1) + "\t" + str(oracle_S2) + "\n")
    


t_oracle.close()

genotype_mat = np.asarray(genotype_mat)
expr_mat = np.asarray(expr_mat)

# Write simulated data in the same format as the PA_h2.py input (plus a variant-gene pairs file)
# Each test is one gene paired with one variant (placed at the gene's TSS)
sample_names = ['sample_' + str(sample_idx) for sample_idx in range(eqtl_sample_size)]
gene_ids = ['gene_' + str(test_idx) for test_idx in range(n_tests)]
variant_ids = ['variant_' + str(test_idx) for test_idx in range(n_tests)]
gene_positions = [1 + test_idx*10000 for test_idx in range(n_tests)]

write_expression_bed(simulated_data_stem + '_expression.bed.gz', gene_ids, gene_positions, sample_names, expr_mat)
write_plink2_genotype(simulated_data_stem + '_genotype_chr', variant_ids, gene_positions, sample_names, genotype_mat)
write_binary_E_file(simulated_data_stem + '_E_binary.txt', sample_names, e_variable_binary)
write_variant_gene_pairs_file(simulated_data_stem + '_variant_gene_pairs.txt', gene_ids, variant_ids)
