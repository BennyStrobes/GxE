import numpy as np
import os
import sys
import pdb


def group_assoc_stats(g, y):
    # Center genotype and expression within group (equivalent to including an intercept)
    g_c = g - np.mean(g)
    y_c = y - np.mean(y)
    beta = np.sum(g_c * y_c) / np.sum(np.square(g_c))
    resid = y_c - beta * g_c
    resid_var = np.sum(np.square(resid)) / (len(y) - 2)  # df: intercept + slope
    se = np.sqrt(resid_var / np.sum(np.square(g_c)))
    residual_sd = np.sqrt(resid_var)
    return beta, se, residual_sd










simulation_iter = sys.argv[1]
eqtl_sample_size = int(sys.argv[2])
n_tests = int(sys.argv[3])
e_variable_ratio = float(sys.argv[4])
prop_PA_tests = float(sys.argv[5])
simulation_oracle_file = sys.argv[6]
pa_mixture_model_input_file = sys.argv[7]
pa_mixture_model_geno_resid_S_input_file = sys.argv[8]

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
t_mixture_model_input = open(pa_mixture_model_input_file,'w')
t_mixture_model_geno_resid_S_input = open(pa_mixture_model_geno_resid_S_input_file,'w')
t_oracle.write("test_idx\ttest_type\tS1\tS2\n")
t_mixture_model_input.write("test_idx\tS1\tS2\tbeta_1\tbeta_2\tbeta_1_se\tbeta_2_se\te1_val\te2_val\n")
t_mixture_model_geno_resid_S_input.write("test_idx\tS1\tS2\tbeta_1\tbeta_2\tbeta_1_se\tbeta_2_se\te1_val\te2_val\n")



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

    # Run association test seperately in each E-variable group
    group_0_indices = np.where(e_variable_binary == 0)[0]
    group_1_indices = np.where(e_variable_binary == 1)[0]
    group_0_genotype = std_genotype_data[group_0_indices]
    group_1_genotype = std_genotype_data[group_1_indices]
    group_0_expr = expr_data[group_0_indices]
    group_1_expr = expr_data[group_1_indices]

    # Calculate association statistics for each group
    group_0_beta, group_0_se, group_0_residual_sd = group_assoc_stats(group_0_genotype, group_0_expr)
    group_1_beta, group_1_se, group_1_residual_sd = group_assoc_stats(group_1_genotype, group_1_expr)

    # Estimate bin-specific scales (S1, S2) from the data
    S1_est = np.sqrt(np.var(group_0_expr))
    S2_est = np.sqrt(np.var(group_1_expr))

    # Write data required for input to Mixture model inference
    t_mixture_model_input.write(str(test_idx) + "\t" + str(S1_est) + "\t" + str(S2_est) + "\t" + str(group_0_beta) + "\t" + str(group_1_beta) + "\t" + str(group_0_se) + "\t" + str(group_1_se) + "\t" + str(np.mean(e_variable_binary_std[group_0_indices])) + "\t" + str(np.mean(e_variable_binary_std[group_1_indices])) + "\n")
    t_mixture_model_geno_resid_S_input.write(str(test_idx) + "\t" + str(group_0_residual_sd) + "\t" + str(group_1_residual_sd) + "\t" + str(group_0_beta) + "\t" + str(group_1_beta) + "\t" + str(group_0_se) + "\t" + str(group_1_se) + "\t" + str(np.mean(e_variable_binary_std[group_0_indices])) + "\t" + str(np.mean(e_variable_binary_std[group_1_indices])) + "\n")

    # Write oracle information to file
    t_oracle.write(str(test_idx) + "\t" + str(test_assignment[test_idx]) + "\t" + str(oracle_S1) + "\t" + str(oracle_S2) + "\n")
    


t_oracle.close()
t_mixture_model_input.close()
t_mixture_model_geno_resid_S_input.close()