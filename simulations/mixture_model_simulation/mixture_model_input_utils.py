import numpy as np


# Shared by simulate -> make_pa_mixture_model_input_from_data.py -> inference, so that
# simulated and real data produce mixture model input in exactly the same way.

MIXTURE_MODEL_INPUT_HEADER = "test_idx\tS1\tS2\tbeta_1\tbeta_2\tbeta_1_se\tbeta_2_se\te1_val\te2_val\n"


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


def compute_test_summary_stats(std_genotype, expr, e_variable_binary):
    # std_genotype: genotype standardized globally (across both E bins)
    # expr: expression for the same samples
    # e_variable_binary: 0/1 E variable. E bin 1 is E == 0, E bin 2 is E == 1
    # Returns (raw_S_row, geno_resid_S_row): values of each mixture model input line (excluding test_idx)
    e_variable_binary_std = (e_variable_binary - np.mean(e_variable_binary))/np.std(e_variable_binary)

    # Run association test seperately in each E-variable group
    group_0_indices = np.where(e_variable_binary == 0)[0]
    group_1_indices = np.where(e_variable_binary == 1)[0]
    group_0_genotype = std_genotype[group_0_indices]
    group_1_genotype = std_genotype[group_1_indices]
    group_0_expr = expr[group_0_indices]
    group_1_expr = expr[group_1_indices]

    # Calculate association statistics for each group
    group_0_beta, group_0_se, group_0_residual_sd = group_assoc_stats(group_0_genotype, group_0_expr)
    group_1_beta, group_1_se, group_1_residual_sd = group_assoc_stats(group_1_genotype, group_1_expr)

    # Estimate bin-specific scales (S1, S2) from the data
    S1_est = np.sqrt(np.var(group_0_expr))
    S2_est = np.sqrt(np.var(group_1_expr))

    e1_val = np.mean(e_variable_binary_std[group_0_indices])
    e2_val = np.mean(e_variable_binary_std[group_1_indices])

    raw_S_row = [S1_est, S2_est, group_0_beta, group_1_beta, group_0_se, group_1_se, e1_val, e2_val]
    geno_resid_S_row = [group_0_residual_sd, group_1_residual_sd, group_0_beta, group_1_beta, group_0_se, group_1_se, e1_val, e2_val]
    return raw_S_row, geno_resid_S_row


def format_input_line(test_idx, row):
    return str(test_idx) + "\t" + "\t".join([str(val) for val in row]) + "\n"
