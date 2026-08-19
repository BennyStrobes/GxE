import numpy as np
import os
import sys
import pdb
import argparse



def load_input_data(mixture_model_input_file):
    """
    Load in the input data for the PA mixture model inference
    """
    data = {}
    S1 = []
    S2 = []
    beta_1 = []
    beta_2 = []
    beta_1_se = []
    beta_2_se = []
    e1_val = []
    e2_val = []
    f = open(mixture_model_input_file, 'r')
    header = f.readline().strip().split("\t")
    for line in f:
        line = line.strip().split("\t")
        test_idx = int(line[0])
        S1.append(float(line[1]))
        S2.append(float(line[2]))
        beta_1.append(float(line[3]))
        beta_2.append(float(line[4]))
        beta_1_se.append(float(line[5]))
        beta_2_se.append(float(line[6]))
        e1_val.append(float(line[7]))
        e2_val.append(float(line[8]))
    f.close()
    data['S1'] = np.array(S1)
    data['S2'] = np.array(S2)
    data['beta_1'] = np.array(beta_1)
    data['beta_2'] = np.array(beta_2)
    data['beta_1_se'] = np.array(beta_1_se)
    data['beta_2_se'] = np.array(beta_2_se)
    data['e1_val'] = np.array(e1_val)
    data['e2_val'] = np.array(e2_val)
    return data


def univariate_normal_log_pdf(x, mean, variance):
    return (-0.5 * np.log(2.0 * np.pi * variance)) - (0.5 * np.square(x - mean) / variance)


def sample_inverse_gamma(shape_param, rate_param):
    # If X ~ Gamma(shape, rate) then 1/X ~ IG(shape, rate). Note numpy's gamma takes scale = 1/rate
    return 1.0 / np.random.gamma(shape=shape_param, scale=1.0 / rate_param)


def run_pa_mixture_model_inference(mixture_model_input_data, output_stem, max_iter=5000, n_burn_in_iters=4000, ig_a=1e-16, ig_b=1e-16, pi_prior_alpha=.1, pi_prior_beta=.1):
    # Standard (non-collapsed) Gibbs sampler for the PA mixture model.
    # Every test carries latent effects for BOTH components (eta_i for PA; beta_i, gamma_i for INT).
    # The unassigned component's latents have no likelihood term, so their full conditional is the
    # prior and they are refreshed from it each sweep. All updates are univariate:
    #   Z given latents (univariate normal likelihoods), eta given Z, beta given (gamma, Z),
    #   gamma given (beta, Z), pi given Z, variances given all latents.

    # Unpack data
    S1 = mixture_model_input_data['S1']
    S2 = mixture_model_input_data['S2']
    beta_1 = mixture_model_input_data['beta_1']
    beta_2 = mixture_model_input_data['beta_2']
    var_1 = np.square(mixture_model_input_data['beta_1_se'])
    var_2 = np.square(mixture_model_input_data['beta_2_se'])
    e1 = mixture_model_input_data['e1_val'][0]  # shared across tests
    e2 = mixture_model_input_data['e2_val'][0]
    n_tests = len(S1)

    # Initialize model parameters
    pi = 0.5
    sig_sq_beta = 0.1
    sig_sq_gamma = 0.1
    sig_sq_eta = 0.1
    etas = np.zeros(n_tests)
    betas = np.zeros(n_tests)
    gammas = np.zeros(n_tests)

    # Accumulators for posterior summaries (post burn-in)
    pa_membership_prob_sum = np.zeros(n_tests)
    eta_post_mean_sum = np.zeros(n_tests)
    beta_post_mean_sum = np.zeros(n_tests)
    gamma_post_mean_sum = np.zeros(n_tests)
    n_posterior_samples = 0
    pi_samples = []
    sig_sq_beta_samples = []
    sig_sq_gamma_samples = []
    sig_sq_eta_samples = []

    for itera in range(max_iter):
        ###################################################
        # Step 1: Sample component memberships Z given current per-test effects
        ###################################################
        log_like_pa = univariate_normal_log_pdf(beta_1, etas * S1, var_1) + univariate_normal_log_pdf(beta_2, etas * S2, var_2)
        log_like_int = univariate_normal_log_pdf(beta_1, betas + (e1 * gammas), var_1) + univariate_normal_log_pdf(beta_2, betas + (e2 * gammas), var_2)
        # Membership probabilities in log space (logsumexp) to avoid underflow
        log_weight_pa = np.log(pi) + log_like_pa
        log_weight_int = np.log(1.0 - pi) + log_like_int
        max_log_weight = np.maximum(log_weight_pa, log_weight_int)
        weight_pa = np.exp(log_weight_pa - max_log_weight)
        weight_int = np.exp(log_weight_int - max_log_weight)
        pa_membership_prob = weight_pa / (weight_pa + weight_int)
        Z = 1.0 * (np.random.rand(n_tests) < pa_membership_prob)

        ###################################################
        # Step 2: Sample per-test effects given Z
        # Assigned component: univariate conjugate normal update. Unassigned: refresh from the prior.
        ###################################################
        # PA effects (eta_i)
        eta_post_var = 1.0 / ((1.0 / sig_sq_eta) + (np.square(S1) / var_1) + (np.square(S2) / var_2))
        eta_post_mean = eta_post_var * ((S1 * beta_1 / var_1) + (S2 * beta_2 / var_2))
        eta_posterior_draw = eta_post_mean + (np.sqrt(eta_post_var) * np.random.normal(size=n_tests))
        eta_prior_draw = np.sqrt(sig_sq_eta) * np.random.normal(size=n_tests)
        etas = np.where(Z == 1, eta_posterior_draw, eta_prior_draw)

        # INT intercepts (beta_i) given current gamma_i
        beta_resid_1 = beta_1 - (e1 * gammas)
        beta_resid_2 = beta_2 - (e2 * gammas)
        beta_post_var = 1.0 / ((1.0 / sig_sq_beta) + (1.0 / var_1) + (1.0 / var_2))
        beta_post_mean = beta_post_var * ((beta_resid_1 / var_1) + (beta_resid_2 / var_2))
        beta_posterior_draw = beta_post_mean + (np.sqrt(beta_post_var) * np.random.normal(size=n_tests))
        beta_prior_draw = np.sqrt(sig_sq_beta) * np.random.normal(size=n_tests)
        betas = np.where(Z == 0, beta_posterior_draw, beta_prior_draw)

        # INT interaction effects (gamma_i) given current beta_i
        gamma_resid_1 = beta_1 - betas
        gamma_resid_2 = beta_2 - betas
        gamma_post_var = 1.0 / ((1.0 / sig_sq_gamma) + (np.square(e1) / var_1) + (np.square(e2) / var_2))
        gamma_post_mean = gamma_post_var * ((e1 * gamma_resid_1 / var_1) + (e2 * gamma_resid_2 / var_2))
        gamma_posterior_draw = gamma_post_mean + (np.sqrt(gamma_post_var) * np.random.normal(size=n_tests))
        gamma_prior_draw = np.sqrt(sig_sq_gamma) * np.random.normal(size=n_tests)
        gammas = np.where(Z == 0, gamma_posterior_draw, gamma_prior_draw)

        ###################################################
        # Step 3: Sample mixture proportion pi
        ###################################################
        n_pa = np.sum(Z == 1)
        n_int = n_tests - n_pa
        pi = np.random.beta(pi_prior_alpha + n_pa, pi_prior_beta + n_int)

        ###################################################
        # Step 4: Sample effect variances (inverse-gamma conjugate updates)
        # All tests' latents are in the state (unassigned ones are prior draws), so sums run over all tests
        ###################################################
        sig_sq_eta = sample_inverse_gamma(ig_a + (0.5 * n_tests), ig_b + (0.5 * np.sum(np.square(etas))))
        sig_sq_beta = sample_inverse_gamma(ig_a + (0.5 * n_tests), ig_b + (0.5 * np.sum(np.square(betas))))
        sig_sq_gamma = sample_inverse_gamma(ig_a + (0.5 * n_tests), ig_b + (0.5 * np.sum(np.square(gammas))))

        ###################################################
        # Record samples
        ###################################################
        if itera >= n_burn_in_iters:
            # Rao-Blackwellized posterior PA membership probability (average of membership probs, not binary draws)
            pa_membership_prob_sum = pa_membership_prob_sum + pa_membership_prob
            # Conditional posterior means of the per-test effects (as-if assigned to the relevant component)
            eta_post_mean_sum = eta_post_mean_sum + eta_post_mean
            beta_post_mean_sum = beta_post_mean_sum + beta_post_mean
            gamma_post_mean_sum = gamma_post_mean_sum + gamma_post_mean
            n_posterior_samples = n_posterior_samples + 1
            pi_samples.append(pi)
            sig_sq_beta_samples.append(sig_sq_beta)
            sig_sq_gamma_samples.append(sig_sq_gamma)
            sig_sq_eta_samples.append(sig_sq_eta)

    # Write per-test posterior summaries
    posterior_pa_prob = pa_membership_prob_sum / n_posterior_samples
    posterior_eta = eta_post_mean_sum / n_posterior_samples
    posterior_beta = beta_post_mean_sum / n_posterior_samples
    posterior_gamma = gamma_post_mean_sum / n_posterior_samples
    t_per_test = open(output_stem + '_per_test_posterior_summary.txt', 'w')
    t_per_test.write("test_idx\tPA_posterior_prob\tposterior_mean_eta\tposterior_mean_beta\tposterior_mean_gamma\n")
    for test_idx in range(n_tests):
        t_per_test.write(str(test_idx) + "\t" + str(posterior_pa_prob[test_idx]) + "\t" + str(posterior_eta[test_idx]) + "\t" + str(posterior_beta[test_idx]) + "\t" + str(posterior_gamma[test_idx]) + "\n")
    t_per_test.close()

    # Write posterior mean and 95% credible interval for the global parameters
    t_global = open(output_stem + '_global_parameter_posterior_summary.txt', 'w')
    t_global.write("parameter\tposterior_mean\tci_95_lower\tci_95_upper\n")
    for param_name, param_samples in [('pi', pi_samples), ('sig_sq_beta', sig_sq_beta_samples), ('sig_sq_gamma', sig_sq_gamma_samples), ('sig_sq_eta', sig_sq_eta_samples)]:
        param_samples = np.array(param_samples)
        ci_lower, ci_upper = np.percentile(param_samples, [2.5, 97.5])
        t_global.write(param_name + "\t" + str(np.mean(param_samples)) + "\t" + str(ci_lower) + "\t" + str(ci_upper) + "\n")
    t_global.close()

    return


parser = argparse.ArgumentParser()
parser.add_argument('--mixture-model-input-file', type=str, required=True, help='File containing per-test summary statistics (S1, S2, betas, ses, e-values) used as input to the PA mixture model')
parser.add_argument('--output-stem', type=str, required=True, help='Stem of output files')
parser.add_argument('--seed', type=int, default=1, help='Random seed for the Gibbs sampler')
parser.add_argument('--ig-a', type=float, default=1e-5, help='Shape parameter of the inverse-gamma hyperprior on the effect variances')
parser.add_argument('--ig-b', type=float, default=1e-5, help='Rate parameter of the inverse-gamma hyperprior on the effect variances')
parser.add_argument('--pi-prior-alpha', type=float, default=.1, help='Alpha parameter of the Beta hyperprior on the mixture proportion pi')
parser.add_argument('--pi-prior-beta', type=float, default=.1, help='Beta parameter of the Beta hyperprior on the mixture proportion pi')
args = parser.parse_args()

mixture_model_input_file = args.mixture_model_input_file
output_stem = args.output_stem

np.random.seed(args.seed)

# Load in the input data
mixture_model_input_data = load_input_data(mixture_model_input_file)


run_pa_mixture_model_inference(mixture_model_input_data, output_stem, ig_a=args.ig_a, ig_b=args.ig_b, pi_prior_alpha=args.pi_prior_alpha, pi_prior_beta=args.pi_prior_beta)
