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


def initialize_model_params(mixture_model_input_data):
    model_params = {}
    model_params['n_tests'] = len(mixture_model_input_data['S1'])
    model_params['pi'] = 0.5  # Initial guess for the proportion of PA tests
    model_params['Z'] = np.random.choice([0, 1], size=len(mixture_model_input_data['S1']), p=[1-model_params['pi'], model_params['pi']])  # Random initial assignment of tests to PA or non-PA
    model_params['betas'] = np.zeros(model_params['n_tests'])
    model_params['etas'] = np.zeros(model_params['n_tests'])
    model_params['gammas'] = np.zeros(model_params['n_tests'])
    model_params['sig_sq_beta'] = 0.1
    model_params['sig_sq_eta'] = 0.1
    model_params['sig_sq_gamma'] = 0.1
    return model_params


def bivariate_zero_mean_normal_log_pdf(x1, x2, c11, c12, c22):
    """
    Log density of N(0, [[c11, c12], [c12, c22]]) evaluated at (x1, x2). Vectorized over tests.
    """
    det = (c11 * c22) - np.square(c12)
    quad_form = ((c22 * np.square(x1)) - (2.0 * c12 * x1 * x2) + (c11 * np.square(x2))) / det
    return -np.log(2.0 * np.pi) - (0.5 * np.log(det)) - (0.5 * quad_form)


def sample_inverse_gamma(shape_param, rate_param):
    # If X ~ Gamma(shape, rate) then 1/X ~ IG(shape, rate). Note numpy's gamma takes scale = 1/rate
    return 1.0 / np.random.gamma(shape=shape_param, scale=1.0 / rate_param)


def run_pa_mixture_model_inference(mixture_model_input_data, output_stem, max_iter=5000, n_burn_in_iters=4000, ig_a=1e-16, ig_b=1e-16, pi_prior_alpha=.1, pi_prior_beta=.1):
    # ig_a, ig_b: inverse-gamma hyperprior on the effect variances (conjugate to the Gaussian effect priors)
    #             default IG(1e-3, 1e-3) is weak/near-flat; Gibbs update: IG(a + n_k/2, b + sum(effects^2)/2)
    # pi_prior_alpha, pi_prior_beta: Beta hyperprior on the mixture proportion pi (default Beta(1,1) = uniform)
    model_params = initialize_model_params(mixture_model_input_data)

    # Unpack data
    S1 = mixture_model_input_data['S1']
    S2 = mixture_model_input_data['S2']
    beta_1 = mixture_model_input_data['beta_1']
    beta_2 = mixture_model_input_data['beta_2']
    var_1 = np.square(mixture_model_input_data['beta_1_se'])
    var_2 = np.square(mixture_model_input_data['beta_2_se'])
    e1 = mixture_model_input_data['e1_val'][0]  # shared across tests
    e2 = mixture_model_input_data['e2_val'][0]
    n_tests = model_params['n_tests']

    # Unpack model parameters
    pi = model_params['pi']
    sig_sq_beta = model_params['sig_sq_beta']
    sig_sq_gamma = model_params['sig_sq_gamma']
    sig_sq_eta = model_params['sig_sq_eta']

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
        # Step 1: Sample component memberships Z (collapsed: per-test effects integrated out)
        ###################################################
        # PA marginal: beta_hat_i ~ N(0, sig_sq_eta * s_i s_i^T + diag(se^2))
        log_like_pa = bivariate_zero_mean_normal_log_pdf(
            beta_1, beta_2,
            (sig_sq_eta * np.square(S1)) + var_1,
            sig_sq_eta * S1 * S2,
            (sig_sq_eta * np.square(S2)) + var_2)
        # INT marginal: beta_hat_i ~ N(0, X diag(sig_sq_beta, sig_sq_gamma) X^T + diag(se^2)) with X = [[1, e1], [1, e2]]
        log_like_int = bivariate_zero_mean_normal_log_pdf(
            beta_1, beta_2,
            sig_sq_beta + (sig_sq_gamma * e1 * e1) + var_1,
            sig_sq_beta + (sig_sq_gamma * e1 * e2),
            sig_sq_beta + (sig_sq_gamma * e2 * e2) + var_2)
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
        # (sampled for all tests vectorized; draws for tests assigned to the other component are never used)
        ###################################################
        # PA effects (eta_i): conjugate normal update
        eta_post_var = 1.0 / ((1.0 / sig_sq_eta) + (np.square(S1) / var_1) + (np.square(S2) / var_2))
        eta_post_mean = eta_post_var * ((S1 * beta_1 / var_1) + (S2 * beta_2 / var_2))
        etas = eta_post_mean + (np.sqrt(eta_post_var) * np.random.normal(size=n_tests))

        # INT effects (beta_i, gamma_i): bivariate conjugate normal update
        # Posterior precision A = prior precision + X^T diag(1/se^2) X
        A11 = (1.0 / sig_sq_beta) + (1.0 / var_1) + (1.0 / var_2)
        A12 = (e1 / var_1) + (e2 / var_2)
        A22 = (1.0 / sig_sq_gamma) + (np.square(e1) / var_1) + (np.square(e2) / var_2)
        rhs1 = (beta_1 / var_1) + (beta_2 / var_2)
        rhs2 = (e1 * beta_1 / var_1) + (e2 * beta_2 / var_2)
        det_A = (A11 * A22) - np.square(A12)
        # Posterior covariance C = A^{-1}
        C11 = A22 / det_A
        C12 = -A12 / det_A
        C22 = A11 / det_A
        beta_post_mean = (C11 * rhs1) + (C12 * rhs2)
        gamma_post_mean = (C12 * rhs1) + (C22 * rhs2)
        # Sample via cholesky of C: [[L11, 0], [L21, L22]]
        L11 = np.sqrt(C11)
        L21 = C12 / L11
        L22 = np.sqrt(C22 - (np.square(C12) / C11))
        std_normal_1 = np.random.normal(size=n_tests)
        std_normal_2 = np.random.normal(size=n_tests)
        betas = beta_post_mean + (L11 * std_normal_1)
        gammas = gamma_post_mean + (L21 * std_normal_1) + (L22 * std_normal_2)

        ###################################################
        # Step 3: Sample mixture proportion pi
        ###################################################
        n_pa = np.sum(Z == 1)
        n_int = n_tests - n_pa
        pi = np.random.beta(pi_prior_alpha + n_pa, pi_prior_beta + n_int)

        ###################################################
        # Step 4: Sample effect variances (inverse-gamma conjugate updates; only tests assigned to the component)
        ###################################################
        sig_sq_eta = sample_inverse_gamma(ig_a + (0.5 * n_pa), ig_b + (0.5 * np.sum(np.square(etas[Z == 1]))))
        sig_sq_beta = sample_inverse_gamma(ig_a + (0.5 * n_int), ig_b + (0.5 * np.sum(np.square(betas[Z == 0]))))
        sig_sq_gamma = sample_inverse_gamma(ig_a + (0.5 * n_int), ig_b + (0.5 * np.sum(np.square(gammas[Z == 0]))))

        ###################################################
        # Record samples
        ###################################################
        if itera >= n_burn_in_iters:
            # Rao-Blackwellized posterior PA membership probability (average of membership probs, not binary draws)
            pa_membership_prob_sum = pa_membership_prob_sum + pa_membership_prob
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
parser.add_argument('--ig-a', type=float, default=1e-5, help='Shape parameter of the inverse-gamma hyperprior on the effect variances (default 1e-3: weak/near-flat)')
parser.add_argument('--ig-b', type=float, default=1e-5, help='Rate parameter of the inverse-gamma hyperprior on the effect variances (default 1e-3: weak/near-flat)')
parser.add_argument('--pi-prior-alpha', type=float, default=.1, help='Alpha parameter of the Beta hyperprior on the mixture proportion pi')
parser.add_argument('--pi-prior-beta', type=float, default=.1, help='Beta parameter of the Beta hyperprior on the mixture proportion pi')
args = parser.parse_args()

mixture_model_input_file = args.mixture_model_input_file
output_stem = args.output_stem

np.random.seed(args.seed)

# Load in the input data
mixture_model_input_data = load_input_data(mixture_model_input_file)


run_pa_mixture_model_inference(mixture_model_input_data, output_stem, ig_a=args.ig_a, ig_b=args.ig_b, pi_prior_alpha=args.pi_prior_alpha, pi_prior_beta=args.pi_prior_beta)
