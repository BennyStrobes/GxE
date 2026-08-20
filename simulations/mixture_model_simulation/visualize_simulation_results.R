args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_pi_estimate_barplot <- function(pi_df, plot_title) {
	# True simulated pi for each bar (x values are the simulated fractions)
	pi_df$true_pi <- as.numeric(as.character(pi_df$simulated_prop_PA_tests))
	pi_df$simulated_prop_PA_tests <- factor(pi_df$simulated_prop_PA_tests, levels=pi_df$simulated_prop_PA_tests)

	p <- ggplot(pi_df, aes(x=simulated_prop_PA_tests, y=avg_posterior_mean_pi)) +
		geom_bar(stat="identity", fill="#5598e7", color="#1c5cab", linewidth=0.3, width=0.62) +
		geom_errorbar(aes(ymin=ci_95_lower, ymax=ci_95_upper), width=0.18, linewidth=0.4, color="gray20") +
		# Dashed line across each bar marking the simulated (true) pi
		geom_errorbar(aes(ymin=true_pi, ymax=true_pi), width=0.62, linetype="dashed", linewidth=0.45, color="gray10") +
		scale_y_continuous(limits=c(0, 1.0), expand=c(0, 0)) +
		labs(x="Simulated fraction of PA tests", y="Estimated fraction of PA tests", title=plot_title) +
		figure_theme() +
		theme(plot.subtitle = element_text(size=9, color="gray30"))
	return(p)
}


load_in_calibration_data <- function(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio, n_bins=10) {
	# Pool per-test posterior PA probabilities (and true PA labels from the oracle) across all simulation iterations
	all_probs <- c()
	all_labels <- c()
	for (simulation_iter in 1:n_simulations) {
		simulation_name_stem <- paste0("simulation_iter_", simulation_iter, "_eqtl_sample_size_", eqtl_sample_size, "_n_tests_", n_tests, "_e_variable_ratio_", e_variable_ratio, "_prop_PA_tests_", prop_PA_tests)
		posterior_file <- paste0(simulation_results_dir, simulation_name_stem, inference_method_suffix, "_per_test_posterior_summary.txt")
		oracle_file <- paste0(simulated_data_dir, simulation_name_stem, "_oracle.txt")
		posterior_df <- read.table(posterior_file, header=TRUE, sep="\t")
		oracle_df <- read.table(oracle_file, header=TRUE, sep="\t")
		# Align oracle rows to posterior rows via test_idx
		matched_labels <- oracle_df$test_type[match(posterior_df$test_idx, oracle_df$test_idx)]
		all_probs <- c(all_probs, posterior_df$PA_posterior_prob)
		all_labels <- c(all_labels, matched_labels)
	}
	# Bin tests by posterior PA probability
	bin_assignments <- cut(all_probs, breaks=seq(0.0, 1.0, length.out=(n_bins+1)), include.lowest=TRUE)
	expected_arr <- c()
	observed_arr <- c()
	ci_lower_arr <- c()
	ci_upper_arr <- c()
	n_bin_tests_arr <- c()
	for (bin_level in levels(bin_assignments)) {
		bin_indices <- bin_assignments == bin_level
		n_bin_tests <- sum(bin_indices)
		if (n_bin_tests == 0) {
			next
		}
		# Expected: average posterior probability in the bin. Observed: fraction of tests in bin that are truly PA
		expected_prob <- mean(all_probs[bin_indices])
		observed_frac <- mean(all_labels[bin_indices])
		se_observed <- sqrt((observed_frac*(1.0-observed_frac))/n_bin_tests)
		expected_arr <- c(expected_arr, expected_prob)
		observed_arr <- c(observed_arr, observed_frac)
		ci_lower_arr <- c(ci_lower_arr, max(0.0, observed_frac - (1.96*se_observed)))
		ci_upper_arr <- c(ci_upper_arr, min(1.0, observed_frac + (1.96*se_observed)))
		n_bin_tests_arr <- c(n_bin_tests_arr, n_bin_tests)
	}
	calibration_df <- data.frame(expected_pa_prob=expected_arr, observed_pa_fraction=observed_arr, ci_95_lower=ci_lower_arr, ci_95_upper=ci_upper_arr, n_tests_in_bin=n_bin_tests_arr)
	return(calibration_df)
}


make_calibration_plot <- function(calibration_df, plot_title) {
	p <- ggplot(calibration_df, aes(x=expected_pa_prob, y=observed_pa_fraction)) +
		# y=x line corresponding to perfect calibration
		geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", linewidth=0.4) +
		geom_errorbar(aes(ymin=ci_95_lower, ymax=ci_95_upper), width=0.0, linewidth=0.4, color="gray20") +
		geom_line(color="#2a78d6", linewidth=0.4) +
		geom_point(color="#2a78d6", size=2.0) +
		scale_x_continuous(limits=c(0.0, 1.0)) +
		scale_y_continuous(limits=c(0.0, 1.0)) +
		labs(x="Expected PA probability\n(bin average posterior probability)", y="Observed fraction of PA tests", title=plot_title) +
		figure_theme()
	return(p)
}


load_in_posterior_prob_distribution_data <- function(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio) {
	# Pool per-test posterior PA probabilities (and true PA labels from the oracle) across all simulation iterations
	all_probs <- c()
	all_labels <- c()
	for (simulation_iter in 1:n_simulations) {
		simulation_name_stem <- paste0("simulation_iter_", simulation_iter, "_eqtl_sample_size_", eqtl_sample_size, "_n_tests_", n_tests, "_e_variable_ratio_", e_variable_ratio, "_prop_PA_tests_", prop_PA_tests)
		posterior_file <- paste0(simulation_results_dir, simulation_name_stem, inference_method_suffix, "_per_test_posterior_summary.txt")
		oracle_file <- paste0(simulated_data_dir, simulation_name_stem, "_oracle.txt")
		posterior_df <- read.table(posterior_file, header=TRUE, sep="\t")
		oracle_df <- read.table(oracle_file, header=TRUE, sep="\t")
		# Align oracle rows to posterior rows via test_idx
		matched_labels <- oracle_df$test_type[match(posterior_df$test_idx, oracle_df$test_idx)]
		all_probs <- c(all_probs, posterior_df$PA_posterior_prob)
		all_labels <- c(all_labels, matched_labels)
	}
	distribution_df <- data.frame(pa_posterior_prob=all_probs, true_test_type=factor(ifelse(all_labels == 1, "PA", "non-PA"), levels=c("non-PA", "PA")))
	return(distribution_df)
}


make_posterior_prob_distribution_plot <- function(distribution_df, plot_title) {
	p <- ggplot(distribution_df, aes(x=pa_posterior_prob, fill=true_test_type)) +
		geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.55, breaks=seq(0.0, 1.0, length.out=31), color=NA) +
		scale_fill_manual(values=c("non-PA"="#2a78d6", "PA"="#eb6834"), name="True test type") +
		scale_x_continuous(limits=c(0.0, 1.0)) +
		labs(x="Posterior PA probability", y="Density", title=plot_title) +
		figure_theme()
	return(p)
}


load_in_avg_pi_estimates <- function(simulation_results_dir, inference_method_suffix, prop_PA_tests_arr, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio) {
	simulated_prop_arr <- c()
	avg_pi_arr <- c()
	ci_lower_arr <- c()
	ci_upper_arr <- c()
	for (prop_PA_tests in prop_PA_tests_arr) {
		# Extract posterior mean pi from each simulation iteration
		pi_estimates <- c()
		for (simulation_iter in 1:n_simulations) {
			file_name <- paste0(simulation_results_dir, "simulation_iter_", simulation_iter, "_eqtl_sample_size_", eqtl_sample_size, "_n_tests_", n_tests, "_e_variable_ratio_", e_variable_ratio, "_prop_PA_tests_", prop_PA_tests, inference_method_suffix, "_global_parameter_posterior_summary.txt")
			df <- read.table(file_name, header=TRUE, sep="\t")
			pi_estimates <- c(pi_estimates, df$posterior_mean[as.character(df$parameter) == "pi"])
		}
		# Average across simulations and 95% CI of the average (normal approximation)
		avg_pi <- mean(pi_estimates)
		se_pi <- sd(pi_estimates)/sqrt(length(pi_estimates))
		simulated_prop_arr <- c(simulated_prop_arr, prop_PA_tests)
		avg_pi_arr <- c(avg_pi_arr, avg_pi)
		ci_lower_arr <- c(ci_lower_arr, avg_pi - (1.96*se_pi))
		ci_upper_arr <- c(ci_upper_arr, avg_pi + (1.96*se_pi))
	}
	pi_df <- data.frame(simulated_prop_PA_tests=simulated_prop_arr, avg_posterior_mean_pi=avg_pi_arr, ci_95_lower=ci_lower_arr, ci_95_upper=ci_upper_arr)
	return(pi_df)
}










simulation_results_dir <- args[1]
visualization_results_dir <- args[2]
simulated_data_dir <- args[3]

# Simulation parameters (must match driver_key.sh / run_single_simulation.sh)
n_simulations <- 50
eqtl_sample_size <- "500"
n_tests <- "2000"
e_variable_ratio <- "0.4"
prop_PA_tests_arr <- c("0.2", "0.4", "0.6", "0.8")


######################
# Plot ability to estimate pi (fraction of PA tests) across simulations
######################
# Load in average posterior mean pi (and 95% CI of the average) across simulations, per simulated pi
pi_df <- load_in_avg_pi_estimates(simulation_results_dir, "_mixture_model_results", prop_PA_tests_arr, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
pi_geno_resid_df <- load_in_avg_pi_estimates(simulation_results_dir, "_geno_resid_mixture_model_results", prop_PA_tests_arr, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)

# Make bar plots of estimated pi (with 95% CIs) vs simulated pi and save to pdf
pi_barplot <- make_pi_estimate_barplot(pi_df, "Mixture model inference (total S)")
output_file <- paste0(visualization_results_dir, "estimated_pi_barplot_mixture_model_results.pdf")
ggsave(pi_barplot, file=output_file, width=7.2, height=4.0, units="in")

pi_geno_resid_barplot <- make_pi_estimate_barplot(pi_geno_resid_df, "Mixture model inference (genotype-residualized S)")
output_file <- paste0(visualization_results_dir, "estimated_pi_barplot_geno_resid_mixture_model_results.pdf")
ggsave(pi_geno_resid_barplot, file=output_file, width=7.2, height=4.0, units="in")


######################
# Plot calibration of posterior PA membership probabilities across simulations
######################
inference_method_suffixes <- c("_mixture_model_results", "_geno_resid_mixture_model_results")
for (inference_method_suffix in inference_method_suffixes) {
	prop_PA_tests <- "0.2"
	calibration_df <- load_in_calibration_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	calibration_plot2 <- make_calibration_plot(calibration_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.4"
	calibration_df <- load_in_calibration_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	calibration_plot4 <- make_calibration_plot(calibration_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.6"
	calibration_df <- load_in_calibration_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	calibration_plot6 <- make_calibration_plot(calibration_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.8"
	calibration_df <- load_in_calibration_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	calibration_plot8 <- make_calibration_plot(calibration_df, paste0("Simulated pi = ", prop_PA_tests))

	joint_calibration_plot <- plot_grid(calibration_plot2, calibration_plot4, calibration_plot6, calibration_plot8, nrow=2, ncol=2)
	output_file <- paste0(visualization_results_dir, "posterior_PA_prob", inference_method_suffix, "_calibration_plots.pdf")
	ggsave(joint_calibration_plot, file=output_file, width=7.2, height=6.0, units="in")
}

######################
# Plot the distribution of posterior PA probabilities across simulations, for each simulated pi
######################
for (inference_method_suffix in inference_method_suffixes) {
	prop_PA_tests <- "0.2"
	distribution_df <- load_in_posterior_prob_distribution_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	distribution_plot2 <- make_posterior_prob_distribution_plot(distribution_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.4"
	distribution_df <- load_in_posterior_prob_distribution_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	distribution_plot4 <- make_posterior_prob_distribution_plot(distribution_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.6"
	distribution_df <- load_in_posterior_prob_distribution_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	distribution_plot6 <- make_posterior_prob_distribution_plot(distribution_df, paste0("Simulated pi = ", prop_PA_tests))
	prop_PA_tests <- "0.8"
	distribution_df <- load_in_posterior_prob_distribution_data(simulation_results_dir, simulated_data_dir, inference_method_suffix, prop_PA_tests, n_simulations, eqtl_sample_size, n_tests, e_variable_ratio)
	distribution_plot8 <- make_posterior_prob_distribution_plot(distribution_df, paste0("Simulated pi = ", prop_PA_tests))

    legend <- get_legend(distribution_plot2 + theme(legend.position="bottom"))
	joint_distribution_plot <- plot_grid(distribution_plot2+theme(legend.position="none"), distribution_plot4 +theme(legend.position="none"), distribution_plot6 +theme(legend.position="none"), distribution_plot8 +theme(legend.position="none"), nrow=2, ncol=2)
	
    joint_distribution_plot <- plot_grid(joint_distribution_plot, legend, nrow=2, ncol=1, rel_heights=c(1, 0.1))
    output_file <- paste0(visualization_results_dir, "posterior_PA_prob", inference_method_suffix, "_distribution_plots.pdf")
	ggsave(joint_distribution_plot, file=output_file, width=7.2, height=6.5, units="in")
}
