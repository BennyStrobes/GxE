args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
library(hash)
library(RColorBrewer)
options(warn=1)








figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



load_in_interaction_variance_estimates_for_fixed_sample_size_inference_method_but_vary_scaling_factors <- function(simulation_results_dir, simulation_version, sample_size, inference_method, pa_scaling_factors) {
	all_dfs <- list()
	for (scaling_factor_iter in 1:length(pa_scaling_factors)) {
		scaling_factor = pa_scaling_factors[scaling_factor_iter]

		file_name <- paste0(simulation_results_dir, simulation_version, "_simulation_propAmpScale_", scaling_factor, "_sampleSize_", sample_size, "_organized_simulation_results.txt")

		df <- read.table(file_name, header=TRUE, sep="\t")
		df = df[as.character(df$inference_method) == inference_method,]
		df = df[as.character(df$variance_component) %in% c("interaction_genetic_variance","PA_interaction_genetic_variance"),]

    	# Add scaling factor column
    	df$scaling_factor <- rep(scaling_factor, dim(df)[1])


    	# Store in list
    	all_dfs[[as.character(scaling_factor)]] <- df
    }
    # Concatenate all data frames vertically
  	combined_df <- do.call(rbind, all_dfs)
  	combined_df$scaling_factor = factor(combined_df$scaling_factor, levels=pa_scaling_factors)
  	return(combined_df)
}


make_interaction_variance_estimate_plot_old <- function(df, pa_scaling_factors, simulation_version, sample_size, inference_method) {
	pp<-ggplot(data=df, aes(x=scaling_factor, y=estimated_variance, fill=variance_component)) +
  		geom_bar(stat="identity", position=position_dodge()) +
  		geom_errorbar(aes(ymin=estimated_variance_lb, ymax=estimated_variance_ub), width=.4, position=position_dodge(.9))  +
  		#scale_fill_manual(values=c("skyblue", "grey", "grey"))+
  		figure_theme() +
  		labs(x="Proportional amplification scaling factor", y="h2", fill="") #+
  		#geom_hline(yintercept=sim_h2, linetype=2) 

  	return(pp)

}

make_interaction_variance_estimate_plot <- function(df, pa_scaling_factors, simulation_version, sample_size, inference_method) {
dodge <- position_dodge(width = 0.9)

pp <- ggplot(df, aes(x = scaling_factor,
                     y = estimated_variance,
                     fill = variance_component)) +
  geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(aes(ymin = estimated_variance_lb,
                    ymax = estimated_variance_ub),
                width = 0.25,
                position = dodge) +
  geom_point(aes(y = true_simulated_variance,
                 color = variance_component),
             position = dodge,
             shape = 95,   # short horizontal line symbol
             size = 10,    # controls the line length
             stroke = 1, show.legend=FALSE) +
  scale_fill_manual(values=c("grey", "steelblue2")) + 
  scale_color_manual(values = c("interaction_genetic_variance" = "purple3",
                                "PA_interaction_genetic_variance" = "firebrick3")) +
  figure_theme() +
  labs(x = "Proportional amplification scaling factor",
       y = "h2",
       fill = "",
       color = "true simulated",
       title=paste0(inference_method, "   (", simulation_version, " / N=", sample_size, ")"))

  return(pp)
}


simulation_results_dir = args[1]
visualization_results_dir = args[2]


#################################
#################################
# Bar plot showing estimate interaction variance (standard interaction variance vs PA-interaction)
simulation_version = "total_PA"
sample_size = "450"
pa_scaling_factors = c("1.5", "2.0", "4.0", "6.0", "8.0", "10.0")

# Inference using no proportional amp
inference_method = "interaction-VCA"
# Load in data
df <- load_in_interaction_variance_estimates_for_fixed_sample_size_inference_method_but_vary_scaling_factors(simulation_results_dir, simulation_version, sample_size, inference_method, pa_scaling_factors)
pp1 <- make_interaction_variance_estimate_plot(df, pa_scaling_factors, simulation_version, sample_size, inference_method)
# Inference using prop amp
inference_method = "PA+interaction-VCA"
# Load in data
df <- load_in_interaction_variance_estimates_for_fixed_sample_size_inference_method_but_vary_scaling_factors(simulation_results_dir, simulation_version, sample_size, inference_method, pa_scaling_factors)
pp2 <- make_interaction_variance_estimate_plot(df, pa_scaling_factors, simulation_version, sample_size, inference_method)
# Plot
legender=get_legend(pp2)
joint_pp <- plot_grid(pp1 + theme(legend.position="none"), pp2+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.34))
# save
output_file <- paste0(visualization_results_dir, "interaction_variance_estimates_bar_plot_", simulation_version, "_", sample_size, ".pdf")
ggsave(joint_pp, file=output_file, width=7.2, height=5.5, units="in")

#################################
#################################
# Bar plot showing estimate interaction variance (standard interaction variance vs PA-interaction)
simulation_version = "mixed_variance_component"
sample_size = "450"
pa_scaling_factors = c("1.5", "2.0", "4.0", "6.0", "8.0", "10.0")

# Inference using no proportional amp
inference_method = "interaction-VCA"
# Load in data
df <- load_in_interaction_variance_estimates_for_fixed_sample_size_inference_method_but_vary_scaling_factors(simulation_results_dir, simulation_version, sample_size, inference_method, pa_scaling_factors)
pp1 <- make_interaction_variance_estimate_plot(df, pa_scaling_factors, simulation_version, sample_size, inference_method)
# Inference using prop amp
inference_method = "PA+interaction-VCA"
# Load in data
df <- load_in_interaction_variance_estimates_for_fixed_sample_size_inference_method_but_vary_scaling_factors(simulation_results_dir, simulation_version, sample_size, inference_method, pa_scaling_factors)
pp2 <- make_interaction_variance_estimate_plot(df, pa_scaling_factors, simulation_version, sample_size, inference_method)
# Plot
legender=get_legend(pp2)
joint_pp <- plot_grid(pp1 + theme(legend.position="none"), pp2+ theme(legend.position="none"), legender, ncol=1, rel_heights=c(1,1,.34))
# save
output_file <- paste0(visualization_results_dir, "interaction_variance_estimates_bar_plot_", simulation_version, "_", sample_size, ".pdf")
ggsave(joint_pp, file=output_file, width=7.2, height=5.5, units="in")

