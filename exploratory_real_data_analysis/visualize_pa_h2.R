#############################################
# Visualize PA-H2 results across gene category bins.
# For each statistic that genes were binned on, plot PA_interaction_fraction across the bins,
# with the un-stratified (all genes) estimate drawn as a dashed reference line.
# Uses the covariate-adjusted run, which run_pa_h2.sh writes at the bare output stem.
#############################################
args = commandArgs(trailingOnly=TRUE)
library(cowplot)
library(ggplot2)
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


pa_h2_summary_file <- function(pa_h2_results_dir, tissue_name, normalization_method, cell_type, stem_suffix) {
	# run_pa_h2.sh writes the covariate-adjusted run at the bare output stem; the permuted
	# null is the same path with "_permuted" before "_PA_h2_summary.txt"
	return(paste0(pa_h2_results_dir, "pa_h2_results_", tissue_name, "_", normalization_method, "_", cell_type, "_", stem_suffix, "_PA_h2_summary.txt"))
}


pa_h2_downsample_summary_file <- function(pa_h2_results_dir, tissue_name, normalization_method, cell_type, downsampling_percentage, gene_filter_tag) {
	# The downsampled stems put the downsampling percentage (and the gene-filter tag) ahead of
	# the normalization method, so they do not share a shape with pa_h2_summary_file
	return(paste0(pa_h2_results_dir, "pa_h2_results_", tissue_name, "_downsample_", downsampling_percentage, gene_filter_tag, "_", normalization_method, "_", cell_type, "_all_genes_PA_h2_summary.txt"))
}


load_variance_class <- function(file_name, variance_class) {
	# Pull a single variance class out of one PA-H2 summary file. Returns NULL if the file is
	# not there yet (eg. that job is still queued) so a partial sweep still plots.
	if (file.exists(file_name) == FALSE) {
		print(paste0("WARNING: missing ", file_name))
		return(NULL)
	}
	df <- read.table(file_name, header=TRUE, sep="\t")
	df <- df[as.character(df$variance_class) == variance_class, ]
	if (nrow(df) != 1) {
		print(paste0("WARNING: expected exactly one ", variance_class, " row in ", file_name))
		return(NULL)
	}
	return(df)
}


load_pa_h2_across_bins <- function(pa_h2_results_dir, tissue_name, normalization_method, cell_type, gene_category_statistic, gene_category_bins, variance_class) {
	all_dfs <- list()
	for (bin_iter in 1:length(gene_category_bins)) {
		gene_category_bin <- gene_category_bins[bin_iter]
		file_name <- pa_h2_summary_file(pa_h2_results_dir, tissue_name, normalization_method, cell_type, paste0(gene_category_statistic, "_", gene_category_bin))
		df <- load_variance_class(file_name, variance_class)
		if (is.null(df)) {
			next
		}
		df$gene_category_bin <- gene_category_bin
		all_dfs[[as.character(gene_category_bin)]] <- df
	}
	if (length(all_dfs) == 0) {
		return(NULL)
	}
	combined_df <- do.call(rbind, all_dfs)
	combined_df$gene_category_bin <- factor(combined_df$gene_category_bin, levels=gene_category_bins)
	return(combined_df)
}


make_pa_h2_bin_plot <- function(df, gene_category_statistic, all_gene_estimate, y_label) {
	pp <- ggplot(df, aes(x=gene_category_bin, y=variance_estimate)) +
		geom_bar(stat="identity", fill="steelblue2") +
		geom_errorbar(aes(ymin=variance_estimate - standard_error, ymax=variance_estimate + standard_error), width=0.25) +
		# Keep every bin on the axis, so a job that has not finished shows up as a gap
		# rather than silently closing ranks
		scale_x_discrete(drop=FALSE) +
		# Fixed y range so panels are comparable across statistics. coord_cartesian zooms
		# rather than filtering, so a bar or error bar running past the range is clipped
		# at the edge instead of being dropped from the plot entirely.
		coord_cartesian(ylim=c(-0.1, 0.4)) +
		figure_theme() +
		labs(x=paste0(gene_category_statistic, " bin  (0 = lowest)"), y=y_label, title=gene_category_statistic)

	# Un-stratified estimate, for reference
	if (is.null(all_gene_estimate) == FALSE) {
		pp <- pp + geom_hline(yintercept=all_gene_estimate, linetype=2)
	}
	return(pp)
}


load_pa_h2_across_downsampling <- function(pa_h2_results_dir, tissue_name, normalization_method, cell_type, downsampling_percentages, variance_class) {
	# Both gene sets are run at every downsampling percentage, so load them together and keep
	# the gene set as a column to group on
	gene_set_tags <- c("all genes"="", "lowly expressed genes removed"="_filtered_genes")
	all_dfs <- list()
	for (gene_set_iter in 1:length(gene_set_tags)) {
		gene_set_label <- names(gene_set_tags)[gene_set_iter]
		gene_filter_tag <- gene_set_tags[[gene_set_iter]]
		for (downsampling_iter in 1:length(downsampling_percentages)) {
			downsampling_percentage <- downsampling_percentages[downsampling_iter]
			file_name <- pa_h2_downsample_summary_file(pa_h2_results_dir, tissue_name, normalization_method, cell_type, downsampling_percentage, gene_filter_tag)
			df <- load_variance_class(file_name, variance_class)
			if (is.null(df)) {
				next
			}
			df$downsampling_percentage <- downsampling_percentage
			df$gene_set <- gene_set_label
			all_dfs[[paste0(gene_set_label, "_", downsampling_percentage)]] <- df
		}
	}
	if (length(all_dfs) == 0) {
		return(NULL)
	}
	combined_df <- do.call(rbind, all_dfs)
	combined_df$downsampling_percentage <- factor(combined_df$downsampling_percentage, levels=downsampling_percentages)
	combined_df$gene_set <- factor(combined_df$gene_set, levels=names(gene_set_tags))
	return(combined_df)
}


load_pa_h2_across_normalization_methods <- function(pa_h2_results_dir, tissue_name, normalization_methods, cell_type, variance_class) {
	# The un-stratified (all genes) run for each normalization method
	all_dfs <- list()
	for (normalization_iter in 1:length(normalization_methods)) {
		normalization_method <- normalization_methods[normalization_iter]
		file_name <- pa_h2_summary_file(pa_h2_results_dir, tissue_name, normalization_method, cell_type, "all_genes")
		df <- load_variance_class(file_name, variance_class)
		if (is.null(df)) {
			next
		}
		df$normalization_method <- normalization_method
		all_dfs[[normalization_method]] <- df
	}
	if (length(all_dfs) == 0) {
		return(NULL)
	}
	combined_df <- do.call(rbind, all_dfs)
	combined_df$normalization_method <- factor(combined_df$normalization_method, levels=normalization_methods)
	return(combined_df)
}


make_pa_h2_normalization_plot <- function(df, y_label) {
	pp <- ggplot(df, aes(x=normalization_method, y=variance_estimate)) +
		geom_bar(stat="identity", fill="steelblue2", width=0.6) +
		geom_errorbar(aes(ymin=variance_estimate - standard_error, ymax=variance_estimate + standard_error), width=0.25) +
		scale_x_discrete(drop=FALSE) +
		figure_theme() +
		labs(x="Normalization method", y=y_label, title="PA-H2 across normalization methods")
	return(pp)
}


make_pa_h2_downsample_plot <- function(df, all_gene_estimate, y_label) {
	dodge <- position_dodge(width=0.9)

	pp <- ggplot(df, aes(x=downsampling_percentage, y=variance_estimate, fill=gene_set)) +
		geom_bar(stat="identity", position=dodge) +
		geom_errorbar(aes(ymin=variance_estimate - standard_error, ymax=variance_estimate + standard_error), width=0.25, position=dodge) +
		scale_x_discrete(drop=FALSE) +
		scale_fill_manual(values=c("grey", "steelblue2")) +
		figure_theme() +
		labs(x="Downsampling percentage of reads", y=y_label, fill="", title="PA-H2 across sequencing depth")

	if (is.null(all_gene_estimate) == FALSE) {
		pp <- pp + geom_hline(yintercept=all_gene_estimate, linetype=2)
	}
	return(pp)
}


#####################
# Command line args
#####################
pa_h2_results_dir        = args[1]
tissue_name              = args[2]
normalization_method     = args[3]
cell_type                = args[4]
gene_category_statistics = strsplit(args[5], ",")[[1]]   # comma separated
num_gene_category_bins   = as.numeric(args[6])
downsampling_percentages = strsplit(args[7], ",")[[1]]   # comma separated
normalization_methods    = strsplit(args[8], ",")[[1]]   # comma separated
visualization_output_stem = args[9]

variance_class <- "PA_interaction_fraction"
y_label <- "PA interaction fraction"
gene_category_bins <- 0:(num_gene_category_bins - 1)

# Un-stratified (all genes) estimate, drawn as a dashed line on every panel
all_gene_df <- load_variance_class(pa_h2_summary_file(pa_h2_results_dir, tissue_name, normalization_method, cell_type, "all_genes"), variance_class)
all_gene_estimate <- NULL
if (is.null(all_gene_df) == FALSE) {
	all_gene_estimate <- all_gene_df$variance_estimate
}

plot_list <- list()
for (statistic_iter in 1:length(gene_category_statistics)) {
	gene_category_statistic <- gene_category_statistics[statistic_iter]

	df <- load_pa_h2_across_bins(pa_h2_results_dir, tissue_name, normalization_method, cell_type, gene_category_statistic, gene_category_bins, variance_class)
	if (is.null(df)) {
		print(paste0("WARNING: no results found for ", gene_category_statistic, "; skipping"))
		next
	}

	pp <- make_pa_h2_bin_plot(df, gene_category_statistic, all_gene_estimate, y_label)

	output_file <- paste0(visualization_output_stem, "pa_h2_", variance_class, "_by_", gene_category_statistic, ".pdf")
	ggsave(pp, file=output_file, width=7.2, height=3.3, units="in")
	print(output_file)

	plot_list[[gene_category_statistic]] <- pp
}

# Joint figure with one panel per binning statistic
if (length(plot_list) > 0) {
	joint_pp <- plot_grid(plotlist=plot_list, ncol=2)
	joint_output_file <- paste0(visualization_output_stem, "pa_h2_", variance_class, "_by_gene_category_statistic.pdf")
	ggsave(joint_pp, file=joint_output_file, width=11.0, height=3.3*ceiling(length(plot_list)/2), units="in")
	print(joint_output_file)
}


#####################
# Separate figure: across sequencing depth rather than gene category bins
#####################
downsample_df <- load_pa_h2_across_downsampling(pa_h2_results_dir, tissue_name, normalization_method, cell_type, downsampling_percentages, variance_class)
if (is.null(downsample_df)) {
	print("WARNING: no downsampled results found; skipping the downsampling figure")
} else {
	downsample_pp <- make_pa_h2_downsample_plot(downsample_df, all_gene_estimate, y_label)
	downsample_output_file <- paste0(visualization_output_stem, "pa_h2_", variance_class, "_by_downsampling_percentage.pdf")
	ggsave(downsample_pp, file=downsample_output_file, width=7.2, height=3.3, units="in")
	print(downsample_output_file)
}


#####################
# Separate figure: across normalization methods (all genes)
#####################
normalization_df <- load_pa_h2_across_normalization_methods(pa_h2_results_dir, tissue_name, normalization_methods, cell_type, variance_class)
if (is.null(normalization_df)) {
	print("WARNING: no all-genes results found; skipping the normalization method figure")
} else {
	normalization_pp <- make_pa_h2_normalization_plot(normalization_df, y_label)
	normalization_output_file <- paste0(visualization_output_stem, "pa_h2_", variance_class, "_by_normalization_method.pdf")
	ggsave(normalization_pp, file=normalization_output_file, width=7.2, height=3.3, units="in")
	print(normalization_output_file)
}
