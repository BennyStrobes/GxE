#############################################
# Visualize per-gene expression variance stratified by a binary E-variable,
# comparing across normalization methods.
# Produces:
#   (1) var_E1 vs var_E0 scatter (log-log), faceted by normalization method
#   (2) log2 variance ratio  log2(var_E1 / var_E0), density overlaid by method
#   (3) log2 ratio vs mean stratum variance (MA-style artifact check), faceted by method
#   (4) method-agreement scatter: log2 ratio of each method vs the first method
# Requires: ggplot2, gridExtra
#############################################

library(ggplot2)
library(gridExtra)

# ---- Command line args ----
args <- commandArgs(trailingOnly = TRUE)
processed_expression_dir <- args[1]
tissue_name              <- args[2]
cell_type                <- args[3]
normalization_methods    <- strsplit(args[4], ",")[[1]]  # comma-separated
output_stem              <- args[5]

# ---- Load and combine per-method variance files ----
load_method <- function(method) {
    f <- file.path(processed_expression_dir,
                   paste0(tissue_name, ".", method, ".", cell_type,
                          ".per_gene_variance.txt"))
    d <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Drop genes with a non-finite or non-positive variance in either stratum
    d <- d[is.finite(d$var_E0) & is.finite(d$var_E1) &
           d$var_E0 > 0 & d$var_E1 > 0, ]
    d$log2_ratio    <- log2(d$var_E1 / d$var_E0)
    d$mean_variance <- (d$var_E0 + d$var_E1) / 2
    d$method        <- method
    d
}

if (length(normalization_methods) == 0 || all(normalization_methods == "")) {
    stop("No normalization methods provided (arg 4 was empty).")
}

df <- do.call(rbind, lapply(normalization_methods, load_method))
# Keep method ordering stable in facets/legends
df$method <- factor(df$method, levels = normalization_methods)

subtitle_str <- paste0(tissue_name, "  |  E = ", cell_type,
                       "  |  methods: ", paste(normalization_methods, collapse = ", "))

# ---- Plot 1: var_E1 vs var_E0 scatter, faceted by method ----
p1 <- ggplot(df, aes(x = var_E0, y = var_E1)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~ method) +
    labs(x = "var (E = 0, low)", y = "var (E = 1, high)",
         title = "Per-gene variance: E=1 vs E=0") +
    theme_bw()

# ---- Plot 2: density of log2 variance ratio, overlaid by method ----
medians <- aggregate(log2_ratio ~ method, data = df, FUN = median)
p2 <- ggplot(df, aes(x = log2_ratio, color = method, fill = method)) +
    geom_density(alpha = 0.15) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_vline(data = medians, aes(xintercept = log2_ratio, color = method),
               linetype = "dotted") +
    labs(x = "log2(var_E1 / var_E0)", y = "density",
         title = "Variance ratio by method (dotted = median)") +
    theme_bw()

# ---- Plot 3: log2 ratio vs mean stratum variance, faceted by method ----
p3 <- ggplot(df, aes(x = mean_variance, y = log2_ratio)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    scale_x_log10() +
    facet_wrap(~ method) +
    labs(x = "mean stratum variance  (var_E0 + var_E1)/2",
         y = "log2(var_E1 / var_E0)",
         title = "Variance ratio vs expression level") +
    theme_bw()

plot_list <- list(p1, p2, p3)

# ---- Plot 4: method-agreement scatter (only if >= 2 methods) ----
if (length(normalization_methods) >= 2) {
    ref_method <- normalization_methods[1]
    ref <- df[df$method == ref_method, c("gene_name", "log2_ratio")]
    colnames(ref) <- c("gene_name", "log2_ratio_ref")
    other <- df[df$method != ref_method, ]
    merged <- merge(other, ref, by = "gene_name")

    p4 <- ggplot(merged, aes(x = log2_ratio_ref, y = log2_ratio)) +
        geom_point(alpha = 0.3, size = 0.6) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        facet_wrap(~ method) +
        labs(x = paste0("log2 ratio  (", ref_method, ")"),
             y = "log2 ratio  (other method)",
             title = paste0("Method agreement vs ", ref_method)) +
        theme_bw()
    plot_list <- c(plot_list, list(p4))
}

# ---- Save combined figure ----
output_file <- paste0(output_stem, ".png")
g <- arrangeGrob(grobs = plot_list, ncol = 2, top = subtitle_str)
ggsave(output_file, g, width = 14, height = 5 * ceiling(length(plot_list) / 2),
       dpi = 200)
cat("Wrote:", output_file, "\n")

# ---- Separate figure: mean-variance relationship, one panel per (stratum x method) ----
# Reshape to long format: one row per (gene, method, stratum)
long_E0 <- data.frame(method = df$method, stratum = "E0 (low)",
                      mean = df$mean_E0, var = df$var_E0)
long_E1 <- data.frame(method = df$method, stratum = "E1 (high)",
                      mean = df$mean_E1, var = df$var_E1)
mv_df <- rbind(long_E0, long_E1)

p_mv <- ggplot(mv_df, aes(x = mean, y = var)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    scale_y_log10() +
    facet_grid(stratum ~ method) +   # rows = stratum, cols = method -> 2 x n_methods grid
    labs(x = "per-gene mean (within stratum)",
         y = "per-gene variance (within stratum)",
         title = "Mean-variance relationship", subtitle = subtitle_str) +
    theme_bw()

mv_output_file <- paste0(output_stem, ".mean_variance.png")
ggsave(mv_output_file, p_mv, width = 5 * length(normalization_methods), height = 9, dpi = 200)
cat("Wrote:", mv_output_file, "\n")

# ---- Separate figure: delta(var) vs delta(mean), one panel per normalization method ----
df$delta_var  <- df$var_E1 - df$var_E0    # E1 - E0
df$delta_mean <- df$mean_E1 - df$mean_E0  # E1 - E0

p_delta <- ggplot(df, aes(x = delta_mean, y = delta_var)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    facet_wrap(~ method) +
    labs(x = "delta mean  (mean_E1 - mean_E0)",
         y = "delta var  (var_E1 - var_E0)",
         title = "Change in variance vs change in mean across strata",
         subtitle = subtitle_str) +
    theme_bw()

delta_output_file <- paste0(output_stem, ".delta_var_vs_delta_mean.png")
ggsave(delta_output_file, p_delta, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", delta_output_file, "\n")

# ---- Separate figure: mean_E0 vs mean_E1 colored by delta(var), one panel per method ----
p_mm <- ggplot(df, aes(x = mean_E0, y = mean_E1, color = delta_var)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_color_gradient2(low = "blue", mid = "grey85", high = "red", midpoint = 0,
                          name = "delta var\n(E1 - E0)") +
    facet_wrap(~ method) +
    labs(x = "mean (E0, low)", y = "mean (E1, high)",
         title = "Stratum means colored by change in variance",
         subtitle = subtitle_str) +
    theme_bw()

mm_output_file <- paste0(output_stem, ".mean_E0_vs_mean_E1_by_delta_var.png")
ggsave(mm_output_file, p_mm, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", mm_output_file, "\n")

# ---- Separate figure: delta(var) vs mean_E0, one panel per method ----
p_dm0 <- ggplot(df, aes(x = mean_E0, y = delta_var)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    facet_wrap(~ method) +
    labs(x = "mean (E0, low)", y = "delta var  (var_E1 - var_E0)",
         title = "Change in variance vs baseline (E0) mean",
         subtitle = subtitle_str) +
    theme_bw()

dm0_output_file <- paste0(output_stem, ".delta_var_vs_mean_E0.png")
ggsave(dm0_output_file, p_dm0, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", dm0_output_file, "\n")

# ---- Separate figure: log2 variance ratio vs mean_E0, one panel per method ----
p_rm0 <- ggplot(df, aes(x = mean_E0, y = log2_ratio)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    facet_wrap(~ method) +
    labs(x = "mean (E0, low)", y = "log2(var_E1 / var_E0)",
         title = "Variance ratio vs baseline (E0) mean",
         subtitle = subtitle_str) +
    theme_bw()

rm0_output_file <- paste0(output_stem, ".log2_var_ratio_vs_mean_E0.png")
ggsave(rm0_output_file, p_rm0, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", rm0_output_file, "\n")

# ---- Separate figure: delta mean vs average mean, colored by log2 variance ratio, one panel per method ----
df$avg_mean <- (df$mean_E0 + df$mean_E1) / 2

p_dm_avg <- ggplot(df, aes(x = delta_mean, y = avg_mean, color = log2_ratio)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    scale_color_gradient2(low = "blue", mid = "grey85", high = "red", midpoint = 0,
                          name = "log2(var_E1\n/ var_E0)") +
    facet_wrap(~ method) +
    labs(x = "delta mean  (mean_E1 - mean_E0)",
         y = "average mean  (mean_E0 + mean_E1)/2",
         title = "Delta mean vs average mean, colored by variance ratio",
         subtitle = subtitle_str) +
    theme_bw()

dm_avg_output_file <- paste0(output_stem, ".delta_mean_vs_avg_mean_by_log2_var_ratio.png")
ggsave(dm_avg_output_file, p_dm_avg, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", dm_avg_output_file, "\n")

# ---- Separate figure: log2 variance ratio vs average mean, colored by delta mean, one panel per method ----
p_ratio_avg <- ggplot(df, aes(x = avg_mean, y = log2_ratio, color = delta_mean)) +
    geom_point(alpha = 0.6, size = 0.7) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_color_gradient2(low = "blue", mid = "grey85", high = "red", midpoint = 0,
                          name = "delta mean\n(E1 - E0)") +
    facet_wrap(~ method) +
    labs(x = "average mean  (mean_E0 + mean_E1)/2",
         y = "log2(var_E1 / var_E0)",
         title = "Variance ratio vs average mean, colored by delta mean",
         subtitle = subtitle_str) +
    theme_bw()

ratio_avg_output_file <- paste0(output_stem, ".log2_var_ratio_vs_avg_mean_by_delta_mean.png")
ggsave(ratio_avg_output_file, p_ratio_avg, width = 6 * length(normalization_methods), height = 5.5, dpi = 200)
cat("Wrote:", ratio_avg_output_file, "\n")
