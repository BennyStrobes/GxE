#############################################
# Visualize per-gene features (from the 1K1K single-cell gene feature file)
# against the per-gene variance stratified by a binary E-variable.
#
# Produces four figures:
#   (1) violin   : GTEx (full depth) vs 1K1K
#   (2) mean+CI  : GTEx (full depth) vs 1K1K
#   (3) violin   : the asymmetric downsampling bins
#   (4) mean+CI  : the asymmetric downsampling bins
#############################################
library(ggplot2)

extract_gtex_abs_log_var_ratio <- function(per_gene_variance_file) {
    aa <- read.table(per_gene_variance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    ratios <- aa$var_E0 / aa$var_E1
    log_ratios <- log(ratios)
    return(abs(log_ratios))

}

extract_sc_1k1k_abs_log_var_ratio <- function(sc_1k1k_gene_features_file) {
    aa <- read.table(sc_1k1k_gene_features_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # The log_var_* columns are log2, so the variance ratio is 2^(log_var_B - log_var_CD4).
    # Multiplying the log2 difference by ln(2) puts this on the natural-log scale used by
    # extract_gtex_abs_log_var_ratio, so the two datasets are directly comparable. Taking the
    # difference directly also avoids exponentiating, which overflows to Inf (and then NaN)
    # for large-magnitude log variances.
    log_ratios <- (aa$log_var_B - aa$log_var_CD4) * log(2)
    return(abs(log_ratios))
}

# Drop non-finite ratios (a zero stratum variance makes the log ratio infinite), order the
# ticks by sort_key, and fold the per-group gene count into the tick label.
finalize_ratio_df <- function(df) {
    n_dropped <- sum(!is.finite(df$abs_log_var_ratio))
    if (n_dropped > 0) {
        print(paste0("WARNING: dropping ", n_dropped, " genes with a non-finite abs log variance ratio"))
        df <- df[is.finite(df$abs_log_var_ratio), ]
    }

    # Gene counts are per (label, gene_set): the lowly-expressed-gene filter is re-applied
    # at each downsampling depth, so a bin's n differs between the two gene sets.
    group_key <- if ("gene_set" %in% names(df)) {
        paste(df$gene_set, df$label, sep = "\r")
    } else {
        df$label
    }
    counts <- table(group_key)
    df$tick <- paste0(df$label, "\n(n = ", counts[group_key], ")")
    df$tick <- factor(df$tick, levels = unique(df$tick[order(df$sort_key)]))
    df$series <- factor(df$series, levels = unique(df$series[order(df$sort_key)]))
    return(df)
}

# GTEx at full depth vs 1K1K
build_comparison_df <- function(gtex_ratios, sc_1k1k_ratios) {
    finalize_ratio_df(rbind(
        data.frame(label = "GTEx", series = "GTEx", sort_key = 1,
                   abs_log_var_ratio = gtex_ratios),
        data.frame(label = "1K1K", series = "1K1K", sort_key = 2,
                   abs_log_var_ratio = sc_1k1k_ratios)
    ))
}

# Discover the asymmetric-downsampling per-gene variance files under a file stem and read
# the abs log variance ratio out of each. Files are named
#   <stem>.downsample_<pct>[.filtered_genes].<normalization>.<cell_type>.per_gene_variance.txt
# so the downsampling percentage is recovered from the file name itself rather than being
# re-specified here. Only the ".filtered_genes" variant is plotted -- the all-genes matrices
# are left out because the lowly-expressed-gene filter is what makes the bins comparable.
build_downsampling_df <- function(asymmetric_per_gene_variance_stem) {
    files <- Sys.glob(paste0(asymmetric_per_gene_variance_stem,
                             ".downsample_*.filtered_genes.*per_gene_variance.txt"))
    if (length(files) == 0) {
        print(paste0("WARNING: no asymmetric filtered-gene per-gene variance files found under stem ",
                     asymmetric_per_gene_variance_stem))
        return(NULL)
    }

    df <- do.call(rbind, lapply(files, function(f) {
        pct <- sub(".*\\.downsample_([0-9]+\\.[0-9]+).*", "\\1", basename(f))
        data.frame(label = pct,
                   series = "GTEx asymmetric downsample",
                   sort_key = as.numeric(pct),
                   abs_log_var_ratio = extract_gtex_abs_log_var_ratio(f))
    }))
    return(finalize_ratio_df(df))
}

# Shared styling so all four figures read as one set. Colour follows the series, never the
# downsampling percentage -- the percentage is ordered, so the x-axis carries it instead.
SERIES_COLOURS <- c("GTEx" = "#2a78d6",
                    "1K1K" = "#eb6834",
                    "GTEx asymmetric downsample" = "#2a78d6")

ratio_plot_theme <- function(aesthetic = "fill") {
    list(
        if (aesthetic == "fill") scale_fill_manual(values = SERIES_COLOURS, guide = "none")
        else scale_colour_manual(values = SERIES_COLOURS, guide = "none"),
        theme_classic(base_size = 11),
        theme(legend.position = "none",
              strip.background = element_blank(),
              strip.text = element_text(size = 9, colour = "#0b0b0b"),
              plot.title = element_text(size = 11, colour = "#0b0b0b"),
              plot.subtitle = element_text(size = 9, colour = "#52514e"),
              axis.text = element_text(colour = "#0b0b0b"),
              axis.text.x = element_text(size = 8),
              axis.title = element_text(colour = "#52514e"),
              axis.line = element_line(colour = "#b8b7b1", linewidth = 0.3),
              axis.ticks = element_line(colour = "#b8b7b1", linewidth = 0.3))
    )
}

# Facet by gene set when the frame carries one (the downsampling figures). Scales are free
# on x so each panel shows only its own ticks, which carry that panel's gene counts.
maybe_facet <- function(df) {
    if ("gene_set" %in% names(df)) facet_wrap(~gene_set, nrow = 1, scales = "free_x") else NULL
}

# Violin of the abs log variance-ratio distributions. Groups differ in gene count, so the
# violins are width-scaled and the tick labels carry n.
make_abs_log_var_ratio_violin <- function(df, x_label, title, subtitle) {
    ggplot(df, aes(x = tick, y = abs_log_var_ratio, fill = series)) +
        geom_violin(scale = "width", width = 0.8, colour = "#52514e", linewidth = 0.3,
                    alpha = 0.75, trim = TRUE) +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "#fcfcfb",
                     colour = "#52514e", linewidth = 0.3) +
        maybe_facet(df) +
        labs(x = x_label, y = "|log variance ratio|", title = title, subtitle = subtitle) +
        ratio_plot_theme(aesthetic = "fill")
}

# Mean abs log variance ratio with a 95% CI built from the standard error of the mean
# (mean +/- 1.96 * sd / sqrt(n)). Point-and-interval rather than bars: the interval is what
# carries the information, and bars anchored at zero would compress it to nothing.
make_abs_log_var_ratio_mean_ci <- function(df, x_label, title, subtitle) {
    group_cols <- intersect(c("tick", "series", "gene_set"), names(df))
    groups <- unique(df[, group_cols, drop = FALSE])
    summary_df <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
        g <- groups[i, , drop = FALSE]
        keep <- rep(TRUE, nrow(df))
        for (col in group_cols) keep <- keep & df[[col]] == g[[col]]
        vals <- df$abs_log_var_ratio[keep]
        se <- sd(vals) / sqrt(length(vals))
        cbind(g, data.frame(mean = mean(vals),
                            lower = mean(vals) - 1.96 * se,
                            upper = mean(vals) + 1.96 * se))
    }))

    ggplot(summary_df, aes(x = tick, y = mean, colour = series)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08, linewidth = 0.6) +
        geom_point(size = 3) +
        geom_text(aes(label = sprintf("%.3f", mean)), hjust = -0.35, size = 2.6,
                  colour = "#0b0b0b") +
        maybe_facet(df) +
        labs(x = x_label, y = "mean |log variance ratio|", title = title, subtitle = subtitle) +
        ratio_plot_theme(aesthetic = "colour")
}

# ---- Command line args ----
args <- commandArgs(trailingOnly = TRUE)
sc_1k1k_gene_features_file <- args[1]
per_gene_variance_file     <- args[2]
visualize_gene_features_dir <- args[3]
# Optional 4th arg: file stem for the asymmetrically-downsampled per-gene variance files.
# When supplied, the two downsampling figures are produced as well.
asymmetric_per_gene_variance_stem <- if (length(args) >= 4) args[4] else NA

# Extract the absolute log variance ratios from the per-gene variance file
gtex_abs_log_var_ratios <- extract_gtex_abs_log_var_ratio(per_gene_variance_file)
sc_1k1k_log_var_ratios <- extract_sc_1k1k_abs_log_var_ratio(sc_1k1k_gene_features_file)

comparison_df <- build_comparison_df(gtex_abs_log_var_ratios, sc_1k1k_log_var_ratios)
comparison_subtitle <- "GTEx: high vs low neutrophil proportion  |  1K1K: B cell vs CD4 T cell"

# ---- (1) GTEx vs 1K1K, violin ----
output_file <- file.path(visualize_gene_features_dir,
                         "gtex_vs_sc_1k1k_abs_log_var_ratio_violin.pdf")
pp <- make_abs_log_var_ratio_violin(comparison_df, NULL,
                                    "Per-gene variance heterogeneity between E strata",
                                    comparison_subtitle)
ggsave(pp, file = output_file, width = 5.0, height = 4.0, units = "in")

# ---- (2) GTEx vs 1K1K, mean with 95% CI ----
output_file <- file.path(visualize_gene_features_dir,
                         "gtex_vs_sc_1k1k_abs_log_var_ratio_mean_ci.pdf")
pp <- make_abs_log_var_ratio_mean_ci(comparison_df, NULL,
                                     "Mean per-gene variance heterogeneity between E strata",
                                     "Error bars are 95% CI on the mean (1.96 x SEM)")
ggsave(pp, file = output_file, width = 5.0, height = 4.0, units = "in")

# ---- (3) and (4): the asymmetric downsampling bins ----
if (!is.na(asymmetric_per_gene_variance_stem)) {
    downsampling_df <- build_downsampling_df(asymmetric_per_gene_variance_stem)
    if (!is.null(downsampling_df)) {
        downsampling_subtitle <- paste0("GTEx, reads thinned in the high-neutrophil stratum only  |  ",
                                        "lowly-expressed genes removed")

        output_file <- file.path(visualize_gene_features_dir,
                                 "asymmetric_downsampling_abs_log_var_ratio_violin.pdf")
        pp <- make_abs_log_var_ratio_violin(downsampling_df, "downsampling percentage",
                                            "Per-gene variance heterogeneity by downsampling depth",
                                            downsampling_subtitle)
        ggsave(pp, file = output_file, width = 6.0, height = 4.0, units = "in")

        output_file <- file.path(visualize_gene_features_dir,
                                 "asymmetric_downsampling_abs_log_var_ratio_mean_ci.pdf")
        pp <- make_abs_log_var_ratio_mean_ci(downsampling_df, "downsampling percentage",
                                             "Mean per-gene variance heterogeneity by downsampling depth",
                                             paste0("Error bars are 95% CI on the mean (1.96 x SEM)  |  ",
                                                    "lowly-expressed genes removed"))
        ggsave(pp, file = output_file, width = 6.0, height = 4.0, units = "in")
    }
}
