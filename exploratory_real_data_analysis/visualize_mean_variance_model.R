#############################################
# Visualize the fitted mean-variance relationship, separately within each E stratum.
# Input is the output of fit_mean_variance_model.py (per-gene variance file with the
# fitted variances appended).
# Produces:
#   (1) observed log2 variance vs mean, per stratum, with the fitted curve overlaid
#   (2) the two fitted curves alone, so the between-stratum difference is readable
#   (3) predicted vs observed log2 variance ratio (pred_abs_var_diff vs abs_var_diff)
# Requires: ggplot2, gridExtra
#############################################

library(ggplot2)
library(gridExtra)

# ---- Command line args ----
args <- commandArgs(trailingOnly = TRUE)
mean_variance_model_file <- args[1]
tissue_name              <- args[2]
cell_type                <- args[3]
output_stem              <- args[4]

d <- read.table(mean_variance_model_file, header = TRUE, sep = "\t",
                stringsAsFactors = FALSE, na.strings = c("nan", "NA"))

subtitle_str <- paste0(tissue_name, "  |  E = ", cell_type)

# Two strata, in fixed order. Color carries identity, linetype repeats it so the
# fitted curves stay distinguishable without color.
stratum_levels <- c("E0 (low)", "E1 (high)")
stratum_colors <- c("E0 (low)" = "#377EB8", "E1 (high)" = "#E41A1C")
stratum_lines  <- c("E0 (low)" = "solid",   "E1 (high)" = "22")

# ---- Reshape to one row per (gene, stratum) ----
to_long <- function(stratum_label, mean_col, var_col, pred_col) {
    data.frame(gene_name = d$gene_name,
               stratum   = stratum_label,
               mean      = d[[mean_col]],
               log2_var  = log2(d[[var_col]]),
               pred_log2_var = d[[pred_col]],
               stringsAsFactors = FALSE)
}
mv <- rbind(to_long("E0 (low)",  "mean_E0", "var_E0", "pred_log2_var_E0"),
            to_long("E1 (high)", "mean_E1", "var_E1", "pred_log2_var_E1"))
mv <- mv[is.finite(mv$mean) & is.finite(mv$log2_var) & is.finite(mv$pred_log2_var), ]
mv$stratum <- factor(mv$stratum, levels = stratum_levels)

# The fitted curve is a function of the mean, so ordering by mean draws it as a line
curve_df <- mv[order(mv$stratum, mv$mean), ]

# ---- Plot 1: observed points + fitted curve, one panel per stratum ----
p1 <- ggplot(mv, aes(x = mean, y = log2_var)) +
    geom_point(alpha = 0.20, size = 0.5, color = "grey40") +
    geom_line(data = curve_df,
              aes(y = pred_log2_var, color = stratum, linetype = stratum),
              linewidth = 1) +
    scale_color_manual(values = stratum_colors, name = "E stratum") +
    scale_linetype_manual(values = stratum_lines, name = "E stratum") +
    facet_wrap(~ stratum) +
    labs(x = "per-gene mean (within stratum)",
         y = "log2 per-gene variance (within stratum)",
         title = "Mean-variance relationship with fitted polynomial",
         subtitle = subtitle_str) +
    theme_bw()

# ---- Plot 2: the two fitted curves overlaid ----
p2 <- ggplot(curve_df, aes(x = mean, y = pred_log2_var,
                           color = stratum, linetype = stratum)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = stratum_colors, name = "E stratum") +
    scale_linetype_manual(values = stratum_lines, name = "E stratum") +
    labs(x = "per-gene mean (within stratum)",
         y = "fitted log2 variance",
         title = "Fitted mean-variance curves, both strata",
         subtitle = subtitle_str) +
    theme_bw()

# ---- Plot 3: predicted vs observed variance ratio ----
ratio_df <- data.frame(
    abs_var_diff      = abs(log2(d$var_E0 / d$var_E1)),
    pred_abs_var_diff = abs(d$pred_log2_var_E0 - d$pred_log2_var_E1))
ratio_df <- ratio_df[is.finite(ratio_df$abs_var_diff) &
                     is.finite(ratio_df$pred_abs_var_diff), ]

p3 <- ggplot(ratio_df, aes(x = abs_var_diff, y = pred_abs_var_diff)) +
    geom_point(alpha = 0.25, size = 0.5, color = "grey40") +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(x = "abs_var_diff   |log2(var_E0 / var_E1)|  (observed)",
         y = "pred_abs_var_diff  (fitted)",
         title = "Fitted vs observed variance contrast",
         subtitle = subtitle_str) +
    theme_bw()

# ---- Save ----
output_file <- paste0(output_stem, ".mean_variance_model.png")
g <- arrangeGrob(grobs = list(p1, p2, p3), ncol = 2, top = subtitle_str)
ggsave(output_file, g, width = 14, height = 10, dpi = 200)
cat("Wrote:", output_file, "\n")
