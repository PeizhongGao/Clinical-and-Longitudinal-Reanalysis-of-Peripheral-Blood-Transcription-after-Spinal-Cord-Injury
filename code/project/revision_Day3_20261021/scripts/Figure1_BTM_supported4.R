# Public-package copy: see metadata/code_changes.csv for provenance.
# Figure contract: one quantitative forest plot highlighting the four longitudinal
# BTM contrasts with BH16 q < 0.05. The remaining eight longitudinal contrasts and four
# clinical contrasts are not discarded: the full family remains in Supplementary
# Table S4 and Figure S9. Input is the frozen 19 October 2026 result table.

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else getwd()
input_file <- file.path(root, "revision_numeric_freeze_20261019", "tables", "BTM_primary16.csv")
output_dir <- file.path(root, "revision_Day3_20261021", "main_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8-BOM")
required <- c("estimate", "lower_CL", "upper_CL", "p_value", "n_observations",
              "n_participants", "spec", "variant", "module_id", "q_BH_primary16",
              "source_branch", "contrast", "scope")
stopifnot(all(required %in% names(results)))

# Apply the revision-stage display threshold within the designated BH16
# family. The full 16-test family remains available in the supplement.
family <- subset(
  results,
  source_branch == "longitudinal" & variant == "published" &
    spec == "disjoint_both" & scope == "all37"
)
shown <- subset(family, !is.na(q_BH_primary16) & q_BH_primary16 < 0.05)
expected <- c("M7.0|3m - acute", "M7.0|6m - acute", "M7.0|12m - acute",
              "M223|6m - acute")
observed <- paste(shown$module_id, shown$contrast, sep = "|")
stopifnot(nrow(family) == 12L, nrow(shown) == 4L, setequal(observed, expected))
stopifnot(all(shown$lower_CL <= shown$estimate), all(shown$estimate <= shown$upper_CL))

month_label <- c("3m - acute" = "3 months vs acute",
                 "6m - acute" = "6 months vs acute",
                 "12m - acute" = "12 months vs acute")
shown$contrast_label <- paste(
  shown$module_id,
  unname(month_label[shown$contrast]),
  sep = "   "
)
order_bottom_to_top <- c("M223   6 months vs acute",
                         "M7.0   12 months vs acute",
                         "M7.0   6 months vs acute",
                         "M7.0   3 months vs acute")
shown$contrast_label <- factor(shown$contrast_label, levels = order_bottom_to_top)
shown$module_id <- factor(shown$module_id, levels = c("M7.0", "M223"))

p <- ggplot(shown, aes(y = contrast_label, colour = module_id)) +
  geom_vline(xintercept = 0, linewidth = 0.45, linetype = "22", colour = "#777777") +
  geom_segment(aes(x = lower_CL, xend = upper_CL,
                   yend = contrast_label), linewidth = 0.8, lineend = "round") +
  geom_point(aes(x = estimate), size = 2.6) +
  scale_colour_manual(values = c("M7.0" = "#2C6E8F", "M223" = "#B6763D")) +
  scale_x_continuous(breaks = seq(0, 0.6, by = 0.2),
                     limits = c(-0.025, 0.75),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Estimate (source-specific log-expression units)",
    y = NULL,
    subtitle = "Four longitudinal contrasts meeting the designated BH16 q < 0.05 threshold"
  ) +
  theme_classic(base_size = 8.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(colour = "#222222"),
    axis.text.y = element_text(size = 8.5, margin = margin(r = 5)),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 7)),
    plot.subtitle = element_text(size = 8, colour = "#444444", margin = margin(b = 8)),
    legend.position = "none",
    plot.margin = margin(8, 9, 6, 6)
  )

base <- file.path(output_dir, "Figure1_BTM_supported4")
width_in <- 6.35
height_in <- 3.05
svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(p)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in,
                     family = "Arial")
print(p)
grDevices::dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_in, height = height_in,
              units = "in", res = 600, background = "white")
print(p)
grDevices::dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, background = "white", compression = "lzw")
print(p)
grDevices::dev.off()

source_data <- shown[, c("module_id", "contrast", "estimate", "lower_CL", "upper_CL",
                         "p_value", "q_BH_primary16", "n_participants",
                         "n_observations", "spec", "variant", "scope")]
write.csv(source_data, file.path(output_dir, "Figure1_BTM_supported4_source_data.csv"),
          row.names = FALSE, na = "")
cat(sprintf("BH16 family contrasts: %d\nDisplayed q<0.05 contrasts: %d\n", nrow(family), nrow(shown)))
print(source_data, row.names = FALSE)
