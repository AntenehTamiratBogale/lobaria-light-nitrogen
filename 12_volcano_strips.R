# Wider, shorter volcano plots with comparison strips above the plotting area.
# Run after script 10, in the same R session. No statistical tests are rerun.
# Coloured strips show direction of abundance, not statistical significance.

# Replace the old labels inside each plot with three strips above it.
volcano_strips <- function(plot, left, right, left_colour, right_colour, arrow = TRUE) {
  old_labels <- vapply(plot$layers, function(layer) {
    labels <- layer$aes_params$label
    if (is.null(labels) && is.data.frame(layer$data)) labels <- layer$data$label
    inherits(layer$geom, "GeomText") && !is.null(labels) &&
      any(grepl("LLS|HLS|N0\\.32|N1\\.28|RAT|DTP|Negative interaction|Positive interaction", labels))
  }, logical(1))
  plot$layers <- plot$layers[!old_labels]
  panel <- ggplot_build(plot)$layout$panel_params[[1]]
  xr <- panel$x.range
  yr <- panel$y.range
  xr <- c(min(xr[1], -1.2), max(xr[2], 1.2))
  height <- diff(yr)
  strip_bottom <- yr[2] + 0.02 * height
  strip_top <- yr[2] + 0.10 * height
  strip_middle <- (strip_bottom + strip_top) / 2

  # Plotmath arrows also work in the PDF export.
  left_label <- if (arrow) paste0('"" %up% bold("', left, '")') else left
  right_label <- if (arrow) paste0('"" %up% bold("', right, '")') else right

  plot +
    annotate("rect", xmin = xr[1], xmax = -1, ymin = strip_bottom,
             ymax = strip_top, fill = left_colour, colour = NA) +
    annotate("rect", xmin = -1, xmax = 1, ymin = strip_bottom,
             ymax = strip_top, fill = "grey85", colour = NA) +
    annotate("rect", xmin = 1, xmax = xr[2], ymin = strip_bottom,
             ymax = strip_top, fill = right_colour, colour = NA) +
    annotate("text", x = (xr[1] - 1) / 2, y = strip_middle,
             label = left_label, parse = arrow, size = 9 / .pt,
             colour = "white", fontface = "bold") +
    annotate("text", x = (1 + xr[2]) / 2, y = strip_middle,
             label = right_label, parse = arrow, size = 9 / .pt,
             colour = "white", fontface = "bold") +
    coord_cartesian(xlim = xr, ylim = yr, expand = FALSE, clip = "off") +
    labs(subtitle = " ") +
    theme(plot.subtitle = element_text(size = 1, margin = margin(b = 30)),
          plot.margin = margin(6, 6, 6, 6), legend.position = "right")
}

# Amplicon panels. LLS, N0.32 and RAT remain on the left.
v16_light <- volcano_strips(p_16L, "LLS", "HLS", col_LLS, col_HLS)
v16_nitrogen <- volcano_strips(p_16N, "N0.32", "N1.28", col_N0.32, col_N1.28)
v16_transplant <- volcano_strips(p_16T, "RAT", "DTP", col_RAT, col_DTP)
v16_interaction <- volcano_strips(p_16LN, "Negative interaction", "Positive interaction",
                                  "grey80", "grey65", arrow = FALSE)

# 18S results were calculated in script 10.
p18_light <- ampl_volcano(res_18S_L, vol_colors_18S, "18S — Light") +
  labs(x = "Log2 fold change (HLS / LLS)", y = "−log10 adjusted P")
p18_nitrogen <- ampl_volcano(res_18S_N, vol_colors_18S, "18S — Nitrogen") +
  labs(x = "Log2 fold change (N1.28 / N0.32)", y = "−log10 adjusted P")
p18_interaction <- ampl_volcano(res_18S_LN, vol_colors_18S, "18S — Light × Nitrogen") +
  labs(x = "Interaction log2 fold change", y = "−log10 adjusted P")
v18_light <- volcano_strips(p18_light, "LLS", "HLS", col_LLS, col_HLS)
v18_nitrogen <- volcano_strips(p18_nitrogen, "N0.32", "N1.28", col_N0.32, col_N1.28)
v18_interaction <- volcano_strips(p18_interaction, "Negative interaction", "Positive interaction",
                                  "grey80", "grey65", arrow = FALSE)

# Global metaproteome panels; no ten-functional-group comparisons are changed.
vmeta_light <- volcano_strips(p_vol_light, "LLS", "HLS", col_LLS, col_HLS)
vmeta_nitrogen <- volcano_strips(p_vol_nitrogen, "N0.32", "N1.28", col_N0.32, col_N1.28)
vmeta_transplant <- volcano_strips(p_vol_transplant, "RAT", "DTP", col_RAT, col_DTP)
vmeta_interaction <- volcano_strips(p_vol_interaction, "Negative interaction", "Positive interaction",
                                   "grey80", "grey65", arrow = FALSE)

# Supplementary Fig. S5: three RAT/DTP panels, with common axis titles.
s5_A <- volcano_strips(vol_T16, "RAT", "DTP", col_RAT, col_DTP)
s5_B <- volcano_strips(vol_T18, "RAT", "DTP", col_RAT, col_DTP)
s5_C <- volcano_strips(vol_Tmeta, "RAT", "DTP", col_RAT, col_DTP)
fig_s5_core <- (s5_A | s5_B | s5_C) + plot_annotation(tag_levels = "A")
fig_s5 <- shared_axes(fig_s5_core, "Log2 fold change (DTP / RAT)", "−log10 adjusted P")
save_supp("FigureS5_RAT_DTP_volcano", fig_s5, 14.8, 4.8)

# Save JPG, PDF and SVG versions of the main figures.
save_volcano <- function(path, plot, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (extension in c("jpg", "pdf", "svg")) {
    ggsave(paste0(path, ".", extension), plot,
           width = width, height = height, dpi = 300, bg = "white")
  }
}
# Keep the class legend shorter than the plotting panel.
v16_light <- v16_light +
  labs(x = expression(Log[2]~"fold change (HLS / LLS)")) +
  theme(plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10), axis.text = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.height = unit(8, "pt"), legend.key.width = unit(7, "pt"),
    legend.spacing.y = unit(0, "pt"), legend.margin = margin(2, 2, 2, 2),
    legend.box.spacing = unit(3, "pt")) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.5)))
save_volcano(file.path(OUTDIR, "Figure6_16S_light_volcano"), v16_light, 4.2, 3.5)

# Figure 7 has single-line axis titles and one grey legend below both panels.
meta_A <- vmeta_light +
  labs(x = expression(Log[2]~"fold change (HLS / LLS)"))
meta_B <- vmeta_nitrogen +
  labs(x = expression(Log[2]~"fold change (N1.28 / N0.32)")) +
  theme(axis.title.y = element_blank())
# Use the same partner scale in both panels to collect one legend.
partner_scale <- scale_colour_manual(values = partner_colors,
  limits = partner_levels, drop = FALSE, name = "Partner")
meta_A <- meta_A + partner_scale + guides(colour = guide_legend(nrow = 1))
meta_B <- meta_B + partner_scale + guides(colour = guide_legend(nrow = 1))
fig7 <- ((meta_A | meta_B) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")) &
  theme(legend.position = "bottom", legend.direction = "horizontal",
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.tag = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 9), axis.text = element_text(size = 8),
    legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
    legend.key.width = unit(6, "pt"), legend.key.height = unit(7, "pt"),
    legend.spacing.x = unit(0, "pt"), legend.margin = margin(2, 2, 2, 2),
    legend.box.background = element_rect(fill = "grey90", colour = "grey65"),
    legend.background = element_rect(fill = "grey90", colour = NA),
    legend.key = element_rect(fill = "grey90", colour = NA),
    legend.box.margin = margin(2, 3, 2, 3),
    legend.box.spacing = unit(8, "pt"))
save_volcano(file.path(OUTDIR, "Figure7_metaproteome_volcano"), fig7, 4.8, 3.5)

# Individual volcano panels for each dataset and comparison.
# These are separate exports, not additional numbered supplementary figures.
volcano_panels <- list(
  "16S_Light" = v16_light, "16S_Nitrogen" = v16_nitrogen,
  "16S_RAT_DTP" = v16_transplant, "16S_Interaction" = v16_interaction,
  "18S_Light" = v18_light, "18S_Nitrogen" = v18_nitrogen,
  "18S_RAT_DTP" = s5_B + labs(x = "Log2 fold change (DTP / RAT)", y = "−log10 adjusted P"),
  "18S_Interaction" = v18_interaction,
  "Metaproteome_Light" = vmeta_light, "Metaproteome_Nitrogen" = vmeta_nitrogen,
  "Metaproteome_RAT_DTP" = vmeta_transplant, "Metaproteome_Interaction" = vmeta_interaction)
for (name in names(volcano_panels)) {
  save_volcano(file.path(OUTDIR, "volcano_panels", name), volcano_panels[[name]], 5.5, 3.5)
}
# Keep the existing Figure 7 individual-panel filenames up to date.
save_volcano(file.path(INDIVIDUAL_DIR, "Figure7", "Figure7A_metaproteome_Light"), vmeta_light, 4.2, 3.1)
save_volcano(file.path(INDIVIDUAL_DIR, "Figure7", "Figure7B_metaproteome_Nitrogen"), vmeta_nitrogen, 4.2, 3.1)
cat("Updated volcano plots with upward-arrow strips. Statistical results are unchanged.\n")
