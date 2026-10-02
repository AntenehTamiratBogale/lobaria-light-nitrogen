# Refresh the four supplementary figures without rerunning statistical analyses.
# Run in the same R session after 10_supplementary.R has finished.
# Taxonomy and volcano legends are on the right.
# NMDS width is increased slightly. Alpha-diversity P values remain inside.

# Figure style: readable black text and bold panel titles.
supp_theme <- theme_bw(base_size = 18) + theme(
  panel.grid = element_blank(), axis.text = element_text(colour = "black"),
  axis.title = element_text(size = 18),
  plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
  plot.tag = element_text(face = "bold", size = 18),
  legend.title = element_text(face = "bold", size = 13),
  legend.text = element_text(size = 14), plot.margin = margin(6, 6, 6, 6))

# Widen the specimen bars and keep compact legends on the right.
taxonomy_panel <- function(plot, title) {
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, "GeomBar") ||
        inherits(plot$layers[[i]]$geom, "GeomCol")) {
      plot$layers[[i]]$geom_params$width <- 0.75
    }
  }
  plot + labs(title = title, x = NULL, y = NULL) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE, reverse = TRUE, title.position = "top")) +
    theme(legend.position = "right", legend.box = "vertical",
          legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.35, "cm"),
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
          axis.text = element_text(colour = "black", size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(face = "bold", size = 13),
          plot.margin = margin(6, 6, 6, 6))
}
tax_T16 <- taxonomy_panel(p_tax_16S_T, "16S amplicons")
tax_T18 <- taxonomy_panel(p_tax_18S_T, "18S amplicons")
tax_T18_noA <- taxonomy_panel(p_tax_18S_noA_T, "18S without Ascomycota")

# Shared axis titles; each NMDS keeps the scale of its own ordination.
shared_axes <- function(plot, x, y) {
  plot <- add_common_y_title(plot, y, text_size = 18, width = 0.025)
  x_title <- grid::textGrob(x, gp = grid::gpar(fontsize = 18, col = "black"))
  patchwork::wrap_plots(plot, patchwork::wrap_elements(full = x_title),
                       ncol = 1, heights = c(1, 0.075))
}
# 18S retains the all-specimen NMDS coordinates from script 03.
# Its PERMANOVA uses only the six RAT/DTP specimens.
nmds_T16 <- p_nmds_T16 + labs(title = "16S amplicons", x = NULL, y = NULL) + supp_theme + theme(legend.position = "none")
nmds_T18 <- p_nmds_T18 + labs(title = "18S amplicons", x = NULL, y = NULL) + supp_theme + theme(legend.position = "none")
nmds_Tmeta <- p_nmds_Tp + labs(title = "Metaproteome", x = NULL, y = NULL) + supp_theme + theme(legend.position = "none")

# RAT/DTP box colours match those used in the other plots.
alpha_plot <- function(values, index, p, dataset) {
  values <- values[!is.na(values$Transplant) & values$Transplant %in% c("RAT", "DTP"), ]
  values$Transplant <- factor(values$Transplant, levels = c("RAT", "DTP"))
  p_text <- if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
  ggplot(values, aes(x = Transplant, y = .data[[index]], fill = Transplant)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, colour = "black") +
    geom_point(size = 2.5, colour = "black") +
    scale_fill_manual(values = c("RAT" = col_RAT, "DTP" = col_DTP), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    annotate("text", x = 1.5, y = Inf, label = p_text,
             vjust = 1.3, size = 18 / .pt, colour = "black") +
    labs(x = NULL, y = NULL, title = dataset) + supp_theme
}
# Panel order matches main Fig. 5: richness A/B, alpha diversity C/D.
alpha_A <- alpha_plot(alpha_16S, "Observed", wt_T_obs_16S, "16S amplicons")
alpha_B <- alpha_plot(alpha_18S, "Observed", wt_T_obs_18S, "18S amplicons")
alpha_C <- alpha_plot(alpha_16S, "Shannon", wt_T_sha_16S, "16S amplicons")
alpha_D <- alpha_plot(alpha_18S, "Shannon", wt_T_sha_18S, "18S amplicons")
# Show the dataset titles once, and RAT/DTP labels only on the bottom row.
alpha_A <- alpha_A + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
alpha_B <- alpha_B + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
alpha_C <- alpha_C + labs(title = NULL)
alpha_D <- alpha_D + labs(title = NULL)


# Colour significant amplicon ASVs by class; use grey for other ASVs.
# The 18S class colours match the supplementary 18S taxonomy panel.
top_classes_18S <- ps_18S_T %>% tax_top(n = 15, rank = "Class") %>% sort()
vol_colors_18S <- setNames(colors_18S[seq_along(top_classes_18S)], top_classes_18S)
ampl_volcano <- function(result, palette, title) {
  significant <- result[result$significant & is.finite(result$log2FoldChange) &
                          is.finite(result$neg_logp), ]
  nonsignificant <- result[!result$significant & is.finite(result$log2FoldChange) &
                             is.finite(result$neg_logp), ]
  missing_classes <- setdiff(unique(as.character(significant$Class)), names(palette))
  missing_classes <- missing_classes[!is.na(missing_classes)]
  palette <- c(palette, setNames(rep("grey35", length(missing_classes)), missing_classes))
  ggplot() +
    geom_point(data = nonsignificant, aes(log2FoldChange, neg_logp),
               colour = "grey65", size = 1.4, alpha = 0.65) +
    geom_point(data = significant, aes(log2FoldChange, neg_logp, colour = Class),
               size = 1.4, alpha = 0.8) +
    scale_colour_manual(values = palette, na.value = "grey35", name = "Class") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey45") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    labs(title = title, x = NULL, y = NULL) + supp_theme
}
vol_T16 <- ampl_volcano(res_16S_T, vol_colors_16S, "16S amplicons")
vol_T18 <- ampl_volcano(res_18S_T, vol_colors_18S, "18S amplicons")
vol_Tmeta <- p_vol_transplant + labs(title = "Metaproteome", x = NULL, y = NULL) + supp_theme

# Keep compact legends on the right and use wider figure exports.
wide_volcano <- function(plot) {
  plot + guides(colour = guide_legend(ncol = 1, byrow = TRUE, title.position = "top")) +
    theme(legend.position = "right", legend.box = "vertical",
          legend.text = element_text(size = 12),
          legend.key.height = unit(0.35, "cm"), legend.key.width = unit(0.35, "cm"))
}
vol_T16 <- wide_volcano(vol_T16)
vol_T18 <- wide_volcano(vol_T18)
vol_Tmeta <- wide_volcano(vol_Tmeta)


# Save exactly four combined RAT/DTP supplementary figures.
# S1 remains the existing canopy figure. The new figures are S2–S5.
fig_s2_core <- (tax_T16 | tax_T18 | tax_T18_noA) + plot_annotation(tag_levels = "A")
fig_s2_core <- fig_s2_core & theme(plot.tag = element_text(face = "bold", size = 18))
fig_s2 <- shared_axes(fig_s2_core, "Specimen", "Relative abundance")
fig_s3_core <- (nmds_T16 | nmds_T18 | nmds_Tmeta) + plot_annotation(tag_levels = "A")
fig_s3 <- shared_axes(fig_s3_core, "NMDS1", "NMDS2")
# Separate y-axis titles for richness and alpha diversity, shared within rows.
alpha_top <- (alpha_A | alpha_B) + plot_annotation(tag_levels = list(c("A", "B")))
alpha_bottom <- (alpha_C | alpha_D) + plot_annotation(tag_levels = list(c("C", "D")))
alpha_top <- add_common_y_title(alpha_top, "Observed ASV richness", text_size = 18, width = 0.05)
alpha_bottom <- add_common_y_title(alpha_bottom, "Alpha diversity", text_size = 18, width = 0.05)
fig_s4 <- alpha_top / alpha_bottom
fig_s5_core <- (vol_T16 | vol_T18 | vol_Tmeta) + plot_annotation(tag_levels = "A")
fig_s5 <- shared_axes(fig_s5_core, "Log2 fold change (DTP / RAT)", "−log10 adjusted P")
save_supp("FigureS2_RAT_DTP_taxonomy", fig_s2, 14, 5.5)
save_supp("FigureS3_RAT_DTP_NMDS", fig_s3, 10.5, 4.5)
save_supp("FigureS4_RAT_DTP_alpha_diversity", fig_s4, 9, 7.5)
save_supp("FigureS5_RAT_DTP_volcano", fig_s5, 14, 5.5)


cat("Updated all four supplementary figures. Statistical results are unchanged.\n")
