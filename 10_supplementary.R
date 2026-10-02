# Supplementary figures and tables.
# Run after scripts 01–08, or use source("run_all.R").
# RAT and DTP are analysed separately from the factorial experiment.

cat("\nSupplementary figures and tables\n")
library(openxlsx)
SUPP_DIR <- file.path(OUTDIR, "supplementary_figures")
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
# Remove the earlier plot exports that have been replaced by four figures.
old_plots <- c("RAT_DTP_taxonomy_16S", "RAT_DTP_taxonomy_18S",
  "RAT_DTP_taxonomy_18S_without_Ascomycota", "RAT_DTP_NMDS_16S",
  "RAT_DTP_NMDS_18S", "RAT_DTP_NMDS_metaproteome",
  "RAT_DTP_alpha_diversity_16S", "RAT_DTP_alpha_diversity_18S",
  "RAT_DTP_volcano_16S", "RAT_DTP_volcano_18S", "RAT_DTP_volcano_metaproteome",
  "RAT_DTP_partner_composition_metaproteome", "FigureS2_RAT_DTP_composition")
for (name in old_plots) {
  paths <- file.path(SUPP_DIR, paste0(name, c(".jpg", ".pdf", ".svg")))
  unlink(paths[file.exists(paths)])
}


# Save each plot as JPG, PDF and SVG.
save_supp <- function(name, plot, width = 6, height = 5) {
  for (extension in c("jpg", "pdf", "svg")) {
    ggsave(file.path(SUPP_DIR, paste0(name, ".", extension)), plot,
           width = width, height = height, dpi = 300, bg = "white")
  }
}

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
# Prepare the three supplementary volcano panels.
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
save_supp("FigureS2_RAT_DTP_taxonomy", fig_s2, 14, 5.5)
save_supp("FigureS3_RAT_DTP_NMDS", fig_s3, 10.5, 4.5)
save_supp("FigureS4_RAT_DTP_alpha_diversity", fig_s4, 9, 7.5)

# Supplementary Fig. S5: three RAT/DTP panels, with common axis titles.
s5_A <- volcano_strips(vol_T16, "RAT", "DTP", col_RAT, col_DTP)
s5_B <- volcano_strips(vol_T18, "RAT", "DTP", col_RAT, col_DTP)
s5_C <- volcano_strips(vol_Tmeta, "RAT", "DTP", col_RAT, col_DTP)
fig_s5_core <- (s5_A | s5_B | s5_C) + plot_annotation(tag_levels = "A")
fig_s5 <- shared_axes(fig_s5_core, "Log2 fold change (DTP / RAT)", "−log10 adjusted P")
save_supp("FigureS5_RAT_DTP_volcano", fig_s5, 14.8, 4.8)
# Save the separate 18S RAT/DTP volcano panel.
save_volcano(file.path(OUTDIR, "volcano_panels", "18S_RAT_DTP"),
  s5_B + labs(x = "Log2 fold change (DTP / RAT)", y = "−log10 adjusted P"), 5.5, 3.5)


# Class-level relative abundances and Wilcoxon tests.
# BH correction is applied across classes within each contrast and dataset.
class_values <- function(ps) {
  aggregated <- tax_glom(ps, taxrank = "Class", NArm = TRUE)
  relative <- transform_sample_counts(aggregated, function(x) x / sum(x))
  values <- as(otu_table(relative), "matrix")
  if (!taxa_are_rows(relative)) values <- t(values)
  values <- t(values)
  colnames(values) <- as(tax_table(relative), "matrix")[, "Class"]
  metadata <- data.frame(sample_data(relative))[rownames(values), , drop = FALSE]
  list(values = values, metadata = metadata)
}
class_tests <- function(ps, variable, levels) {
  data <- class_values(ps)
  group <- factor(data$metadata[[variable]], levels = levels)
  result <- lapply(seq_len(ncol(data$values)), function(i) {
    x <- data$values[group == levels[1], i]
    y <- data$values[group == levels[2], i]
    data.frame(contrast = paste(levels, collapse = " vs "), Class = colnames(data$values)[i],
               Group1 = levels[1], Group2 = levels[2],
               mean_group1_pct = mean(x) * 100, mean_group2_pct = mean(y) * 100,
               W = unname(wilcox.test(x, y, exact = FALSE)$statistic),
               P = wilcox.test(x, y, exact = FALSE)$p.value)
  })
  result <- bind_rows(result)
  result$FDR <- p.adjust(result$P, method = "BH")
  result
}
classes_16S <- bind_rows(class_tests(ps_treated_16S, "Light", c("LLS", "HLS")),
                        class_tests(ps_treated_16S, "Nitrogen", c("N0.32", "N1.28")),
                        class_tests(ps_transplant_16S, "Transplant", c("RAT", "DTP")))
classes_18S <- bind_rows(class_tests(ps_treated_18S, "Light", c("LLS", "HLS")),
                        class_tests(ps_treated_18S, "Nitrogen", c("N0.32", "N1.28")),
                        class_tests(ps_transplant_18S, "Transplant", c("RAT", "DTP")))

# Preserve the existing class-level factorial interaction analysis.
class_interaction <- function(ps) {
  data <- class_values(ps)
  result <- lapply(seq_len(ncol(data$values)), function(i) {
    y <- data$values[, i]
    if (sum(y > 0, na.rm = TRUE) < 4) return(NULL)
    fit <- anova(lm(y ~ Light * Nitrogen, data = data$metadata))
    data.frame(Class = colnames(data$values)[i], Light_P = fit["Light", "Pr(>F)"],
               Nitrogen_P = fit["Nitrogen", "Pr(>F)"],
               Interaction_P = fit["Light:Nitrogen", "Pr(>F)"])
  })
  bind_rows(result)
}

# Dispersion tests accompanying PERMANOVA.
# Test each grouping separately, using its corresponding distance matrix.
dispersion_test <- function(distance, metadata, variable, dataset, comparison) {
  group <- factor(metadata[[variable]])
  fit <- betadisper(distance, group, type = "median", bias.adjust = TRUE)
  set.seed(345)
  test <- permutest(fit, permutations = 9999)
  data.frame(Dataset = dataset, Comparison = comparison, F = test$tab[1, "F"],
             P = test$tab[1, "Pr(>F)"], Method = "PERMDISP; median; bias adjustment")
}
dispersion_results <- bind_rows(
  dispersion_test(dist_transplant_16S, metadata_transplant_16S, "Transplant", "16S", "RAT vs DTP"),
  dispersion_test(dist_transplant_18S, metadata_transplant_18S, "Transplant", "18S", "RAT vs DTP"),
  dispersion_test(dist_transplant_prot, meta_prot_transplant, "Transplant", "Metaproteome", "RAT vs DTP"),
  dispersion_test(dist_treated_16S, metadata_treated_16S, "Light", "16S", "Light"),
  dispersion_test(dist_treated_16S, metadata_treated_16S, "Nitrogen", "16S", "Nitrogen"),
  dispersion_test(dist_treated_18S, metadata_treated_18S, "Light", "18S", "Light"),
  dispersion_test(dist_treated_18S, metadata_treated_18S, "Nitrogen", "18S", "Nitrogen"),
  dispersion_test(dist_treated_prot, meta_prot_treat, "Light", "Metaproteome", "Light"),
  dispersion_test(dist_treated_prot, meta_prot_treat, "Nitrogen", "Metaproteome", "Nitrogen"))

# NMDS diagnostics and coordinates for all displayed ordinations.
nmds_stress <- data.frame(
  Dataset = c("16S", "16S", "18S", "Metaproteome", "Metaproteome"),
  Ordination = c("Factorial", "RAT vs DTP", "All 18 specimens", "Factorial", "RAT vs DTP"),
  Stress = c(nmds_treated_16S$stress, nmds_transplant_16S$stress, nmds_full_18S$stress,
             nmds_prot_LN$stress, nmds_prot_T$stress), Dimensions = 2,
  Distance = c("Bray-Curtis", "Bray-Curtis", "Bray-Curtis", "Euclidean", "Euclidean"))
coordinate_table <- function(points, dataset, analysis) {
  data.frame(Dataset = dataset, Analysis = analysis, Specimen = rownames(points),
             NMDS1 = points$NMDS1, NMDS2 = points$NMDS2)
}
nmds_coordinates <- bind_rows(
  coordinate_table(pts_treated_16S, "16S", "Factorial"),
  coordinate_table(pts_transplant_16S, "16S", "RAT vs DTP"),
  coordinate_table(pts_treated_18S, "18S", "Factorial; all-specimen fit"),
  coordinate_table(pts_transplant_18S, "18S", "RAT vs DTP; all-specimen fit"),
  coordinate_table(pts_prot_LN, "Metaproteome", "Factorial"),
  coordinate_table(pts_prot_T, "Metaproteome", "RAT vs DTP"))

# Write plain worksheets with bold headings and no coloured table styles.
write_sheet <- function(workbook, name, values, title) {
  if (name %in% names(workbook)) removeWorksheet(workbook, name)
  addWorksheet(workbook, name)
  writeData(workbook, name, title, startRow = 1)
  writeData(workbook, name, values, startRow = 3, withFilter = TRUE)
  addStyle(workbook, name, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  addStyle(workbook, name, createStyle(textDecoration = "bold"),
           rows = 3, cols = seq_len(ncol(values)), gridExpand = TRUE)
  numeric_columns <- which(vapply(values, is.numeric, logical(1)))
  if (length(numeric_columns) > 0 && nrow(values) > 0) {
    addStyle(workbook, name, createStyle(numFmt = "0.###############"),
             rows = seq_len(nrow(values)) + 3, cols = numeric_columns, gridExpand = TRUE)
  }
  # Use scientific notation for very small values so they do not display as zero.
  for (column in numeric_columns) {
    small_rows <- which(!is.na(values[[column]]) & values[[column]] != 0 &
                        abs(values[[column]]) < 1e-12)
    if (length(small_rows) > 0) {
      addStyle(workbook, name, createStyle(numFmt = "0.00000000000000E+00"),
               rows = small_rows + 3, cols = column)
    }
  }
  freezePane(workbook, name, firstActiveRow = 4)
  setColWidths(workbook, name, cols = seq_len(ncol(values)), widths = 20)
}

# Table S1 and S6 metadata come from the supplied six-table workbook.
# All statistical/result worksheets below are rebuilt from this R run.
wb <- loadWorkbook("supplementary_table_template.xlsx")
write_sheet(wb, "S2A_16S_PERMANOVA", filter(perm_full_table, Dataset == "16S"), "Table S2A — 16S PERMANOVA")
write_sheet(wb, "S2B_18S_PERMANOVA", filter(perm_full_table, Dataset == "18S"), "Table S2B — 18S PERMANOVA")
write_sheet(wb, "S2C_Meta_PERMANOVA", filter(perm_full_table, Dataset == "Metaproteome"), "Table S2C — metaproteome PERMANOVA")
alpha_stats <- alpha_all_stats
alpha_stats$Metric <- ifelse(alpha_stats$Metric == "Observed", "Observed ASV richness", "Alpha diversity")
alpha_values <- bind_rows(mutate(alpha_16S, Dataset = "16S"), mutate(alpha_18S, Dataset = "18S"))
names(alpha_values)[names(alpha_values) == "Shannon"] <- "Alpha_diversity"
write_sheet(wb, "S2D_Alpha_Statistics", alpha_stats, "Table S2D — alpha-diversity tests, including RAT vs DTP")
write_sheet(wb, "S2E_Alpha_Values", alpha_values, "Table S2E — specimen alpha diversity and observed richness")
write_sheet(wb, "S2F_Dispersion", dispersion_results, "Table S2F — multivariate dispersion tests")
curated <- bind_rows(lapply(names(fam_data), function(name) mutate(fam_data[[name]], Functional_group = name)))
write_sheet(wb, "S3B_Curated_Proteins", curated, "Table S3B — final assignments to ten functional protein groups")
write_sheet(wb, "S4A_16S_Classes", classes_16S, "Table S4A — 16S class-level comparisons")
write_sheet(wb, "S4B_18S_Classes", classes_18S, "Table S4B — 18S class-level comparisons")
write_sheet(wb, "S4C_Partner_Pairs", pc_tbl, "Table S4C — partner-pair comparisons; not RAT/DTP tests")
write_sheet(wb, "S4D_16S_Interaction", class_interaction(ps_treated_16S), "Table S4D — 16S class-level factorial models; raw P values")
write_sheet(wb, "S4E_18S_Interaction", class_interaction(ps_treated_18S), "Table S4E — 18S class-level factorial models; raw P values")
write_sheet(wb, "S5A_NMDS_Stress", nmds_stress, "Table S5A — NMDS diagnostics")
write_sheet(wb, "S5B_NMDS_Coordinates", nmds_coordinates, "Table S5B — NMDS coordinates")
contents <- data.frame(Table = paste0("S", 1:6), Contents = c(
  "Specimen and tree metadata", "PERMANOVA, alpha diversity and dispersion",
  "Search terms and final functional-group assignments",
  "Class-level tests and partner-pair comparisons",
  "NMDS stress and specimen coordinates", "Lopu2020nr database taxa"))
write_sheet(wb, "README", contents, "Supplementary Tables S1–S6 — worksheets grouped by table number")
# Use current terminology in the search-term worksheet.
keywords <- read.xlsx("supplementary_table_template.xlsx", sheet = "S3A_Keywords", startRow = 3)
keywords[] <- lapply(keywords, function(x) if (is.character(x)) {
  x <- gsub("subsetted|subsetting", "functional-group selection", x, ignore.case = TRUE)
  gsub("functional subset", "functional group", x, ignore.case = TRUE)
} else x)
# Record the search and exclusion rules used in script 08.
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("RuBisCO", "RuBisCO_large"), na.rm = TRUE)))
keywords[group_row, 5] <- RuBisCO_large_search
keywords$Curation[group_row] <- paste("RuBisCO large and small subunits; exclude activase, carboxysome shell proteins, RuBisCO-like proteins and subunit-modifying methyltransferases.", "Exclusion pattern:", RuBisCO_large_exclude)
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("Nitrogenase", "Nitrogenase"), na.rm = TRUE)))
keywords[group_row, 5] <- Nitrogenase_search
keywords$Curation[group_row] <- paste("Explicit nitrogenase annotations; exclude accessory proteins and shared family annotations that do not distinguish nitrogenase from other functions.", "Exclusion pattern:", Nitrogenase_exclude)
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("GS-GOGAT", "GS_GOGAT"), na.rm = TRUE)))
keywords[group_row, 5] <- GS_GOGAT_search
keywords$Curation[group_row] <- paste("Explicit glutamine synthetase or glutamate synthase annotations; exclude regulators, unrelated biosynthetic enzymes and generic amidotransferases.", "Exclusion pattern:", GS_GOGAT_exclude)
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("Peroxiredoxin", "Peroxiredoxin"), na.rm = TRUE)))
keywords[group_row, 5] <- Peroxiredoxin_search
keywords$Curation[group_row] <- paste("Peroxiredoxin/AhpC proteins; exclude AhpF and proteins described only as regenerating AhpC.", "Exclusion pattern:", Peroxiredoxin_exclude)
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("GroEL/Hsp60", "GroEL_Hsp60"), na.rm = TRUE)))
keywords[group_row, 5] <- GroEL_Hsp60_search
keywords$Curation[group_row] <- paste("GroEL/Hsp60 chaperonins and explicitly annotated GroES/Cpn60 cochaperonins; exclude alcohol dehydrogenases, negative regulators, TCP-1 family and generic GroES-like domains.", "Exclusion pattern:", GroEL_Hsp60_exclude)
group_row <- which(apply(keywords, 1, function(x) any(x %in% c("Catalase", "Catalase"), na.rm = TRUE)))
keywords[group_row, 5] <- Catalase_search
keywords$Curation[group_row] <- paste("Catalase and catalase-peroxidase annotations; hydroperoxidase alone is insufficient.", "Exclusion pattern:", Catalase_exclude)
psii_row <- which(apply(keywords, 1, function(x) any(x == "PSII", na.rm = TRUE)))
keywords$Curation[psii_row] <- "PSII components and associated proteins; cytochrome b6-f components excluded."
group_row <- which(apply(keywords, 1, function(x) any(x == "Urease", na.rm = TRUE)))
keywords$Curation[group_row] <- "Urease catalytic subunits and explicitly annotated urease accessory proteins."
group_row <- which(apply(keywords, 1, function(x) any(x == "Phycobilisome", na.rm = TRUE)))
keywords$Curation[group_row] <- "Phycobiliproteins, linker proteins and explicitly annotated assembly-related proteins."
write_sheet(wb, "S3A_Keywords", keywords, "Table S3A — search terms and manual functional-group selection")
saveWorkbook(wb, file.path(OUTDIR, "Supplementary_Tables_S1-S6.xlsx"), overwrite = TRUE)

# Supplemental Data S1 contains complete differential-abundance results.
# Functional_Groups contains only light, nitrogen and their interaction.
wb_data <- createWorkbook()
results_tables <- list(
  "16S_Light" = res_16S_L, "16S_Nitrogen" = res_16S_N,
  "16S_Interaction" = res_16S_LN, "16S_RAT_DTP" = res_16S_T,
  "18S_Light" = res_18S_L, "18S_Nitrogen" = res_18S_N,
  "18S_Interaction" = res_18S_LN, "18S_RAT_DTP" = res_18S_T,
  "Protein_All" = all_global_raw, "Protein_RAT_DTP" = res_transplant_prot,
  "Functional_Groups" = fam_df)
for (name in names(results_tables)) {
  write_sheet(wb_data, name, results_tables[[name]], paste("Supplemental Data S1 —", name))
  write.csv(results_tables[[name]], file.path(OUTDIR, paste0(name, ".csv")), row.names = FALSE)
}
write_sheet(wb_data, "README", data.frame(Worksheet = names(results_tables)),
            "Full differential-abundance results.")
saveWorkbook(wb_data, file.path(OUTDIR, "Supplemental_Data_S1.xlsx"), overwrite = TRUE)
cat("Saved supplementary figures, Supplementary_Tables_S1-S6.xlsx and Supplemental_Data_S1.xlsx.\n")

# RAT/DTP PERMANOVA table for the supplementary PDF.
rat_permanova <- perm_full_table[perm_full_table$Contrast == "RAT vs DTP", ]
write.csv(rat_permanova, file.path(OUTDIR, "RAT_DTP_PERMANOVA.csv"), row.names = FALSE)
latex_rows <- vapply(seq_len(nrow(rat_permanova)), function(i) {
  p_text <- if (rat_permanova$p.value[i] < 0.001) "$<0.001$" else
              sprintf("%.3f", rat_permanova$p.value[i])
  sprintf("%s & %d & %.2f & %.3f & %s \\\\", rat_permanova$Dataset[i],
          rat_permanova$Df[i], rat_permanova$R2[i], rat_permanova$F[i], p_text)
}, character(1))
writeLines(c(
  "\\renewcommand{\\thetable}{S\\arabic{table}}",
  "\\setcounter{table}{1}",
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{RAT--DTP PERMANOVA results. Bray--Curtis dissimilarities were used for amplicons and Euclidean distances of log$_2$-transformed, imputed NSAF values for the metaproteome. Tests used only the six RAT/DTP specimens. Complete results are supplied in the Table S2 worksheets.}",
  "\\label{tab:rat-dtp-permanova}",
  "\\begin{tabular}{lrrrr}",
  "\\hline", "Dataset & Df & $R^2$ & Pseudo-$F$ & \\textit{P} \\\\",
  "\\hline", latex_rows, "\\hline", "\\end{tabular}", "\\end{table}"),
  file.path(OUTDIR, "TableS2_RAT_DTP_PERMANOVA.tex"))

# Keep Supplementary Tables S1–S6 in six worksheets.
# Sections A, B and C stay within their respective table worksheet.
cat("\nCombine supplementary table sections into six worksheets\n")
table_path <- file.path(OUTDIR, "Supplementary_Tables_S1-S6.xlsx")
table_names <- getSheetNames(table_path)
table_workbook <- createWorkbook()
for (table_number in 1:6) {
  name <- paste0("S", table_number)
  sections <- table_names[grepl(paste0("^", name, "[A-Z_]"), table_names)]
  addWorksheet(table_workbook, name)
  writeData(table_workbook, name, paste("Supplementary Table", name), startRow = 1)
  addStyle(table_workbook, name, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  next_row <- 3
  for (section in sections) {
    values <- read.xlsx(table_path, sheet = section, startRow = 3, check.names = FALSE)
    names(values) <- sub("Functional.subset", "Protein.group", names(values), fixed = TRUE)
    for (column in names(values)) {
      if (is.character(values[[column]])) {
        values[[column]] <- gsub("Photobiont", "Phycobiont", values[[column]], fixed = TRUE)
        if (section == "S3A_Keywords" && !grepl("keyword|search", column, ignore.case = TRUE))
          values[[column]] <- gsub("final subset", "final functional group", values[[column]], fixed = TRUE)
      }
    }
    writeData(table_workbook, name, section, startRow = next_row)
    writeData(table_workbook, name, values, startRow = next_row + 1)
    addStyle(table_workbook, name, createStyle(textDecoration = "bold"), rows = next_row, cols = 1)
    addStyle(table_workbook, name, createStyle(textDecoration = "bold", wrapText = TRUE),
             rows = next_row + 1, cols = seq_len(ncol(values)), gridExpand = TRUE)
    numeric_columns <- which(vapply(values, is.numeric, logical(1)))
    if (length(numeric_columns) && nrow(values))
      addStyle(table_workbook, name, createStyle(numFmt = "General"),
               rows = next_row + 1 + seq_len(nrow(values)), cols = numeric_columns, gridExpand = TRUE)
    next_row <- next_row + nrow(values) + 4
  }
  freezePane(table_workbook, name, firstActiveRow = 5)
  setColWidths(table_workbook, name, cols = 1:30, widths = 24)
}
saveWorkbook(table_workbook, table_path, overwrite = TRUE)
cat("Saved Supplementary Tables S1–S6 with six worksheets.\n")
