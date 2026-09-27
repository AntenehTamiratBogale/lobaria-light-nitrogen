# Figure 9: focused analysis of ten functional protein groups.
# Search terms in the protein descriptions assign proteins to groups.
# The script then compares group NSAF values across treatments and partners.

cat("\nFunctional protein groups\n")

desc_filt <- tolower(meta_nsaf_filt$desc)

# Match protein descriptions to the ten groups and save the matches.
make_family_df <- function(idx, fam_name) {
  df_meta <- meta_nsaf_filt[idx, , drop = FALSE]
  df_nsaf <- as.data.frame(mat_nsaf_filt[idx, , drop = FALSE], optional = TRUE)
  df      <- cbind(df_meta, df_nsaf)
  df$family  <- fam_name
  df$partner <- if ("symbionts" %in% colnames(df)) assign_partner(df$symbionts) else "Other"
  df
}

cat("  Selecting PSII...\n")
# Photosynthesis: PSII, RuBisCO, and phycobilisome proteins.
idx_PSII <- grepl(
  paste("photosystem ii", "psii", "oxygen evolving",
        "psba", "psbb", "psbc", "psbd", "psbe", "psbf", "psbh",
        "psbi", "psbj", "psbk", "psbl", "psbm", "psbn", "psbo",
        "psbp", "psbq", "psbr", "psbs", "psbt", "psbu", "psbv",
        "psbw", "psbx", "psby", "psbz",
        "d1 protein", "d2 protein", "cp47", "cp43", sep = "|"),
  desc_filt, ignore.case = TRUE)
PSII_df   <- make_family_df(idx_PSII, "PSII")
PSII_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "PSII_nsaf.csv")
write.csv(PSII_df, PSII_path, row.names = FALSE)
imp_PSII  <- read.csv(PSII_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_PSII),
            paste(sort(unique(imp_PSII$partner)), collapse = ", ")))

cat("  Selecting RuBisCO_large...\n")
idx_RuBisCO_large <- grepl(
  paste("ribulose bisphosphate carboxylase",
        "ribulose-1,5-bisphosphate",
        "rubisco large", "rbcl", "rbcs", "carboxysome", sep = "|"),
  desc_filt, ignore.case = TRUE)
RuBisCO_large_df   <- make_family_df(idx_RuBisCO_large, "RuBisCO_large")
RuBisCO_large_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "RuBisCO_large_nsaf.csv")
write.csv(RuBisCO_large_df, RuBisCO_large_path, row.names = FALSE)
imp_RuBisCO_large  <- read.csv(RuBisCO_large_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_RuBisCO_large),
            paste(sort(unique(imp_RuBisCO_large$partner)), collapse = ", ")))

cat("  Selecting Phycobilisome...\n")
idx_Phycobilisome <- grepl(
  paste("phycobili", "phycocyanin", "phycoerythrin",
        "allophycocyanin", "phycobilisome linker",
        "phycoerythrocyanin", "phycobiliprotein", sep = "|"),
  desc_filt, ignore.case = TRUE)
Phycobilisome_df   <- make_family_df(idx_Phycobilisome, "Phycobilisome")
Phycobilisome_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Phycobilisome_nsaf.csv")
write.csv(Phycobilisome_df, Phycobilisome_path, row.names = FALSE)
imp_Phycobilisome  <- read.csv(Phycobilisome_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Phycobilisome),
            paste(sort(unique(imp_Phycobilisome$partner)), collapse = ", ")))

cat("  Selecting Nitrogenase...\n")
# Nitrogen metabolism: nitrogenase, urease, GS-GOGAT, arginine.
idx_Nitrogenase <- grepl(
  paste("nitrogenase", "nifh", "nifd", "nifk",
        "nitrogen fixat", "dinitrogenase",
        "nifb", "nifu", "nifs", "nifn", sep = "|"),
  desc_filt, ignore.case = TRUE)
Nitrogenase_df   <- make_family_df(idx_Nitrogenase, "Nitrogenase")
Nitrogenase_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Nitrogenase_nsaf.csv")
write.csv(Nitrogenase_df, Nitrogenase_path, row.names = FALSE)
imp_Nitrogenase  <- read.csv(Nitrogenase_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Nitrogenase),
            paste(sort(unique(imp_Nitrogenase$partner)), collapse = ", ")))

cat("  Selecting Urease...\n")
idx_Urease <- grepl(
  paste("urease", "urea amidohydrolase",
        "ureab", "ureac", "uread",
        "urease alpha", "urease beta", "urease accessory", sep = "|"),
  desc_filt, ignore.case = TRUE)
Urease_df   <- make_family_df(idx_Urease, "Urease")
Urease_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Urease_nsaf.csv")
write.csv(Urease_df, Urease_path, row.names = FALSE)
imp_Urease  <- read.csv(Urease_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Urease),
            paste(sort(unique(imp_Urease$partner)), collapse = ", ")))

cat("  Selecting GS_GOGAT...\n")
idx_GS_GOGAT <- grepl(
  paste("glutamine synthetase", "glutamate synthase", "gogat",
        "glna", "gltb", "gltd", "glutamine amidotransferase",
        "ferredoxin-dependent glutamate synthase",
        "nadph.*glutamate synthase", sep = "|"),
  desc_filt, ignore.case = TRUE)
GS_GOGAT_df   <- make_family_df(idx_GS_GOGAT, "GS_GOGAT")
GS_GOGAT_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "GS_GOGAT_nsaf.csv")
write.csv(GS_GOGAT_df, GS_GOGAT_path, row.names = FALSE)
imp_GS_GOGAT  <- read.csv(GS_GOGAT_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_GS_GOGAT),
            paste(sort(unique(imp_GS_GOGAT$partner)), collapse = ", ")))

cat("  Selecting Arginine...\n")
arginine_patterns <- c(

  "argininosuccinate synthase",
  "argininosuccinate synthetase",
  "arginosuccinate synthase",
  "argininosuccinate lyase",

  "\\barginase\\b",
  "arginine deiminase",
  "arginine decarboxylase",

  "n-acetylglutamate synthase",
  "acetylglutamate kinase",
  "n-acetyl-5-\\s*glutamyl phosphate",
  "acetylornithine aminotransferase",
  "acetylornithine deacetylase",
  "ornithine acetyltransferase",
  "glutamate n-acetyltransferase",
  "cyclic version of arginine biosynthesis",

  "ornithine carbamoyltransferase",
  "ornithine transcarbamylase",
  "produce l-citrulline"
)
idx_Arginine <- grepl(
  paste(arginine_patterns, collapse = "|"),
  desc_filt,
  ignore.case = TRUE,
  perl = TRUE
)
Arginine_df   <- make_family_df(idx_Arginine, "Arginine")
Arginine_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Arginine_nsaf.csv")
write.csv(Arginine_df, Arginine_path, row.names=FALSE)
imp_Arginine  <- read.csv(Arginine_path, stringsAsFactors=FALSE, check.names=FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Arginine),
            paste(sort(unique(imp_Arginine$partner)), collapse=", ")))

cat("  Selecting Peroxiredoxin...\n")
# Cellular stress: peroxiredoxin, GroEL/Hsp60, and catalase.
idx_Peroxiredoxin <- grepl(
  paste("peroxiredoxin", "thioredoxin peroxidase",
        "1-cys peroxiredoxin", "2-cys peroxiredoxin",
        "alkyl hydroperoxide reductase",
        "ahpc", "ahpf", "bcp", "thiol peroxidase", sep = "|"),
  desc_filt, ignore.case = TRUE)
Peroxiredoxin_df   <- make_family_df(idx_Peroxiredoxin, "Peroxiredoxin")
Peroxiredoxin_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Peroxiredoxin_nsaf.csv")
write.csv(Peroxiredoxin_df, Peroxiredoxin_path, row.names = FALSE)
imp_Peroxiredoxin  <- read.csv(Peroxiredoxin_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Peroxiredoxin),
            paste(sort(unique(imp_Peroxiredoxin$partner)), collapse = ", ")))

cat("  Selecting GroEL_Hsp60...\n")
idx_GroEL_Hsp60 <- grepl(
  paste("tcp-1/cpn60", "chaperonin 60", "groel", "hsp60", "cpn60",
        "groes", "chaperonin, 60 kd", "60 kda chaperonin", sep = "|"),
  desc_filt, ignore.case = TRUE)
GroEL_Hsp60_df   <- make_family_df(idx_GroEL_Hsp60, "GroEL_Hsp60")
GroEL_Hsp60_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "GroEL_Hsp60_nsaf.csv")
write.csv(GroEL_Hsp60_df, GroEL_Hsp60_path, row.names = FALSE)
imp_GroEL_Hsp60  <- read.csv(GroEL_Hsp60_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_GroEL_Hsp60),
            paste(sort(unique(imp_GroEL_Hsp60$partner)), collapse = ", ")))

cat("  Selecting Catalase...\n")
idx_Catalase <- grepl(
  paste("catalase", "catalase-peroxidase", "katg",
        "manganese containing catalase", "hydroperoxidase", sep = "|"),
  desc_filt, ignore.case = TRUE)
Catalase_df   <- make_family_df(idx_Catalase, "Catalase")
Catalase_path <- file.path(OUTDIR, "protein_group_analysis/functional_groups", "Catalase_nsaf.csv")
write.csv(Catalase_df, Catalase_path, row.names = FALSE)
imp_Catalase  <- read.csv(Catalase_path, stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("    %d proteins | partners: %s\n\n",
            nrow(imp_Catalase),
            paste(sort(unique(imp_Catalase$partner)), collapse = ", ")))

# Put the ten groups into one list for the following comparisons.
fam_data <- list(
  PSII          = imp_PSII,
  RuBisCO_large = imp_RuBisCO_large,
  Phycobilisome = imp_Phycobilisome,
  Nitrogenase   = imp_Nitrogenase,
  Urease        = imp_Urease,
  GS_GOGAT      = imp_GS_GOGAT,
  Arginine      = imp_Arginine,
  Peroxiredoxin = imp_Peroxiredoxin,
  GroEL_Hsp60   = imp_GroEL_Hsp60,
  Catalase      = imp_Catalase
)

cat(sprintf("  Comparison matrix assembled: %d protein groups\n", length(fam_data)))
cat(sprintf("  Protein groups: %s\n", paste(names(fam_data), collapse = " | ")))
cat("  N-metabolism groups: Nitrogenase | Urease | GS_GOGAT | Arginine\n\n")
cat("Section 15 complete.\n\n")

cat("\nCompare functional groups\n")

con_tar <- file(
  file.path(OUTDIR, "protein_group_analysis/summaries/protein_group_analysis_results.txt"), "wt")
writeLines(c(
  strrep("=", 70),
  "NSAF PROTEIN-GROUP ANALYSIS — FULL RESULTS",
  paste("Generated:", Sys.time()),
  "Method: mean of nonzero NSAF values per functional group per specimen",
  "        log2 transformed → lmFit + eBayes",
  "        BH correction per comparison across all 10 curated protein groups",
  strrep("=", 70), ""), con_tar)

cat("  Computing log2 functional-group means per specimen:\n")

# For each specimen, average nonzero NSAF values within each group.
# Transform those means to log2 before fitting the models.
nsaf_family_means <- list()

for (fam in names(fam_data)) {
  df      <- fam_data[[fam]]
  nsaf_av <- nsaf_avail[nsaf_avail %in% colnames(df)]
  mat_f   <- as.matrix(df[, nsaf_av])
  mat_f[mat_f == 0] <- NA
  mn      <- colMeans(mat_f, na.rm = TRUE)
  mn[is.nan(mn)] <- 0
  nsaf_family_means[[fam]] <- mn
  cat(sprintf("    %-18s: %d proteins | mean range [%.2e, %.2e]\n",
              fam, nrow(df), min(mn, na.rm = TRUE), max(mn, na.rm = TRUE)))
}

fam_log2_mat <- do.call(rbind, lapply(nsaf_family_means, function(x) log2(x + 1e-15)))
colnames(fam_log2_mat) <- nsaf_avail
cat(sprintf("\n  log2 functional-group matrix: %d groups × %d specimens\n\n",
            nrow(fam_log2_mat), ncol(fam_log2_mat)))

n0_32_lls_cols <- c("nsaf_N0.32-LLS1", "nsaf_N0.32-LLS2", "nsaf_N0.32-LLS3")[
  c("nsaf_N0.32-LLS1", "nsaf_N0.32-LLS2", "nsaf_N0.32-LLS3") %in% colnames(fam_log2_mat)]
n0_32_hls_cols <- c("nsaf_N0.32-HLS1", "nsaf_N0.32-HLS2", "nsaf_N0.32-HLS3")[
  c("nsaf_N0.32-HLS1", "nsaf_N0.32-HLS2", "nsaf_N0.32-HLS3") %in% colnames(fam_log2_mat)]
n1_28_lls_cols <- c("nsaf_N1.28-LLS1", "nsaf_N1.28-LLS2", "nsaf_N1.28-LLS3")[
  c("nsaf_N1.28-LLS1", "nsaf_N1.28-LLS2", "nsaf_N1.28-LLS3") %in% colnames(fam_log2_mat)]
n1_28_hls_cols <- c("nsaf_N1.28-HLS1", "nsaf_N1.28-HLS2", "nsaf_N1.28-HLS3")[
  c("nsaf_N1.28-HLS1", "nsaf_N1.28-HLS2", "nsaf_N1.28-HLS3") %in% colnames(fam_log2_mat)]
all_treat <- c(n0_32_lls_cols, n0_32_hls_cols, n1_28_lls_cols, n1_28_hls_cols)

Light_grp <- factor(c(rep("LLS", 3), rep("HLS", 3), rep("LLS", 3), rep("HLS", 3)))
N_grp     <- factor(c(rep("N0.32", 3), rep("N0.32", 3), rep("N1.28", 3), rep("N1.28", 3)))
des_int   <- model.matrix(~ Light_grp * N_grp)

cat(sprintf("  Interaction design: %d samples × %d terms\n",
            nrow(des_int), ncol(des_int)))
cat(sprintf("  Interaction term  : %s\n\n", colnames(des_int)[ncol(des_int)]))

cat("  Running lmFit per comparison per functional group...\n\n")

# Fit each group contrast with limma; adjust P values across ten
# groups separately for each comparison.
all_raw_results <- list()

for (cmp in all_cmps) {

  g1c      <- cmp$g1[cmp$g1 %in% colnames(fam_log2_mat)]
  g2c      <- cmp$g2[cmp$g2 %in% colnames(fam_log2_mat)]
  if (!length(g1c) || !length(g2c)) next

  cmp_clean <- gsub("\n", " ", cmp$name)
  ctx        <- get_ctx(cmp$name)

  fams_in_ctx <- switch(ctx,
                        Light    = names(fam_data),
                        Nitrogen = names(fam_data),
                        Stress   = names(fam_data),
                        names(fam_data))

  cat(sprintf("  [%-9s] %s\n", ctx, cmp_clean))
  cat(sprintf("             families tested: %s\n",
              paste(fams_in_ctx, collapse = ", ")))

  for (fam in fams_in_ctx) {

    sub_mat <- fam_log2_mat[fam, c(g1c, g2c), drop = FALSE]
    grp     <- factor(c(rep(cmp$g1_lbl, length(g1c)), rep(cmp$g2_lbl, length(g2c))),
                      levels = c(cmp$g1_lbl, cmp$g2_lbl))
    des     <- model.matrix(~ grp)

    fit <- tryCatch(
      eBayes(lmFit(sub_mat, des), trend = FALSE, robust = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) next

    tt <- as.data.frame(
      topTable(fit, coef = 2, number = Inf, sort.by = "p", adjust.method = "none"))

    all_raw_results[[paste(fam, cmp_clean)]] <- data.frame(
      family     = fam,
      comparison = cmp_clean,
      context    = ctx,
      g1_lbl     = cmp$g1_lbl,
      g2_lbl     = cmp$g2_lbl,
      log2FC     = tt$logFC[1],
      AveExpr    = tt$AveExpr[1],
      P.Value    = tt$P.Value[1],
      stringsAsFactors = FALSE)

    direction <- if (tt$logFC[1] > 0)
      paste0("right=", cmp$g2_lbl, " > left=", cmp$g1_lbl)
    else
      paste0("left=",  cmp$g1_lbl, " > right=", cmp$g2_lbl)

    cat(sprintf("             %-18s log2FC = %+.3f  p = %.4f  [%s]\n",
                fam, tt$logFC[1], tt$P.Value[1], direction))
  }
  cat("\n")
}

cat("  Computing light x nitrogen interaction per functional group...\n\n")

for (fam in names(fam_data)) {

  sub_int <- fam_log2_mat[fam, all_treat, drop = FALSE]

  fit_int <- tryCatch(
    eBayes(lmFit(sub_int, des_int), trend = FALSE, robust = FALSE),
    error = function(e) NULL)
  if (is.null(fit_int)) next

  tt_int <- as.data.frame(
    topTable(fit_int, coef = ncol(des_int), number = Inf,
             sort.by = "p", adjust.method = "none"))

  n0_32_lls_mean <- mean(fam_log2_mat[fam, n0_32_lls_cols], na.rm = TRUE)
  n0_32_hls_mean <- mean(fam_log2_mat[fam, n0_32_hls_cols], na.rm = TRUE)
  n1_28_lls_mean <- mean(fam_log2_mat[fam, n1_28_lls_cols], na.rm = TRUE)
  n1_28_hls_mean <- mean(fam_log2_mat[fam, n1_28_hls_cols], na.rm = TRUE)
  int_lfc <- (n0_32_hls_mean - n0_32_lls_mean) - (n1_28_hls_mean - n1_28_lls_mean)

  all_raw_results[[paste(fam, "Light x Nitrogen (Interaction)")]] <- data.frame(
    family     = fam,
    comparison = "Light x Nitrogen (Interaction)",
    context    = "Interaction",
    g1_lbl     = "LLS (N0.32)",
    g2_lbl     = "HLS (N0.32)",
    log2FC     = int_lfc,
    AveExpr    = tt_int$AveExpr[1],
    P.Value    = tt_int$P.Value[1],
    stringsAsFactors = FALSE)

  int_direction <- if (int_lfc > 0)
    "right=HLS > left=LLS (effect larger in left=N0.32)"
  else
    "left=LLS > right=HLS (effect larger in right=N1.28)"

  cat(sprintf("    %-18s interaction log2FC = %+.3f  p = %.4f  [%s]\n",
              fam, int_lfc, tt_int$P.Value[1], int_direction))
}
cat("\n")

cat("  Applying BH correction per comparison across families...\n\n")

FC_THRESH_FAMILY <- 0.5

fam_df_raw          <- dplyr::bind_rows(all_raw_results)
fam_df_raw$adj.P.Val <- NA_real_

for (cmp_name in unique(fam_df_raw$comparison)) {
  idx <- fam_df_raw$comparison == cmp_name
  fam_df_raw$adj.P.Val[idx] <- p.adjust(
    fam_df_raw$P.Value[idx], method = "BH")
  cat(sprintf("    %-45s : %d families, BH applied\n",
              cmp_name, sum(idx)))
}

fam_df_raw$stars <- p_stars(fam_df_raw$adj.P.Val)
fam_df_raw$sig   <- !is.na(fam_df_raw$adj.P.Val) &
  fam_df_raw$adj.P.Val < 0.05    &
  abs(fam_df_raw$log2FC) >= FC_THRESH_FAMILY
fam_df           <- fam_df_raw
sig_df_fam       <- fam_df %>%
  dplyr::filter(sig) %>%
  dplyr::arrange(context, family)

cat(sprintf("\n  --- NSAF protein-group Results (BH per comparison, |log2FC|>=%.1f) ---\n",
            FC_THRESH_FAMILY))
cat(sprintf("  Total tests                          : %d\n", nrow(fam_df)))
cat(sprintf("  Significant (FDR<0.05 & |FC|>=%.1f) : %d\n\n",
            FC_THRESH_FAMILY, nrow(sig_df_fam)))

print(
  fam_df %>%
    mutate(adj.P.Val = signif(adj.P.Val, 3), log2FC = round(log2FC, 3),
           left  = g1_lbl, right = g2_lbl,
           direction = ifelse(log2FC > 0,
                              paste0("right=", g2_lbl, " higher"),
                              paste0("left=",  g1_lbl, " higher"))) %>%
    select(context, family, comparison, left, right,
           log2FC, adj.P.Val, stars, direction) %>%
    arrange(context, family, comparison) %>%
    as.data.frame(),
  row.names = FALSE)

cat("\n  --- Significant families (FDR < 0.05 & |log2FC| >= 0.5) ---\n")
if (nrow(sig_df_fam) > 0) {

  print(
    sig_df_fam %>%
      select(context, family, comparison, g1_lbl, g2_lbl,
             log2FC, adj.P.Val, stars) %>%
      mutate(adj.P.Val = signif(adj.P.Val, 3), log2FC = round(log2FC, 3),
             left  = g1_lbl, right = g2_lbl,
             direction = ifelse(log2FC > 0,
                                paste0("right=", g2_lbl, " higher"),
                                paste0("left=",  g1_lbl, " higher"))) %>%
      select(context, family, comparison, left, right,
             log2FC, adj.P.Val, stars, direction) %>%
      as.data.frame(),
    row.names = FALSE)

  int_sig <- sig_df_fam %>%
    filter(comparison == "Light x Nitrogen (Interaction)")

  if (nrow(int_sig) > 0) {
    cat("\n  Families with significant Light x Nitrogen INTERACTION:\n")
    print(
      int_sig %>%
        select(family, log2FC, adj.P.Val, stars) %>%
        mutate(adj.P.Val = signif(adj.P.Val, 3), log2FC = round(log2FC, 3)),
      row.names = FALSE)
    cat("\n  NOTE: Significant interaction supports ATP-Nitrogenase rescue hypothesis\n")
  } else {
    cat("\n  No significant Light x Nitrogen interaction detected.\n")
    cat("  NOTE: Light and Nitrogen act independently on these protein families.\n")
  }

} else {
  cat("  None\n")
}

writeLines(c(
  strrep("-", 70),
  paste0("BH per comparison across all 10 curated protein groups"),
  paste0("Significance: FDR < 0.05 AND |log2FC| >= 0.5 (protein-group threshold)"),
  paste0("Global volcano threshold: |log2FC| >= 1 (set separately in the global metaproteome analysis)"),
  paste0("Total tests        : ", nrow(fam_df)),
  paste0("Significant        : ", nrow(sig_df_fam)), ""), con_tar)

for (i in seq_len(nrow(fam_df))) {
  writeLines(
    sprintf("  %-20s %-45s log2FC=%+.3f  p.adj=%.4f  %s",
            fam_df$family[i],
            fam_df$comparison[i],
            fam_df$log2FC[i],
            fam_df$adj.P.Val[i],
            fam_df$stars[i]),
    con_tar)
}
close(con_tar)

write.csv(
  fam_df,
  file.path(OUTDIR, "protein_group_analysis/tables/nsaf_family_results.csv"),
  row.names = FALSE)

cat("\nSaved: protein_group_analysis_results.txt | nsaf_family_results.csv\n")
cat("Section 16 complete.\n\n")

cat("\nCompare partner contributions\n")
cat("  Test          : paired Wilcoxon signed-rank (paired by sample)\n")
cat(sprintf("  Samples       : %d (all six specimen groups pooled)\n", ncol(mat_nsaf_log2)))
cat("  Per-sample value per partner = mean NSAF across the group's proteins\n")
cat("                                 assigned to that partner (zeros as NA)\n")
cat("  Correction    : Benjamini-Hochberg, applied separately within each of the\n")
cat("                  6 partner-pair comparisons, across the 10 families\n\n")

# Compare the partner contributions within a protein group using
# paired specimen values, where both partners have proteins.
partner_pair_compare <- function(fam_name, df, partner_a, partner_b) {
  nsaf_av <- nsaf_avail[nsaf_avail %in% colnames(df)]
  if (!length(nsaf_av)) return(NULL)
  mat_all <- as.matrix(df[, nsaf_av, drop = FALSE])
  storage.mode(mat_all) <- "numeric"

  rows_a <- df$partner == partner_a
  rows_b <- df$partner == partner_b
  n_a <- sum(rows_a, na.rm = TRUE); n_b <- sum(rows_b, na.rm = TRUE)
  if (n_a == 0 || n_b == 0) {
    cat(sprintf("  [%s vs %s] %-16s : SKIPPED — %s n=%d protein(s), %s n=%d protein(s) (need >=1 in each)\n",
                partner_a, partner_b, fam_name, partner_a, n_a, partner_b, n_b))
    return(NULL)
  }

  mat_a <- mat_all[rows_a, , drop = FALSE]; mat_a[mat_a == 0] <- NA
  mat_b <- mat_all[rows_b, , drop = FALSE]; mat_b[mat_b == 0] <- NA
  mean_a <- colMeans(mat_a, na.rm = TRUE); mean_a[is.nan(mean_a)] <- 0
  mean_b <- colMeans(mat_b, na.rm = TRUE); mean_b[is.nan(mean_b)] <- 0

  common_samples <- intersect(names(mean_a), names(mean_b))
  x_a <- mean_a[common_samples]; x_b <- mean_b[common_samples]

  wt <- tryCatch(wilcox.test(x_a, x_b, paired = TRUE, exact = FALSE),
                 error = function(e) list(statistic = NA_real_, p.value = NA_real_))

  data.frame(
    comparison    = paste0(partner_a, " vs ", partner_b),
    family        = fam_name,
    partner_A     = partner_a, n_proteins_A = n_a,
    mean_A        = mean(x_a, na.rm = TRUE), sd_A = sd(x_a, na.rm = TRUE),
    partner_B     = partner_b, n_proteins_B = n_b,
    mean_B        = mean(x_b, na.rm = TRUE), sd_B = sd(x_b, na.rm = TRUE),
    n_samples     = length(common_samples),
    V_statistic   = unname(wt$statistic),
    P.Value       = unname(wt$p.value),
    stringsAsFactors = FALSE)
}

partner_pairs <- combn(target_symbionts, 2, simplify = FALSE)

pc_results <- list()
for (pr in partner_pairs) {
  pa <- pr[1]; pb <- pr[2]
  cat(sprintf("  --- %s vs %s ---\n", pa, pb))
  for (fam in names(fam_data)) {
    res <- partner_pair_compare(fam, fam_data[[fam]], pa, pb)
    if (!is.null(res)) pc_results[[paste(pa, pb, fam)]] <- res
  }
  cat("\n")
}
pc_tbl <- dplyr::bind_rows(pc_results)

if (nrow(pc_tbl) > 0) {

  pc_tbl$adj.P.Val <- NA_real_
  for (cmp_nm in unique(pc_tbl$comparison)) {
    idx <- pc_tbl$comparison == cmp_nm
    pc_tbl$adj.P.Val[idx] <- p.adjust(pc_tbl$P.Value[idx], method = "BH")
  }
  pc_tbl$stars     <- p_stars(pc_tbl$adj.P.Val)
  pc_tbl$direction <- ifelse(pc_tbl$mean_A > pc_tbl$mean_B,
                             paste0(pc_tbl$partner_A, " > ", pc_tbl$partner_B),
                             paste0(pc_tbl$partner_B, " > ", pc_tbl$partner_A))
  pc_tbl$log2FC_A_vs_B <- log2((pc_tbl$mean_A + 1e-15) / (pc_tbl$mean_B + 1e-15))

  cat("  --- Detailed results, per partner-pair comparison, per functional group ---\n\n")
  for (cmp_nm in unique(pc_tbl$comparison)) {
    sub <- pc_tbl[pc_tbl$comparison == cmp_nm, ]
    cat(sprintf("  [%s]\n", cmp_nm))
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      cat(sprintf("    %-16s : %s (n=%d proteins) mean=%.4e +/- %.4e | %s (n=%d proteins) mean=%.4e +/- %.4e\n",
                  r$family, r$partner_A, r$n_proteins_A, r$mean_A, r$sd_A,
                  r$partner_B, r$n_proteins_B, r$mean_B, r$sd_B))
      cat(sprintf("    %-16s   n paired samples=%d  V=%s  p=%.4e  p.adj(BH)=%.4e  %s  log2FC(A/B)=%+.3f  [%s]\n\n",
                  "", r$n_samples,
                  ifelse(is.na(r$V_statistic), "NA", sprintf("%.1f", r$V_statistic)),
                  r$P.Value, r$adj.P.Val, r$stars, r$log2FC_A_vs_B, r$direction))
    }
  }

  n_sig_pc <- sum(pc_tbl$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  --- Summary: %d / %d functional-group x partner-pair tests are significant\n",
              n_sig_pc, nrow(pc_tbl)))
  cat("      (FDR<0.05, paired Wilcoxon, corrected within each partner-pair) ---\n\n")
  if (n_sig_pc > 0) {
    sig_pc <- pc_tbl[!is.na(pc_tbl$adj.P.Val) & pc_tbl$adj.P.Val < 0.05, ]
    for (cmp_nm in unique(sig_pc$comparison)) {
      sub <- sig_pc[sig_pc$comparison == cmp_nm, ]
      cat(sprintf("    [%s]\n", cmp_nm))
      for (i in seq_len(nrow(sub))) {
        cat(sprintf("      * %-16s : %-30s  p.adj=%.4e  %s\n",
                    sub$family[i], sub$direction[i], sub$adj.P.Val[i], sub$stars[i]))
      }
    }
    cat("\n")
  } else {
    cat("    None of the 60 functional-group x partner-pair tests are significant.\n\n")
  }

  cat("  --- Full results table (all 6 partner-pairs x 10 families) ---\n")
  print(
    pc_tbl %>%
      dplyr::mutate(mean_A = signif(mean_A, 4), mean_B = signif(mean_B, 4),
                    P.Value = signif(P.Value, 4), adj.P.Val = signif(adj.P.Val, 4),
                    log2FC_A_vs_B = round(log2FC_A_vs_B, 3)) %>%
      dplyr::select(comparison, family, partner_A, n_proteins_A, mean_A,
                    partner_B, n_proteins_B, mean_B,
                    n_samples, P.Value, adj.P.Val, stars, log2FC_A_vs_B, direction) %>%
      dplyr::arrange(comparison, family) %>%
      as.data.frame(),
    row.names = FALSE)
  cat("\n")
} else {
  cat("  No protein families had enough proteins in any partner pair —\n")
  cat("  all comparisons skipped.\n\n")
}

write.csv(pc_tbl,
          file.path(OUTDIR, "protein_group_analysis/tables/partner_pairwise_family_stats.csv"),
          row.names = FALSE)
cat("Saved: partner_pairwise_family_stats.csv (all 6 pairwise partner comparisons x 10 families)\n")
cat("Section 16B complete.\n\n")

cat("\nPrepare Figure 9\n")

.sub_base <- paste0(
  "y=mean NSAF +/- SE | stacked by symbiotic partner | ",
  "stars from family-mean lmFit (BH across 10 groups) | ",
  "*p<0.05 **p<0.01 ***p<0.001 | ")

# Prepare the stacked protein bars for the treatment panels.
build_stacked_data <- function(fam_list, cmp_list) {

  rows <- list()
  for (cmp in cmp_list) {
    cmp_name <- cmp$name
    g1s      <- cmp$g1[cmp$g1 %in% nsaf_avail]
    g2s      <- cmp$g2[cmp$g2 %in% nsaf_avail]
    if (!length(g1s) || !length(g2s)) next

    for (fam in fam_list) {
      if (!fam %in% names(fam_data)) next
      df      <- fam_data[[fam]]
      nsaf_av <- nsaf_avail[nsaf_avail %in% colnames(df)]
      if (!length(nsaf_av)) next

      mat_raw <- as.matrix(df[, nsaf_av])
      mat_raw[mat_raw < 0 | is.na(mat_raw)] <- 0
      mat_na  <- mat_raw
      mat_na[mat_na == 0] <- NA
      fam_mean_per_samp <- colMeans(mat_na, na.rm = TRUE)
      fam_mean_per_samp[is.nan(fam_mean_per_samp)] <- 0
      fam_sum_per_samp  <- colSums(mat_raw, na.rm = TRUE)

      for (grp_idx in 1:2) {
        smp     <- if (grp_idx == 1L) g1s else g2s
        grp_lbl <- if (grp_idx == 1L) cmp$g1_lbl else cmp$g2_lbl
        smp_av  <- smp[smp %in% nsaf_av]
        if (!length(smp_av)) next

        fam_mean_grp <- fam_mean_per_samp[smp_av]
        fam_sum_grp  <- fam_sum_per_samp[smp_av]

        for (pt in partner_levels) {
          pt_rows <- df$partner == pt
          n_pt    <- sum(pt_rows, na.rm = TRUE)

          if (n_pt == 0) {
            rows[[length(rows) + 1]] <- data.frame(
              comparison = cmp_name,
              grp_order  = as.integer(grp_idx),
              group      = grp_lbl,
              partner    = pt,
              family     = fam,
              mean_nsaf  = 0,
              se_nsaf    = 0,
              stringsAsFactors = FALSE)
            next
          }

          pt_sum   <- colSums(mat_raw[pt_rows, smp_av, drop = FALSE], na.rm = TRUE)
          pt_prop  <- ifelse(fam_sum_grp > 0, pt_sum / fam_sum_grp, 0)
          pt_contr <- fam_mean_grp * pt_prop

          rows[[length(rows) + 1]] <- data.frame(
            comparison = cmp_name,
            grp_order  = as.integer(grp_idx),
            group      = grp_lbl,
            partner    = pt,
            family     = fam,
            mean_nsaf  = mean(pt_contr, na.rm = TRUE),
            se_nsaf    = ifelse(length(pt_contr) > 1,
                                sd(pt_contr, na.rm = TRUE) / sqrt(length(pt_contr)),
                                0),
            stringsAsFactors = FALSE)
        }
      }
    }
  }

  summ <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(partner = factor(partner, levels = partner_levels))

  total_se_rows <- list()
  for (cmp in cmp_list) {
    g1s <- cmp$g1[cmp$g1 %in% nsaf_avail]
    g2s <- cmp$g2[cmp$g2 %in% nsaf_avail]
    for (fam in fam_list) {
      if (!fam %in% names(nsaf_family_means)) next
      mn <- nsaf_family_means[[fam]]
      for (g in list(list(s = g1s, o = 1L, l = cmp$g1_lbl),
                     list(s = g2s, o = 2L, l = cmp$g2_lbl))) {
        sv <- mn[g$s[g$s %in% names(mn)]]
        if (!length(sv)) next
        total_se_rows[[length(total_se_rows) + 1]] <- data.frame(
          comparison = cmp$name,
          grp_order  = g$o,
          group      = g$l,
          family     = fam,
          total_nsaf = mean(sv, na.rm = TRUE),
          total_se   = ifelse(length(sv) > 1,
                              sd(sv, na.rm = TRUE) / sqrt(length(sv)),
                              0),
          stringsAsFactors = FALSE)
      }
    }
  }
  total_se <- dplyr::bind_rows(total_se_rows)

  cmp_names_clean <- gsub("\n", " ", sapply(cmp_list, `[[`, "name"))
  wt <- fam_df %>%
    dplyr::filter(
      family     %in% fam_list,
      comparison %in% cmp_names_clean) %>%
    dplyr::select(comparison, family, adj.P.Val, stars, sig, log2FC)

  cat(sprintf("    Annotations: %d rows matched (%d significant)\n",
              nrow(wt), sum(wt$sig, na.rm = TRUE)))

  list(summ = summ, total = total_se, wt = wt)
}

.build_stacked_ggplot <- function(pd, fam_order, cmp_list,
                                  ctx_title, ctx_subtitle) {

  n_fam          <- length(fam_order)
  cmp_names      <- sapply(cmp_list, `[[`, "name")
  cmp_names_clean <- gsub("\n", " ", cmp_names)

  summ <- pd$summ %>%
    dplyr::filter(
      family    %in% fam_order,
      comparison %in% cmp_names,
      mean_nsaf > 0) %>%
    dplyr::mutate(
      comparison = gsub("\n", " ", comparison),
      fam_idx    = match(family, fam_order),
      xpos       = ifelse(grp_order == 1L, xpos_g1(fam_idx), xpos_g2(fam_idx)),
      cmp_title  = factor(comparison, levels = cmp_names_clean),
      partner    = factor(partner, levels = partner_levels))

  total_se <- pd$total %>%
    dplyr::filter(
      family    %in% fam_order,
      comparison %in% cmp_names) %>%
    dplyr::mutate(
      comparison = gsub("\n", " ", comparison),
      fam_idx    = match(family, fam_order),
      xpos       = ifelse(grp_order == 1L, xpos_g1(fam_idx), xpos_g2(fam_idx)),
      cmp_title  = factor(comparison, levels = cmp_names_clean))

  x_breaks <- sapply(seq_len(n_fam), xpos_ctr)
  n_map    <- sapply(fam_order, function(f)
    if (f %in% names(fam_data)) nrow(fam_data[[f]]) else "?")
  display_fams <- gsub("RuBisCO_large", "RuBisCO", fam_order)
  display_fams <- gsub("GS_GOGAT", "GS-GOGAT", display_fams)
  x_labels <- paste0(display_fams, "\n(n=", n_map, ")")
  vlines   <- vline_at(n_fam)

  top_fam <- total_se %>%
    dplyr::group_by(cmp_title, family, fam_idx) %>%
    dplyr::summarise(
      top_y = max(total_nsaf + total_se, na.rm = TRUE),
      .groups = "drop")

  star_df <- pd$wt %>%
    dplyr::filter(
      family    %in% fam_order,
      comparison %in% cmp_names_clean,
      !is.na(stars),
      stars != "") %>%
    dplyr::mutate(
      fam_idx   = match(family, fam_order),
      xpos_star = xpos_ctr(fam_idx),
      cmp_title = factor(comparison, levels = cmp_names_clean)) %>%
    dplyr::left_join(top_fam, by = c("cmp_title", "family", "fam_idx"))

  ggplot(summ, aes(x = xpos, y = mean_nsaf, fill = partner)) +
    geom_vline(
      xintercept = vlines,
      linetype   = "dashed",
      color      = "gray40",
      linewidth  = 0.4) +
    geom_col(
      position  = "stack",
      width     = BAR_W,
      colour    = "white",
      linewidth = 0.2,
      alpha     = 0.95) +
    geom_errorbar(
      data = total_se %>% dplyr::filter(total_nsaf > 0),
      aes(x    = xpos,
          ymin = pmax(0, total_nsaf - total_se),
          ymax = total_nsaf + total_se),
      inherit.aes = FALSE,
      width     = 0.18,
      linewidth = 0.5,
      colour    = "grey25") +
    {if (nrow(star_df) > 0)
      list(
        geom_segment(
          data = star_df,
          aes(x    = xpos_star - BAR_W/2,
              xend = xpos_star + BAR_W/2,
              y    = top_y * 1.12,
              yend = top_y * 1.12),
          inherit.aes = FALSE,
          colour    = "grey20",
          linewidth = 0.6,
          lineend   = "round"),
        geom_text(
          data = star_df,
          aes(x     = xpos_star,
              y     = top_y * 1.18,
              label = stars),
          inherit.aes = FALSE,
          size     = 7,
          fontface = "bold",
          colour   = "grey10",
          hjust    = 0.5))
      else
        list(geom_blank())} +
    scale_fill_manual(
      values = partner_colors,
      name   = "Symbiotic partner",
      guide  = guide_legend(
        reverse    = FALSE,
        override.aes = list(colour = "grey40", linewidth = 0.4))) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(add = 0.02)) +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.18))) +
    facet_wrap(
      ~ cmp_title,
      ncol           = 1,
      strip.position = "right",
      scales         = "fixed") +
    nsaf_bar_theme +
    theme(
      panel.spacing.y = unit(0.6, "lines"),
      strip.placement = "outside") +
    labs(
      title    = ctx_title,
      subtitle = ctx_subtitle,
      x        = NULL,
      y        = "Mean NSAF ± SE (stacked by partner)")
}

photo_fams <- c("PSII", "RuBisCO_large", "Phycobilisome")
nmet_fams  <- c("Nitrogenase", "Urease", "GS_GOGAT", "Arginine")

.add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag          = element_text(size = 23, face = "bold", colour = "black"),
          plot.tag.position = "topleft")
}

.col_header <- function(label, fill_col, mapping_text=NULL) {
  p <- ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
             fill=fill_col, colour=NA) +
    annotate("text", x=0.5, y=0.68, label=label,
             colour="white", fontface="bold", size=7.6) +
    theme_void() +
    theme(plot.margin = margin(3, 2, 3, 2, "pt"))
  if (!is.null(mapping_text)) {
    p <- p + annotate("text", x=0.5, y=0.25, label=mapping_text,
                      colour="white", fontface="bold", size=5.5)
  }
  p
}

.row_label <- function(label) {
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
             fill="grey82", colour=NA) +
    annotate("text", x=0.5, y=0.5, label=label, angle=90,
             colour="black", fontface="bold", size=18/.pt, hjust=0.5, vjust=0.5,
             lineheight=0.92) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
}

.na_panel <- function() {
  ggplot() +
    annotate("text", x=0.5, y=0.5, label="not applicable",
             colour="grey55", size=4.5, fontface="italic") +
    theme_void() +
    theme(panel.border = element_rect(colour="grey80", fill=NA, linewidth=0.3),
          plot.margin  = margin(4, 4, 4, 4, "pt"))
}

.strip_y <- function(p) {
  p + theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
}

.no_legend <- function(p) {
  p + theme(legend.position = "none")
}

cat("  Stress-response and RAT-DTP plots omitted; statistics retained in tables.\n")
cat("Section 17 complete.\n\n")

cat("\nArrange Figure 9\n")

light_main21    <- light_cmps[3]
nitrogen_main21 <- nitrogen_cmps[3]

photo_light_main_data21 <- build_stacked_data(photo_fams, light_main21)
photo_nitro_main_data21 <- build_stacked_data(photo_fams, nitrogen_main21)
photo_light_main21 <- .build_stacked_ggplot(
  pd=photo_light_main_data21, fam_order=photo_fams, cmp_list=light_main21,
  ctx_title=NULL, ctx_subtitle=NULL) + guides(fill="none")
photo_nitro_main21 <- .build_stacked_ggplot(
  pd=photo_nitro_main_data21, fam_order=photo_fams, cmp_list=nitrogen_main21,
  ctx_title=NULL, ctx_subtitle=NULL) + guides(fill="none")
cat("  Built photosynthesis main-effect panels (3 protein groups, no legend).\n")

nmet_light_data21 <- build_stacked_data(nmet_fams, light_main21)
nmet_nitro_data21 <- build_stacked_data(nmet_fams, nitrogen_main21)
nmet_light_main21 <- .build_stacked_ggplot(
  pd=nmet_light_data21, fam_order=nmet_fams, cmp_list=light_main21,
  ctx_title=NULL, ctx_subtitle=NULL)
nmet_nitro_main21 <- .build_stacked_ggplot(
  pd=nmet_nitro_data21, fam_order=nmet_fams, cmp_list=nitrogen_main21,
  ctx_title=NULL, ctx_subtitle=NULL)
cat("  Built N-metabolism main-effect panels (4 protein groups: Nitrogenase, Urease, GS_GOGAT, Arginine).\n\n")

.font_fix21 <- function(p, has_legend=FALSE) {
  p <- p + theme(
    axis.text.x  = element_text(size=24, colour="black", lineheight=1.06),
    axis.text.y  = element_text(size=15, colour="black"),
    axis.title   = element_text(size=24),
    plot.title   = element_text(size=20, face="bold", hjust=0.5, lineheight=1.15),
    plot.subtitle= element_text(size=20, colour="grey35", hjust=0.5, lineheight=1.2),
    strip.text   = element_text(size=24, colour="black"),
    strip.text.y.right = element_text(size=24, face="bold", angle=90, lineheight=1.10),
    plot.margin=margin(0.5,0.5,0.5,0.5,"pt"))
  if (has_legend) p <- p + theme(legend.text=element_text(size=22), legend.title=element_text(size=23,face="bold"),legend.margin=margin(0),legend.box.margin=margin(0),legend.box.spacing=unit(0,"pt"))
  p
}
photo_light_main21 <- .font_fix21(photo_light_main21)
photo_nitro_main21 <- .font_fix21(photo_nitro_main21)
nmet_light_main21  <- .font_fix21(nmet_light_main21,  has_legend=TRUE)
nmet_nitro_main21  <- .font_fix21(nmet_nitro_main21,  has_legend=TRUE)

# Save photosynthesis and nitrogen-metabolism panels separately.
fig9_dir <- file.path(INDIVIDUAL_DIR,"Figure9")

.no_cmp_strip21 <- function(p) {
  p + theme(
    strip.text = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    strip.text.y.left = element_blank(),
    strip.text.y.right = element_blank(),
    strip.background = element_blank()
  )
}

photo_light_ind21 <- .no_cmp_strip21(photo_light_main21) +
  labs(title="Light", subtitle="Left bar = LLS   |   Right bar = HLS")
photo_nitro_ind21 <- .no_cmp_strip21(photo_nitro_main21) +
  labs(title="Nitrogen", subtitle="Left bar = N0.32   |   Right bar = N1.28")
nmet_light_ind21 <- .no_cmp_strip21(nmet_light_main21) +
  labs(title="Light", subtitle="Left bar = LLS   |   Right bar = HLS")
nmet_nitro_ind21 <- .no_cmp_strip21(nmet_nitro_main21) +
  labs(title="Nitrogen", subtitle="Left bar = N0.32   |   Right bar = N1.28")

save_svg(file.path(fig9_dir,"Figure9A_photosynthesis_Light"), photo_light_ind21, width=9.2, height=5.8)
save_svg(file.path(fig9_dir,"Figure9B_N_metabolism_Light"), nmet_light_ind21, width=11.4, height=5.8)
save_svg(file.path(fig9_dir,"Figure9C_photosynthesis_Nitrogen"), photo_nitro_ind21, width=9.2, height=5.8)
save_svg(file.path(fig9_dir,"Figure9D_N_metabolism_Nitrogen"), nmet_nitro_ind21, width=11.4, height=5.8)

.group_col_header21 <- function(title) {
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
             fill="grey82", colour=NA) +
    annotate("text", x=0.5, y=0.5, label=title,
             colour="black", fontface="bold", size=27/.pt,
             hjust=0.5, vjust=0.5) +

    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme_void() +
    theme(plot.margin=margin(0,0,0,0,"pt"))
}

.effect_row_header21 <- function(label, fill_col) {
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
             fill=fill_col, colour=NA) +
    annotate("text", x=0.5, y=0.5, label=label, angle=90,
             colour="black", fontface="bold", size=26/.pt,
             hjust=0.5, vjust=0.5) +
    theme_void() +
    theme(plot.margin=margin(0,0,0,0,"pt"))
}

.comparison_row_strip21 <- function(label) {
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf,
             fill="grey82", colour=NA) +
    annotate("text", x=0.5, y=0.5, label=label, angle=90,
             colour="black", fontface="bold", size=24/.pt,
             hjust=0.5, vjust=0.5) +
    theme_void() +
    theme(plot.margin=margin(0,0,0,0,"pt"))
}

.hide_x21 <- function(p) {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

pA21 <- .add_tag(.hide_x21(.no_cmp_strip21(photo_light_main21)) +
                   theme(axis.title.y=element_blank()), "A")
pB21 <- .add_tag(.hide_x21(.no_cmp_strip21(nmet_light_main21)) |>
                   .strip_y(), "B")
pC21 <- .add_tag(.no_cmp_strip21(photo_nitro_main21) +
                   theme(axis.title.y=element_blank()), "C")
pD21 <- .add_tag(.no_cmp_strip21(nmet_nitro_main21) |>
                   .strip_y(), "D")

left_strip_w21  <- 0.34
right_strip_w21 <- 0.46
col_header_h21  <- 1.02
cell_w21_photo  <- length(photo_fams)*1.85 + 2.35
cell_w21_nmet   <- length(nmet_fams)*1.85 + 2.35
cell_h21        <- 5.55

hdr_photo21 <- .group_col_header21("Photosynthesis protein groups")
hdr_nmet21 <- .group_col_header21("N metabolism protein groups")
row_light21 <- .effect_row_header21("Light", "grey82")
row_nitro21 <- .effect_row_header21("Nitrogen", "grey82")
cmp_light21 <- .comparison_row_strip21("LLS vs HLS")
cmp_nitro21 <- .comparison_row_strip21("N0.32 vs N1.28")

# Combine the four main treatment panels into Figure 9.
fig9_design21 <- "
#AB#
CDEF
GHIJ
"

fig_21_core <- wrap_plots(
  A=hdr_photo21, B=hdr_nmet21,
  C=row_light21, D=pA21, E=pB21, F=cmp_light21,
  G=row_nitro21, H=pC21, I=pD21, J=cmp_nitro21,
  design=fig9_design21,
  widths=c(left_strip_w21, cell_w21_photo, cell_w21_nmet, right_strip_w21),
  heights=c(col_header_h21, cell_h21, cell_h21),
  guides="collect") &
  theme(legend.position="bottom",plot.margin=margin(0),legend.margin=margin(0),legend.box.margin=margin(0),legend.box.spacing=unit(0,"pt"))

fig_21 <- add_common_y_title(
  fig_21_core, "Mean NSAF ± SE (stacked by partner)", text_size=25, width=0.022)
fig_21_w <- left_strip_w21 + cell_w21_photo + cell_w21_nmet + right_strip_w21 + 0.06
fig_21_h <- col_header_h21 + cell_h21*2 + 0.06

save_svg(file.path(OUTDIR,"Figure9_protein_groups"), fig_21,
         width=fig_21_w, height=fig_21_h)
cat("Saved: Figure9_protein_groups.svg\n")
cat("       Left strips: Light | Nitrogen (90-degree text, thin boxes)\n")
cat("       Column headers show the two sets of functional groups\n")
cat("       Far-right strips: LLS vs HLS | N0.32 vs N1.28 (one strip per row; no internal duplicates)\n")
cat("       Shared x-axis labels: lower row supplies labels for upper row\n\n")
cat("Section 21 complete.\n\n")
