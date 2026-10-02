# Run from the extracted folder containing the scripts and input files.
# Save printed results, messages and warnings in a console text file.
dir.create("output", recursive = TRUE, showWarnings = FALSE)

run_analysis <- function() {
  console_log <- file("output/console_output.txt", open = "wt")
  sink(console_log, split = TRUE)
  sink(console_log, type = "message")
  previous_warn <- getOption("warn")
  options(warn = 1)

  # Close the log and save package versions even if a script stops.
  on.exit({
    writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
    sink(type = "message")
    sink()
    close(console_log)
    options(warn = previous_warn)
  }, add = TRUE)

  tryCatch({
    cat("Analysis started:", format(Sys.time()), "\n")

    # Load packages and input data.
    source("01_setup.R")
    # Taxonomic composition, Figure 2.
    source("02_taxonomy.R")
    # NMDS and marginal PERMANOVA, Figures 3 and 4.
    source("03_ordination.R")
    # Observed richness and alpha diversity, Figure 5.
    source("04_alpha_diversity.R")
    # ASV differential abundance, Figure 6.
    source("05_asv_abundance.R")
    # Global protein differential abundance, Figure 7.
    source("06_protein_abundance.R")
    # Partner contributions, Figure 8.
    source("07_partner_contributions.R")
    # Functional groups with corrected annotation selection, Figure 9.
    source("08_functional_groups.R")

    # Canopy calculation and Figure S1, when photographs are present.
    if (file.exists("lower_light_tree.png") && file.exists("higher_light_tree.png")) {
      source("09_canopy_openness.R")
    } else {
      cat("Canopy photographs are absent; Figure S1 was not regenerated.\n")
    }

    # Supplementary figures and datasets.
    source("10_supplementary.R")

    cat("\nAnalysis finished:", format(Sys.time()), "\n")
    cat("\nR and package versions\n")
    print(sessionInfo())
  }, error = function(e) {
    cat("\nAnalysis stopped:", conditionMessage(e), "\n")
    stop(e)
  })
}

run_analysis()
