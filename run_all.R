# Run from the folder with the scripts, six CSV files, and two canopy PNGs.
# The scripts run in order because later files use objects from earlier files.
# Figures and tables are written to output/.

# Load packages, specimen labels, and input data.
source("01_setup.R")
# Figure 2: class-level taxonomic composition.
source("02_taxonomy.R")
# Figures 3 and 4: NMDS and PERMANOVA.
source("03_ordination.R")
# Figure 5: observed richness and alpha diversity.
source("04_alpha_diversity.R")
# Figure 6: bacterial ASV differential abundance.
source("05_asv_abundance.R")
# Figure 7: global protein differential abundance.
source("06_protein_abundance.R")
# Figure 8: partner and functional-role contributions.
source("07_partner_contributions.R")
# Figure 9: the ten functional protein groups.
source("08_functional_groups.R")
# Fig. S1: calculate openness and draw the two-tree figure.
source("09_canopy_openness.R")

cat("\nAnalysis finished. Package versions are listed below.\n")
print(sessionInfo())
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
