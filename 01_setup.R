# Shared setup for Figures 2–9 and the analysis tables.
# Run run_all.R from the folder containing these scripts and CSV files.
# This file loads the inputs and defines names and plot settings used later.

# Input CSV files are in the working folder. Results go into output/.
PROJECT_DIR <- normalizePath(getwd())
DATA_DIR <- PROJECT_DIR
OUTDIR <- file.path(PROJECT_DIR, "output")
INDIVIDUAL_DIR <- file.path(OUTDIR, "individual_panels")
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
dir.create(INDIVIDUAL_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR,"protein_group_analysis","tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR,"protein_group_analysis","functional_groups"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUTDIR,"protein_group_analysis","summaries"), recursive=TRUE, showWarnings=FALSE)
set.seed(42)

# Install missing packages and load them
# phyloseq holds ASV data; DESeq2 tests ASVs; limma tests proteins.
bioc_packages <- c("phyloseq", "DESeq2", "limma")

# vegan: diversity and PERMANOVA; ggordiplots: NMDS plots.
# ggh4x and patchwork: figure panels; ggplot2: plotting.
# dplyr, tidyr, tibble: data tables; RColorBrewer: plot palettes.
# png: canopy photographs; svglite: SVG figures; openxlsx: supplementary tables.
cran_packages <- c("vegan", "ggordiplots", "ggh4x", "ggplot2",
                   "dplyr", "tidyr", "tibble", "patchwork",
                   "RColorBrewer", "png", "svglite", "openxlsx")

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# microViz: taxonomy bar plots; install from the author's R-universe.
if (!requireNamespace("microViz", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  microviz_bioc <- c("microbiome", "ComplexHeatmap")
  microviz_missing <- microviz_bioc[!vapply(microviz_bioc, requireNamespace,
                                          logical(1), quietly = TRUE)]
  if (length(microviz_missing)) {
    BiocManager::install(microviz_missing, ask = FALSE, update = FALSE)
  }
  install.packages("microViz", repos = c(
    "https://david-barnett.r-universe.dev", "https://cloud.r-project.org"))
}

# Load all packages used in scripts 01 to 09.
for (package in c(bioc_packages, cran_packages, "microViz")) {
  library(package, character.only = TRUE)
}
library(grid)  # R package for arranging the canopy photographs.
cat("Packages loaded.\n\n")

# Labels used throughout the study.
# LLS and HLS are lower- and higher-light specimens.
# N0.32 and N1.28 are the two NH4NO3 concentrations in mM.
# RAT remained attached; DTP was detached and transplanted back.
# These colours are for figures only and do not affect any tests.
col_LLS<-"#1A75FF"; col_HLS<-"darkblue"
col_N0.32<-"#e31a1c"; col_N1.28<-"darkred"
col_RAT<-"#7AB84F"; col_DTP<-"darkgreen"
color_values<-c("LLS"=col_LLS,"HLS"=col_HLS,"N0.32"=col_N0.32,"N1.28"=col_N1.28,"RAT"=col_RAT,"DTP"=col_DTP)

disp_label <- c("LLS"="LLS", "HLS"="HLS", "N0.32"="N0.32", "N1.28"="N1.28", "RAT"="RAT", "DTP"="DTP")

# Plot colours for the lichen partners and other associated organisms.
partner_colors<-c(Mycobiont="#A66A3F",Phycobiont="#7EC84F",
                  Cyanobiont="#2A9D8F",Bacteria="#F5D46B",Other="#C77CFF")
partner_levels<-names(partner_colors)
sym_colors_ext<-partner_colors

# Class colours were chosen after inspecting the data for clear figures.
# Readers may change these colours without changing the results.
# The next two palettes are used for the class-level taxonomy bars.
class_colors_16S<-c(
  "Alphaproteobacteria"="darkblue","Verrucomicrobiae"="darkred",
  "Acidobacteriae"="green","Actinobacteria"="black",
  "Gammaproteobacteria"="deeppink","Chloroflexia"="cyan",
  "Bacteroidia"="magenta","Thermoleophilia"="red",
  "Cyanobacteriia"="darkgreen","Planctomycetes"="#FB9A99",
  "Phycisphaerae"="yellow4","Fimbriimonadia"="dodgerblue2",
  "Blastocatellia"="steelblue4","Clostridia"="darkorange")
class_colors_18S<-c(
  "Ascomycota"="darkred","Basidiomycota"="#FF7F00","Chytridiomycota"="#FDBF6F",
  "Microsporidiomycota"="#CAB2D6","Chlorophyceae"="green","Trebouxiophyceae"="darkgreen",
  "Embryophyceae"="#33A02C","Arthropoda"="#1F78B4","Colpodea"="#A6CEE3",
  "Oligohymenophorea"="#B2DF8A","Gregarinomorphea"="#E31A1C",
  "Filosa-Sarcomonadea"="#6A3D9A","Filosa-Thecofilosea"="#B15928")

# Class colours for the significant ASVs in the volcano plots.
vol_colors_16S<-c(
  "Alphaproteobacteria"="#1F78B4","Verrucomicrobiae"="#E31A1C",
  "Acidobacteriae"="#33A02C","Actinobacteria"="#FF7F00",
  "Gammaproteobacteria"="#6A3D9A","Chloroflexia"="#A6CEE3",
  "Bacteroidia"="#B2DF8A","Thermoleophilia"="#FDBF6F",
  "Cyanobacteriia"="#2C6E49","Planctomycetes"="#FB9A99",
  "Phycisphaerae"="#CAB2D6","Fimbriimonadia"="#B15928",
  "Blastocatellia"="#4D9221","Clostridia"="#762A83")
vol_colors_18S<-c(
  "Ascomycota"="#D73027","Basidiomycota"="#FF7F00","Chytridiomycota"="#FDBF6F",
  "Microsporidiomycota"="#CAB2D6","Chlorophyceae"="#1A9850","Trebouxiophyceae"="#4575B4",
  "Embryophyceae"="#33A02C","Arthropoda"="#1F78B4","Colpodea"="#A6CEE3",
  "Oligohymenophorea"="#B2DF8A","Gregarinomorphea"="#E31A1C",
  "Filosa-Sarcomonadea"="#6A3D9A","Filosa-Thecofilosea"="#B15928")

# Font sizes for text, titles, and treatment labels.
TXT<-12; TTL<-14; STRIP<-11

# The common plotting style is reused across the manuscript figures.
base_theme<-theme_bw()+theme(
  plot.background=element_rect(fill="transparent",colour="transparent"),
  panel.background=element_rect(fill="transparent",colour="gray"),
  panel.border=element_rect(fill="transparent",colour="gray50"),
  panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
  axis.text=element_text(size=TXT,color="black"),
  axis.title=element_text(size=TTL,color="black"),
  strip.background=element_rect(fill="gray",colour="black"),
  strip.text=element_text(size=STRIP,face="bold",colour="black"),
  legend.text=element_text(size=TXT),
  legend.title=element_text(size=TXT,face="bold"),
  legend.key.size=unit(0.4,"cm"),
  legend.background=element_rect(fill="transparent",colour="transparent"),
  plot.title=element_text(size=TTL,face="bold",hjust=0.5))
alpha_theme<-base_theme+theme(axis.title.x=element_blank(),axis.text.y=element_text(size=12),legend.position="none")

# Save each ggplot as an SVG file in the output folder.
save_svg <- function(base_path, plot_obj, width, height) {
  svg_path <- paste0(base_path, ".svg")
  dir.create(dirname(svg_path), recursive=TRUE, showWarnings=FALSE)
  tryCatch(
    ggsave(svg_path, plot_obj, width=width, height=height, bg="white"),
    error=function(e) cat("  WARNING SVG:", conditionMessage(e), "\n")
  )
  invisible(NULL)
}

# Put one vertical y-axis title beside a set of figure panels.
add_common_y_title <- function(plot_obj, label, text_size=14, width=0.055) {
  y_title_plot <- ggplot() +
    annotate("text", x=0.5, y=0.5, label=label, angle=90,
             size=text_size/.pt, colour="black", hjust=0.5, vjust=0.5) +
    xlim(0,1) + ylim(0,1) + theme_void() +
    theme(plot.margin=margin(0,0,0,0,"pt"))
  patchwork::wrap_plots(
    patchwork::wrap_elements(full=y_title_plot),
    patchwork::wrap_elements(full=plot_obj),
    ncol=2, widths=c(width,1)
  )
}

# Plot style for the protein abundance bars in Figure 9.
nsaf_bar_theme<-theme_bw()+theme(
  panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
  panel.grid.major.y=element_line(color="grey88",linewidth=0.35),
  strip.background=element_rect(fill="grey93",colour="grey65"),
  strip.text.y.right=element_text(size=14,face="bold",angle=90,lineheight=1.2),
  axis.text.x=element_text(size=14,colour="black",lineheight=1.1),
  axis.text.y=element_text(size=12,colour="black"),
  axis.title=element_text(size=11),
  plot.title=element_text(size=13,face="bold",hjust=0.5),
  plot.subtitle=element_text(size=8.5,colour="grey35",hjust=0.5,lineheight=1.3),
  plot.background=element_rect(fill="white",colour="transparent"),
  legend.position="bottom",legend.key.size=unit(0.45,"cm"),
  legend.text=element_text(size=9))

# Treatment labels at the top of the plot panels.
strip_light<-strip_themed(
  text_x=elem_list_text(colour=c(col_LLS,col_HLS),size=c(STRIP,STRIP)),
  background_x=list(element_rect(color=col_LLS,linewidth=1,fill="transparent"),
                    element_rect(color=col_HLS,linewidth=1,fill="transparent")))
strip_nitrogen<-strip_themed(
  text_x=elem_list_text(colour=c(col_N0.32,col_N1.28),size=c(STRIP,STRIP)),
  background_x=list(element_rect(color=col_N0.32,linewidth=1,fill="transparent"),
                    element_rect(color=col_N1.28,linewidth=1,fill="transparent")))
strip_transplant<-strip_themed(
  text_x=elem_list_text(colour=c(col_RAT,col_DTP),size=c(STRIP,STRIP)),
  background_x=list(element_rect(color=col_RAT,linewidth=1,fill="transparent"),
                    element_rect(color=col_DTP,linewidth=1,fill="transparent")))

# Use the manuscript specimen IDs already present in the input CSV files.
# The order is RAT, DTP, then the four light × nitrogen groups.
prot_ids <- c(paste0("RAT", 1:3), paste0("DTP", 1:3),
              paste0("N0.32-LLS", 1:3), paste0("N0.32-HLS", 1:3),
              paste0("N1.28-LLS", 1:3), paste0("N1.28-HLS", 1:3))
all_nsaf_cols <- paste0("nsaf_", prot_ids)
cnt_cols <- paste0("count_", prot_ids)
donor_cols <- paste0("count_RAT", 1:3)
trans_cols <- paste0("count_DTP", 1:3)

target_symbionts_raw <- c("Bacteria","Cyanobiont","Photobiont","Mycobiont")
target_symbionts <- c("Bacteria","Cyanobiont","Phycobiont","Mycobiont")

# Spacing for paired bars in the protein group figure.
BAR_W    <- 0.30
STEP     <- 0.80
HALF_GAP <- 0.17

xpos_g1<-function(i)(i-1)*STEP+(1-HALF_GAP)
xpos_g2<-function(i)(i-1)*STEP+(1+HALF_GAP)
xpos_ctr<-function(i)(i-1)*STEP+1
vline_at<-function(n){if(n<2)return(numeric(0));sapply(seq_len(n-1),function(i)(xpos_g2(i)+xpos_g1(i+1))/2)}

# Ten functional protein groups used in the focused analysis.
photo_fams    <- c("PSII", "RuBisCO_large", "Phycobilisome")

nmet_fams     <- c("Nitrogenase", "Urease", "GS_GOGAT", "Arginine")

stress_fams   <- c("Catalase", "Peroxiredoxin", "GroEL_Hsp60")

light_fams    <- c(photo_fams, stress_fams)

nitrogen_fams <- c(photo_fams, nmet_fams, stress_fams)

all_fam_order <- c(photo_fams, nmet_fams, stress_fams)

# Specimens in each light contrast, within nitrogen groups and pooled.
light_cmps<-list(
  list(name="LLS vs HLS\n(N0.32)",g1=c("nsaf_N0.32-LLS1","nsaf_N0.32-LLS2","nsaf_N0.32-LLS3"),g2=c("nsaf_N0.32-HLS1","nsaf_N0.32-HLS2","nsaf_N0.32-HLS3"),g1_lbl="LLS (N0.32)",g2_lbl="HLS (N0.32)"),
  list(name="LLS vs HLS\n(N1.28)",g1=c("nsaf_N1.28-LLS1","nsaf_N1.28-LLS2","nsaf_N1.28-LLS3"),g2=c("nsaf_N1.28-HLS1","nsaf_N1.28-HLS2","nsaf_N1.28-HLS3"),g1_lbl="LLS (N1.28)",g2_lbl="HLS (N1.28)"),
  list(name="LLS vs HLS",g1=c("nsaf_N0.32-LLS1","nsaf_N0.32-LLS2","nsaf_N0.32-LLS3","nsaf_N1.28-LLS1","nsaf_N1.28-LLS2","nsaf_N1.28-LLS3"),g2=c("nsaf_N0.32-HLS1","nsaf_N0.32-HLS2","nsaf_N0.32-HLS3","nsaf_N1.28-HLS1","nsaf_N1.28-HLS2","nsaf_N1.28-HLS3"),g1_lbl="LLS",g2_lbl="HLS"))
# Specimens in each nitrogen contrast, within light groups and pooled.
nitrogen_cmps<-list(
  list(name="N0.32 vs N1.28\n(LLS)",g1=c("nsaf_N0.32-LLS1","nsaf_N0.32-LLS2","nsaf_N0.32-LLS3"),g2=c("nsaf_N1.28-LLS1","nsaf_N1.28-LLS2","nsaf_N1.28-LLS3"),g1_lbl="N0.32 (LLS)",g2_lbl="N1.28 (LLS)"),
  list(name="N0.32 vs N1.28\n(HLS)",g1=c("nsaf_N0.32-HLS1","nsaf_N0.32-HLS2","nsaf_N0.32-HLS3"),g2=c("nsaf_N1.28-HLS1","nsaf_N1.28-HLS2","nsaf_N1.28-HLS3"),g1_lbl="N0.32 (HLS)",g2_lbl="N1.28 (HLS)"),
  list(name="N0.32 vs N1.28",g1=c("nsaf_N0.32-LLS1","nsaf_N0.32-LLS2","nsaf_N0.32-LLS3","nsaf_N0.32-HLS1","nsaf_N0.32-HLS2","nsaf_N0.32-HLS3"),g2=c("nsaf_N1.28-LLS1","nsaf_N1.28-LLS2","nsaf_N1.28-LLS3","nsaf_N1.28-HLS1","nsaf_N1.28-HLS2","nsaf_N1.28-HLS3"),g1_lbl="N0.32",g2_lbl="N1.28"))
# Separate comparison of the two donor-tree groups.
# The ten functional groups are compared for light and nitrogen only.
all_cmps<-c(light_cmps,nitrogen_cmps)
get_ctx<-function(n)dplyr::case_when(n%in%sapply(light_cmps,`[[`,"name")~"Light",n%in%sapply(nitrogen_cmps,`[[`,"name")~"Nitrogen",TRUE~"Stress")

# Convert P values into the stars shown in figures and tables.
sig_stars <- function(p) dplyr::case_when(
  is.na(p)   ~ "",
  p < 0.001  ~ "***",
  p < 0.01   ~ "**",
  p < 0.05   ~ "*",
  TRUE       ~ "")

p_stars  <- sig_stars
sig_fmt  <- function(p) dplyr::case_when(
  is.na(p)   ~ "ns",
  p < 0.001  ~ "***",
  p < 0.01   ~ "**",
  p < 0.05   ~ "*",
  TRUE       ~ "ns")
assign_partner<-function(sv){sv<-trimws(as.character(sv));dplyr::case_when(sv=="Mycobiont"~"Mycobiont",sv=="Photobiont"~"Phycobiont",sv=="Cyanobiont"~"Cyanobiont",sv=="Bacteria"~"Bacteria",TRUE~"Other")}
make_strip_box<-function(label){ggplot()+annotate("rect",xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill="gray70",color="black",linewidth=0.5)+annotate("text",x=0,y=0,label=label,fontface="bold",size=2.8,color="black")+theme_void()+theme(plot.margin=unit(c(0,0,0,0),"pt"))}

# The metadata and data tables already use the manuscript specimen names.
# Read metadata and set the treatment labels for all specimens.
meta_raw<-read.csv(file.path(DATA_DIR,"meta_data.csv"),header=TRUE,check.names=FALSE,
                   stringsAsFactors=FALSE)

required_meta_cols<-c("ID","Treatment_group","Light","Nitrogen")
missing_meta_cols<-setdiff(required_meta_cols,colnames(meta_raw))
if(length(missing_meta_cols)>0)
  stop("meta_data.csv is missing required columns: ",paste(missing_meta_cols,collapse=", "))

if(!all(meta_raw$ID%in%prot_ids))
  stop("Unknown specimen IDs: ",paste(setdiff(meta_raw$ID,prot_ids),collapse=", "))
if(anyDuplicated(meta_raw$ID))
  stop("Duplicate specimen IDs found in meta_data.csv")
if(nrow(meta_raw)!=18)
  stop("Expected 18 specimens in meta_data.csv, found ",nrow(meta_raw))

# Keep the input ID for matching specimen columns in the other CSV files.
meta_raw$OriginalID<-meta_raw$ID
meta_raw$Treatment_group<-factor(meta_raw$Treatment_group,
  levels=c("RAT","DTP","N0.32-LLS","N0.32-HLS","N1.28-LLS","N1.28-HLS"))
meta_raw$Light<-factor(meta_raw$Light,levels=c("LLS","HLS"))
meta_raw$Nitrogen<-factor(meta_raw$Nitrogen,levels=c("N0.32","N1.28"))

meta_raw$RAT_DTP<-ifelse(as.character(meta_raw$Treatment_group)%in%c("RAT","DTP"),
                         as.character(meta_raw$Treatment_group),NA_character_)
meta_raw$RAT_DTP<-factor(meta_raw$RAT_DTP,levels=c("RAT","DTP"))

meta_raw$Transplant<-meta_raw$RAT_DTP

.expected_groups<-c("RAT","DTP","N0.32-LLS","N0.32-HLS","N1.28-LLS","N1.28-HLS")
.group_counts<-table(meta_raw$Treatment_group,useNA="ifany")
if(any(.group_counts[.expected_groups]!=3L))
  stop("Unexpected Treatment_group counts: ",paste(names(.group_counts),.group_counts,collapse="; "))
if(sum(!is.na(meta_raw$Light)&!is.na(meta_raw$Nitrogen))!=12L)
  stop("Expected 12 Light x Nitrogen factorial specimens")
if(sum(!is.na(meta_raw$RAT_DTP))!=6L)
  stop("Expected 6 specimens in the RAT-DTP assessment")

row.names(meta_raw)<-meta_raw$ID
samp_data<-sample_data(data.frame(meta_raw,row.names=meta_raw$ID))
cat(sprintf("  Samples loaded: %d\n\n",nrow(meta_raw)))

# Read 16S counts and taxonomy; retain classified bacterial ASVs.
asv_16S_raw<-read.csv(file.path(DATA_DIR,"amplicon_16S_asv.csv"),header=TRUE,row.names=1,check.names=FALSE)
tax_16S_raw<-read.csv(file.path(DATA_DIR,"amplicon_16S_tax.csv"),header=TRUE,row.names=1,check.names=FALSE)
tax_16S_raw<-tax_16S_raw[,-6]
asv_16S_raw<-asv_16S_raw[,meta_raw$OriginalID]
colnames(asv_16S_raw)<-meta_raw$ID
phyloseq_16S<-phyloseq(otu_table(as.matrix(asv_16S_raw),taxa_are_rows=TRUE),tax_table(as.matrix(tax_16S_raw)),samp_data)
phyloseq_16S_filt<-subset_taxa(phyloseq_16S,Kingdom=="Bacteria"&!is.na(Kingdom)&!is.na(Phylum)&!is.na(Class)&!is.na(Order)&!is.na(Family)&Family!="Mitochondria"&Order!="Chloroplast")
ps_treated_16S<-prune_samples(!is.na(sample_data(phyloseq_16S_filt)$Light),phyloseq_16S_filt)
ps_transplant_16S<-prune_samples(!is.na(sample_data(phyloseq_16S_filt)$Transplant),phyloseq_16S_filt)
tax_df_16S<-as.data.frame(tax_table(phyloseq_16S_filt))
cat(sprintf("  16S filtered taxa: %d\n",ntaxa(phyloseq_16S_filt)))

# Read 18S counts and taxonomy; retain classified eukaryotic ASVs.
# A second object excludes Ascomycota for the Figure 2 display.
asv_18S_raw<-read.csv(file.path(DATA_DIR,"amplicon_18S_asv.csv"),header=TRUE,row.names=1,check.names=FALSE)
tax_18S_raw<-read.csv(file.path(DATA_DIR,"amplicon_18S_tax.csv"),header=TRUE,row.names=1,check.names=FALSE)
tax_18S_raw<-tax_18S_raw[,-c(7:8)]
asv_18S_raw<-asv_18S_raw[,meta_raw$OriginalID]
colnames(asv_18S_raw)<-meta_raw$ID
phyloseq_18S<-phyloseq(otu_table(as.matrix(asv_18S_raw),taxa_are_rows=TRUE),tax_table(as.matrix(tax_18S_raw)),samp_data)
phyloseq_18S_filt<-subset_taxa(phyloseq_18S,Kingdom=="Eukaryota"&!is.na(Kingdom)&!is.na(Supergroup)&!is.na(Division)&!is.na(Class)&!is.na(Order)&!is.na(Family))

phyloseq_18S_noAsco<-subset_taxa(phyloseq_18S_filt,Class!="Ascomycota")
ps_treated_18S<-prune_samples(!is.na(sample_data(phyloseq_18S_filt)$Light),phyloseq_18S_filt)
ps_transplant_18S<-prune_samples(!is.na(sample_data(phyloseq_18S_filt)$Transplant),phyloseq_18S_filt)
ps_treated_18S_noAsco<-prune_samples(!is.na(sample_data(phyloseq_18S_noAsco)$Light),phyloseq_18S_noAsco)
ps_transplant_18S_noAsco<-prune_samples(!is.na(sample_data(phyloseq_18S_noAsco)$Transplant),phyloseq_18S_noAsco)
tax_df_18S<-as.data.frame(tax_table(phyloseq_18S_filt))
cat(sprintf("  18S filtered taxa: %d\n",ntaxa(phyloseq_18S_filt)))

# Read protein counts and NSAF values from proteome.csv.
# NSAF means normalized spectral abundance factor.
prot<-read.csv(file.path(DATA_DIR,"proteome.csv"),stringsAsFactors=FALSE,check.names=FALSE)

# The proteome columns already carry these specimen IDs.
prot$d<-tolower(prot$desc)
cat(sprintf("  Protein groups loaded: %d\n",nrow(prot)))
meta_cols<-intersect(c("members_identifier","desc","symbol","og","symbionts","subRole"),colnames(prot))

nsaf_avail<-all_nsaf_cols[all_nsaf_cols%in%colnames(prot)]
mat_nsaf_raw<-as.matrix(prot[,nsaf_avail]); rownames(mat_nsaf_raw)<-prot$members_identifier
storage.mode(mat_nsaf_raw)<-"numeric"

prot_nsaf_all<-prot%>%select(all_of(meta_cols),any_of(nsaf_avail))
prot_nsaf_tgt<-prot_nsaf_all%>%filter(symbionts%in%target_symbionts_raw)

nr_raw<-prot%>%mutate(Partner_cat=case_when(symbionts=="Mycobiont"~"Mycobiont",symbionts=="Photobiont"~"Phycobiont",symbionts=="Cyanobiont"~"Cyanobiont",symbionts=="Bacteria"~"Bacteria",TRUE~NA_character_))%>%filter(!is.na(Partner_cat))%>%mutate(Partner_cat=factor(Partner_cat,levels=names(partner_colors)[1:4]))
m_prot<-as.matrix(nr_raw[,cnt_cols]); rownames(m_prot)<-nr_raw$members_identifier; storage.mode(m_prot)<-"integer"

# Summarize protein counts by orthologous group for ordination.
prot_agg<-prot%>%mutate(across(all_of(cnt_cols),as.numeric))%>%filter(rowSums(select(.,all_of(cnt_cols))==0)/18<0.9)%>%group_by(og)%>%summarise(across(all_of(cnt_cols),sum),.groups="drop")%>%as.data.frame()
rownames(prot_agg)<-prot_agg$og; prot_agg<-prot_agg[,-1]
colnames(prot_agg)<-gsub("^count_","",colnames(prot_agg)); prot_agg<-prot_agg[,prot_ids]
prot_matrix<-t(as.matrix(prot_agg)); meta_prot<-meta_raw[prot_ids,]
treated_ids_prot<-prot_ids[!is.na(meta_prot$Light)]; transplant_ids_prot<-prot_ids[!is.na(meta_prot$Transplant)]
metadata_treated_prot<-meta_prot[treated_ids_prot,]; metadata_transplant_prot<-meta_prot[transplant_ids_prot,]
cat("\nAll data loaded.\n\n")

cat("\nPrepare NSAF values\n")
keep_nsaf<-rowSums(mat_nsaf_raw>0,na.rm=TRUE)>=1
# Remove proteins with only zero NSAF values.
# Replace remaining zero or missing entries with one fifth of the
# smallest nonzero NSAF for that protein before the log2 step.
mat_nsaf_filt<-mat_nsaf_raw[keep_nsaf,,drop=FALSE]
meta_nsaf_filt<-prot[keep_nsaf,meta_cols,drop=FALSE]
prot_filt<-prot[keep_nsaf,]
cat(sprintf("  Before filtering : %d proteins\n",nrow(mat_nsaf_raw)))
cat(sprintf("  After  filtering : %d proteins\n",nrow(mat_nsaf_filt)))
cat(sprintf("  Removed          : %d proteins (all-zero NSAF)\n\n",nrow(mat_nsaf_raw)-nrow(mat_nsaf_filt)))
n_imp_nsaf<-sum(mat_nsaf_filt==0|is.na(mat_nsaf_filt))
mat_nsaf_imp<-mat_nsaf_filt
for(i in seq_len(nrow(mat_nsaf_filt))){row_vals<-mat_nsaf_filt[i,];nonzero_vals<-row_vals[!is.na(row_vals)&row_vals>0];if(length(nonzero_vals)==0)next;row_imp_val<-min(nonzero_vals)/5;zeros<-is.na(row_vals)|row_vals==0;mat_nsaf_imp[i,zeros]<-row_imp_val}
row_imp_vals<-apply(mat_nsaf_filt,1,function(r){nz<-r[!is.na(r)&r>0];if(length(nz)==0)return(NA);min(nz)/5})
cat(sprintf("  Imputation: per-row minimum non-zero / 5\n"))
cat(sprintf("  Cells imputed    : %d\n",n_imp_nsaf))
cat(sprintf("  Imputation range : %.3e - %.3e\n\n",min(row_imp_vals,na.rm=TRUE),max(row_imp_vals,na.rm=TRUE)))
mat_nsaf_log2<-log2(mat_nsaf_imp); mat_nsaf_log2[is.infinite(mat_nsaf_log2)]<-NA
write.csv(cbind(meta_nsaf_filt,as.data.frame(mat_nsaf_imp,optional=TRUE)),file.path(OUTDIR,"protein_group_analysis/tables/nsaf_filtered_imputed.csv"),row.names=FALSE)
write.csv(cbind(meta_nsaf_filt,as.data.frame(mat_nsaf_log2,optional=TRUE)),file.path(OUTDIR,"protein_group_analysis/tables/nsaf_filtered_log2.csv"),row.names=FALSE)
cat("Saved: nsaf_filtered_imputed.csv | nsaf_filtered_log2.csv\n")
cat("Section 5 complete.\n\n")

cat("\nFit global protein models\n")
# Compare the six specimen groups using limma. Positive log2 fold
# changes mean more protein in the group named first in the contrast.
# The model needs syntactic column names; the specimen IDs stay unchanged.
model_groups <- c("RAT", "DTP", "N0_32_LLS", "N0_32_HLS",
                  "N1_28_LLS", "N1_28_HLS")
groups_prot <- factor(rep(model_groups, each = 3), levels = model_groups)
design_prot <- model.matrix(~ 0 + groups_prot)
colnames(design_prot) <- model_groups
contrasts_global <- makeContrasts(
  HLS_vs_LLS = (N0_32_HLS + N1_28_HLS)/2 - (N0_32_LLS + N1_28_LLS)/2,
  N1.28_vs_N0.32 = (N1_28_LLS + N1_28_HLS)/2 - (N0_32_LLS + N0_32_HLS)/2,
  DTP_vs_RAT = DTP - RAT,
  Light_x_Nitrogen = (N1_28_HLS - N1_28_LLS) - (N0_32_HLS - N0_32_LLS),
  levels = design_prot)
cat(sprintf("  Proteins in global NSAF matrix : %d\n",nrow(mat_nsaf_log2)))
cat(sprintf("  Samples                        : %d\n",ncol(mat_nsaf_log2)))
cat("  Contrasts: HLS_vs_LLS | N1.28_vs_N0.32 | DTP_vs_RAT | Light_x_Nitrogen\n\n")
fit_global<-lmFit(mat_nsaf_log2,design_prot)
fit_global<-contrasts.fit(fit_global,contrasts_global)
fit_global<-eBayes(fit_global,trend=FALSE,robust=FALSE)
tt_light_raw<-topTable(fit_global,coef="HLS_vs_LLS",number=Inf,sort.by="P",adjust.method="none")
tt_nitrogen_raw<-topTable(fit_global,coef="N1.28_vs_N0.32",number=Inf,sort.by="P",adjust.method="none")
tt_transplant_raw<-topTable(fit_global,coef="DTP_vs_RAT",number=Inf,sort.by="P",adjust.method="none")
tt_interact_raw<-topTable(fit_global,coef="Light_x_Nitrogen",number=Inf,sort.by="P",adjust.method="none")
tt_light_raw$members_identifier<-rownames(tt_light_raw); tt_light_raw$contrast<-"HLS_vs_LLS"
tt_nitrogen_raw$members_identifier<-rownames(tt_nitrogen_raw); tt_nitrogen_raw$contrast<-"N1.28_vs_N0.32"
tt_transplant_raw$members_identifier<-rownames(tt_transplant_raw); tt_transplant_raw$contrast<-"DTP_vs_RAT"
tt_interact_raw$members_identifier<-rownames(tt_interact_raw); tt_interact_raw$contrast<-"Light_x_Nitrogen"
all_global_raw<-rbind(tt_light_raw,tt_nitrogen_raw,tt_transplant_raw,tt_interact_raw)
cat(sprintf("  Total tests for joint BH correction: %d (= %d proteins x 4 contrasts)\n\n",nrow(all_global_raw),nrow(mat_nsaf_log2)))
# Adjust P values jointly across proteins and the four contrasts.
all_global_raw$adj.P.Val<-p.adjust(all_global_raw$P.Value,method="BH")
all_global_raw$significant<-abs(all_global_raw$logFC)>=1&all_global_raw$adj.P.Val<0.05
res_light_prot<-all_global_raw[all_global_raw$contrast=="HLS_vs_LLS",]%>%left_join(nr_raw%>%select(members_identifier,Partner_cat),by="members_identifier")%>%mutate(padj=adj.P.Val,log2FoldChange=logFC,abs_logFC=abs(logFC))
res_nitrogen_prot<-all_global_raw[all_global_raw$contrast=="N1.28_vs_N0.32",]%>%left_join(nr_raw%>%select(members_identifier,Partner_cat),by="members_identifier")%>%mutate(padj=adj.P.Val,log2FoldChange=logFC,abs_logFC=abs(logFC))
res_transplant_prot<-all_global_raw[all_global_raw$contrast=="DTP_vs_RAT",]%>%left_join(nr_raw%>%select(members_identifier,Partner_cat),by="members_identifier")%>%mutate(padj=adj.P.Val,log2FoldChange=logFC,abs_logFC=abs(logFC))
res_interact_prot<-all_global_raw[all_global_raw$contrast=="Light_x_Nitrogen",]%>%left_join(nr_raw%>%select(members_identifier,Partner_cat),by="members_identifier")%>%mutate(padj=adj.P.Val,log2FoldChange=logFC,abs_logFC=abs(logFC))
cat("  --- Global NSAF limma results (joint BH FDR<0.05 & |log2FC|>=1) ---\n")
cat(sprintf("  HLS vs LLS           : %d significant proteins\n",sum(res_light_prot$significant,na.rm=TRUE)))
cat(sprintf("  N1.28 vs N0.32      : %d significant proteins\n",sum(res_nitrogen_prot$significant,na.rm=TRUE)))
cat(sprintf("  DTP vs RAT           : %d significant proteins\n",sum(res_transplant_prot$significant,na.rm=TRUE)))
cat(sprintf("  Light x Nitrogen    : %d significant proteins\n\n",sum(res_interact_prot$significant,na.rm=TRUE)))
write.csv(all_global_raw,file.path(OUTDIR,"protein_group_analysis/tables/global_nsaf_limma_all.csv"),row.names=FALSE)
write.csv(all_global_raw%>%filter(significant==TRUE),file.path(OUTDIR,"protein_group_analysis/tables/global_nsaf_limma_significant.csv"),row.names=FALSE)
cat("Saved: global_nsaf_limma_all.csv | global_nsaf_limma_significant.csv\n")
cat("Section 6 complete.\n\n")

# Figure style: readable black text and bold panel titles.
supp_theme <- theme_bw(base_size = 18) + theme(
  panel.grid = element_blank(), axis.text = element_text(colour = "black"),
  axis.title = element_text(size = 18),
  plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
  plot.tag = element_text(face = "bold", size = 18),
  legend.title = element_text(face = "bold", size = 13),
  legend.text = element_text(size = 14), plot.margin = margin(6, 6, 6, 6))

# Shared axis titles; each NMDS keeps the scale of its own ordination.
shared_axes <- function(plot, x, y) {
  plot <- add_common_y_title(plot, y, text_size = 18, width = 0.025)
  x_title <- grid::textGrob(x, gp = grid::gpar(fontsize = 18, col = "black"))
  patchwork::wrap_plots(plot, patchwork::wrap_elements(full = x_title),
                       ncol = 1, heights = c(1, 0.075))
}

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

# Save JPG, PDF and SVG versions of the main figures.
save_volcano <- function(path, plot, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (extension in c("jpg", "pdf", "svg")) {
    ggsave(paste0(path, ".", extension), plot,
           width = width, height = height, dpi = 300, bg = "white")
  }
}
