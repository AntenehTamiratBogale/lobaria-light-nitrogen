# Figure 2: class-level composition of the 16S and 18S datasets.
# The light and nitrogen panels use the specimens from the factorial experiment.
# The separate RAT/DTP panels are prepared but not in the main Figure 2.

cat("\nTaxonomic composition\n")
# Set colours for the displayed classes. These are plotting choices.
colors_16S<-c(brewer.pal(8,"Set1"),brewer.pal(7,"Pastel2"))
colors_18S<-c(brewer.pal(7,"Dark2"),brewer.pal(8,"Pastel1"))

# Keep the 18S class colours consistent when Ascomycota is hidden.
ref_taxa_18S<-phyloseq_18S_filt%>%prune_samples(!is.na(sample_data(.)$Light),.)%>%tax_top(n=15,rank="Class")%>%sort()
ref_pal_18S<-colors_18S[seq_along(ref_taxa_18S)]; names(ref_pal_18S)<-ref_taxa_18S; ref_pal_18S["Other"]<-"white"

ref_pal_18S_noAsco<-ref_pal_18S[names(ref_pal_18S)!="Ascomycota"]

# The function keeps up to 15 classes and places the rest in Other.
.tax_barplot<-function(ps_obj,grp_var,grp_lvls,pal,strip_obj,title_str,show_y=TRUE,show_leg=TRUE,show_x=TRUE){
  sample_data(ps_obj)[[grp_var]]<-factor(sample_data(ps_obj)[[grp_var]],levels=grp_lvls)
  top_taxa<-ps_obj%>%tax_top(n=15,rank="Class")%>%sort()

  if(any(top_taxa%in%names(pal))){
    p_pal<-pal[top_taxa]; names(p_pal)<-top_taxa
    p_pal[is.na(p_pal)]<-"gray70"
  }else{
    p_pal<-pal[seq_along(top_taxa)]; names(p_pal)<-top_taxa
  }
  p_pal["Other"]<-if("Other"%in%names(pal))pal["Other"] else "white"
  p<-ps_obj%>%tax_sort(by=sum,at="Class")%>%comp_barplot(tax_order=rev(top_taxa),tax_level="Class",n_taxa=15,merge_other=TRUE,other_name="Other",sample_order="asis",palette=p_pal)+
    facet_grid2(as.formula(paste0("~",grp_var)),space="free",scales="free",strip=strip_obj,labeller=as_labeller(disp_label))+
    scale_x_discrete()+
    scale_y_continuous(expand=c(0.01,0))+labs(fill="Class",y="Relative abundance",title=NULL)+
    base_theme+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=12),axis.text.y=element_text(size=12),legend.position="right",legend.text=element_text(size=12),panel.spacing=unit(0,"pt"),legend.spacing.x=unit(0,"pt"),legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0,0,0,0),plot.margin=margin(0,0,0,0,"pt"))+
    guides(fill=guide_legend(reverse=TRUE))
  if(!show_y)p<-p+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  if(!show_leg)p<-p+guides(fill="none")
  if(!show_x)p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
  p}

# Make 16S panels for light, nitrogen, and RAT/DTP.
ps_16S_L<-prune_samples(!is.na(sample_data(phyloseq_16S_filt)$Light),phyloseq_16S_filt)
ps_16S_N<-prune_samples(!is.na(sample_data(phyloseq_16S_filt)$Nitrogen),phyloseq_16S_filt)
ps_16S_T<-prune_samples(!is.na(sample_data(phyloseq_16S_filt)$Transplant),phyloseq_16S_filt)
p_tax_16S_L<-.tax_barplot(ps_16S_L,"Light",c("LLS","HLS"),colors_16S,strip_light,"16S — Light",show_leg=TRUE)
p_tax_16S_N<-.tax_barplot(ps_16S_N,"Nitrogen",c("N0.32","N1.28"),colors_16S,strip_nitrogen,"16S — Nitrogen",show_y=FALSE,show_leg=FALSE)
p_tax_16S_T<-.tax_barplot(ps_16S_T,"Transplant",c("RAT","DTP"),colors_16S,strip_transplant,"16S — Transplant")
p_tax_16S_LN<-(p_tax_16S_L|p_tax_16S_N)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))

p_tax_16S_L_noX<-.tax_barplot(ps_16S_L,"Light",c("LLS","HLS"),colors_16S,strip_light,"16S — Light",show_leg=TRUE,show_x=FALSE)
p_tax_16S_N_noX<-.tax_barplot(ps_16S_N,"Nitrogen",c("N0.32","N1.28"),colors_16S,strip_nitrogen,"16S — Nitrogen",show_y=FALSE,show_leg=FALSE,show_x=FALSE)
p_tax_16S_T_noX<-.tax_barplot(ps_16S_T,"Transplant",c("RAT","DTP"),colors_16S,strip_transplant,"16S — Transplant",show_x=FALSE)
p_tax_16S_LN_noX<-(p_tax_16S_L_noX|p_tax_16S_N_noX)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))

# Make 18S panels for the same comparisons.
ps_18S_L<-prune_samples(!is.na(sample_data(phyloseq_18S_filt)$Light),phyloseq_18S_filt)
ps_18S_N<-prune_samples(!is.na(sample_data(phyloseq_18S_filt)$Nitrogen),phyloseq_18S_filt)
ps_18S_T<-prune_samples(!is.na(sample_data(phyloseq_18S_filt)$Transplant),phyloseq_18S_filt)
p_tax_18S_L<-.tax_barplot(ps_18S_L,"Light",c("LLS","HLS"),colors_18S,strip_light,"18S — Light")
p_tax_18S_N<-.tax_barplot(ps_18S_N,"Nitrogen",c("N0.32","N1.28"),colors_18S,strip_nitrogen,"18S — Nitrogen",show_y=FALSE,show_leg=FALSE)
p_tax_18S_T<-.tax_barplot(ps_18S_T,"Transplant",c("RAT","DTP"),colors_18S,strip_transplant,"18S — Transplant")
p_tax_18S_LN<-(p_tax_18S_L|p_tax_18S_N)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))

p_tax_18S_L_noX<-.tax_barplot(ps_18S_L,"Light",c("LLS","HLS"),colors_18S,strip_light,"18S — Light",show_x=FALSE)
p_tax_18S_N_noX<-.tax_barplot(ps_18S_N,"Nitrogen",c("N0.32","N1.28"),colors_18S,strip_nitrogen,"18S — Nitrogen",show_y=FALSE,show_leg=FALSE,show_x=FALSE)
p_tax_18S_T_noX<-.tax_barplot(ps_18S_T,"Transplant",c("RAT","DTP"),colors_18S,strip_transplant,"18S — Transplant",show_x=FALSE)
p_tax_18S_LN_noX<-(p_tax_18S_L_noX|p_tax_18S_N_noX)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))

# Repeat the 18S plots without Ascomycota to show less abundant classes.
ps_18S_noA_L<-prune_samples(!is.na(sample_data(phyloseq_18S_noAsco)$Light),phyloseq_18S_noAsco)
ps_18S_noA_N<-prune_samples(!is.na(sample_data(phyloseq_18S_noAsco)$Nitrogen),phyloseq_18S_noAsco)
ps_18S_noA_T<-prune_samples(!is.na(sample_data(phyloseq_18S_noAsco)$Transplant),phyloseq_18S_noAsco)
p_tax_18S_noA_L<-.tax_barplot(ps_18S_noA_L,"Light",c("LLS","HLS"),ref_pal_18S_noAsco,strip_light,"18S (no Ascomycota) — Light")
p_tax_18S_noA_N<-.tax_barplot(ps_18S_noA_N,"Nitrogen",c("N0.32","N1.28"),ref_pal_18S_noAsco,strip_nitrogen,"18S (no Ascomycota) — Nitrogen",show_y=FALSE,show_leg=FALSE)
p_tax_18S_noA_T<-.tax_barplot(ps_18S_noA_T,"Transplant",c("RAT","DTP"),ref_pal_18S_noAsco,strip_transplant,"18S (no Ascomycota) — Transplant")
p_tax_18S_noA_LN<-(p_tax_18S_noA_L|p_tax_18S_noA_N)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))

p_tax_18S_noA_L_noX<-.tax_barplot(ps_18S_noA_L,"Light",c("LLS","HLS"),ref_pal_18S_noAsco,strip_light,"18S (no Ascomycota) — Light",show_x=FALSE)
p_tax_18S_noA_N_noX<-.tax_barplot(ps_18S_noA_N,"Nitrogen",c("N0.32","N1.28"),ref_pal_18S_noAsco,strip_nitrogen,"18S (no Ascomycota) — Nitrogen",show_y=FALSE,show_leg=FALSE,show_x=FALSE)
p_tax_18S_noA_T_noX<-.tax_barplot(ps_18S_noA_T,"Transplant",c("RAT","DTP"),ref_pal_18S_noAsco,strip_transplant,"18S (no Ascomycota) — Transplant",show_x=FALSE)
p_tax_18S_noA_LN_noX<-(p_tax_18S_noA_L_noX|p_tax_18S_noA_N_noX)+plot_layout(widths=c(1,1),guides="collect")&theme(legend.position="right",legend.box.spacing=unit(3,"pt"),legend.margin=margin(0,0,0,2,"pt"),legend.box.margin=margin(0),plot.margin=margin(0))
cat("Prepared Light/Nitrogen taxonomy panels for manuscript Figure 2.\n")
cat("Section 7 complete.\n\n")

# Save the three manuscript panels and the combined Figure 2 SVG.
fig2_dir <- file.path(INDIVIDUAL_DIR,"Figure2")
save_svg(file.path(fig2_dir,"Figure2A_16S_taxonomy"), p_tax_16S_LN_noX, width=6.1, height=3.82)
save_svg(file.path(fig2_dir,"Figure2B_18S_taxonomy"), p_tax_18S_LN_noX, width=6.1, height=3.82)
save_svg(file.path(fig2_dir,"Figure2C_18S_no_Ascomycota_taxonomy"), p_tax_18S_noA_LN_noX, width=6.1, height=3.82)

fig2_A_comb <- p_tax_16S_LN_noX & theme(axis.title.y=element_blank())
fig2_B_comb <- p_tax_18S_LN_noX & theme(axis.title.y=element_blank())
fig2_C_comb <- p_tax_18S_noA_LN_noX & theme(axis.title.y=element_blank())
Figure2_core <- patchwork::wrap_plots(
  patchwork::wrap_elements(full = fig2_A_comb),
  patchwork::wrap_elements(full = fig2_B_comb),
  patchwork::wrap_elements(full = fig2_C_comb),
  nrow = 1,
  widths = c(1,1,1)
) + patchwork::plot_annotation(
  tag_levels="A",
  theme=theme(plot.tag=element_text(size=22,face="bold"),plot.tag.position=c(0.005,1.055),plot.margin=margin(3,0,0,0,"pt"))
)
Figure2 <- add_common_y_title(Figure2_core, "Relative abundance", text_size=14, width=0.022)
save_svg(file.path(OUTDIR,"Figure2_taxonomy"), Figure2, width=15.2, height=3.62)
