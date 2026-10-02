# Figure 7: global metaproteome differential abundance.
# The limma results were calculated in 01_setup.R.
# Volcano plots show the light and nitrogen contrasts.

# Set plot and significance thresholds.
VOLCANO_POINT_SIZE <- 1.15
VOLCANO_SINGLE_WIDTH <- 4.5
VOLCANO_PAIR_WIDTH <- 8.4
VOLCANO_HEIGHT <- 5.0
METAPROTEOME_SINGLE_WIDTH <- 4.0
METAPROTEOME_PAIR_WIDTH <- 7.4
METAPROTEOME_HEIGHT <- 6.2
VOLCANO_DPI <- 300
VOLCANO_FDR <- 0.05
VOLCANO_LOG2FC <- 1

volcano_theme_compact <- theme_bw(base_size=11)+theme(
  panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
  panel.border=element_rect(colour="grey75",fill=NA,linewidth=0.55),
  axis.text=element_text(size=12,colour="black"),
  axis.title=element_text(size=11.5,colour="black"),
  axis.title.x=element_text(size=10.5,colour="black",lineheight=0.95,margin=margin(t=5)),
  plot.title=element_text(size=12.5,face="bold",hjust=0.5,margin=margin(b=4)),
  legend.position="right",legend.direction="vertical",
  legend.text=element_text(size=8.5),legend.title=element_text(size=9,face="bold"),
  legend.key.size=unit(0.34,"cm"),
  plot.margin=margin(1,1,1,1,"pt"))

ps_16S_da <- ps_treated_16S
ps_16S_T_da <- ps_transplant_16S
cat(sprintf("  16S retained ASVs: factorial=%d | transplant=%d\n",ntaxa(ps_16S_da),ntaxa(ps_16S_T_da)))

cat("\nProtein differential abundance\n")

# Attach each protein group to its lichen partner and mark those
# passing FDR < 0.05 and absolute log2 fold change >= 1.
tt_light_vol <- res_light_prot%>%left_join(meta_nsaf_filt%>%select(members_identifier,symbionts)%>%distinct(),by="members_identifier")%>%mutate(Partner=factor(assign_partner(symbionts),levels=partner_levels),neg_logp=-log10(pmax(adj.P.Val,1e-300)),sig=!is.na(adj.P.Val)&adj.P.Val<VOLCANO_FDR&abs(logFC)>=VOLCANO_LOG2FC)
tt_nitrogen_vol <- res_nitrogen_prot%>%left_join(meta_nsaf_filt%>%select(members_identifier,symbionts)%>%distinct(),by="members_identifier")%>%mutate(Partner=factor(assign_partner(symbionts),levels=partner_levels),neg_logp=-log10(pmax(adj.P.Val,1e-300)),sig=!is.na(adj.P.Val)&adj.P.Val<VOLCANO_FDR&abs(logFC)>=VOLCANO_LOG2FC)
tt_transplant_vol <- res_transplant_prot%>%left_join(meta_nsaf_filt%>%select(members_identifier,symbionts)%>%distinct(),by="members_identifier")%>%mutate(Partner=factor(assign_partner(symbionts),levels=partner_levels),neg_logp=-log10(pmax(adj.P.Val,1e-300)),sig=!is.na(adj.P.Val)&adj.P.Val<VOLCANO_FDR&abs(logFC)>=VOLCANO_LOG2FC)
tt_interact_vol <- res_interact_prot%>%left_join(meta_nsaf_filt%>%select(members_identifier,symbionts)%>%distinct(),by="members_identifier")%>%mutate(Partner=factor(assign_partner(symbionts),levels=partner_levels),neg_logp=-log10(pmax(adj.P.Val,1e-300)),sig=!is.na(adj.P.Val)&adj.P.Val<VOLCANO_FDR&abs(logFC)>=VOLCANO_LOG2FC)

# Draw light, nitrogen, RAT/DTP, and interaction plots.
p_vol_light <- ggplot()+geom_point(data=tt_light_vol[!tt_light_vol$sig,],aes(logFC,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+geom_point(data=tt_light_vol[tt_light_vol$sig,],aes(logFC,neg_logp,colour=Partner),size=VOLCANO_POINT_SIZE,alpha=0.6)+scale_colour_manual(values=partner_colors,na.value="grey35",name="Partner")+geom_vline(xintercept=c(-1,1),linetype="dashed",colour="grey45")+geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey45")+annotate("text",x=-Inf,y=Inf,label="↑ LLS",hjust=-.15,vjust=1.25,colour=col_LLS,fontface="bold",size=4)+annotate("text",x=Inf,y=Inf,label="↑ HLS",hjust=1.15,vjust=1.25,colour=col_HLS,fontface="bold",size=4)+labs(title="Metaproteome — Light",x=expression(atop(Log[2]~Fold~Change,"(HLS / LLS)")),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact
p_vol_nitrogen <- ggplot()+geom_point(data=tt_nitrogen_vol[!tt_nitrogen_vol$sig,],aes(logFC,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+geom_point(data=tt_nitrogen_vol[tt_nitrogen_vol$sig,],aes(logFC,neg_logp,colour=Partner),size=VOLCANO_POINT_SIZE,alpha=0.6)+scale_colour_manual(values=partner_colors,na.value="grey35",name="Partner")+geom_vline(xintercept=c(-1,1),linetype="dashed",colour="grey45")+geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey45")+annotate("text",x=-Inf,y=Inf,label="↑ N0.32",hjust=-.15,vjust=1.25,colour=col_N0.32,fontface="bold",size=4)+annotate("text",x=Inf,y=Inf,label="↑ N1.28",hjust=1.15,vjust=1.25,colour=col_N1.28,fontface="bold",size=4)+labs(title="Metaproteome — Nitrogen",x=expression(atop(Log[2]~Fold~Change,"(N1.28 / N0.32)")),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact
p_vol_light_comb <- p_vol_light + labs(tag="A") +
  theme(plot.tag=element_text(size=21,face="bold",margin=margin(0,0,1,1,"pt")),plot.tag.position=c(0.012,0.94),plot.margin=margin(1.5,0,0.5,0.5,"pt"))
p_vol_nitrogen_comb <- p_vol_nitrogen + labs(tag="B") +
  theme(axis.title.y=element_blank(),plot.tag=element_text(size=21,face="bold",margin=margin(0,0,1,1,"pt")),plot.tag.position=c(0.012,0.94),plot.margin=margin(1.5,0.5,0.5,0,"pt"))
p_vol_light_nitrogen <- p_vol_light_comb+p_vol_nitrogen_comb+plot_layout(nrow=1,guides="collect")&
  theme(legend.position="right",legend.box.spacing=unit(0,"pt"),legend.margin=margin(0),legend.box.margin=margin(0),plot.margin=margin(0))
p_vol_transplant <- ggplot()+geom_point(data=tt_transplant_vol[!tt_transplant_vol$sig,],aes(logFC,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+geom_point(data=tt_transplant_vol[tt_transplant_vol$sig,],aes(logFC,neg_logp,colour=Partner),size=VOLCANO_POINT_SIZE,alpha=0.6)+scale_colour_manual(values=partner_colors,na.value="grey35",name="Partner")+geom_vline(xintercept=c(-1,1),linetype="dashed",colour="grey45")+geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey45")+annotate("text",x=-Inf,y=Inf,label="↑ RAT",hjust=-.15,vjust=1.25,colour=col_RAT,fontface="bold",size=4)+annotate("text",x=Inf,y=Inf,label="↑ DTP",hjust=1.15,vjust=1.25,colour=col_DTP,fontface="bold",size=4)+labs(title="Metaproteome — RAT vs DTP",x=expression(Log[2]~Fold~Change~"(DTP / RAT)"),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact
p_vol_interaction <- ggplot()+geom_point(data=tt_interact_vol[!tt_interact_vol$sig,],aes(logFC,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+geom_point(data=tt_interact_vol[tt_interact_vol$sig,],aes(logFC,neg_logp,colour=Partner),size=VOLCANO_POINT_SIZE,alpha=0.6)+scale_colour_manual(values=partner_colors,na.value="grey35",name="Partner")+geom_vline(xintercept=c(-1,1),linetype="dashed",colour="grey45")+geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey45")+annotate("text",x=-Inf,y=Inf,label="Negative interaction",hjust=-.05,vjust=1.25,fontface="bold",size=3.5)+annotate("text",x=Inf,y=Inf,label="Positive interaction",hjust=1.05,vjust=1.25,fontface="bold",size=3.5)+labs(title="Metaproteome — Light × Nitrogen interaction",x="Interaction log2 fold change",y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact

# Save the two main treatment panels separately.
fig7_dir <- file.path(INDIVIDUAL_DIR,"Figure7")
save_svg(file.path(fig7_dir,"Figure7A_metaproteome_Light"), p_vol_light, width=3.7, height=3.40)
save_svg(file.path(fig7_dir,"Figure7B_metaproteome_Nitrogen"), p_vol_nitrogen, width=3.7, height=3.40)

global_volcano_all <- bind_rows(mutate(tt_light_vol,Contrast="Light"),mutate(tt_nitrogen_vol,Contrast="Nitrogen"),mutate(tt_transplant_vol,Contrast="RAT vs DTP"),mutate(tt_interact_vol,Contrast="Light × Nitrogen"))
global_vol_summary <- data.frame(
  Contrast=c("Light: LLS vs HLS","Nitrogen: N0.32 vs N1.28","Transplant: RAT vs DTP","Light × Nitrogen interaction"),
  n_total=c(nrow(tt_light_vol),nrow(tt_nitrogen_vol),nrow(tt_transplant_vol),nrow(tt_interact_vol)),
  n_sig=c(sum(tt_light_vol$sig,na.rm=TRUE),sum(tt_nitrogen_vol$sig,na.rm=TRUE),sum(tt_transplant_vol$sig,na.rm=TRUE),sum(tt_interact_vol$sig,na.rm=TRUE)),
  n_up=c(sum(tt_light_vol$sig&tt_light_vol$logFC>0,na.rm=TRUE),sum(tt_nitrogen_vol$sig&tt_nitrogen_vol$logFC>0,na.rm=TRUE),sum(tt_transplant_vol$sig&tt_transplant_vol$logFC>0,na.rm=TRUE),sum(tt_interact_vol$sig&tt_interact_vol$logFC>0,na.rm=TRUE)),
  n_down=c(sum(tt_light_vol$sig&tt_light_vol$logFC<0,na.rm=TRUE),sum(tt_nitrogen_vol$sig&tt_nitrogen_vol$logFC<0,na.rm=TRUE),sum(tt_transplant_vol$sig&tt_transplant_vol$logFC<0,na.rm=TRUE),sum(tt_interact_vol$sig&tt_interact_vol$logFC<0,na.rm=TRUE)),
  stringsAsFactors=FALSE)
global_vol_summary$pct_sig <- round(100*global_vol_summary$n_sig/global_vol_summary$n_total,2)
# Export all protein results and a count of significant changes.
write.csv(global_volcano_all,file.path(OUTDIR,"protein_group_analysis/tables/nsaf_global_volcano_all.csv"),row.names=FALSE)
write.csv(global_volcano_all%>%filter(sig),file.path(OUTDIR,"protein_group_analysis/tables/nsaf_global_volcano_significant.csv"),row.names=FALSE)
write.csv(global_vol_summary,file.path(OUTDIR,"protein_group_analysis/tables/nsaf_global_volcano_summary.csv"),row.names=FALSE)
cat("Section 11 complete.\n\n")

save_svg(file.path(OUTDIR,"Figure7_metaproteome_volcano"), p_vol_light_nitrogen, width=6.75, height=3.40)

# Global metaproteome panels; no ten-functional-group comparisons are changed.
vmeta_light <- volcano_strips(p_vol_light, "LLS", "HLS", col_LLS, col_HLS)
vmeta_nitrogen <- volcano_strips(p_vol_nitrogen, "N0.32", "N1.28", col_N0.32, col_N1.28)
vmeta_transplant <- volcano_strips(p_vol_transplant, "RAT", "DTP", col_RAT, col_DTP)
vmeta_interaction <- volcano_strips(p_vol_interaction, "Negative interaction", "Positive interaction",
                                   "grey80", "grey65", arrow = FALSE)

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

# Save individual metaproteome volcano panels.
protein_volcano_panels <- list(
  "Metaproteome_Light" = vmeta_light, "Metaproteome_Nitrogen" = vmeta_nitrogen,
  "Metaproteome_RAT_DTP" = vmeta_transplant, "Metaproteome_Interaction" = vmeta_interaction)
for (name in names(protein_volcano_panels)) {
  save_volcano(file.path(OUTDIR, "volcano_panels", name), protein_volcano_panels[[name]], 5.5, 3.5)
}

# Keep the existing Figure 7 individual-panel filenames up to date.
save_volcano(file.path(INDIVIDUAL_DIR, "Figure7", "Figure7A_metaproteome_Light"), vmeta_light, 4.2, 3.1)
save_volcano(file.path(INDIVIDUAL_DIR, "Figure7", "Figure7B_metaproteome_Nitrogen"), vmeta_nitrogen, 4.2, 3.1)
