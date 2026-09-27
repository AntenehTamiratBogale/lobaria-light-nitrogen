# Figure 6: differential abundance of bacterial ASVs.
# DESeq2 compares HLS with LLS and N1.28 with N0.32.
# RAT/DTP and light × nitrogen are also saved as separate results.

cat("\nASV differential abundance\n")

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
  plot.title=element_text(size=12.5,face="bold",hjust=0.5,margin=margin(b=4)),
  legend.position="right",legend.direction="vertical",
  legend.text=element_text(size=8.5),legend.title=element_text(size=9,face="bold"),
  legend.key.size=unit(0.34,"cm"),
  legend.box.spacing=unit(0,"pt"),legend.margin=margin(0),legend.box.margin=margin(0),
  plot.margin=margin(0,0,0,0,"pt"))

# Use the unrarefied 16S counts for DESeq2.
ps_16S_da <- ps_treated_16S
ps_16S_T_da <- ps_transplant_16S
cat(sprintf("  16S retained ASVs: factorial=%d | transplant=%d\n",ntaxa(ps_16S_da),ntaxa(ps_16S_T_da)))

# Fit the main-effects model and mark ASVs passing FDR and fold change.
dds_16S_main <- DESeq(phyloseq_to_deseq2(ps_16S_da,~Nitrogen+Light),
                      test="Wald",fitType="parametric")

res_16S_L <- as.data.frame(results(dds_16S_main,name="Light_HLS_vs_LLS",alpha=VOLCANO_FDR))
res_16S_L$ASV <- rownames(res_16S_L)
res_16S_L$Class <- tax_df_16S[res_16S_L$ASV,"Class"]
res_16S_L$neg_logp <- -log10(pmax(res_16S_L$padj,1e-300))
res_16S_L$significant <- !is.na(res_16S_L$padj)&res_16S_L$padj<VOLCANO_FDR&abs(res_16S_L$log2FoldChange)>=VOLCANO_LOG2FC

res_16S_N <- as.data.frame(results(dds_16S_main,name="Nitrogen_N1.28_vs_N0.32",alpha=VOLCANO_FDR))
res_16S_N$ASV <- rownames(res_16S_N)
res_16S_N$Class <- tax_df_16S[res_16S_N$ASV,"Class"]
res_16S_N$neg_logp <- -log10(pmax(res_16S_N$padj,1e-300))
res_16S_N$significant <- !is.na(res_16S_N$padj)&res_16S_N$padj<VOLCANO_FDR&abs(res_16S_N$log2FoldChange)>=VOLCANO_LOG2FC

dds_16S_T <- DESeq(phyloseq_to_deseq2(ps_16S_T_da,~Transplant),
                   test="Wald",fitType="parametric")
res_16S_T <- as.data.frame(results(dds_16S_T,name="Transplant_DTP_vs_RAT",alpha=VOLCANO_FDR))
res_16S_T$ASV <- rownames(res_16S_T)
res_16S_T$Class <- tax_df_16S[res_16S_T$ASV,"Class"]
res_16S_T$neg_logp <- -log10(pmax(res_16S_T$padj,1e-300))
res_16S_T$significant <- !is.na(res_16S_T$padj)&res_16S_T$padj<VOLCANO_FDR&abs(res_16S_T$log2FoldChange)>=VOLCANO_LOG2FC

# Test the light-by-nitrogen interaction.
dds_16S_LN <- DESeq(phyloseq_to_deseq2(ps_16S_da,~Nitrogen*Light),
                    test="Wald",fitType="parametric")
int_name_16S <- grep("Nitrogen.*Light|Light.*Nitrogen",resultsNames(dds_16S_LN),value=TRUE)
if(length(int_name_16S)!=1)stop("Cannot uniquely identify the 16S interaction coefficient: ",paste(resultsNames(dds_16S_LN),collapse=", "))
res_16S_LN <- as.data.frame(results(dds_16S_LN,name=int_name_16S,alpha=VOLCANO_FDR))
res_16S_LN$ASV <- rownames(res_16S_LN)
res_16S_LN$Class <- tax_df_16S[res_16S_LN$ASV,"Class"]
res_16S_LN$neg_logp <- -log10(pmax(res_16S_LN$padj,1e-300))
res_16S_LN$significant <- !is.na(res_16S_LN$padj)&res_16S_LN$padj<VOLCANO_FDR&abs(res_16S_LN$log2FoldChange)>=VOLCANO_LOG2FC

# Keep exact DESeq2 results for each contrast in CSV files.
write.csv(res_16S_L,file.path(OUTDIR,"deseq2_16S_light_adjusted.csv"),row.names=FALSE)
write.csv(res_16S_N,file.path(OUTDIR,"deseq2_16S_nitrogen_adjusted.csv"),row.names=FALSE)
write.csv(res_16S_T,file.path(OUTDIR,"deseq2_16S_transplant.csv"),row.names=FALSE)
write.csv(res_16S_LN,file.path(OUTDIR,"deseq2_16S_interaction.csv"),row.names=FALSE)

# Colour significant ASVs by their assigned bacterial class.
.class_scale_16S <- function(df) {
  lev <- unique(as.character(df$Class[df$significant & !is.na(df$Class)]))
  lev <- lev[nzchar(lev)]
  if (!length(lev)) return(NULL)
  pal <- vol_colors_16S
  miss <- setdiff(lev, names(pal))
  if (length(miss)) pal <- c(pal, stats::setNames(rep("grey35", length(miss)), miss))
  scale_colour_manual(values=pal, na.value="grey35", name="Class")
}

p_16L <- ggplot()+
  geom_point(data=res_16S_L[!res_16S_L$significant&is.finite(res_16S_L$log2FoldChange)&is.finite(res_16S_L$neg_logp),],aes(log2FoldChange,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+
  geom_point(data=res_16S_L[res_16S_L$significant&is.finite(res_16S_L$log2FoldChange)&is.finite(res_16S_L$neg_logp),],aes(log2FoldChange,neg_logp,colour=Class),size=VOLCANO_POINT_SIZE,alpha=0.6)+
  .class_scale_16S(res_16S_L)+
  geom_vline(xintercept=c(-VOLCANO_LOG2FC,VOLCANO_LOG2FC),linetype="dashed",colour="grey45")+
  geom_hline(yintercept=-log10(VOLCANO_FDR),linetype="dashed",colour="grey45")+
  annotate("text",x=-Inf,y=Inf,label="↑ LLS",hjust=-0.15,vjust=1.25,colour=col_LLS,fontface="bold",size=4)+
  annotate("text",x=Inf,y=Inf,label="↑ HLS",hjust=1.15,vjust=1.25,colour=col_HLS,fontface="bold",size=4)+
  labs(title="16S — Light",x=expression(atop(Log[2]~Fold~Change,"(HLS / LLS)")),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact

p_16N <- ggplot()+
  geom_point(data=res_16S_N[!res_16S_N$significant&is.finite(res_16S_N$log2FoldChange)&is.finite(res_16S_N$neg_logp),],aes(log2FoldChange,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+
  geom_point(data=res_16S_N[res_16S_N$significant&is.finite(res_16S_N$log2FoldChange)&is.finite(res_16S_N$neg_logp),],aes(log2FoldChange,neg_logp,colour=Class),size=VOLCANO_POINT_SIZE,alpha=0.6)+
  .class_scale_16S(res_16S_N)+
  geom_vline(xintercept=c(-VOLCANO_LOG2FC,VOLCANO_LOG2FC),linetype="dashed",colour="grey45")+
  geom_hline(yintercept=-log10(VOLCANO_FDR),linetype="dashed",colour="grey45")+
  annotate("text",x=-Inf,y=Inf,label="↑ N0.32",hjust=-0.15,vjust=1.25,colour=col_N0.32,fontface="bold",size=4)+
  annotate("text",x=Inf,y=Inf,label="↑ N1.28",hjust=1.15,vjust=1.25,colour=col_N1.28,fontface="bold",size=4)+
  labs(title="16S — Nitrogen",x=expression(Log[2]~Fold~Change~"(N1.28 / N0.32)"),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact

p_16_light_nitrogen <- p_16L+p_16N+plot_layout(nrow=1,guides="collect")&theme(legend.position="right")

p_16T <- ggplot()+
  geom_point(data=res_16S_T[!res_16S_T$significant&is.finite(res_16S_T$log2FoldChange)&is.finite(res_16S_T$neg_logp),],aes(log2FoldChange,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+
  geom_point(data=res_16S_T[res_16S_T$significant&is.finite(res_16S_T$log2FoldChange)&is.finite(res_16S_T$neg_logp),],aes(log2FoldChange,neg_logp,colour=Class),size=VOLCANO_POINT_SIZE,alpha=0.6)+
  .class_scale_16S(res_16S_T)+
  geom_vline(xintercept=c(-VOLCANO_LOG2FC,VOLCANO_LOG2FC),linetype="dashed",colour="grey45")+
  geom_hline(yintercept=-log10(VOLCANO_FDR),linetype="dashed",colour="grey45")+
  annotate("text",x=-Inf,y=Inf,label="↑ RAT",hjust=-0.15,vjust=1.25,colour=col_RAT,fontface="bold",size=4)+
  annotate("text",x=Inf,y=Inf,label="↑ DTP",hjust=1.15,vjust=1.25,colour=col_DTP,fontface="bold",size=4)+
  labs(title="16S — RAT vs DTP",x=expression(Log[2]~Fold~Change~"(DTP / RAT)"),y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact

p_16LN <- ggplot()+
  geom_point(data=res_16S_LN[!res_16S_LN$significant&is.finite(res_16S_LN$log2FoldChange)&is.finite(res_16S_LN$neg_logp),],aes(log2FoldChange,neg_logp),colour="grey55",size=VOLCANO_POINT_SIZE,alpha=0.6)+
  geom_point(data=res_16S_LN[res_16S_LN$significant&is.finite(res_16S_LN$log2FoldChange)&is.finite(res_16S_LN$neg_logp),],aes(log2FoldChange,neg_logp,colour=Class),size=VOLCANO_POINT_SIZE,alpha=0.6)+
  .class_scale_16S(res_16S_LN)+
  geom_vline(xintercept=c(-VOLCANO_LOG2FC,VOLCANO_LOG2FC),linetype="dashed",colour="grey45")+
  geom_hline(yintercept=-log10(VOLCANO_FDR),linetype="dashed",colour="grey45")+
  annotate("text",x=-Inf,y=Inf,label="Negative interaction",hjust=-0.05,vjust=1.25,fontface="bold",size=3.5)+
  annotate("text",x=Inf,y=Inf,label="Positive interaction",hjust=1.05,vjust=1.25,fontface="bold",size=3.5)+
  labs(title="16S — Light × Nitrogen interaction",x="Interaction log2 fold change",y=expression(-log[10](adjusted~italic(p))))+volcano_theme_compact

save_svg(file.path(OUTDIR,"Figure6_16S_light_volcano"), p_16L, width=4.00, height=2.80)
