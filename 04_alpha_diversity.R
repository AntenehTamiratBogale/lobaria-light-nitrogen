# Figure 5: within-specimen diversity of 16S and 18S ASVs.
# Observed is ASV richness. The other index includes relative abundances.
# The code calls that index Shannon so its calculation remains explicit.

cat("\nAlpha diversity\n")
# Rarefy each dataset to its smallest specimen read count.
set.seed(1); phyloseq_rarefied_16S<-rarefy_even_depth(phyloseq_16S_filt,sample.size=min(sample_sums(phyloseq_16S_filt)),rngseed=1,replace=TRUE,trimOTUs=TRUE,verbose=FALSE)
set.seed(1); phyloseq_rarefied_18S<-rarefy_even_depth(phyloseq_18S_filt,sample.size=min(sample_sums(phyloseq_18S_filt)),rngseed=1,replace=TRUE,trimOTUs=TRUE,verbose=FALSE)
cat(sprintf("  16S rarefaction depth: %d reads\n",min(sample_sums(phyloseq_16S_filt))))
cat(sprintf("  18S rarefaction depth: %d reads\n\n",min(sample_sums(phyloseq_18S_filt))))
otu_16S<-as.matrix(otu_table(phyloseq_rarefied_16S)); if(taxa_are_rows(phyloseq_rarefied_16S))otu_16S<-t(otu_16S)
otu_18S<-as.matrix(otu_table(phyloseq_rarefied_18S)); if(taxa_are_rows(phyloseq_rarefied_18S))otu_18S<-t(otu_18S)
sample_data_16S<-as.data.frame(sample_data(phyloseq_rarefied_16S)); sample_data_18S<-as.data.frame(sample_data(phyloseq_rarefied_18S))
# Calculate observed richness and the alpha-diversity index.
alpha_16S<-data.frame(SampleID=rownames(sample_data_16S),Observed=rowSums(otu_16S>0),Shannon=vegan::diversity(otu_16S,index="shannon"),Light=sample_data_16S$Light,Nitrogen=sample_data_16S$Nitrogen,Transplant=sample_data_16S$Transplant,stringsAsFactors=FALSE)
alpha_18S<-data.frame(SampleID=rownames(sample_data_18S),Observed=rowSums(otu_18S>0),Shannon=vegan::diversity(otu_18S,index="shannon"),Light=sample_data_18S$Light,Nitrogen=sample_data_18S$Nitrogen,Transplant=sample_data_18S$Transplant,stringsAsFactors=FALSE)
# Print mean and SD for each light, nitrogen, and RAT/DTP group.
cat("  16S alpha diversity (mean +/- SD per group):\n")
for(grp in c("LLS","HLS","N0.32","N1.28","RAT","DTP")){col<-if(grp%in%c("LLS","HLS"))"Light" else if(grp%in%c("N0.32","N1.28"))"Nitrogen" else "Transplant";sub<-alpha_16S[!is.na(alpha_16S[[col]])&alpha_16S[[col]]==grp,];cat(sprintf("    %s  n=%d  Observed=%.1f+/-%.1f  Shannon=%.3f+/-%.3f\n",grp,nrow(sub),mean(sub$Observed,na.rm=TRUE),sd(sub$Observed,na.rm=TRUE),mean(sub$Shannon,na.rm=TRUE),sd(sub$Shannon,na.rm=TRUE)))}
cat("  18S alpha diversity (mean +/- SD per group):\n")
for(grp in c("LLS","HLS","N0.32","N1.28","RAT","DTP")){col<-if(grp%in%c("LLS","HLS"))"Light" else if(grp%in%c("N0.32","N1.28"))"Nitrogen" else "Transplant";sub<-alpha_18S[!is.na(alpha_18S[[col]])&alpha_18S[[col]]==grp,];cat(sprintf("    %s  n=%d  Observed=%.1f+/-%.1f  Shannon=%.3f+/-%.3f\n",grp,nrow(sub),mean(sub$Observed,na.rm=TRUE),sd(sub$Observed,na.rm=TRUE),mean(sub$Shannon,na.rm=TRUE),sd(sub$Shannon,na.rm=TRUE)))}
cat("\n")
# Gather observed richness and the index for each treatment group.
get_alpha<-function(otu_mat,ids,index="shannon"){if(index=="observed")rowSums(otu_mat[ids,,drop=FALSE]>0) else vegan::diversity(otu_mat[ids,,drop=FALSE],index=index)}
.run_alpha<-function(ps_rar,sd_var_list){otu_m<-as.matrix(otu_table(ps_rar));if(taxa_are_rows(ps_rar))otu_m<-t(otu_m);sd_d<-as.data.frame(sample_data(ps_rar));res<-list();for(grp in sd_var_list){ids<-rownames(sd_d)[!is.na(sd_d[[grp$var]])&sd_d[[grp$var]]==grp$lbl];res[[paste0("obs_",grp$lbl)]]<-get_alpha(otu_m,ids,"observed");res[[paste0("sha_",grp$lbl)]]<-get_alpha(otu_m,ids,"shannon")};res}
res16<-.run_alpha(phyloseq_rarefied_16S,list(list(var="Light",lbl="LLS"),list(var="Light",lbl="HLS"),list(var="Nitrogen",lbl="N0.32"),list(var="Nitrogen",lbl="N1.28"),list(var="Transplant",lbl="RAT"),list(var="Transplant",lbl="DTP")))
res18<-.run_alpha(phyloseq_rarefied_18S,list(list(var="Light",lbl="LLS"),list(var="Light",lbl="HLS"),list(var="Nitrogen",lbl="N0.32"),list(var="Nitrogen",lbl="N1.28"),list(var="Transplant",lbl="RAT"),list(var="Transplant",lbl="DTP")))
# Compare RAT and DTP specimens separately.
wt_T_obs_16S<-wilcox.test(res16$obs_RAT,res16$obs_DTP,exact=FALSE)$p.value; wt_T_sha_16S<-wilcox.test(res16$sha_RAT,res16$sha_DTP,exact=FALSE)$p.value
wt_T_obs_18S<-wilcox.test(res18$obs_RAT,res18$obs_DTP,exact=FALSE)$p.value; wt_T_sha_18S<-wilcox.test(res18$sha_RAT,res18$sha_DTP,exact=FALSE)$p.value
cat("  --- Separate RAT–DTP Wilcoxon comparisons ---\n")
cat(sprintf("  16S Observed: P=%.4f | alpha diversity: P=%.4f\n",wt_T_obs_16S,wt_T_sha_16S))
cat(sprintf("  18S Observed: P=%.4f | alpha diversity: P=%.4f\n\n",wt_T_obs_18S,wt_T_sha_18S))
# Fit light × nitrogen models using only the twelve factorial specimens.
alpha_16S_LN<-alpha_16S[!is.na(alpha_16S$Light)&!is.na(alpha_16S$Nitrogen),]
alpha_18S_LN<-alpha_18S[!is.na(alpha_18S$Light)&!is.na(alpha_18S$Nitrogen),]
alpha_16S_LN$Light<-factor(alpha_16S_LN$Light,levels=c("LLS","HLS")); alpha_16S_LN$Nitrogen<-factor(alpha_16S_LN$Nitrogen,levels=c("N0.32","N1.28"))
alpha_18S_LN$Light<-factor(alpha_18S_LN$Light,levels=c("LLS","HLS")); alpha_18S_LN$Nitrogen<-factor(alpha_18S_LN$Nitrogen,levels=c("N0.32","N1.28"))
lm_16S_obs<-lm(Observed~Light*Nitrogen,data=alpha_16S_LN); lm_16S_sha<-lm(Shannon~Light*Nitrogen,data=alpha_16S_LN)
lm_18S_obs<-lm(Observed~Light*Nitrogen,data=alpha_18S_LN); lm_18S_sha<-lm(Shannon~Light*Nitrogen,data=alpha_18S_LN)
# Print the ANOVA tables for richness and the abundance-based index.
aov_16S_obs<-anova(lm_16S_obs); aov_16S_sha<-anova(lm_16S_sha)
aov_18S_obs<-anova(lm_18S_obs); aov_18S_sha<-anova(lm_18S_sha)

lm_L_obs_16S<-aov_16S_obs["Light","Pr(>F)"]; lm_N_obs_16S<-aov_16S_obs["Nitrogen","Pr(>F)"]
lm_L_sha_16S<-aov_16S_sha["Light","Pr(>F)"]; lm_N_sha_16S<-aov_16S_sha["Nitrogen","Pr(>F)"]
lm_L_obs_18S<-aov_18S_obs["Light","Pr(>F)"]; lm_N_obs_18S<-aov_18S_obs["Nitrogen","Pr(>F)"]
lm_L_sha_18S<-aov_18S_sha["Light","Pr(>F)"]; lm_N_sha_18S<-aov_18S_sha["Nitrogen","Pr(>F)"]
cat("  --- Linear model: Observed ~ Light * Nitrogen ---\n  16S:\n"); print(aov_16S_obs); cat("  18S:\n"); print(aov_18S_obs)
cat("\n  --- Linear model: Shannon ~ Light * Nitrogen ---\n  16S:\n"); print(aov_16S_sha); cat("  18S:\n"); print(aov_18S_sha)
p_int_16S_obs<-aov_16S_obs["Light:Nitrogen","Pr(>F)"]; p_int_16S_sha<-aov_16S_sha["Light:Nitrogen","Pr(>F)"]
p_int_18S_obs<-aov_18S_obs["Light:Nitrogen","Pr(>F)"]; p_int_18S_sha<-aov_18S_sha["Light:Nitrogen","Pr(>F)"]

# Apply BH correction to the four light × nitrogen interaction tests.
int_p_raw <- c(p_int_16S_obs, p_int_16S_sha, p_int_18S_obs, p_int_18S_sha)
int_p_adj <- p.adjust(int_p_raw, method = "BH")
p_int_16S_obs_adj <- int_p_adj[1]
p_int_16S_sha_adj <- int_p_adj[2]
p_int_18S_obs_adj <- int_p_adj[3]
p_int_18S_sha_adj <- int_p_adj[4]

cat("\n  --- Interaction p-values (Light x Nitrogen) — BH corrected across 4 tests ---\n")
cat(sprintf("  16S Observed: raw p = %.4f  adj p = %.4f %s\n",
            p_int_16S_obs, p_int_16S_obs_adj, sig_fmt(p_int_16S_obs_adj)))
cat(sprintf("  16S Shannon:  raw p = %.4f  adj p = %.4f %s\n",
            p_int_16S_sha, p_int_16S_sha_adj, sig_fmt(p_int_16S_sha_adj)))
cat(sprintf("  18S Observed: raw p = %.4f  adj p = %.4f %s\n",
            p_int_18S_obs, p_int_18S_obs_adj, sig_fmt(p_int_18S_obs_adj)))
cat(sprintf("  18S Shannon:  raw p = %.4f  adj p = %.4f %s\n\n",
            p_int_18S_sha, p_int_18S_sha_adj, sig_fmt(p_int_18S_sha_adj)))

# Write the exact treatment P values to CSV files.
alpha_interact_tbl <- data.frame(
  Dataset          = rep(c("16S","18S"), each = 2),
  Metric           = rep(c("Observed","Shannon"), 2),
  LM_Light_p       = c(aov_16S_obs["Light","Pr(>F)"],    aov_16S_sha["Light","Pr(>F)"],
                       aov_18S_obs["Light","Pr(>F)"],    aov_18S_sha["Light","Pr(>F)"]),
  LM_Nitrogen_p    = c(aov_16S_obs["Nitrogen","Pr(>F)"], aov_16S_sha["Nitrogen","Pr(>F)"],
                       aov_18S_obs["Nitrogen","Pr(>F)"], aov_18S_sha["Nitrogen","Pr(>F)"]),
  LM_Interaction_p_raw = int_p_raw,
  LM_Interaction_p_adj = int_p_adj,
  stringsAsFactors = FALSE) %>%
  mutate(
    LM_Light_sig       = sig_fmt(LM_Light_p),
    LM_Nitrogen_sig    = sig_fmt(LM_Nitrogen_p),
    LM_Interaction_sig = sig_fmt(LM_Interaction_p_adj))

write.csv(alpha_interact_tbl,
          file.path(OUTDIR, "alpha_interaction_lm_results.csv"),
          row.names = FALSE)
cat("Saved: alpha_interaction_lm_results.csv\n")

alpha_all_stats <- data.frame(
  Dataset = rep(c("16S","18S"), each=2),
  Metric = rep(c("Observed","Shannon"), 2),
  Light_p = c(lm_L_obs_16S,lm_L_sha_16S,lm_L_obs_18S,lm_L_sha_18S),
  Nitrogen_p = c(lm_N_obs_16S,lm_N_sha_16S,lm_N_obs_18S,lm_N_sha_18S),
  Light_x_Nitrogen_p_raw = int_p_raw,
  Light_x_Nitrogen_p_BH = int_p_adj,
  RAT_DTP_p = c(wt_T_obs_16S,wt_T_sha_16S,wt_T_obs_18S,wt_T_sha_18S),
  stringsAsFactors = FALSE
)
write.csv(alpha_all_stats, file.path(OUTDIR,"alpha_all_effects_statistics.csv"), row.names=FALSE)
cat("Saved: alpha_all_effects_statistics.csv (Light, Nitrogen, interaction, RAT-DTP).\n")

any_interact_sig <- any(int_p_adj < 0.05, na.rm = TRUE)
cat("\n  Light x Nitrogen interaction results retained in tables only (no interaction figures exported).\n")
# Draw four boxplots with individual specimen values and significance.
.alpha_LN_plot<-function(obs_D,obs_B,obs_M,obs_S,sha_D,sha_B,sha_M,sha_S,wL_o,wN_o,wL_s,wN_s,pfx=""){
  dc<-rbind(data.frame(id=names(obs_D),Observed=obs_D,Shannon=sha_D,group="LLS"),data.frame(id=names(obs_B),Observed=obs_B,Shannon=sha_B,group="HLS"),data.frame(id=names(obs_M),Observed=obs_M,Shannon=sha_M,group="N0.32"),data.frame(id=names(obs_S),Observed=obs_S,Shannon=sha_S,group="N1.28"))
  dc$group<-factor(dc$group,levels=c("LLS","HLS","N0.32","N1.28")); om<-max(c(obs_D,obs_B,obs_M,obs_S))*1.1; sm<-max(c(sha_D,sha_B,sha_M,sha_S))*1.1
  .bp<-function(dc,y_col,y_max,y_lab,wL,wN,show_x=FALSE){ggplot(dc,aes(x=group,y=.data[[y_col]],fill=group,color=group))+geom_boxplot(alpha=0.6,width=0.5,outlier.shape=NA)+geom_jitter(width=0.12,size=3,alpha=0.8)+geom_vline(xintercept=2.5,linetype="dashed",color="gray40",linewidth=0.4)+annotate("segment",x=1,xend=2,y=y_max*0.95,yend=y_max*0.95,color="black",linewidth=0.6,lineend="round")+annotate("text",x=1.5,y=y_max,label=sig_fmt(wL),size=5,fontface="bold")+annotate("segment",x=3,xend=4,y=y_max*0.95,yend=y_max*0.95,color="black",linewidth=0.6,lineend="round")+annotate("text",x=3.5,y=y_max,label=sig_fmt(wN),size=5,fontface="bold")+scale_fill_manual(values=color_values)+scale_color_manual(values=color_values)+alpha_theme+labs(y=y_lab)+theme(plot.margin=margin(0,0,0,0,"pt"))+{if(show_x)theme(axis.text.x=element_text(size=12,colour="black",angle=30,hjust=0.5,vjust=1,margin=margin(t=2))) else theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())}}
  list(obs=.bp(dc,"Observed",om,"Observed diversity",wL_o,wN_o),sha=.bp(dc,"Shannon",sm,"Shannon diversity",wL_s,wN_s,show_x=TRUE))}
al16<-.alpha_LN_plot(res16$obs_LLS,res16$obs_HLS,res16$obs_N0.32,res16$obs_N1.28,res16$sha_LLS,res16$sha_HLS,res16$sha_N0.32,res16$sha_N1.28,lm_L_obs_16S,lm_N_obs_16S,lm_L_sha_16S,lm_N_sha_16S)
al18<-.alpha_LN_plot(res18$obs_LLS,res18$obs_HLS,res18$obs_N0.32,res18$obs_N1.28,res18$sha_LLS,res18$sha_HLS,res18$sha_N0.32,res18$sha_N1.28,lm_L_obs_18S,lm_N_obs_18S,lm_L_sha_18S,lm_N_sha_18S)

cat("Section 9 complete.\n\n")

# Save separate panels and the combined Figure 5 SVG.
fig5_dir <- file.path(INDIVIDUAL_DIR,"Figure5")
save_svg(file.path(fig5_dir,"Figure5A_16S_observed_diversity"), al16$obs, width=3.0, height=3.0)
save_svg(file.path(fig5_dir,"Figure5B_18S_observed_diversity"), al18$obs, width=3.0, height=3.0)
save_svg(file.path(fig5_dir,"Figure5C_16S_Shannon_diversity"), al16$sha, width=3.0, height=3.0)
save_svg(file.path(fig5_dir,"Figure5D_18S_Shannon_diversity"), al18$sha, width=3.0, height=3.0)

.tag_alpha <- function(p, tag) p + labs(tag=tag) +
  theme(plot.tag=element_text(size=22,face="bold",hjust=0,vjust=1),
        plot.tag.position=c(0.012,0.985),
        plot.margin=margin(0.5,2.5,0.5,2.5,"pt"))
# Combine and label the four panels.
obs_left  <- .tag_alpha(al16$obs + theme(axis.title.y=element_blank()), "A")
obs_right <- .tag_alpha(al18$obs + theme(axis.title.y=element_blank()), "B") + theme(plot.margin=margin(0.5,4.5,0.5,2.5,"pt"))
sha_left  <- .tag_alpha(al16$sha + theme(axis.title.y=element_blank()), "C")
sha_right <- .tag_alpha(al18$sha + theme(axis.title.y=element_blank()), "D") + theme(plot.margin=margin(0.5,4.5,0.5,2.5,"pt"))
obs_row <- patchwork::wrap_plots(obs_left, obs_right, ncol=2)
sha_row <- patchwork::wrap_plots(sha_left, sha_right, ncol=2)
obs_row <- add_common_y_title(obs_row, "Observed diversity", text_size=14, width=0.048)
sha_row <- add_common_y_title(sha_row, "Shannon diversity", text_size=14, width=0.048)
Figure5 <- patchwork::wrap_plots(obs_row, sha_row, ncol=1, heights=c(1,1.14)) +
  patchwork::plot_annotation(theme=theme(plot.margin=margin(0)))
save_svg(file.path(OUTDIR,"Figure5_alpha_diversity"), Figure5, width=5.55, height=4.98)
