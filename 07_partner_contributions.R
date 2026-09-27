# Figure 8: contributions of mycobiont, phycobiont, cyanobiont,
# bacteria, and other organisms to proteins and NSAF abundance.

cat("\nPartner contributions\n")
# Use the manuscript term Phycobiont for entries named Photobiont.
sym5<-function(sv){sv<-trimws(as.character(sv));dplyr::case_when(sv=="Photobiont"~"Phycobiont",sv%in%target_symbionts~sv,TRUE~"Other")}

# Count protein groups and sum NSAF values for each partner.
pg_sym<-prot_filt%>%dplyr::mutate(symbiont5=sym5(symbionts))%>%dplyr::group_by(symbiont5)%>%dplyr::summarise(n=dplyr::n(),.groups="drop")%>%dplyr::mutate(pct=n/sum(n)*100,metric="Protein groups (%)")

filt_ids<-prot_filt$members_identifier; mat_ov_raw<-mat_nsaf_raw[rownames(mat_nsaf_raw)%in%filt_ids,,drop=FALSE]
nsaf_row_sums<-rowSums(mat_ov_raw,na.rm=TRUE)
nsaf_sym<-prot_filt%>%dplyr::mutate(symbiont5=sym5(symbionts),nsaf_total=nsaf_row_sums[match(members_identifier,names(nsaf_row_sums))],nsaf_total=tidyr::replace_na(nsaf_total,0))%>%dplyr::group_by(symbiont5)%>%dplyr::summarise(total_nsaf=sum(nsaf_total,na.rm=TRUE),.groups="drop")%>%dplyr::mutate(pct=total_nsaf/sum(total_nsaf)*100,metric="NSAF abundance (%)")

sym_stack_order<-nsaf_sym%>%arrange(pct)%>%pull(symbiont5)%>%as.character()

.expected_syms21<-c("Mycobiont","Phycobiont","Cyanobiont","Bacteria","Other")
if(!setequal(sym_stack_order,.expected_syms21)){
  cat("  WARNING: sym_stack_order does not match the 5 expected categories.\n")
  cat("           expected:",paste(.expected_syms21,collapse=", "),"\n")
  cat("           got:     ",paste(sym_stack_order,collapse=", "),"\n")
  missing_syms21<-setdiff(.expected_syms21,sym_stack_order)
  cat("           missing: ",paste(missing_syms21,collapse=", "),"- appending to preserve them.\n\n")
  sym_stack_order<-c(sym_stack_order,missing_syms21)
}

pg_sym2<-pg_sym%>%dplyr::mutate(symbiont5=factor(symbiont5,levels=sym_stack_order),value=pct,label=ifelse(pct>2,paste0(format(n,big.mark=",")," (",round(pct,1),"%)"),""),metric="Protein groups")
nsaf_sym2<-nsaf_sym%>%dplyr::mutate(symbiont5=factor(symbiont5,levels=sym_stack_order),value=pct,label=ifelse(pct>2,paste0("(",round(pct,1),"%)"),""),n=NA_integer_,metric="NSAF abundance")
overview_left<-dplyr::bind_rows(pg_sym2,nsaf_sym2)%>%dplyr::mutate(metric=factor(metric,levels=c("Protein groups","NSAF abundance")))

# Plot protein and NSAF proportions by partner.
p_left<-ggplot(overview_left,aes(x=metric,y=value,fill=symbiont5))+
  geom_col(position="stack",color="grey25",linewidth=0.4,width=0.80)+
  geom_text(aes(label=label),position=position_stack(vjust=0.5),size=3.25,color="white",fontface="bold",lineheight=0.9)+
  geom_vline(xintercept=1.5,linetype="dashed",color="gray40",linewidth=0.6)+
  scale_fill_manual(values=sym_colors_ext,name="Symbiont")+
  scale_x_discrete(labels=c("Protein groups"="Protein\ngroups","NSAF abundance"="NSAF\nabundance"),expand=expansion(add=0.62))+
  scale_y_continuous(expand=expansion(mult=c(0,0.01)),limits=c(0,100),labels=function(x)paste0(x,"%"))+
  labs(x=NULL,y="Proportion (%)")+base_theme+
  theme(legend.position="none",
        axis.text.x=element_text(size=12.5,lineheight=1,margin=margin(t=2)),
        plot.margin=margin(1,4,1,4,"pt"))

nsaf_imp_df<-cbind(meta_nsaf_filt,as.data.frame(mat_nsaf_imp,optional=TRUE))
# Select the ten most abundant annotated functional roles.
sr_nsaf_long_ov<-nsaf_imp_df%>%mutate(subRole_clean=dplyr::case_when(is.na(subRole)|subRole==""~NA_character_,subRole=="various"~NA_character_,grepl("unknown|uncharacterized|hypothetical|DUF|unclassified",subRole,ignore.case=TRUE)~NA_character_,TRUE~subRole),symbiont5=sym5(symbionts))%>%filter(!is.na(subRole_clean))%>%select(members_identifier,symbiont5,subRole_clean,all_of(nsaf_avail))%>%pivot_longer(cols=all_of(nsaf_avail),names_to="sample",values_to="nsaf_val")
top10_sr_ov<-sr_nsaf_long_ov%>%group_by(subRole_clean)%>%summarise(total=sum(nsaf_val,na.rm=TRUE),.groups="drop")%>%slice_max(total,n=10,with_ties=FALSE)%>%arrange(total)%>%pull(subRole_clean)

sr_sym_summ_ov<-sr_nsaf_long_ov%>%filter(subRole_clean%in%top10_sr_ov)%>%group_by(symbiont5,subRole_clean,sample)%>%summarise(sample_sum=sum(nsaf_val,na.rm=TRUE),.groups="drop")%>%group_by(symbiont5,subRole_clean)%>%summarise(mean_nsaf=mean(sample_sum,na.rm=TRUE),.groups="drop")

subrole_order<-tibble(subRole_clean=top10_sr_ov)%>%
  left_join(sr_sym_summ_ov%>%filter(symbiont5=="Mycobiont")%>%select(subRole_clean,myco_mean=mean_nsaf),by="subRole_clean")%>%
  mutate(myco_mean=tidyr::replace_na(myco_mean,0))%>%
  arrange(myco_mean)%>%
  pull(subRole_clean)
sr_sym_summ_ov<-sr_sym_summ_ov%>%mutate(symbiont5=factor(symbiont5,levels=sym_stack_order),subRole_clean=factor(subRole_clean,levels=subrole_order))

# Plot the most abundant functional roles.
p_right<-ggplot(sr_sym_summ_ov,aes(x=subRole_clean,y=mean_nsaf,fill=symbiont5))+
  geom_col(position="stack",color="grey25",linewidth=0.3,width=0.90)+
  coord_flip(clip="off")+
  scale_fill_manual(values=sym_colors_ext,name="Symbiont")+
  scale_y_continuous(breaks=c(0,0.05,0.10),expand=expansion(mult=c(0,0.12)))+
  labs(x=NULL,y="Mean NSAF per specimen\n(filtered + imputed)")+base_theme+
  theme(
    axis.text.y=element_text(size=12,hjust=1,margin=margin(r=10)),
    panel.grid.major.y=element_blank(),
    legend.position="inside",
    legend.position.inside=c(0.98,0.04),
    legend.justification=c(1,0),
    legend.background=element_rect(fill="white",colour="grey70",linewidth=0.3),
    legend.key.size=unit(0.34,"cm"),
    legend.text=element_text(size=9.5),
    legend.title=element_text(size=10,face="bold"),
    plot.margin=margin(1,8,1,2,"pt")
  )

p_left_comb<-p_left+labs(tag="A")+
  theme(plot.tag=element_text(size=21,face="bold"),plot.tag.position="topleft")

p_right_comb<-p_right+labs(tag="B")+
  theme(plot.tag=element_text(size=21,face="bold",hjust=1),plot.tag.position=c(0.57,0.965))
p_sym_overview<-(p_left_comb|p_right_comb)+plot_layout(widths=c(2.20,3.80)) & theme(plot.margin=margin(0))
cat("Saved: symbiont_overview.*\n")

# Keep the full functional-role table, including unclassified roles.
sr_all_nsaf<-nsaf_imp_df%>%mutate(subRole_label=dplyr::case_when(is.na(subRole)|subRole==""~"Unclassified",subRole=="various"~"Various",grepl("unknown|uncharacterized|hypothetical|DUF|unclassified",subRole,ignore.case=TRUE)~"Function unknown",TRUE~subRole))%>%select(members_identifier,subRole_label,all_of(nsaf_avail))%>%pivot_longer(cols=all_of(nsaf_avail),names_to="sample",values_to="nsaf_val")%>%group_by(subRole_label)%>%summarise(total_nsaf=sum(nsaf_val,na.rm=TRUE),n_proteins=n_distinct(members_identifier),.groups="drop")%>%mutate(pct_nsaf=total_nsaf/sum(total_nsaf)*100)%>%arrange(desc(total_nsaf))
write.csv(sr_all_nsaf,file.path(OUTDIR,"subrole_full_table_nsaf.csv"),row.names=FALSE)
cat("Section 12 complete.\n\n")

# Save both panels and the combined Figure 8 SVG.
fig8_dir <- file.path(INDIVIDUAL_DIR,"Figure8")
save_svg(file.path(fig8_dir,"Figure8A_partner_proportions"), p_left, width=4.7, height=5.07)

p_right_ind <- p_right +
  theme(plot.margin=margin(1,8,1,2,"pt"))
save_svg(file.path(fig8_dir,"Figure8B_top10_functional_subroles"), p_right_ind, width=9.6, height=5.6)

save_svg(file.path(OUTDIR,"Figure8_partner_contributions"), p_sym_overview, width=11.90, height=5.00)
