# Figures 3 and 4: community and protein composition.
# NMDS shows distances among specimens; PERMANOVA tests treatment terms.
# RAT/DTP is analysed separately from the light × nitrogen experiment.

cat("\nNMDS and PERMANOVA\n")

# Draw NMDS points, group centroids, and the displayed R2 and P values.
make_nmds_plot<-function(pts,spiders_df,ellipse_df,grp_var,grp_levels,col_vec,perm_p,perm_R2,xr=NULL,yr=NULL,label_offset=NULL){

  x_rng<-if(!is.null(xr)) xr else range(c(pts$NMDS1,spiders_df$NMDS1,ellipse_df$NMDS1),na.rm=TRUE)
  y_rng<-if(!is.null(yr)) yr else range(c(pts$NMDS2,spiders_df$NMDS2,ellipse_df$NMDS2),na.rm=TRUE)
  spiders_df$label_x<-spiders_df$cntr.x
  if(!is.null(label_offset)) spiders_df$label_x<-spiders_df$label_x+unname(label_offset[as.character(spiders_df[[grp_var]])])*diff(x_rng)
  x_ann<-x_rng[1]+0.68*diff(x_rng)
  y_r2 <-y_rng[1]+0.96*diff(y_rng)
  y_p  <-y_rng[1]+0.88*diff(y_rng)
  # Figure labels use two decimal places; exact values remain in the results CSV.
  p_label<-if(is.na(perm_p)) "P=NA" else if(perm_p<0.005) "P<0.01" else sprintf("P=%.2f",perm_p)
  p<-ggplot()+
    geom_point(data=pts,aes(x=NMDS1,y=NMDS2,fill=.data[[grp_var]],color=.data[[grp_var]]),size=2.5,stroke=0.8,alpha=0.7)+
    geom_segment(data=spiders_df,aes(x=cntr.x,y=cntr.y,xend=NMDS1,yend=NMDS2,color=.data[[grp_var]]),linewidth=0.5,alpha=0.7)+
    geom_point(data=spiders_df,aes(x=NMDS1,y=NMDS2,fill=.data[[grp_var]],color=.data[[grp_var]]),size=2.5,stroke=0.9,alpha=0.7)+
    geom_label(data=spiders_df,aes(x=label_x,y=cntr.y,label=unname(disp_label[as.character(.data[[grp_var]])]),color=.data[[grp_var]]),fontface="bold",size=14/.pt,fill="white",alpha=1,label.padding=unit(0.15,"lines"))+
    geom_path(data=ellipse_df,aes(x=NMDS1,y=NMDS2,color=.data[[grp_var]]),linewidth=0.75,linetype=2,alpha=0.5)+
    scale_fill_manual(values=col_vec)+scale_color_manual(values=col_vec)+
    annotate("text",x=x_ann,y=y_r2,label=sprintf("R\u00b2=%.2f",perm_R2),hjust=0,vjust=1,size=11/.pt,color="black",fontface="bold.italic")+
    annotate("text",x=x_ann,y=y_p,label=p_label,hjust=0,vjust=1,size=11/.pt,color="black",fontface="bold.italic")+
    labs(x="NMDS1",y="NMDS2")+base_theme+theme(legend.position="none",plot.margin=margin(0,0,0,0,"pt"))
  if(!is.null(xr)||!is.null(yr)) p<-p+coord_cartesian(xlim=xr,ylim=yr,clip="off")
  p}

.common_range<-function(...,pad=0.06){
  dfs<-list(...)
  xs<-unlist(lapply(dfs,function(d) if(!is.null(d)&&nrow(d)>0) d$NMDS1 else NULL))
  ys<-unlist(lapply(dfs,function(d) if(!is.null(d)&&nrow(d)>0) d$NMDS2 else NULL))
  xr<-range(xs,na.rm=TRUE); yr<-range(ys,na.rm=TRUE)
  xpad<-diff(xr)*pad; ypad<-diff(yr)*pad
  list(x=c(xr[1]-xpad,xr[2]+xpad), y=c(yr[1]-ypad,yr[2]+ypad))}

# Use Bray–Curtis distances for the 16S community analysis.
dist_full_16S<-phyloseq::distance(phyloseq_16S_filt,method="bray")
metadata_full_16S<-data.frame(sample_data(phyloseq_16S_filt))
dist_treated_16S<-phyloseq::distance(ps_treated_16S,method="bray"); metadata_treated_16S<-data.frame(sample_data(ps_treated_16S))
set.seed(345); nmds_treated_16S<-metaMDS(dist_treated_16S,k=2,try=50,trymax=200,autotransform=FALSE,trace=FALSE)
pts_treated_16S_xy<-as.data.frame(nmds_treated_16S$points)%>%setNames(c("NMDS1","NMDS2"))
pts_treated_16S<-cbind(pts_treated_16S_xy,metadata_treated_16S[rownames(pts_treated_16S_xy),,drop=FALSE])
dist_transplant_16S<-phyloseq::distance(ps_transplant_16S,method="bray"); metadata_transplant_16S<-data.frame(sample_data(ps_transplant_16S))
set.seed(345); nmds_transplant_16S<-metaMDS(dist_transplant_16S,k=2,try=50,trymax=200,autotransform=FALSE,trace=FALSE)
pts_transplant_16S_xy<-as.data.frame(nmds_transplant_16S$points)%>%setNames(c("NMDS1","NMDS2"))
pts_transplant_16S<-cbind(pts_transplant_16S_xy,metadata_transplant_16S[rownames(pts_transplant_16S_xy),,drop=FALSE])
cat(sprintf("  NMDS stress (16S factorial): %.4f\n",nmds_treated_16S$stress))
cat(sprintf("  NMDS stress (16S RAT-DTP):   %.4f\n",nmds_transplant_16S$stress))

set.seed(345); perm_treated_16S<-adonis2(dist_treated_16S~Light+Nitrogen+Light:Nitrogen,data=metadata_treated_16S,permutations=9999,by="terms")

set.seed(345); perm_treated_16S_margin<-adonis2(dist_treated_16S~Light+Nitrogen,data=metadata_treated_16S,permutations=9999,by="margin")
set.seed(345); perm_transplant_16S<-adonis2(dist_transplant_16S~Transplant,data=metadata_transplant_16S,permutations=9999,by="terms")

perm_16S_R2_light<-perm_treated_16S_margin["Light","R2"]; perm_16S_p_light<-perm_treated_16S_margin["Light","Pr(>F)"]
perm_16S_R2_nitro<-perm_treated_16S_margin["Nitrogen","R2"]; perm_16S_p_nitro<-perm_treated_16S_margin["Nitrogen","Pr(>F)"]
perm_16S_R2_int<-perm_treated_16S["Light:Nitrogen","R2"]; perm_16S_p_int<-perm_treated_16S["Light:Nitrogen","Pr(>F)"]
perm_16S_R2_trans<-perm_transplant_16S["Transplant","R2"]; perm_16S_p_trans<-perm_transplant_16S["Transplant","Pr(>F)"]

.get_ordi<-function(nmds_obj,metadata,grp){
  o<-gg_ordiplot(nmds_obj,groups=metadata[[grp]],plot=FALSE,kind="se",conf=0.95,ellipse=TRUE,hull=FALSE,spiders=TRUE)
  sp<-o$df_spiders%>%dplyr::rename(NMDS1=x,NMDS2=y,!!grp:=Group)%>%as_tibble()%>%filter(!is.na(.data[[grp]]))
  el<-o$df_ellipse%>%dplyr::rename(NMDS1=x,NMDS2=y,!!grp:=Group)%>%filter(!is.na(.data[[grp]]))
  list(sp=sp,el=el)}

o16L<-.get_ordi(nmds_treated_16S,metadata_treated_16S,"Light"); o16N<-.get_ordi(nmds_treated_16S,metadata_treated_16S,"Nitrogen"); o16T<-.get_ordi(nmds_transplant_16S,metadata_transplant_16S,"Transplant")
p_nmds_L16<-make_nmds_plot(pts_treated_16S,o16L$sp%>%filter(Light%in%c("LLS","HLS")),o16L$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_16S_p_light,perm_16S_R2_light)
p_nmds_N16<-make_nmds_plot(pts_treated_16S,o16N$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),o16N$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_16S_p_nitro,perm_16S_R2_nitro)
p_nmds_T16<-make_nmds_plot(pts_transplant_16S,o16T$sp%>%filter(Transplant%in%c("RAT","DTP")),o16T$el%>%filter(Transplant%in%c("RAT","DTP")),"Transplant",c("RAT","DTP"),c("RAT"=col_RAT,"DTP"=col_DTP),perm_16S_p_trans,perm_16S_R2_trans)

dist_full_18S<-phyloseq::distance(phyloseq_18S_filt,method="bray")
metadata_full_18S<-data.frame(sample_data(phyloseq_18S_filt))
# Use Bray–Curtis distances for the 18S community analysis.
dist_treated_18S<-phyloseq::distance(ps_treated_18S,method="bray"); metadata_treated_18S<-data.frame(sample_data(ps_treated_18S))
dist_transplant_18S<-phyloseq::distance(ps_transplant_18S,method="bray"); metadata_transplant_18S<-data.frame(sample_data(ps_transplant_18S))

set.seed(345); nmds_full_18S<-metaMDS(dist_full_18S,k=2,try=50,trymax=200,
                                      autotransform=FALSE,trace=FALSE)
pts_full_18S_xy<-as.data.frame(nmds_full_18S$points)%>%setNames(c("NMDS1","NMDS2"))
pts_full_18S<-cbind(pts_full_18S_xy,
                    metadata_full_18S[rownames(pts_full_18S_xy),,drop=FALSE])
pts_treated_18S<-pts_full_18S[!is.na(pts_full_18S$Light),,drop=FALSE]
pts_transplant_18S<-pts_full_18S[!is.na(pts_full_18S$Transplant),,drop=FALSE]

nmds_treated_18S<-nmds_full_18S
nmds_transplant_18S<-nmds_full_18S
cat(sprintf("  NMDS stress (18S; all 18 specimens): %.4f\n",nmds_full_18S$stress))

set.seed(345); perm_treated_18S<-adonis2(dist_treated_18S~Light+Nitrogen+Light:Nitrogen,data=metadata_treated_18S,permutations=9999,by="terms")
set.seed(345); perm_treated_18S_margin<-adonis2(dist_treated_18S~Light+Nitrogen,data=metadata_treated_18S,permutations=9999,by="margin")
set.seed(345); perm_transplant_18S<-adonis2(dist_transplant_18S~Transplant,data=metadata_transplant_18S,permutations=9999,by="terms")

perm_18S_R2_light<-perm_treated_18S_margin["Light","R2"]; perm_18S_p_light<-perm_treated_18S_margin["Light","Pr(>F)"]
perm_18S_R2_nitro<-perm_treated_18S_margin["Nitrogen","R2"]; perm_18S_p_nitro<-perm_treated_18S_margin["Nitrogen","Pr(>F)"]
perm_18S_R2_int<-perm_treated_18S["Light:Nitrogen","R2"]; perm_18S_p_int<-perm_treated_18S["Light:Nitrogen","Pr(>F)"]
perm_18S_R2_trans<-perm_transplant_18S["Transplant","R2"]; perm_18S_p_trans<-perm_transplant_18S["Transplant","Pr(>F)"]

o18L<-.get_ordi(nmds_full_18S,metadata_full_18S,"Light"); o18N<-.get_ordi(nmds_full_18S,metadata_full_18S,"Nitrogen"); o18T<-.get_ordi(nmds_full_18S,metadata_full_18S,"Transplant")
p_nmds_L18<-make_nmds_plot(pts_treated_18S,o18L$sp%>%filter(Light%in%c("LLS","HLS")),o18L$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_18S_p_light,perm_18S_R2_light)
p_nmds_N18<-make_nmds_plot(pts_treated_18S,o18N$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),o18N$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_18S_p_nitro,perm_18S_R2_nitro,label_offset=c("N0.32"=0.12,"N1.28"=-0.12))
p_nmds_T18<-make_nmds_plot(pts_transplant_18S,o18T$sp%>%filter(Transplant%in%c("RAT","DTP")),o18T$el%>%filter(Transplant%in%c("RAT","DTP")),"Transplant",c("RAT","DTP"),c("RAT"=col_RAT,"DTP"=col_DTP),perm_18S_p_trans,perm_18S_R2_trans)

mat_prot_nmds<-t(mat_nsaf_log2)
mat_prot_nmds[is.na(mat_prot_nmds)]<-min(mat_nsaf_log2,na.rm=TRUE)
treat_ids_nsaf<-c("nsaf_N0.32-LLS1","nsaf_N0.32-LLS2","nsaf_N0.32-LLS3","nsaf_N0.32-HLS1","nsaf_N0.32-HLS2","nsaf_N0.32-HLS3",
                  "nsaf_N1.28-LLS1","nsaf_N1.28-LLS2","nsaf_N1.28-LLS3","nsaf_N1.28-HLS1","nsaf_N1.28-HLS2","nsaf_N1.28-HLS3")
dist_treated_prot<-as.dist(as.matrix(vegan::vegdist(mat_prot_nmds,method="euclidean"))[treat_ids_nsaf,treat_ids_nsaf])
trans_ids_nsaf<-c("nsaf_RAT1","nsaf_RAT2","nsaf_RAT3","nsaf_DTP1","nsaf_DTP2","nsaf_DTP3")
dist_transplant_prot<-as.dist(as.matrix(vegan::vegdist(mat_prot_nmds,method="euclidean"))[trans_ids_nsaf,trans_ids_nsaf])
meta_prot_treat<-data.frame(
  Light=factor(c(rep("LLS",3),rep("HLS",3),rep("LLS",3),rep("HLS",3)),levels=c("LLS","HLS")),
  Nitrogen=factor(c(rep("N0.32",3),rep("N0.32",3),rep("N1.28",3),rep("N1.28",3)),levels=c("N0.32","N1.28")),
  row.names=treat_ids_nsaf)
meta_prot_transplant<-data.frame(
  Transplant=factor(c(rep("RAT",3),rep("DTP",3)),levels=c("RAT","DTP")),
  row.names=trans_ids_nsaf)
set.seed(345); perm_treated_prot<-adonis2(dist_treated_prot~Light+Nitrogen+Light:Nitrogen,data=meta_prot_treat,permutations=9999,by="terms")
set.seed(345); perm_treated_prot_margin<-adonis2(dist_treated_prot~Light+Nitrogen,data=meta_prot_treat,permutations=9999,by="margin")
set.seed(345); perm_trans_prot<-adonis2(dist_transplant_prot~Transplant,data=meta_prot_transplant,permutations=9999,by="terms")

perm_prot_R2_light<-perm_treated_prot_margin["Light","R2"]; perm_prot_p_light<-perm_treated_prot_margin["Light","Pr(>F)"]
perm_prot_R2_nitro<-perm_treated_prot_margin["Nitrogen","R2"]; perm_prot_p_nitro<-perm_treated_prot_margin["Nitrogen","Pr(>F)"]
perm_prot_R2_int<-perm_treated_prot["Light:Nitrogen","R2"]; perm_prot_p_int<-perm_treated_prot["Light:Nitrogen","Pr(>F)"]
perm_prot_R2_trans<-perm_trans_prot["Transplant","R2"]; perm_prot_p_trans<-perm_trans_prot["Transplant","Pr(>F)"]

# Use Euclidean distances of log2 NSAF values for the proteins.
set.seed(345); nmds_prot_LN<-metaMDS(dist_treated_prot,k=2,try=50,trymax=200,autotransform=FALSE,trace=FALSE)
pts_prot_LN_xy<-data.frame(nmds_prot_LN$points)%>%setNames(c("NMDS1","NMDS2"))
pts_prot_LN<-cbind(pts_prot_LN_xy,meta_prot_treat[rownames(pts_prot_LN_xy),,drop=FALSE])
set.seed(345); nmds_prot_T<-metaMDS(dist_transplant_prot,k=2,try=50,trymax=200,autotransform=FALSE,trace=FALSE)
pts_prot_T_xy<-data.frame(nmds_prot_T$points)%>%setNames(c("NMDS1","NMDS2"))
pts_prot_T<-cbind(pts_prot_T_xy,meta_prot_transplant[rownames(pts_prot_T_xy),,drop=FALSE])
cat(sprintf("  NMDS stress (metaproteome factorial): %.4f\n",nmds_prot_LN$stress))
cat(sprintf("  NMDS stress (metaproteome RAT-DTP):   %.4f\n",nmds_prot_T$stress))

oLp<-.get_ordi(nmds_prot_LN,meta_prot_treat,"Light"); oNp<-.get_ordi(nmds_prot_LN,meta_prot_treat,"Nitrogen"); oTp<-.get_ordi(nmds_prot_T,meta_prot_transplant,"Transplant")
p_nmds_Lp<-make_nmds_plot(pts_prot_LN,oLp$sp%>%filter(Light%in%c("LLS","HLS")),oLp$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_prot_p_light,perm_prot_R2_light)+labs(title=NULL)
p_nmds_Np<-make_nmds_plot(pts_prot_LN,oNp$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),oNp$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_prot_p_nitro,perm_prot_R2_nitro)+labs(title=NULL)
p_nmds_Tp<-make_nmds_plot(pts_prot_T,oTp$sp,oTp$el,"Transplant",c("RAT","DTP"),c("RAT"=col_RAT,"DTP"=col_DTP),perm_prot_p_trans,perm_prot_R2_trans)+labs(title=NULL)

make_strip_box_v<-function(label,top_frac=0,bottom_frac=0){
  y0 <- bottom_frac
  y1 <- 1-top_frac
  ggplot()+
    annotate("rect",xmin=0,xmax=1,ymin=y0,ymax=y1,fill="gray70",color="black",linewidth=0.5)+
    annotate("text",x=0.5,y=(y0+y1)/2,label=label,fontface="bold",size=20/.pt,color="black",angle=90)+
    scale_x_continuous(limits=c(0,1),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),expand=c(0,0))+
    coord_cartesian(clip="off",expand=FALSE)+
    theme_void()+
    theme(plot.margin=margin(0))
}

rng_16S  <- .common_range(pts_treated_16S, o16L$sp, o16L$el, o16N$sp, o16N$el)
rng_18S  <- .common_range(pts_treated_18S, o18L$sp, o18L$el, o18N$sp, o18N$el)
rng_prot <- .common_range(pts_prot_LN,     oLp$sp,  oLp$el,  oNp$sp,  oNp$el)

p_nmds_L16<-make_nmds_plot(pts_treated_16S,o16L$sp%>%filter(Light%in%c("LLS","HLS")),o16L$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_16S_p_light,perm_16S_R2_light,xr=rng_16S$x,yr=rng_16S$y)
p_nmds_N16<-make_nmds_plot(pts_treated_16S,o16N$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),o16N$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_16S_p_nitro,perm_16S_R2_nitro,xr=rng_16S$x,yr=rng_16S$y)
p_nmds_L18<-make_nmds_plot(pts_treated_18S,o18L$sp%>%filter(Light%in%c("LLS","HLS")),o18L$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_18S_p_light,perm_18S_R2_light,xr=rng_18S$x,yr=rng_18S$y)
p_nmds_N18<-make_nmds_plot(pts_treated_18S,o18N$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),o18N$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_18S_p_nitro,perm_18S_R2_nitro,xr=rng_18S$x,yr=rng_18S$y,label_offset=c("N0.32"=0.12,"N1.28"=-0.12))
p_nmds_Lp<-make_nmds_plot(pts_prot_LN,oLp$sp%>%filter(Light%in%c("LLS","HLS")),oLp$el%>%filter(Light%in%c("LLS","HLS")),"Light",c("LLS","HLS"),c("LLS"=col_LLS,"HLS"=col_HLS),perm_prot_p_light,perm_prot_R2_light,xr=rng_prot$x,yr=rng_prot$y)+labs(title=NULL)
p_nmds_Np<-make_nmds_plot(pts_prot_LN,oNp$sp%>%filter(Nitrogen%in%c("N0.32","N1.28")),oNp$el%>%filter(Nitrogen%in%c("N0.32","N1.28")),"Nitrogen",c("N0.32","N1.28"),c("N0.32"=col_N0.32,"N1.28"=col_N1.28),perm_prot_p_nitro,perm_prot_R2_nitro,xr=rng_prot$x,yr=rng_prot$y)+labs(title=NULL)

.no_y<-theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

.tag_nmds<-function(p,tag,xr,yr){
  p + annotate("text",
               x=xr[1]+0.035*diff(xr),
               y=yr[2]-0.035*diff(yr),
               label=tag,hjust=0,vjust=1,
               size=22/.pt,fontface="bold",colour="black")
}

.nmds_row_3<-function(pL,pR,lbl,tagL,tagR,xr,yr,strip_top_frac=0,strip_bottom_frac=0){
  edge_col <- "gray50"
  edge_lwd <- 0.5
  pL_fixed<-pL+
    scale_x_continuous(labels=function(b){if(length(b)>0)b[length(b)]<-"";b},expand=expansion(mult=0))+
    scale_y_continuous(expand=expansion(mult=0))+
    coord_cartesian(xlim=xr,ylim=yr,clip="off",expand=FALSE)+
    theme(panel.border=element_blank())+
    annotate("segment",x=xr[1],xend=xr[2],y=yr[2],yend=yr[2],colour=edge_col,linewidth=edge_lwd)+
    annotate("segment",x=xr[1],xend=xr[2],y=yr[1],yend=yr[1],colour=edge_col,linewidth=edge_lwd)+
    annotate("segment",x=xr[1],xend=xr[1],y=yr[1],yend=yr[2],colour=edge_col,linewidth=edge_lwd)+
    annotate("segment",x=xr[2]-0.004*diff(xr),xend=xr[2]-0.004*diff(xr),y=yr[1],yend=yr[2],linetype="dashed",colour="gray40",linewidth=0.4)
  pR_fixed<-pR+
    scale_x_continuous(expand=expansion(mult=0))+
    scale_y_continuous(expand=expansion(mult=0))+
    coord_cartesian(xlim=xr,ylim=yr,clip="off",expand=FALSE)+
    theme(panel.border=element_blank())+
    annotate("segment",x=xr[1],xend=xr[2],y=yr[2],yend=yr[2],colour=edge_col,linewidth=edge_lwd)+
    annotate("segment",x=xr[1],xend=xr[2],y=yr[1],yend=yr[1],colour=edge_col,linewidth=edge_lwd)+
    annotate("segment",x=xr[2],xend=xr[2],y=yr[1],yend=yr[2],colour=edge_col,linewidth=edge_lwd)
  patchwork::wrap_plots(
    .tag_nmds(pL_fixed + labs(title=NULL) + theme(plot.margin=margin(0)),tagL,xr,yr),
    .tag_nmds(pR_fixed + labs(title=NULL) + .no_y + theme(plot.margin=margin(0)),tagR,xr,yr),
    patchwork::wrap_elements(full=make_strip_box_v(lbl,top_frac=strip_top_frac,bottom_frac=strip_bottom_frac)),
    ncol = 3,
    widths = c(1, 1, 0.058)
  )
}

# Arrange the 16S/18S and protein NMDS panels.
fig3_dir <- file.path(INDIVIDUAL_DIR,"Figure3")
fig4_dir <- file.path(INDIVIDUAL_DIR,"Figure4")
save_svg(file.path(fig3_dir,"Figure3A_16S_Light"),    p_nmds_L16, width=4.5, height=4.5)
save_svg(file.path(fig3_dir,"Figure3B_16S_Nitrogen"), p_nmds_N16, width=4.5, height=4.5)
save_svg(file.path(fig3_dir,"Figure3C_18S_Light"),    p_nmds_L18, width=4.5, height=4.5)
save_svg(file.path(fig3_dir,"Figure3D_18S_Nitrogen"), p_nmds_N18, width=4.5, height=4.5)
save_svg(file.path(fig4_dir,"Figure4A_metaproteome_Light"),    p_nmds_Lp, width=4.5, height=4.5)
save_svg(file.path(fig4_dir,"Figure4B_metaproteome_Nitrogen"), p_nmds_Np, width=4.5, height=4.5)

row_16S  <- .nmds_row_3(p_nmds_L16+theme(axis.title.x=element_blank(),axis.title.y=element_blank()), p_nmds_N16+theme(axis.title.x=element_blank(),axis.title.y=element_blank()), "16S amplicon data", "A", "B", rng_16S$x, rng_16S$y, strip_bottom_frac=0.055)
row_18S  <- .nmds_row_3(p_nmds_L18+labs(x="NMDS1")+theme(axis.title.y=element_blank()), p_nmds_N18+labs(x="NMDS1")+theme(axis.title.y=element_blank()), "18S amplicon data", "C", "D", rng_18S$x, rng_18S$y, strip_bottom_frac=0.115)
row_prot <- .nmds_row_3(p_nmds_Lp+theme(axis.title.y=element_blank()),  p_nmds_Np+theme(axis.title.y=element_blank()),  "metaproteome data", "A", "B", rng_prot$x, rng_prot$y, strip_top_frac=0.018, strip_bottom_frac=0.135)

nmds_combined_amplicons_core <- patchwork::wrap_plots(
  row_16S, row_18S,
  ncol = 1,
  heights = c(1, 1)
) & theme(panel.spacing=unit(0,"pt"),plot.margin=margin(0))
nmds_combined_amplicons <- add_common_y_title(nmds_combined_amplicons_core, "NMDS2", text_size=14, width=0.026)
nmds_metaproteome_only <- add_common_y_title(row_prot, "NMDS2", text_size=14, width=0.026)

.strip_label<-function(lbl){ggplot()+annotate("rect",xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill="gray70",color="black",linewidth=0.5)+annotate("text",x=0,y=0,label=lbl,fontface="bold",size=3.5,color="black")+theme_void()}
col_T16<-(wrap_elements(.strip_label("16S amplicon data"))/(p_nmds_T16+theme(plot.margin=margin(2,0,2,2,"pt"))))+plot_layout(heights=c(0.04,1))
col_T18<-(wrap_elements(.strip_label("18S amplicon data"))/((p_nmds_T18&theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()))+theme(plot.margin=margin(2,2,2,0,"pt"))))+plot_layout(heights=c(0.04,1))

.prow<-function(tbl,term,ds,ctr,type){r<-tbl[term,];data.frame(Dataset=ds,Contrast=ctr,Term=term,Type=type,Df=r[["Df"]],SumOfSqs=round(r[["SumOfSqs"]],4),R2=round(r[["R2"]],4),F=round(r[["F"]],3),p.value=r[["Pr(>F)"]],stringsAsFactors=FALSE)}
# Save sequential full-model and marginal additive-model PERMANOVA
# results, plus the separate RAT/DTP comparison.
perm_full_table<-bind_rows(

  .prow(perm_treated_16S,"Light","16S","Light","Sequential"),.prow(perm_treated_16S,"Nitrogen","16S","Nitrogen","Sequential"),.prow(perm_treated_16S,"Light:Nitrogen","16S","Light x Nitrogen","Sequential"),.prow(perm_treated_16S,"Residual","16S","[Residual]","Sequential"),

  .prow(perm_treated_16S_margin,"Light","16S","Light","Marginal"),.prow(perm_treated_16S_margin,"Nitrogen","16S","Nitrogen","Marginal"),.prow(perm_treated_16S_margin,"Residual","16S","[Residual] (additive model)","Marginal"),

  .prow(perm_transplant_16S,"Transplant","16S","RAT vs DTP","Sequential=Marginal"),

  .prow(perm_treated_18S,"Light","18S","Light","Sequential"),.prow(perm_treated_18S,"Nitrogen","18S","Nitrogen","Sequential"),.prow(perm_treated_18S,"Light:Nitrogen","18S","Light x Nitrogen","Sequential"),.prow(perm_treated_18S,"Residual","18S","[Residual]","Sequential"),

  .prow(perm_treated_18S_margin,"Light","18S","Light","Marginal"),.prow(perm_treated_18S_margin,"Nitrogen","18S","Nitrogen","Marginal"),.prow(perm_treated_18S_margin,"Residual","18S","[Residual] (additive model)","Marginal"),

  .prow(perm_transplant_18S,"Transplant","18S","RAT vs DTP","Sequential=Marginal"),

  .prow(perm_treated_prot,"Light","Metaproteome","Light","Sequential"),.prow(perm_treated_prot,"Nitrogen","Metaproteome","Nitrogen","Sequential"),.prow(perm_treated_prot,"Light:Nitrogen","Metaproteome","Light x Nitrogen","Sequential"),.prow(perm_treated_prot,"Residual","Metaproteome","[Residual]","Sequential"),

  .prow(perm_treated_prot_margin,"Light","Metaproteome","Light","Marginal"),.prow(perm_treated_prot_margin,"Nitrogen","Metaproteome","Nitrogen","Marginal"),.prow(perm_treated_prot_margin,"Residual","Metaproteome","[Residual] (additive model)","Marginal"),

  .prow(perm_trans_prot,"Transplant","Metaproteome","RAT vs DTP","Sequential=Marginal"))%>%
  mutate(sig=case_when(is.na(p.value)~"",p.value<0.001~"***",p.value<0.01~"**",p.value<0.05~"*",p.value<0.1~".",TRUE~"ns"))
write.csv(perm_full_table,file.path(OUTDIR,"permanova_full_table.csv"),row.names=FALSE)
cat("Saved: permanova_full_table.csv (factorial sequential/marginal results plus separate RAT-DTP assessment)\n")
cat("\n--- Sequential (full model) vs Marginal (additive-model) R2 comparison ---\n")
cat("    Light/Nitrogen main effects only — Light:Nitrogen has no Marginal row\n")
cat("    (that term is absent from the additive refit; see note above)\n")
print(perm_full_table%>%filter(Type%in%c("Sequential","Marginal"),Term!="Residual")%>%select(Dataset,Term,Type,R2,p.value,sig)%>%as.data.frame(),row.names=FALSE)
cat("\n--- Separate RAT-DTP transplantation assessment ---\n")
print(perm_full_table%>%filter(Contrast=="RAT vs DTP")%>%select(Dataset,Contrast,R2,p.value,sig)%>%as.data.frame(),row.names=FALSE)
cat("Section 8 complete.\n\n")

save_svg(file.path(OUTDIR,"Figure3_NMDS_amplicons"), nmds_combined_amplicons, width=8.55, height=8.20)
save_svg(file.path(OUTDIR,"Figure4_NMDS_metaproteome"), nmds_metaproteome_only, width=8.55, height=4.00)
