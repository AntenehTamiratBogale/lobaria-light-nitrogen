# Update Figure 9 labels and fonts only.
# Run after run_all.R in the same R session.
# This script does not rerun statistical tests.
.add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(size = 32, face = "bold", colour = "black"),
          plot.tag.position = c(0.02, 0.98))
}

.font_fix21 <- function(p, has_legend=FALSE) {
  p <- p + theme(
    axis.text.x  = element_text(size=24, colour="black", lineheight=1.06),
    axis.text.y  = element_text(size=22, colour="black"),
    axis.title   = element_text(size=30),
    plot.title   = element_text(size=20, face="bold", hjust=0.5, lineheight=1.15),
    plot.subtitle= element_text(size=20, colour="grey35", hjust=0.5, lineheight=1.2),
    strip.text   = element_text(size=24, colour="black"),
    strip.text.y.right = element_text(size=24, face="bold", angle=90, lineheight=1.10),
    plot.margin=margin(0.5,0.5,0.5,0.5,"pt"))
  if (has_legend) p <- p + theme(legend.text=element_text(size=28), legend.title=element_text(size=29,face="bold"),legend.margin=margin(0),legend.box.margin=margin(0),legend.box.spacing=unit(0,"pt"))
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
  labs(title=NULL, subtitle=NULL)
photo_nitro_ind21 <- .no_cmp_strip21(photo_nitro_main21) +
  labs(title=NULL, subtitle=NULL)
nmet_light_ind21 <- .no_cmp_strip21(nmet_light_main21) +
  labs(title=NULL, subtitle=NULL)
nmet_nitro_ind21 <- .no_cmp_strip21(nmet_nitro_main21) +
  labs(title=NULL, subtitle=NULL)

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
             colour="black", fontface="bold", size=32/.pt,
             hjust=0.5, vjust=0.5) +
    theme_void() +
    theme(plot.margin=margin(0,0,0,0,"pt"))
}

.hide_x21 <- function(p) {
  p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

pA21 <- .add_tag(.hide_x21(.no_cmp_strip21(photo_light_main21)) +
                   labs(subtitle=NULL) +
                   theme(axis.title.y=element_blank()), "A")
pB21 <- .add_tag(.hide_x21(.no_cmp_strip21(nmet_light_main21)) |>
                   .strip_y() + labs(subtitle=NULL), "B")
pC21 <- .add_tag(.no_cmp_strip21(photo_nitro_main21) +
                   labs(subtitle=NULL) +
                   theme(axis.title.y=element_blank()), "C")
pD21 <- .add_tag(.no_cmp_strip21(nmet_nitro_main21) |>
                   .strip_y() + labs(subtitle=NULL), "D")

left_strip_w21  <- 0.65
col_header_h21  <- 1.02
cell_w21_photo  <- length(photo_fams)*1.85 + 2.35
cell_w21_nmet   <- length(nmet_fams)*1.85 + 2.35
cell_h21        <- 5.55

hdr_photo21 <- .group_col_header21("Photosynthesis protein groups")
hdr_nmet21 <- .group_col_header21("N metabolism protein groups")
cmp_light21 <- .comparison_row_strip21("LLS vs HLS")
cmp_nitro21 <- .comparison_row_strip21("N0.32 vs N1.28")

# Combine the four main treatment panels into Figure 9.
fig9_design21 <- "
#AB
CDE
FGH
"

fig_21_core <- wrap_plots(
  A=hdr_photo21, B=hdr_nmet21,
  C=cmp_light21, D=pA21, E=pB21,
  F=cmp_nitro21, G=pC21, H=pD21,
  design=fig9_design21,
  widths=c(left_strip_w21, cell_w21_photo, cell_w21_nmet),
  heights=c(col_header_h21, cell_h21, cell_h21),
  guides="collect") &
  theme(legend.position="bottom",plot.margin=margin(0),legend.margin=margin(0),legend.box.margin=margin(0),legend.box.spacing=unit(0,"pt"))

fig_21 <- add_common_y_title(
  fig_21_core, "Mean NSAF ± SE (stacked by partner)", text_size=32, width=0.04)
fig_21_w <- left_strip_w21 + cell_w21_photo + cell_w21_nmet + 0.06
fig_21_h <- col_header_h21 + cell_h21*2 + 0.06

save_svg(file.path(OUTDIR,"Figure9_protein_groups"), fig_21,
         width=fig_21_w, height=fig_21_h)
cat("Saved: Figure9_protein_groups.svg\n")
cat("       Left strips: LLS vs HLS | N0.32 vs N1.28\n")
cat("       Column headers show the two sets of functional groups\n")
cat("       Bar order is explained in the figure caption.\n")
cat("       Shared x-axis labels: lower row supplies labels for upper row\n\n")

# Save a JPG and PDF as well as the SVG.
ggsave(file.path(OUTDIR, "Figure9_protein_groups.jpg"), fig_21,
       width = fig_21_w, height = fig_21_h, dpi = 300, bg = "white")
ggsave(file.path(OUTDIR, "Figure9_protein_groups.pdf"), fig_21,
       width = fig_21_w, height = fig_21_h, device = cairo_pdf, bg = "white")
