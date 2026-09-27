# Fig. S1: canopy openness at the lower- and higher-light trees.
# This uses only the two light-treatment photographs, not the donor tree.
# The blue channel marks sky; 89 rings with 720 points each are averaged.
# Run from the repository root. Place the two PNG photographs there.

photos <- c(
  lower_light_tree = "lower_light_tree.png",
  higher_light_tree = "higher_light_tree.png"
)
labels <- c(
  lower_light_tree = "Lower-light tree",
  higher_light_tree = "Higher-light tree"
)
# Pixel values at or above 0.40 in the blue channel are treated as sky.
threshold <- 0.40
dir.create("output", showWarnings = FALSE)

# Sample sky pixels along each ring and average the 89 ring fractions.
analyse_photo <- function(path) {
  image <- png::readPNG(path)
  if (length(dim(image)) != 3L || dim(image)[3] < 3L) {
    stop("A colour PNG is required: ", path)
  }

  sky <- (image[, , 3] >= threshold) * 1
  nr <- nrow(sky)
  nc <- ncol(sky)
  radius <- min(nr, nc) * 0.95 / 2
  angles <- seq(0, 2 * pi, length.out = 721)[-721]
  gap <- numeric(89)

  for (ring in seq_len(89)) {
    x <- round(nc / 2 + cos(angles) * radius * ring / 89)
    y <- round(nr / 2 + sin(angles) * radius * ring / 89)
    x <- pmax(1, pmin(nc, x))
    y <- pmax(1, pmin(nr, y))
    gap[ring] <- mean(sky[cbind(y, x)])
  }

  # Show transparent pixels against black, as in the original script.
  if (dim(image)[3] == 4L) {
    for (channel in 1:3) {
      image[, , channel] <- image[, , channel] * image[, , 4]
    }
  }

  list(photo = image[, , 1:3, drop = FALSE],
       binary = array(rep(sky, 3), dim = c(nr, nc, 3)),
       openness = 100 * mean(gap))
}

# Analyze the two light-treatment trees.
results <- lapply(photos, analyse_photo)
# Save the exact openness percentages; the figure labels show one decimal.
summary <- data.frame(
  Tree = unname(labels),
  Canopy_openness_percent = vapply(results, function(x) x$openness, numeric(1))
)
print(summary, row.names = FALSE)
write.csv(summary, "output/canopy_openness.csv", row.names = FALSE)

# Draw an image without stretching its original proportions.
draw_image <- function(image) {
  aspect <- dim(image)[2] / dim(image)[1]
  width <- min(0.96, 0.96 * aspect)
  height <- min(0.96, 0.96 / aspect)
  grid::grid.raster(image, width = grid::unit(width, "snpc"),
                    height = grid::unit(height, "snpc"), interpolate = FALSE)
}

# Put the original photographs above their thresholded versions.
draw_figure <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(4, 2,
      heights = grid::unit(c(0.105, 0.075, 0.41, 0.41), "npc"))
  ))

  grid::grid.text(
    "Hemispherical canopy photographs\nand binarised images",
    gp = grid::gpar(fontface = "bold", fontsize = 26, lineheight = 1.05),
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
  )

  for (column in seq_along(results)) {
    item <- results[[column]]
    grid::pushViewport(grid::viewport(layout.pos.row = 2,
                                      layout.pos.col = column))
    grid::grid.rect(width = 0.94, height = 0.82,
                    gp = grid::gpar(fill = "white", col = "black", lwd = 1))
    grid::grid.text(
      labels[[column]],
      gp = grid::gpar(fontface = "bold", fontsize = 24)
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 3,
      layout.pos.col = column,
      layout = grid::grid.layout(2, 1,
        heights = grid::unit(c(0.12, 0.88), "npc"))))
    grid::grid.text(sprintf("Canopy openness: %.1f%%", item$openness),
      gp = grid::gpar(fontface = "bold", fontsize = 23),
      vp = grid::viewport(layout.pos.row = 1))
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    draw_image(item$photo)
    grid::popViewport()
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 4,
      layout.pos.col = column,
      layout = grid::grid.layout(2, 1,
        heights = grid::unit(c(0.17, 0.83), "npc"))))
    grid::grid.text("Binarised image\nBlue-channel threshold: 0.40",
      gp = grid::gpar(fontface = "bold", fontsize = 21, lineheight = 1.05),
      vp = grid::viewport(layout.pos.row = 1))
    grid::pushViewport(grid::viewport(layout.pos.row = 2))
    draw_image(item$binary)
    grid::popViewport()
    grid::popViewport()
  }
  grid::popViewport()
}

# Save Fig. S1 as PNG and SVG.
png("output/FigureS1_canopy_openness.png", width = 10, height = 10.5,
    units = "in", res = 300, bg = "white")
draw_figure()
dev.off()

svg("output/FigureS1_canopy_openness.svg", width = 10, height = 10.5,
    bg = "white")
draw_figure()
dev.off()
