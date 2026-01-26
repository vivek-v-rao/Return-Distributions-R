# plot_utils.r
#
# Small utilities for plotting density fits.

#' Build an x-grid for density plots using quantile range and padding.
density_x_grid <- function(x, lower = 0.01, upper = 0.99, pad_frac = 0.05, n = 400) {
  q <- quantile(x, probs = c(lower, upper), na.rm = TRUE, names = FALSE)
  x_min <- q[1]
  x_max <- q[2]
  pad <- pad_frac * (x_max - x_min)
  if (!is.finite(pad) || pad == 0) pad <- 1
  seq(x_min - pad, x_max + pad, length.out = n)
}

#' Plot density curves and save to a PNG file.
plot_density_fits <- function(x_grid, dens_list, outfile, title, xlab = "returns", ylab = "density") {
  y_max <- max(unlist(dens_list), na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) return(FALSE)

  colors <- c(
    Normal = "#1f77b4",
    T = "#ff7f0e",
    SkewNormal = "#2ca02c",
    SkewT = "#d62728",
    FSSkewNormal = "#17becf",
    FSSkewT = "#8c564b",
    NCT = "#bcbd22",
    JFSkewT = "#e377c2",
    Laplace = "#7f7f7f",
    ALaplace = "#aec7e8",
    GED = "#ffbb78",
    SGED = "#98df8a",
    Sech = "#c5b0d5",
    GSH = "#8c6d31",
    SGSH = "#c49c94",
    NEFGHS = "#f7b6d2",
    Logistic = "#9edae5",
    EGB2 = "#dbdb8d",
    NIG = "#c7c7c7",
    Hyperbolic = "#c7b0ad",
    GH = "#c49c94",
    VG = "#ff9896",
    Champernowne = "#b5cf6b",
    NormalLaplace = "#8dd3c7",
    GT = "#ff9896",
    SGT = "#c5b0d5",
    LogCondens = "#ffbb78",
    LogSpline = "#1f77b4",
    Cauchy = "#aec7e8",
    KDE = "#9467bd"
  )
  line_types <- c(
    Normal = "solid",
    T = "solid",
    SkewNormal = "solid",
    SkewT = "solid",
    FSSkewNormal = "solid",
    FSSkewT = "solid",
    NCT = "solid",
    JFSkewT = "solid",
    Laplace = "solid",
    ALaplace = "solid",
    GED = "solid",
    SGED = "solid",
    Sech = "solid",
    GSH = "solid",
    SGSH = "solid",
    NEFGHS = "solid",
    Logistic = "solid",
    EGB2 = "solid",
    NIG = "solid",
    Hyperbolic = "solid",
    GH = "solid",
    VG = "solid",
    Champernowne = "solid",
    NormalLaplace = "solid",
    GT = "solid",
    SGT = "solid",
    LogCondens = "dotted",
    LogSpline = "dotdash",
    Cauchy = "solid",
    KDE = "dashed"
  )

  png(outfile, width = 900, height = 600)
  plot(x_grid, dens_list[[1]], type = "n",
       main = title, xlab = xlab, ylab = ylab, ylim = c(0, y_max))
  for (model_name in names(dens_list)) {
    col <- if (!is.null(colors[[model_name]])) colors[[model_name]] else "black"
    lty <- if (!is.null(line_types[[model_name]])) line_types[[model_name]] else "solid"
    lines(x_grid, dens_list[[model_name]], col = col, lwd = 2, lty = lty)
  }
  legend("topright", legend = names(dens_list),
         col = colors[names(dens_list)], lwd = 2, lty = line_types[names(dens_list)])
  dev.off()
  TRUE
}
