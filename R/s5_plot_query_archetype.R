#' Plot archetype weight projections
#'
#' @description
#' Visualizes archetype mixture weights for query and reference samples
#' using ternary projections and composition plots.
#'
#' @details
#' For each sample folder in `models_dir`, this function reads
#' `<sample_id>_archetype_weights_all_samples.csv` produced by
#' `infer_fuzzy_weights()` and generates a  PDF containing:
#'
#' \itemize{
#'   \item Ternary plots showing sample positions in archetype space
#'   \item Highlighted query sample projection
#'   \item A stacked bar plot showing archetype composition of the query sample
#' }
#'
#' If the number of archetypes exceeds three, multiple ternary panels are
#' generated for all combinations of three archetypes.
#'
#' Output PDFs are saved alongside the input CSV files in the `plots/` directory.
#'
#' @name plot_fuzzy_query_sample
#'
#' @param models_dir        Character. Root directory.
#' @param matrices_subdir Character. Folder name where `.RData` files are stored.
#' @param sample_pattern  Regex to filter sample folder names. Default \code{""} (all).
#' @param prefix          Optional prefix for output PDF filenames.
#'
#' @return Invisibly returns a named logical vector (sample processed = TRUE).
#'
#'
#' @examples
#' plot_fuzzy_query_sample(
#'   models_dir = "query_output",
#'   matrices_subdir = "inputs"
#' )
#'
#' @export
plot_fuzzy_query_sample <- function(
    models_dir,
    matrices_subdir ,
    sample_pattern   = "",
    prefix           = ""
) {

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (!requireNamespace("Ternary", quietly = TRUE))
    stop("Package 'Ternary' is required. Install with: install.packages('Ternary')")

  if (!requireNamespace("scales", quietly = TRUE))
    stop("Package 'scales' is required. Install with: install.packages('scales')")

  all_vcols <- c(
    "#B81330", "#58839D", "#79A960", "#F4A261", "#2A9D8F",
    "#9B59B6", "#E67E22", "#1ABC9C", "#E74C3C", "#2C3E50"
  )

  # ---- helpers -------------------------------------------------------------

  ##  ternary plots
  .base_ternary <- function() {
    Ternary::TernaryPlot(
      point            = "up",
      atip             = "",
      btip             = "",
      ctip             = "",
      lab.cex          = 0.0,
      grid.lines       = 0,
      grid.minor.lines = 0,
      grid.lty         = "blank",
      grid.col         = "white",
      col              = grDevices::rgb(1, 1, 1),
      axis.col         = "black",
      ticks.col        = "black",
      axis.rotate      = FALSE,
      padding          = 0.12
    )
  }


    .draw_vertices <- function(idx_trio, vertex_labels, vcols) {
    i <- idx_trio[1]; j <- idx_trio[2]; k <- idx_trio[3]


    data_points <- list(
      c(0,   0,   255),   # archetype i -> bottom-right
      c(255, 0,   0  ),   # archetype j -> top
      c(0,   255, 0  )    # archetype k -> bottom-left
    )
    xy <- Ternary::CoordinatesToXY(data_points)

    s <- sqrt(0.1^2 / 2)
    disc_x_off <- c(-s,    0,    s  )
    disc_y_off <- c(-s,    0.1, -s  )

    disc_x <- xy[1, ] + disc_x_off
    disc_y <- xy[2, ] + disc_y_off

    graphics::points(
      disc_x, disc_y,
      pch = 19, cex = 3,
      col = scales::alpha(vcols[c(i, j, k)], 0.7)

    )
    graphics::text(
      disc_x, disc_y,
      labels = as.character(c(i, j, k)),
      cex = 0.80, font = 2, col = "white"
    )

    lab_x_off <- c(-sqrt(0.04^2),  0,    sqrt(0.08^2))
    lab_y_off <- c(-sqrt(0.2^2),   0.2, -sqrt(0.2^2) )

    graphics::text(
      xy[1, ] + lab_x_off,
      xy[2, ] + lab_y_off,
      labels = vertex_labels[c(i, j, k)],
      cex = 0.80,
      font = 2,
      xpd = TRUE
    )
  }

  ## One ternary panel for a trio of archetypes (used when n > 3).
  ## Column order: k (3rd) -> top, i (1st) -> bottom-left, j (2nd) -> bottom-right
  .draw_ternary_panel <- function(idx_trio, weights, samples,
                                  vertex_labels, vcols, query_sample) {
    i <- idx_trio[1]; j <- idx_trio[2]; k <- idx_trio[3]

    # column order for AddToTernary: must match a-tip, b-tip, c-tip
    # a-tip = j (top), b-tip = k (bottom-left), c-tip = i (bottom-right)
    tern_order <- c(j, k, i)

    w_sub  <- weights[, tern_order, drop = FALSE]
    rs     <- rowSums(w_sub)
    rs[rs == 0] <- 1
    w_norm <- w_sub / rs

    hi_idx <- which(samples == query_sample)
    others <- setdiff(seq_len(nrow(w_norm)), hi_idx)

    .base_ternary()

    if (length(others) > 0)
      Ternary::AddToTernary(
        graphics::points,
        w_norm[others, , drop = FALSE],
        pch = 21, cex = 2,
        bg  = "grey75",
        col = "grey40"
      )

    if (length(hi_idx) > 0)
      Ternary::AddToTernary(
        graphics::points,
        w_norm[hi_idx, , drop = FALSE],
        pch = 21, cex =3,
        bg  = "red", col = "black"
      )

    .draw_vertices(idx_trio, vertex_labels, vcols)

    mtext(
      paste(vertex_labels[c(i, j, k)], collapse = " x "),
      side = 1, line = 0.2, cex = 0.48, col = "grey50"
    )
  }

  ## Full simplex for n <= 3.
  ## For n == 3: columns reordered as c(3, 1, 2) -> top, bottom-left, bottom-right.
  .draw_simple_panel <- function(weights, samples, vertex_labels,
                                 vcols, query_sample) {
    n      <- ncol(weights)
    hi_idx <- which(samples == query_sample)
    others <- setdiff(seq_len(nrow(weights)), hi_idx)

    if (n == 1) {
      par(mar = c(2, 2, 2, 2))
      plot(0, 0, pch = 19, col = scales::alpha(vcols[1], 0.90),
           cex = 2.8, axes = FALSE, xlab = "", ylab = "")
      graphics::text(0, 0, labels = "1", col = "white", cex = 0.80, font = 2)
      graphics::text(0, -0.22, labels = vertex_labels[1],
                     cex = 0.80, font = 2, xpd = TRUE)
      return(invisible(NULL))
    }

    if (n == 2) {
      par(mar = c(3, 2, 2, 2))
      x_other <- if (length(others)) weights[others, 1] - 0.5 else numeric(0)
      x_hi    <- if (length(hi_idx)) weights[hi_idx, 1] - 0.5 else NULL
      plot(NULL, xlim = c(-0.85, 0.85), ylim = c(-0.35, 0.35),
           asp = 1, axes = FALSE, xlab = "", ylab = "")
      segments(-0.5, 0, 0.5, 0, col = "black", lwd = 1.2)
      if (length(x_other))
        graphics::points(x_other, rep(0, length(x_other)),
                         pch = 21, cex = 2,
                         bg  = "grey75",
                         col = "grey40")
      graphics::points(c(-0.65, 0.65), c(0, 0),
                       pch = 19, cex = 2.4,
                       col = scales::alpha(vcols[1:2], 0.90))
      graphics::text(c(-0.65, 0.65), c(0, 0),
                     labels = c("1", "2"), col = "white", cex = 0.80, font = 2)
      graphics::text(c(-0.65, 0.65), c(-0.14, -0.14),
                     labels = vertex_labels[1:2],
                     cex = 0.80, font = 2, xpd = TRUE,
                     adj = c(0.5, 1))
      if (!is.null(x_hi))
        graphics::points(x_hi, 0, pch = 21, cex = 2.4,
                         bg = "red", col = "black")
      return(invisible(NULL))
    }

    # n == 3: reorder columns -> 3rd to top, 1st to bottom-left, 2nd to bottom-right
    ordered <- c(3L, 1L, 2L)

    # n == 3: a-tip=col2, b-tip=col3, c-tip=col1  (matches reference c(2,3,1))
    .base_ternary()

    if (length(others) > 0)
      Ternary::AddToTernary(
        graphics::points,
        weights[others, c(2, 3, 1), drop = FALSE],
        pch = 21, cex = 2,
        bg  = "grey75",
        col = "grey40"
      )

    if (length(hi_idx) > 0)
      Ternary::AddToTernary(
        graphics::points,
        weights[hi_idx, c(2, 3, 1), drop = FALSE],
        pch = 21, cex = 3,
        bg  = "red", col = "black"
      )
    else
      warning("  query sample not found.")

    .draw_vertices(1:3, vertex_labels, vcols)
  }

  ## Stacked bar - query sample archetype composition
  .draw_stacked_bar <- function(weights, samples, vertex_labels,
                                vcols, query_sample) {
    idx <- which(samples == query_sample)
    if (!length(idx)) {
      warning("query sample not found - skipping bar page.")
      return(invisible(NULL))
    }

    w <- weights[idx, , drop = FALSE]
    n <- ncol(w)

    par(mar = c(6, 3, 3, 3))
    plot(NULL, xlim = c(0, 2), ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "", bty = "n")

    y_bottom <- 0
    for (a in seq_len(n)) {
      val <- w[1, a]
      if (val > 0) {
        rect(0.6, y_bottom, 1.4, y_bottom + val,
             col = vcols[a], border = "white", lwd = 0.5)
        if (val > 0.05)
          graphics::text(1, y_bottom + val / 2,
                         labels = paste0(round(val * 100), "%"),
                         col = "white", cex = 0.8, font = 2)
        y_bottom <- y_bottom + val
      }
    }

    rect(0.6, 0, 1.4, 1, col = NA, border = "red", lwd = 2)

    axis(2, at = seq(0, 1, 0.2),
         labels = paste0(seq(0, 100, 20), "%"),
         las = 2, cex.axis = 0.8, col = "grey60", col.axis = "grey30")

    mtext(query_sample, side = 1, at = 1,
          las = 1, cex = 0.85, col = "red", font = 2, line = 0.5)

    legend(x = 1, y = -0.18,
           legend  = paste0(seq_len(n), ". ", vertex_labels),
           fill    = vcols,
           border  = "white",
           ncol    = min(n, 3L),
           cex     = 0.78,
           bty     = "n",
           xpd     = TRUE,
           xjust   = 0.5,
           yjust   = 1)

    title(
      main      = paste0("Archetype composition - ", query_sample),
      cex.main  = 0.95,
      font.main = 1,
      col.main  = "grey20"
    )
  }

  # ---- sample loop ---------------------------------------------------------

  all_subdirs <- list.dirs(models_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- all_subdirs[grepl(sample_pattern, basename(all_subdirs))]

  if (!length(sample_dirs))
    stop("No subdirectories matching pattern '", sample_pattern,
         "' found in: ", models_dir)

  message("Found ", length(sample_dirs), " sample folder(s).")

  results <- stats::setNames(logical(length(sample_dirs)),
                             basename(sample_dirs))

  for (sample_dir in sample_dirs) {

    sample_id        <- basename(sample_dir)
    query_sample <- sample_id
    plots_dir        <- file.path(sample_dir, matrices_subdir, "plots")

    message("")
    message(strrep("=", 45))
    message("  ", sample_id, " / ", matrices_subdir)
    message(strrep("=", 45))

    if (!dir.exists(plots_dir)) {
      warning("No plots/ dir found: ", plots_dir, " -- skipping.")
      next
    }

    weight_files <- list.files(
      plots_dir,
      pattern    = "_archetype_weights_all_samples\\.csv$",
      full.names = TRUE
    )

    if (!length(weight_files)) {
      warning("No _archetype_weights_all_samples.csv in ", plots_dir, " -- skipping.")
      next
    }

    for (file in weight_files) {

      message("  Processing: ", basename(file))

      df            <- read.csv(file, check.names = FALSE,
                                stringsAsFactors = FALSE)
      samples       <- df[, 1]
      weights       <- as.matrix(df[, -1])
      n             <- ncol(weights)
      vertex_labels <- colnames(weights)
      vcols         <- all_vcols[seq_len(min(n, length(all_vcols)))]

      rs <- rowSums(weights)
      if (any(abs(rs - 1) > 1e-4)) {
        warning("  Rows do not sum to 1 -- normalising.")
        weights <- weights / rs
      }

      csv_id   <- sub("_archetype_weights_all_samples\\.csv$", "", basename(file))
      out_file <- file.path(plots_dir,
                            paste0(prefix, csv_id,
                                   "_archetype_projection.pdf"))

      # ---- page 1: ternary panels --------------------------------
      if (n <= 3) {

        pdf(out_file, width = 6, height = 6)
        par(mar = c(2, 2, 2, 2))
        .draw_simple_panel(weights, samples, vertex_labels,
                           vcols, query_sample)

      } else {

        combos   <- combn(n, 3, simplify = FALSE)
        n_combos <- length(combos)
        n_cols   <- min(3L, n_combos)
        n_rows   <- ceiling(n_combos / n_cols)

        total_cells     <- n_rows * n_cols
        panel_ids       <- c(seq_len(n_combos),
                             rep(NA_integer_, total_cells - n_combos))
        layout_mat      <- matrix(panel_ids, nrow = n_rows,
                                  ncol = n_cols, byrow = TRUE)
        layout_mat_safe <- layout_mat
        layout_mat_safe[is.na(layout_mat_safe)] <- 0L

        pdf(out_file,
            width  = n_cols * 3.5,
            height = n_rows * 3.5 + 0.6)

        layout(layout_mat_safe, heights = rep(3, n_rows))
        par(mar = c(1.8, 0.5, 1.8, 0.5))

        for (trio in combos)
          .draw_ternary_panel(trio, weights, samples,
                              vertex_labels, vcols, query_sample)
      }

      # ---- page 2: stacked bar ------------------------------------
      layout(matrix(1))
      par(mar = c(6, 3, 3, 3))
      .draw_stacked_bar(weights, samples, vertex_labels,
                        vcols, query_sample)

      dev.off()
      message("  Saved: ", basename(out_file))
      results[sample_id] <- TRUE
    }
  }

  n_ok   <- sum(results,  na.rm = TRUE)
  n_fail <- sum(!results, na.rm = TRUE)
  message("")
  message(strrep("=", 45))
  message("  Done. ", n_ok, " sample(s) processed",
          if (n_fail > 0) paste0(", ", n_fail, " skipped.") else ".")
  message(strrep("=", 45))

  invisible(results)
}
