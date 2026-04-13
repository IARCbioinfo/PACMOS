.save_plot <- function(path, plot, width, height) {
  ggplot2::ggsave(path, plot, width = width, height = height)
  svg_path <- sub("\\.pdf$", ".svg", path)
  svg(svg_path, width = width, height = height)
  print(plot)
  dev.off()
  message("Saved: ", basename(path), " + .svg")
}

generate_LOO_test_sample_QC <- function(
    base_dir,
    global_reference,
    model_folder    = "inputs",
    qc_dir          = file.path(base_dir, "QC_summary"),
    reference_axes  = NULL,
    axis_rename     = NULL,
    axis_cols       = NULL,
    id_col          = NULL
){

  # ---- load reference LFs -----------------------------------------

  ref <- read.csv(global_reference, check.names = FALSE)

  if (is.null(id_col))
    id_col <- if ("Sample" %in% colnames(ref)) "Sample" else colnames(ref)[1]

  rownames(ref) <- as.character(ref[[id_col]])

  if (is.null(reference_axes)) {
    reference_axes <- colnames(ref)[
      sapply(ref, is.numeric) & colnames(ref) != id_col
    ]
  }

  if (!all(reference_axes %in% colnames(ref)))
    stop("Some reference_axes not found in global_reference.")

  # ---- loop over LOO folds ----------------------------------------

  sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  results     <- list()

  for (d in sample_dirs) {

    sample_name <- basename(d)
    pred_file   <- file.path(
      d, model_folder, "plots",
      paste0(sample_name, "_query_sample_LFs.csv")
    )

    if (!file.exists(pred_file)) {
      message("Skipping ", sample_name, ": no _aligned_LFs.csv found.")
      next
    }

    if (!(sample_name %in% rownames(ref))) {
      message("Skipping ", sample_name, ": not found in global_reference.")
      next
    }

    pred      <- read.csv(pred_file, check.names = FALSE)
    pred_vals <- pred[pred$Sample == sample_name, reference_axes, drop = FALSE]
    true_vals <- ref[sample_name, reference_axes, drop = FALSE]

    for (ax in reference_axes) {
      results[[length(results) + 1]] <- data.frame(
        Sample           = sample_name,
        Axis             = ax,
        True             = as.numeric(true_vals[[ax]]),
        Predicted        = as.numeric(pred_vals[[ax]]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(results))
    stop("No predictions found. Check base_dir, model_folder and reference_axes.")

  df <- dplyr::bind_rows(results)
  message("Collected predictions for ", length(unique(df$Sample)), " samples.")

  # ---- optional axis renaming -------------------------------------

  if (!is.null(axis_rename)) {
    df$Axis <- dplyr::recode(df$Axis, !!!axis_rename)
  }

  # ---- per-axis metrics -------------------------------------------

  metrics <- df %>%
    dplyr::group_by(Axis) %>%
    dplyr::summarise(
      Correlation = cor(True, Predicted, use = "pairwise.complete.obs"),
      RMSE        = sqrt(mean((True - Predicted)^2, na.rm = TRUE)),
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      label = paste0(
        "r = ",     formatC(Correlation, digits = 3, format = "f"),
        "\nRMSE = ", formatC(RMSE,        digits = 3, format = "f")
      )
    )

  # ---- save outputs -----------------------------------------------

  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  pred_csv <- file.path(
    qc_dir,
    paste0(model_folder, "_test_sample_prediction_table.csv")
  )
  write.csv(df, pred_csv, row.names = FALSE)
  message("Prediction table saved to: ", pred_csv)

  # ---- auto colour palette if not provided ------------------------

  if (is.null(axis_cols)) {
    axes      <- sort(unique(df$Axis))
    cols      <- grDevices::hcl.colors(length(axes), palette = "Dark3")
    axis_cols <- stats::setNames(cols, axes)
  }

  # ---- scatter plot -----------------------------------------------

  #axis_limits <- df %>%
  #  dplyr::group_by(Axis) %>%
  #  dplyr::summarise(
  ##    x_min = min(True,      na.rm = TRUE),
  #    y_max = max(Predicted, na.rm = TRUE),
  #    .groups = "drop"
  #  )

  #metrics <- dplyr::left_join(metrics, axis_limits, by = "Axis")
  global_min <- min(c(df$True, df$Predicted), na.rm = TRUE)
  global_max <- max(c(df$True, df$Predicted), na.rm = TRUE)

  metrics$x_min <- global_min
  metrics$y_max <- global_max

  n_axes  <- length(unique(df$Axis))
  p_width <- max(4 * n_axes, 6)

  n_axes  <- length(unique(df$Axis))
  p_width <- max(4 * n_axes, 6)

  global_min <- min(c(df$True, df$Predicted), na.rm = TRUE)
  global_max <- max(c(df$True, df$Predicted), na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = True, y = Predicted)) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", colour = "grey60", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      shape = 21, size = 2.5, alpha = 0.7,
      ggplot2::aes(fill = Axis),
      colour = "black", stroke = 0.3
    ) +
    ggplot2::geom_smooth(
      method = "lm", se = TRUE,
      colour = "grey30", fill = "grey70",
      alpha = 0.15, linewidth = 0.7
    ) +
    ggplot2::geom_label(
      data = metrics,
      ggplot2::aes(x = x_min, y = y_max, label = label),
      hjust       = 0,      # left-align
      vjust       = 1,      # top-align
      size        = 3.2,
      linewidth   = 0.3,
      fill        = "white",
      alpha       = 0.85,
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = axis_cols) +
    ggplot2::facet_wrap(~Axis) +
    ggplot2::coord_equal(
      xlim = c(global_min, global_max),
      ylim = c(global_min, global_max)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#EEF2F7",
                                               colour = "grey70"),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none"
    ) +
    ggplot2::labs(
      title = "Leave One Out Test Sample Projection Accuracy",
      x     = "True reference value",
      y     = "Predicted value"
    )

  pred_pdf <- file.path(
    qc_dir,
    paste0(model_folder, "_test_sample_prediction_summary.pdf")
  )
  .save_plot(pred_pdf, p, width = p_width, height = 4)

  message("=====================================")
  message("Test sample QC summary saved in: ", qc_dir)
  message("=====================================")

  invisible(df)
}

generate_LOO_test_sample_QC(
  base_dir         = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO",
  global_reference = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/lungNEN_references/lungNEN_LFs.csv",
  model_folder     = "train_test_Meth_only",
  reference_axes   = c(
    "Factor1",
    "Factor2",
    "Factor5"
  ),
  axis_rename = c(
    "Factor1"        = "Factor1",
    "Factor2" = "Factor2",
    "Factor5" = "Factor5"
  ),
  axis_cols = c(
    "Factor1"        = "#B81330",
    "Factor2" =  "#79A960",
    "Factor5" = "#58839D"
  )
)
