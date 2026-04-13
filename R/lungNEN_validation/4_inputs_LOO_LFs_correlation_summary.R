plot_alignment_QC_summary <- function(
    aggregated_csv,
    out_dir      = dirname(aggregated_csv),
    out_prefix   = "allsamples",
    axis_cols    = NULL,
    axis_rename  = NULL
) {

  # ---- read aggregated CSV -------------------------------------------------

  if (!file.exists(aggregated_csv))
    stop("aggregated_csv not found: ", aggregated_csv)

  metrics <- read.csv(aggregated_csv, check.names = FALSE)
  message("Read ", nrow(metrics), " rows from: ", aggregated_csv)

  required_cols <- c("Model", "Sample", "Axis", "Correlation", "RMSE")
  missing_cols  <- setdiff(required_cols, colnames(metrics))
  if (length(missing_cols))
    stop("Missing columns in aggregated_csv: ",
         paste(missing_cols, collapse = ", "))

  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # ---- optional axis renaming ----------------------------------------------

  if (!is.null(axis_rename)) {
    old_names <- names(axis_rename)
    new_names <- unname(axis_rename)
    for (k in seq_along(old_names))
      metrics$Axis[metrics$Axis == old_names[k]] <- new_names[k]

    # propagate rename into axis_cols if provided
    if (!is.null(axis_cols)) {
      idx <- match(old_names, names(axis_cols))
      idx <- idx[!is.na(idx)]
      if (length(idx))
        names(axis_cols)[idx] <- new_names[!is.na(match(old_names, names(axis_cols)))]
    }
  }

  # ---- auto colour palette -------------------------------------------------

  if (is.null(axis_cols)) {
    axes      <- sort(unique(metrics$Axis))
    cols      <- grDevices::hcl.colors(length(axes), palette = "Dark3")
    axis_cols <- stats::setNames(cols, axes)
  }

  # ---- sample order sorted by mean Pearson r -------------------------------

  sample_order <- metrics |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(mean_r = mean(Correlation, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_r) |>
    dplyr::pull(Sample)

  index_map <- data.frame(
    Index  = seq_along(sample_order),
    Sample = sample_order
  )

  index_file <- file.path(out_dir, paste0(out_prefix, "_sample_index_map.csv"))
  write.csv(index_map, index_file, row.names = FALSE)
  message("Sample index map saved to: ", index_file)

  metrics$Sample      <- factor(metrics$Sample, levels = sample_order)
  metrics$SampleIndex <- as.integer(metrics$Sample)

  # ---- 1. Pearson r line + dot plot ----------------------------------------

  y_min <- floor(min(metrics$Correlation, na.rm = TRUE) * 10) / 10

  p_cor <- ggplot2::ggplot(
    metrics,
    ggplot2::aes(
      x      = SampleIndex,
      y      = Correlation,
      colour = Axis,
      fill   = Axis,
      group  = Axis
    )
  ) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.3, linetype = "dashed") +
    ggplot2::geom_point(shape = 21, size = 2, alpha = 0.8,
                        colour = "black", stroke = 0) +
    ggplot2::scale_fill_manual(values   = axis_cols) +
    ggplot2::scale_colour_manual(values = axis_cols) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(metrics$SampleIndex), by = 10),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(y_min, 1.0, by = 0.1),
      limits = c(y_min, 1.02)
    ) +
    ggplot2::guides(
      fill   = ggplot2::guide_legend(override.aes = list(size = 3)),
      colour = "none"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(size = 6),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "top",
      legend.title       = ggplot2::element_blank(),
      legend.key.width   = ggplot2::unit(1.0, "cm")
    ) +
    ggplot2::labs(
      title = "All-Samples MOFA Alignment \u2014 Pearson r per Factor",
      x     = "Sample index (sorted by mean r)",
      y     = "Pearson r"
    )

  cor_pdf <- file.path(out_dir,
                       paste0(out_prefix, "_LF_correlation_all_samples.pdf"))
  ggplot2::ggsave(cor_pdf, p_cor, width = 16, height = 5)
  message("Correlation plot saved to: ", cor_pdf)

  # ---- 2. RMSE violin + boxplot --------------------------------------------

  p_rmse <- ggplot2::ggplot(
    metrics,
    ggplot2::aes(x = Axis, y = RMSE, fill = Axis)
  ) +
    ggplot2::geom_violin(alpha = 0.35, colour = "grey40",
                         linewidth = 0.5, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          colour = "grey30", fill = "white",
                          linewidth = 0.5) +
    ggplot2::geom_jitter(shape = 21, width = 0.07, size = 2.2,
                         ggplot2::aes(fill = Axis),
                         colour = "black", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = axis_cols) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "none"
    ) +
    ggplot2::labs(
      title = "All-Samples MOFA Alignment \u2014 RMSE Distribution per Factor",
      x     = NULL,
      y     = "RMSE"
    )

  rmse_pdf <- file.path(out_dir,
                        paste0(out_prefix, "_LF_rmse_all_samples.pdf"))
  ggplot2::ggsave(rmse_pdf, p_rmse, width = 10, height = 8)
  message("RMSE plot saved to: ", rmse_pdf)

  invisible(metrics)
}

plot_alignment_QC_summary(
  aggregated_csv = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO/all_samples_quality_check_metrics.csv",
  out_dir        = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO/QC_summary",
  out_prefix     = "lungNEN_input_samples_new_ref_",
  axis_rename    = c(
    "Factor1"            = "Factor1",
    "Factor2"        = "Factor2",
    "Factor5" = "Factor5"
  ),
  axis_cols = c(
    "Factor1"            = "#58839D",
    "Factor2"        = "#B81330",
    "Factor5" = "#79A960"  )
)
