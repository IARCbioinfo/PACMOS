generate_LOO_archetype_QC <- function(
    base_dir,
    model_folder,
    global_file,
    sample_pattern   = "",
    archetype_map    = NULL,
    archetype_labels = NULL,
    archetype_cols   = NULL,
    id_col           = "ID",
    prefix           = "",
    qc_dir           = file.path(base_dir, "QC_summary", "fuzzy")
) {

  default_cols <- c(
    "#B81330", "#58839D", "#79A960", "#F4A261", "#2A9D8F",
    "#9B59B6", "#E67E22", "#1ABC9C", "#E74C3C", "#2C3E50"
  )

  # ── load global weights ───────────────────────────────────────────────────
  global <- read.table(
    global_file,
    header           = TRUE,
    sep              = "\t",
    stringsAsFactors = FALSE,
    check.names      = TRUE
  )

  if (id_col %in% colnames(global)) {
    colnames(global)[colnames(global) == id_col] <- "Sample"
  } else if (!"Sample" %in% colnames(global)) {
    stop("global_file must contain a column named '", id_col, "' or 'Sample'.")
  }

  # ──  filter sample directories ─────────────────────────────────
  all_subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- all_subdirs[grepl(sample_pattern, basename(all_subdirs))]

  if (!length(sample_dirs))
    stop("No subdirectories matching '", sample_pattern, "' in: ", base_dir)

  message("Found ", length(sample_dirs), " sample folder(s).")

  # ── auto-build archetype_map if not supplied ──────────────────────────────
  if (is.null(archetype_map)) {

    recomputed_cols <- NULL
    for (d in sample_dirs) {
      wf <- file.path(d, model_folder, "plots",
                      paste0(prefix, basename(d), "_archetype_weights.csv"))
      if (file.exists(wf)) {
        recomputed_cols <- setdiff(
          colnames(read.csv(wf, nrows = 1, check.names = TRUE,
                            stringsAsFactors = FALSE)),
          "Sample"
        )
        break
      }
    }

    if (is.null(recomputed_cols))
      stop("No _archetype_weights.csv found under model_folder: '",
           model_folder, "'. Check prefix parameter.")

    global_arch_cols <- setdiff(colnames(global), "Sample")
    shared           <- intersect(global_arch_cols, recomputed_cols)

    if (length(shared) > 0) {
      archetype_map <- setNames(shared, shared)
    } else {
      stripped_global <- gsub("^.*\\.", "", global_arch_cols)
      matched         <- match(recomputed_cols, stripped_global)
      valid           <- !is.na(matched)
      archetype_map   <- setNames(
        recomputed_cols[valid],
        global_arch_cols[matched[valid]]
      )
    }

    if (!length(archetype_map))
      stop("Cannot automatically match archetype columns. ",
           "Please supply archetype_map explicitly.")

    message("Auto-detected archetype mapping:")
    for (nm in names(archetype_map))
      message("  '", nm, "'  -->  '", archetype_map[[nm]], "'")

  } else {
    message("Archetype mapping (user-supplied):")
    for (nm in names(archetype_map))
      message("  '", nm, "'  -->  '", archetype_map[[nm]], "'")
  }

  # ── loop over LOO folds ───────────────────────────────────────────────────
  all_metrics <- list()

  for (d in sample_dirs) {

    sample_name  <- basename(d)

    # all_samples file (has both train + test rows); fall back to weights only
    weights_file <- file.path(d, model_folder, "plots",
                              paste0(prefix, sample_name,
                                     "_archetype_weights_all_samples.csv"))
    if (!file.exists(weights_file)) {
      weights_file <- file.path(d, model_folder, "plots",
                                paste0(prefix, sample_name,
                                       "_archetype_weights.csv"))
    }

    if (!file.exists(weights_file)) {
      message("  [SKIP] No weights CSV for: ", sample_name)
      next
    }

    # exclude query sample — training samples only
    recomputed <- read.csv(weights_file, check.names = FALSE,   # ← was TRUE
                           stringsAsFactors = FALSE) %>%
      dplyr::filter(Sample != sample_name)

    common_samples <- intersect(global$Sample, recomputed$Sample)
    if (length(common_samples) < 2) {
      message("  [SKIP] Too few common samples for: ", sample_name)
      next
    }

    g <- global[match(common_samples, global$Sample), ]
    r <- recomputed[match(common_samples, recomputed$Sample), ]

    for (g_col in names(archetype_map)) {

      r_col <- archetype_map[[g_col]]
      if (!(g_col %in% colnames(g)) || !(r_col %in% colnames(r))) next

      x    <- g[[g_col]]
      y    <- r[[r_col]]
      keep <- complete.cases(x, y)
      if (sum(keep) < 2) next

      all_metrics[[length(all_metrics) + 1]] <- data.frame(
        Fold        = sample_name,
        Archetype   = r_col,
        Correlation = cor(x[keep], y[keep]),
        RMSE        = sqrt(mean((x[keep] - y[keep])^2)),
        stringsAsFactors = FALSE
      )
    }

    message("  Processed: ", sample_name)
  }

  if (!length(all_metrics))
    stop("No metrics computed — check base_dir, model_folder, and column names.")

  metrics <- dplyr::bind_rows(all_metrics)

  # ── apply archetype labels ────────────────────────────────────────────────
  if (is.null(archetype_labels)) {
    unique_arch      <- unique(metrics$Archetype)
    archetype_labels <- setNames(
      gsub("[._-]", " ", unique_arch),
      unique_arch
    )
  }

  metrics <- metrics %>%
    dplyr::mutate(Archetype = dplyr::recode(Archetype, !!!archetype_labels))

  # ── apply colours ─────────────────────────────────────────────────────────
  final_labels <- unique(metrics$Archetype)

  if (is.null(archetype_cols)) {
    archetype_cols <- setNames(
      default_cols[seq_along(final_labels)],
      final_labels
    )
  } else {
    missing_labs <- setdiff(final_labels, names(archetype_cols))
    if (length(missing_labs) > 0)
      archetype_cols <- c(
        archetype_cols,
        setNames(default_cols[length(archetype_cols) + seq_along(missing_labs)],
                 missing_labs)
      )
  }

  # ── output directory ──────────────────────────────────────────────────────
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  write.csv(
    metrics,
    file.path(qc_dir, paste0(model_folder, "_archetype_QC_summary.csv")),
    row.names = FALSE
  )

  # ── correlation line plot ─────────────────────────────────────────────────
  sample_order <- metrics %>%
    dplyr::group_by(Fold) %>%
    dplyr::summarise(mean_r = mean(Correlation, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_r) %>%
    dplyr::pull(Fold)

  index_map <- data.frame(Index = seq_along(sample_order),
                          Sample = sample_order)
  write.csv(
    index_map,
    file.path(qc_dir, paste0(model_folder, "_archetype_sample_index_map.csv")),
    row.names = FALSE
  )

  metrics$Fold        <- factor(metrics$Fold, levels = sample_order)
  metrics$SampleIndex <- as.integer(metrics$Fold)

  y_min <- floor(min(metrics$Correlation, na.rm = TRUE) * 10) / 10

  p_cor <- ggplot2::ggplot(
    metrics,
    ggplot2::aes(x = SampleIndex, y = Correlation,
                 colour = Archetype, fill = Archetype, group = Archetype)
  ) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.3, linetype = "dashed") +
    ggplot2::geom_point(shape = 21, size = 2, alpha = 0.8,
                        colour = "black", stroke = 0) +
    ggplot2::scale_fill_manual(values   = archetype_cols) +
    ggplot2::scale_colour_manual(values = archetype_cols) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(metrics$SampleIndex), by = 10),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(y_min, 1.0, 0.1),
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
      title    = paste0("LOO Archetype QC - Pearson r  [", model_folder, "]"),
      subtitle = paste0("Samples sorted by mean r  |  ",
                        "index mapped to _archetype_sample_index_map.csv"),
      x        = "Sample index",
      y        = "Pearson r"
    )

  ggplot2::ggsave(
    file.path(qc_dir, paste0(model_folder,
                             "_archetype_correlation_training_sample_only.pdf")),
    p_cor, width = 16, height = 5
  )

  # ── RMSE violin + boxplot ─────────────────────────────────────────────────
  p_rmse <- ggplot2::ggplot(
    metrics,
    ggplot2::aes(x = Archetype, y = RMSE, fill = Archetype)
  ) +
    ggplot2::geom_violin(alpha = 0.35, colour = "grey40", linewidth = 0.5,
                         trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA,
                          colour = "grey30", fill = "white", linewidth = 0.5) +
    ggplot2::geom_jitter(shape = 21, width = 0.07, size = 2.2,
                         ggplot2::aes(fill = Archetype),
                         colour = "black", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = archetype_cols) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "none"
    ) +
    ggplot2::labs(
      title    = paste0("LOO Archetype QC - RMSE  [", model_folder, "]"),
      subtitle = "Violin + boxplot + individual LOO folds",
      x        = NULL,
      y        = "RMSE"
    )

  ggplot2::ggsave(
    file.path(qc_dir, paste0(model_folder,
                             "_archetype_rmse_training_sample_only.pdf")),
    p_rmse, width = 6, height = 5
  )

  message("=====================================")
  message("Archetype QC written to: ", qc_dir)

  invisible(metrics)
}


generate_LOO_archetype_QC(
  base_dir         = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO",
  model_folder     = "train_test_Meth_only",
  global_file      = "inst/extdata/lungNEN_references/lungNEN_LF_K4_without_label.txt",
  prefix           = "",

  archetype_map = c(
    "K4_sce_proportion" = "SC-enriched",
    "K4_A1_proportion"  = "Ca A1",
    "K4_B_proportion"   = "Ca B",
    "K4_A2_proportion"  = "Ca A2"
  ),

  archetype_labels = c(
    "SC-enriched"            = "SC-enriched",
    "Ca A1" = "Ca A1",
    "Ca B"                    = "Ca B",
    "Ca A2" = "Ca A2"
  )
)
