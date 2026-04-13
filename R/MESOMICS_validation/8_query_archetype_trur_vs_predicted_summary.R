generate_LOO_test_sample_archetype_QC <- function(
    base_dir,
    global_file,
    model_folder,
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
    check.names      = FALSE
  )

  if (id_col %in% colnames(global)) {
    colnames(global)[colnames(global) == id_col] <- "Sample"
  } else if (!"Sample" %in% colnames(global)) {
    stop("global_file must contain a column named '", id_col,
         "' or 'Sample'.")
  }

  global_arch_cols <- setdiff(colnames(global), "Sample")

  # ── discover + filter sample directories ─────────────────────────────────
  all_subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- all_subdirs[grepl(sample_pattern, basename(all_subdirs))]

  if (!length(sample_dirs))
    stop("No subdirectories matching '", sample_pattern,
         "' in: ", base_dir)

  message("Found ", length(sample_dirs), " sample folders.")

  # ── peek at first CSV to get predicted column names ───────────────────────
  pred_arch_cols <- NULL
  for (d in sample_dirs) {
    wf <- file.path(d, model_folder, "plots",
                    paste0(prefix, basename(d), "_archetype_weights.csv")
    )
    if (file.exists(wf)) {
      pred_arch_cols <- setdiff(
        colnames(read.csv(wf, nrows = 1, check.names = FALSE,
                          stringsAsFactors = FALSE)),
        "Sample"
      )
      break
    }
  }

  if (is.null(pred_arch_cols))
    stop("No _archetype_weights.csv found under model_folder: '",
         model_folder, "'.")

  # ── build archetype_map ───────────────────────────────────────────────────
  if (is.null(archetype_map)) {

    n_arch <- min(length(global_arch_cols), length(pred_arch_cols))

    if (length(global_arch_cols) != length(pred_arch_cols))
      warning("Different number of archetypes in global (",
              length(global_arch_cols), ") vs recomputed (",
              length(pred_arch_cols), ") — using first ", n_arch, ".")

    archetype_map <- setNames(
      pred_arch_cols[seq_len(n_arch)],
      global_arch_cols[seq_len(n_arch)]
    )
    message("Archetype mapping (positional fallback):")

  } else {
    message("Archetype mapping (user-supplied):")
  }

  for (nm in names(archetype_map))
    message("  global: '", nm, "'  -->  predicted: '",
            archetype_map[[nm]], "'")

  # ── build display label map keyed by predicted col name ──────────────────
  pred_cols_ordered <- unname(archetype_map)

  if (is.null(archetype_labels)) {
    display_map <- setNames(
      gsub("[._-]", " ", pred_cols_ordered),
      pred_cols_ordered
    )
  } else if (is.null(names(archetype_labels))) {
    display_map <- setNames(
      as.character(archetype_labels)[seq_along(pred_cols_ordered)],
      pred_cols_ordered
    )
  } else {
    display_map <- archetype_labels
    missing <- setdiff(pred_cols_ordered, names(display_map))
    if (length(missing) > 0)
      display_map[missing] <- gsub("[._-]", " ", missing)
  }

  message("Display labels:")
  for (nm in names(display_map))
    message("  '", nm, "'  -->  '", display_map[[nm]], "'")

  # ════════════════════════════════════════════════════════════════════
  # MAIN LOOP
  # ════════════════════════════════════════════════════════════════════
  results <- list()

  for (sample_dir in sample_dirs) {

    sample_name  <- basename(sample_dir)
    weights_file <- file.path(sample_dir, model_folder, "plots",
                              paste0(prefix, sample_name, "_archetype_weights.csv"))

    if (!file.exists(weights_file)) {
      warning(sample_name, ": weights CSV not found — skipping.")
      next
    }

    pred     <- read.csv(weights_file, check.names = FALSE,
                         stringsAsFactors = FALSE)
    pred_row <- pred %>% filter(Sample == sample_name)

    if (nrow(pred_row) == 0) {
      warning(sample_name, ": test sample row not found — skipping.")
      next
    }

    if (!(sample_name %in% global$Sample)) {
      warning(sample_name, ": not in global_file — skipping.")
      next
    }

    true_row <- global %>% filter(Sample == sample_name)

    for (g_col in names(archetype_map)) {

      p_col <- archetype_map[[g_col]]

      if (!(g_col %in% colnames(true_row))) next
      if (!(p_col %in% colnames(pred_row))) next

      results[[length(results) + 1]] <- data.frame(
        Sample    = sample_name,
        Archetype = p_col,
        True      = as.numeric(true_row[[g_col]]),
        Predicted = as.numeric(pred_row[[p_col]]),
        stringsAsFactors = FALSE
      )
    }

    message("  Processed: ", sample_name)
  }

  if (!length(results))
    stop("No matched predictions found. Check base_dir, model_folder, ",
         "sample_pattern, and archetype_map.")

  # ── apply display labels ──────────────────────────────────────────────────
  df <- bind_rows(results) %>%
    filter(!is.na(True), !is.na(Predicted)) %>%
    mutate(Archetype = dplyr::recode(Archetype, !!!display_map))

  # ── colours ───────────────────────────────────────────────────────────────
  final_labels <- unique(df$Archetype)

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
        setNames(default_cols[length(archetype_cols) +
                                seq_along(missing_labs)],
                 missing_labs)
      )
  }

  # ── metrics ───────────────────────────────────────────────────────────────
  metrics <- df %>%
    group_by(Archetype) %>%
    summarise(
      N           = n(),
      Correlation = cor(True, Predicted),
      RMSE        = sqrt(mean((True - Predicted)^2)),
      .groups     = "drop"
    ) %>%
    mutate(label = paste0(
      "r = ",      formatC(Correlation, digits = 3, format = "f"),
      "\nRMSE = ", formatC(RMSE,        digits = 3, format = "f")
    ))

  # ── save outputs ──────────────────────────────────────────────────────────
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  # plot-ready data (display labels already applied)
  write.csv(
    df,
    file.path(qc_dir, paste0(model_folder,
                             "_test_sample_archetype_prediction_table.csv")),
    row.names = FALSE
  )

  # metrics without N
  write.csv(
    metrics[, c("Archetype", "Correlation", "RMSE")],
    file.path(qc_dir, paste0(model_folder,
                             "_test_sample_archetype_metrics.csv")),
    row.names = FALSE
  )

  # ── scatter plot ──────────────────────────────────────────────────────────
  plot_width  <- max(4 * length(final_labels), 8)
  plot_height <- 4.5

  p <- ggplot(df, aes(x = True, y = Predicted)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60", linewidth = 0.5) +
    geom_point(shape = 21, size = 2.5, alpha = 0.7,
               aes(fill = Archetype), colour = "black", stroke = 0.3) +
    geom_smooth(method = "lm", se = TRUE,
                colour = "grey30", fill = "grey70",
                alpha = 0.15, linewidth = 0.7) +
    geom_label(
      data        = metrics,
      aes(x = -Inf, y = Inf, label = label),
      hjust       = -0.05, vjust = 1.1,
      size        = 3.2, linewidth = 0.3,
      fill        = "white", alpha = 0.85,
      fontface    = "bold", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = archetype_cols) +
    facet_wrap(~Archetype, scales = "free") +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#EEF2F7", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    ) +
    labs(
      title    = paste0("LOO Test Sample Archetype Prediction  [",
                        model_folder, "]"),
      subtitle = "True vs predicted archetype weights for each left-out test sample",
      x        = "True archetype weight",
      y        = "Predicted archetype weight"
    )

  ggsave(
    file.path(qc_dir, paste0(model_folder,
                             "_test_sample_archetype_prediction_summary.pdf")),
    p, width = plot_width, height = plot_height
  )

  ggsave(
    file.path(qc_dir, paste0(model_folder,
                             "_test_sample_archetype_prediction_summary.svg")),
    p, width = plot_width, height = plot_height,
    device = svg
  )

  message("=====================================")
  message("Test sample archetype QC written to: ", qc_dir)

  invisible(list(data = df, metrics = metrics))
}

generate_LOO_test_sample_archetype_QC(
  base_dir         = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/Mesomics_LOO",
  global_file      = "inst/extdata/MESOMICS_references/archetypes.txt",
  model_folder     = "train_test_all_omics",
  sample_pattern   = "^MESO_",
  prefix           = "",
  archetype_map = c(
    "MOFA.MESOMICS.Cell division"            = "Cell division",
    "MOFA.MESOMICS.Tumor-immune-interaction" = "Tumor-immune-interaction",
    "MOFA.MESOMICS.Acinar"                   = "Acinar"
  ),
  archetype_labels = c(
    "Cell division"            = "Cell division",
    "Tumor-immune-interaction" = "Tumor immune interaction",
    "Acinar"                   = "Acinar"
  )
)
