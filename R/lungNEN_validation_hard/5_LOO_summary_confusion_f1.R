summarise_loo_kmeans <- function(
    models_dir,
    matrices_subdir = "train_test_all_omics",
    ref_labels_path,
    ref_cols,
    out_dir = file.path(models_dir, "QC_summary"),
    prefix  = ""
) {

  # ---- checks --------------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (missing(ref_labels_path) || !file.exists(ref_labels_path))
    stop("ref_labels_path must be a valid path to an existing CSV file.")

  if (missing(ref_cols) || !is.character(ref_cols) || length(ref_cols) != 2)
    stop("ref_cols must be a character vector of length 2: c(sample_col, label_col).")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # ---- helper: save plot as PDF + SVG -------------------------------------

  .save_plot <- function(p, base_path, width, height) {
    pdf_path <- paste0(base_path, ".pdf")
    svg_path <- paste0(base_path, ".svg")
    ggplot2::ggsave(pdf_path, p, width = width, height = height)
    ggplot2::ggsave(svg_path, p, width = width, height = height,
                    device = svg)   # <-- base R svg() function object, not string
    message("  Saved PDF : ", basename(pdf_path))
    message("  Saved SVG : ", basename(svg_path))
  }
  ref_raw <- read.csv(ref_labels_path, check.names = FALSE,
                      stringsAsFactors = FALSE)

  missing_ref <- setdiff(ref_cols, colnames(ref_raw))
  if (length(missing_ref))
    stop("ref_labels_path missing columns: ",
         paste(missing_ref, collapse = ", "))

  ref_labels <- ref_raw[, ref_cols, drop = FALSE]
  colnames(ref_labels) <- c("sample", "true_label")

  message("Ground truth loaded: ", nrow(ref_labels), " samples from ",
          basename(ref_labels_path))

  # ---- discover valid sample dirs ------------------------------------------

  sample_dirs <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[
    file.exists(file.path(sample_dirs, matrices_subdir, "plots",
                          paste0(basename(sample_dirs),
                                 "_sample_clusters.csv"))) &
      file.exists(file.path(sample_dirs, matrices_subdir, "plots",
                            paste0(basename(sample_dirs),
                                   "_cluster_label_map.csv")))
  ]

  if (!length(sample_dirs))
    stop("No valid sample dirs found. Run infer_kmeans_labels() and ",
         "plot_kmeans_query_sample() first.")

  message("Found ", length(sample_dirs), " sample dir(s).")

  # ---- aggregate predicted labels for query sample only -------------------

  pred_list <- vector("list", length(sample_dirs))

  for (i in seq_along(sample_dirs)) {

    sdir      <- sample_dirs[i]
    sample_id <- basename(sdir)
    out_plots <- file.path(sdir, matrices_subdir, "plots")

    cl_df <- read.csv(
      file.path(out_plots,
                paste0(prefix, sample_id, "_sample_clusters.csv")),
      stringsAsFactors = FALSE
    )

    map_df <- read.csv(
      file.path(out_plots,
                paste0(prefix, sample_id, "_cluster_label_map.csv")),
      stringsAsFactors = FALSE
    )

    cl_df$kmeans_cluster  <- as.character(cl_df$kmeans_cluster)
    map_df$kmeans_cluster <- as.character(map_df$kmeans_cluster)

    query_row <- cl_df[cl_df$sample == sample_id, , drop = FALSE]

    if (nrow(query_row) == 0) {
      message("  [SKIP] query sample not found in clusters: ", sample_id)
      next
    }

    query_row <- merge(query_row, map_df, by = "kmeans_cluster", all.x = TRUE)

    if (is.na(query_row$mapped_bio_label[1])) {
      message("  [SKIP] No mapped label for query: ", sample_id)
      next
    }

    pred_list[[i]] <- data.frame(
      sample          = sample_id,
      predicted_label = query_row$mapped_bio_label[1],
      stringsAsFactors = FALSE
    )
  }

  pred_df <- do.call(rbind, Filter(Negate(is.null), pred_list))

  if (is.null(pred_df) || nrow(pred_df) == 0)
    stop("No predicted labels could be aggregated.")

  message("Aggregated predictions for ", nrow(pred_df), " query samples.")

  # ---- join with ground truth ----------------------------------------------

  results_df <- merge(pred_df, ref_labels, by = "sample", all.x = TRUE)

  missing_truth <- sum(is.na(results_df$true_label))
  if (missing_truth > 0)
    message("  [WARN] ", missing_truth,
            " sample(s) have no ground truth — excluded from metrics.")

  results_df <- results_df[!is.na(results_df$true_label), ]

  # ---- save aggregated predictions CSV ------------------------------------

  results_csv <- file.path(out_dir,
                           paste0(prefix, matrices_subdir,
                                  "_loo_predictions.csv"))
  write.csv(results_df, results_csv, row.names = FALSE)
  message("Saved predictions      : ", basename(results_csv))

  # ---- confusion matrix counts --------------------------------------------

  all_labels <- sort(unique(c(results_df$true_label,
                              results_df$predicted_label)))

  conf_mat <- as.data.frame(
    table(
      true_label      = factor(results_df$true_label,
                               levels = all_labels),
      predicted_label = factor(results_df$predicted_label,
                               levels = all_labels)
    )
  )
  colnames(conf_mat)[colnames(conf_mat) == "Freq"] <- "n"

  conf_mat$pct <- ave(
    conf_mat$n, conf_mat$true_label,
    FUN = function(x) round(100 * x / sum(x), 1)
  )

  conf_mat$is_diag    <- conf_mat$true_label == conf_mat$predicted_label
  conf_mat$font_color <- ifelse(conf_mat$pct > 50, "white", "grey20")

  conf_csv <- file.path(out_dir,
                        paste0(prefix, matrices_subdir,
                               "_confusion_matrix.csv"))
  write.csv(conf_mat, conf_csv, row.names = FALSE)
  message("Saved confusion CSV    : ", basename(conf_csv))

  # ---- overall accuracy ---------------------------------------------------

  total    <- nrow(results_df)
  correct  <- sum(results_df$true_label == results_df$predicted_label)
  accuracy <- round(100 * correct / total, 1)
  message("Overall accuracy       : ", correct, " / ", total,
          " (", accuracy, "%)")

  # ---- per-class metrics --------------------------------------------------

  metrics_list <- lapply(all_labels, function(lbl) {

    tp <- sum(results_df$true_label == lbl &
                results_df$predicted_label == lbl)
    fp <- sum(results_df$true_label != lbl &
                results_df$predicted_label == lbl)
    fn <- sum(results_df$true_label == lbl &
                results_df$predicted_label != lbl)

    precision <- if ((tp + fp) > 0) round(tp / (tp + fp), 3) else NA
    recall    <- if ((tp + fn) > 0) round(tp / (tp + fn), 3) else NA
    f1        <- if (!is.na(precision) && !is.na(recall) &&
                     (precision + recall) > 0)
      round(2 * precision * recall / (precision + recall), 3)
    else NA

    data.frame(
      label     = lbl,
      precision = precision,
      recall    = recall,
      f1        = f1,
      n_true    = sum(results_df$true_label == lbl),
      stringsAsFactors = FALSE
    )
  })

  metrics_df <- do.call(rbind, metrics_list)

  metrics_csv <- file.path(out_dir,
                           paste0(prefix, matrices_subdir,
                                  "_per_class_metrics.csv"))
  write.csv(metrics_df, metrics_csv, row.names = FALSE)
  message("Saved per-class metrics: ", basename(metrics_csv))

  # ---- confusion matrix heatmap -------------------------------------------

  p_conf <- ggplot2::ggplot(
    conf_mat,
    ggplot2::aes(x = predicted_label, y = true_label)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = pct),
      color = "white", linewidth = 1
    ) +
    ggplot2::geom_tile(
      data  = conf_mat[conf_mat$is_diag & conf_mat$n > 0, ],
      ggplot2::aes(x = predicted_label, y = true_label),
      fill  = NA, color = "#01696f", linewidth = 2
    ) +
    ggplot2::geom_text(
      data = conf_mat[conf_mat$n > 0, ],
      ggplot2::aes(label = n, color = font_color),
      size = 4.5, fontface = "bold", vjust = -0.3
    ) +
    ggplot2::geom_text(
      data = conf_mat[conf_mat$n > 0, ],
      ggplot2::aes(label = paste0("(", pct, "%)"), color = font_color),
      size = 3, vjust = 1.3
    ) +
    ggplot2::scale_fill_gradient(
      low      = "#f7fbff",
      high     = "#2166AC",
      na.value = "grey95",
      name     = "% of true class"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::labs(
      title    = "LOO K-means \u2014 Confusion Matrix",
      subtitle = paste0("Overall accuracy: ", correct, " / ", total,
                        "  (", accuracy, "%)"),
      x        = "Predicted Label",
      y        = "True Label"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face  = "bold", size = 15,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(color = "grey40", size = 11,
                                              margin = ggplot2::margin(b = 12)),
      axis.text.x     = ggplot2::element_text(angle = 35, hjust = 1, size = 11),
      axis.text.y     = ggplot2::element_text(size = 11),
      axis.title      = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 10),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin     = ggplot2::margin(16, 16, 16, 16)
    )

  message("")
  message("--- Confusion matrix plot ---")
  .save_plot(
    p_conf,
    base_path = file.path(out_dir, paste0(prefix, matrices_subdir,
                                          "_confusion_matrix")),
    width     = 3 + length(all_labels) * 1.2,
    height    = 3 + length(all_labels) * 1.0
  )

  # ---- per-class F1 bar plot ----------------------------------------------

  p_f1 <- ggplot2::ggplot(
    metrics_df,
    ggplot2::aes(x = reorder(label, -f1), y = f1, fill = label)
  ) +
    ggplot2::geom_col(width = 0.55, show.legend = FALSE,
                      color = "white", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed",
                        color = "grey50", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(is.na(f1), "NA", sprintf("%.2f", f1)),
        y     = ifelse(is.na(f1) | f1 < 0.12, 0.04, f1 - 0.04),
        color = ifelse(is.na(f1) | f1 < 0.12, "grey30", "white")
      ),
      size = 4, fontface = "bold", vjust = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y     = ifelse(is.na(f1), 0.08, f1) + 0.04,
        label = paste0("n=", n_true)
      ),
      size = 3, color = "grey40", vjust = 0
    ) +
    ggplot2::annotate("text", x = Inf, y = 0.815,
                      label = "F1 = 0.8", hjust = 1.1,
                      size = 3, color = "grey50") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::scale_color_identity() +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.12),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title    = "Per-class F1 Score \u2014 LOO K-means",
      subtitle = paste0("Overall accuracy: ", accuracy, "%  |  ",
                        "n = ", total, " samples"),
      x        = "Biological Label",
      y        = "F1 Score"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face  = "bold", size = 15,
                                                 margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(color = "grey40", size = 11,
                                                 margin = ggplot2::margin(b = 12)),
      axis.text.x        = ggplot2::element_text(angle = 35, hjust = 1,
                                                 size = 11),
      axis.text.y        = ggplot2::element_text(size = 11),
      axis.title         = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      plot.background    = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin        = ggplot2::margin(16, 16, 16, 16)
    )

  message("")
  message("--- F1 score plot ---")
  .save_plot(
    p_f1,
    base_path = file.path(out_dir, paste0(prefix, matrices_subdir,
                                          "_f1_scores")),
    width  = 6,
    height = 4.5
  )

  message("")
  message(strrep("=", 45))
  message("  Samples evaluated  : ", total)
  message("  Overall accuracy   : ", accuracy, "%")
  message("  Outputs saved in   : ", out_dir)
  message(strrep("=", 45))

  invisible(results_df)
}


summarise_loo_kmeans(
  models_dir      = "Analysis_010426/lungNEN_hard_LOO",
  matrices_subdir = "train_test_Meth_only",
  ref_labels_path = "Analysis_010426/lungNEN_hard_ref/sample_cluster_labels.csv",
  ref_cols        = c("sample", "bio_label")
)
