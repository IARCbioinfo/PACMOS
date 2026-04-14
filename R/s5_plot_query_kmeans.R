#' Plot K-means Clusters
#'
#' For each sample folder in \code{models_dir}, loads
#' \code{<sample_id>_stable_input.csv} and
#' \code{<sample_id>_sample_clusters.csv}, generates all pairwise
#' scatter plots, maps clusters to biological labels via majority, and saves
#' a heatmap and cluster-label mapping CSV.
#'
#' @name plot_kmeans_query_sample
#'
#' @param models_dir      Character. Root directory.
#'
#' @param matrices_subdir Character. Folder name where `.RData` files are stored..
#'
#' @param lf_cols         Character vector. LF columns to use as features for k-means.
#'   All pairwise combinations are plotted.
#'
#' @param ref_labels_path Character. Path to CSV containing reference labels.
#'
#' @param ref_cols        Character vector of length 2:
#'   \code{c(sample_col, label_col)} — column names in the ref_labels CSV.
#'   E.g. \code{c("sample", "bio_label")} or
#'   \code{c("sample", "mapped_bio_label")}.
#'
#' @param query_sample    Character. Sample ID to highlight (red diamond with
#'   black border) on scatter plots. Restricts processing to this folder only.
#'
#' @param prefix          Optional character prefix for output files.
#'
#' @return Invisibly returns NULL.
#'
#' @export
plot_kmeans_query_sample <- function(
    models_dir,
    matrices_subdir,
    lf_cols,
    ref_labels_path,
    ref_cols,
    sample_pattern = "",
    prefix          = ""
) {

  # ---- checks --------------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (missing(lf_cols) || !is.character(lf_cols) || length(lf_cols) < 2)
    stop("lf_cols must be a character vector with at least 2 column names.")

  if (missing(ref_labels_path) || !file.exists(ref_labels_path))
    stop("ref_labels_path must be a valid path to an existing CSV file.")

  if (missing(ref_cols) || !is.character(ref_cols) || length(ref_cols) != 2)
    stop("ref_cols must be a character vector of length 2: c(sample_col, label_col).")

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Please install it.")

  # ---- load ref_labels -----------------------------------------------------

  ref_raw <- read.csv(ref_labels_path, check.names = FALSE,
                      stringsAsFactors = FALSE)

  missing_ref_cols <- setdiff(ref_cols, colnames(ref_raw))
  if (length(missing_ref_cols))
    stop("ref_labels_path CSV missing columns: ",
         paste(missing_ref_cols, collapse = ", "))

  ref_labels <- ref_raw[, ref_cols, drop = FALSE]
  colnames(ref_labels) <- c("sample", "bio_label")

  message("Loaded ref_labels: ", nrow(ref_labels), " samples from ",
          basename(ref_labels_path))

  # ---- pairwise LF combinations --------------------------------------------

  lf_pairs <- combn(lf_cols, 2, simplify = FALSE)
  message("Pairwise LF combinations: ", length(lf_pairs))

  # ---- discover valid sample dirs ------------------------------------------

  all_subdirs <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- all_subdirs[grepl(sample_pattern, basename(all_subdirs))]

  if (!length(sample_dirs))
    stop("No subdirectories matching pattern '", sample_pattern,
         "' found in: ", models_dir)

  message("Found ", length(sample_dirs), " sample folder(s).")

  # ---- sample loop --------------------------------------------

  results <- stats::setNames(logical(length(sample_dirs)),
                             basename(sample_dirs))

  for (sdir in sample_dirs) {

    sample_id <- basename(sdir)
    query_sample <- sample_id
    out_dir   <- file.path(sdir, matrices_subdir, "plots")

    message("")
    message(strrep("=", 45))
    message("  Processing: ", sample_id)
    message(strrep("=", 45))

    if (!dir.exists(out_dir)) {
      warning("  No plots/ dir found: ", out_dir, " - skipping.")
      next
    }


    si_path <- file.path(out_dir, paste0(sample_id, "_stable_input.csv"))
    if (!file.exists(si_path)) {
      warning("  stable_input not found, skipping: ", sample_id); next
    }


    stable_input <- read.csv(si_path, check.names = FALSE,
                             stringsAsFactors = FALSE)

    missing_lf <- setdiff(lf_cols, colnames(stable_input))
    if (length(missing_lf))
      stop("  Missing LF columns in stable_input: ",
           paste(missing_lf, collapse = ", "))


    # ---- load sample_clusters ------------------------------------------------

    cl_path <- file.path(out_dir,
                         paste0(prefix, sample_id, "_sample_clusters.csv"))
    if (!file.exists(cl_path)) {
      warning("  sample_clusters not found, skipping: ", sample_id); next
    }

    cluster_df <- read.csv(cl_path, stringsAsFactors = FALSE)
    cluster_df$kmeans_cluster <- as.factor(cluster_df$kmeans_cluster)


    # ---- merge LF coords + cluster assignments -------------------------------

    plot_df              <- stable_input[, c(colnames(stable_input)[1],
                                             lf_cols), drop = FALSE]
    colnames(plot_df)[1] <- "sample"
    plot_df              <- merge(plot_df, cluster_df, by = "sample")
    plot_df$is_query     <- plot_df$sample == sample_id


    # ---- cluster to bio label mapping  -----------------
    # do mapping first so mapped_bio_label is available for scatter colouring

    labeled_df <- merge(cluster_df, ref_labels, by = "sample",
                        all.x = FALSE)

    has_mapping <- nrow(labeled_df) > 0

    if (!has_mapping) {
      message("  [WARN] No overlap between sample_clusters and ref_labels. ",
              "Scatter will colour by kmeans_cluster only.")
    }

    if (has_mapping) {

      labeled_df$kmeans_cluster <- as.character(labeled_df$kmeans_cluster)
      labeled_df$bio_label      <- as.character(labeled_df$bio_label)

      overlap <- as.data.frame(
        table(kmeans_cluster = labeled_df$kmeans_cluster,
              bio_label      = labeled_df$bio_label)
      )
      colnames(overlap)[colnames(overlap) == "Freq"] <- "n"

      overlap$pct <- ave(
        overlap$n, overlap$kmeans_cluster,
        FUN = function(x) round(100 * x / sum(x), 1)
      )

      # majority
      cluster_map <- do.call(rbind, lapply(
        split(overlap, overlap$kmeans_cluster),
        function(sub) {
          data.frame(
            kmeans_cluster   = unique(sub$kmeans_cluster),
            mapped_bio_label = as.character(sub$bio_label[which.max(sub$n)]),
            stringsAsFactors = FALSE
          )
        }
      ))

      # merge mapped_bio_label onto plot_df for scatter colouring
      plot_df <- merge(
        plot_df,
        cluster_map[, c("kmeans_cluster", "mapped_bio_label")],
        by.x = "kmeans_cluster", by.y = "kmeans_cluster",
        all.x = TRUE
      )

    } else {
      # fallback: use cluster number as label
      plot_df$mapped_bio_label <- as.character(plot_df$kmeans_cluster)
      cluster_map              <- NULL
      overlap                  <- NULL
    }

    # ensure is_query is still correct after merge reorder
    plot_df$is_query <- plot_df$sample == sample_id

    # ---- multi-page scatter PDF ----------------------------------------------

    scatter_pdf <- file.path(out_dir,
                             paste0(prefix, sample_id, "_kmeans_scatter.pdf"))
    pdf(scatter_pdf, width = 8, height = 6)

    for (pair in lf_pairs) {

      x_col <- pair[1]
      y_col <- pair[2]

      bg_df    <- plot_df[!plot_df$is_query, ]
      query_df <- plot_df[ plot_df$is_query, ]

      p <- ggplot2::ggplot(
        bg_df,
        ggplot2::aes_string(
          x     = x_col,
          y     = y_col,
          fill  = "mapped_bio_label",
          shape = "kmeans_cluster"
        )
      ) +
        # background samples: filled shapes with grey border
        ggplot2::geom_point(size = 2.5, alpha = 0.8, color = "grey30") +

        # query sample: red fill + black border, no label
        ggplot2::geom_point(
          data        = query_df,
          ggplot2::aes_string(
            x     = x_col,
            y     = y_col,
            fill  = "mapped_bio_label",
            shape = "kmeans_cluster"
          ),
          color       = "black",
          size        = 6,
          stroke      = 2,
          inherit.aes = FALSE
        ) +

        # filled shapes (21-25) support both fill and border color
        ggplot2::scale_shape_manual(
          values = setNames(
            21:(21 + nlevels(plot_df$kmeans_cluster) - 1),
            levels(plot_df$kmeans_cluster)
          )
        ) +
        ggplot2::scale_fill_brewer(palette = "Set2", na.value = "grey80") +

        # two separate legends: fill = bio label, shape = cluster number
        ggplot2::guides(
          fill  = ggplot2::guide_legend(
            title        = "Reference Label",
            override.aes = list(
              shape  = 21,
              size   = 3,
              color  = "grey30",
              stroke = 0.4
            )
          ),
          shape = ggplot2::guide_legend(
            title        = "K-means Cluster",
            override.aes = list(
              fill   = "grey60",
              size   = 3,
              color  = "grey30",
              stroke = 0.4
            )
          )
        )  +
        ggplot2::labs(
          title    = paste0(sample_id, " - K-means clusters"),
          subtitle = paste0(x_col, " vs ", y_col
          ),
          x        = x_col,
          y        = y_col
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(legend.position = "right")

      print(p)
    }

    dev.off()
    message("  Saved scatter (", length(lf_pairs), " pages): ",
            basename(scatter_pdf))

    # ---- save cluster map CSV ------------------------------------------------

    if (!is.null(cluster_map)) {
      csv_out <- file.path(out_dir,
                           paste0(prefix, sample_id, "_cluster_label_map.csv"))
      write.csv(cluster_map, csv_out, row.names = FALSE)
      message("  Saved cluster map        : ", basename(csv_out))
    }

    # ---- heatmap -------------------------------------------------------------

    if (!is.null(overlap)) {

      cluster_order <- cluster_map$kmeans_cluster[
        order(cluster_map$mapped_bio_label)
      ]
      label_order <- sort(unique(ref_labels$bio_label))

      overlap$kmeans_cluster <- factor(overlap$kmeans_cluster,
                                       levels = cluster_order)
      overlap$bio_label      <- factor(overlap$bio_label,
                                       levels = label_order)

      p_heatmap <- ggplot2::ggplot(
        overlap,
        ggplot2::aes(x = bio_label, y = kmeans_cluster, fill = pct)
      ) +
        ggplot2::geom_tile(color = "white", linewidth = 0.8) +
        ggplot2::geom_text(
          data = overlap[overlap$n > 0, ],
          ggplot2::aes(label = paste0(n, "\n(", pct, "%)")),
          size = 3, fontface = "bold"
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = "#2166AC") +
        ggplot2::labs(
          title = paste0(sample_id, " - K-means (clusters vs. Reference labels)"),
          x     = "Biological Label",
          y     = "K-means Cluster",
          fill  = "% of cluster"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
        )

      heatmap_pdf <- file.path(out_dir,
                               paste0(prefix, sample_id, "_kmeans_heatmap.pdf"))
      ggplot2::ggsave(heatmap_pdf, p_heatmap, width = 6, height = 4.5)
      message("  Saved heatmap            : ", basename(heatmap_pdf))
    }

    message("")
    message(strrep("=", 45))
    message("  Done: ", sample_id)
    message(strrep("=", 45))

    results[sample_id] <- TRUE


  }

  n_ok   <- sum(results,  na.rm = TRUE)
  n_fail <- sum(!results, na.rm = TRUE)

  message("")
  message(strrep("=", 45))
  message("  Done. ", n_ok, " sample(s) processed",
          if (n_fail > 0) paste0(", ", n_fail, " skipped.") else ".")
  message(strrep("=", 45))

  invisible(NULL)

}

