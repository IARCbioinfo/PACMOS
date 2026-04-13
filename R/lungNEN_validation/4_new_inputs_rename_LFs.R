s3_plot_all_samples_mofa <- function(
    models_dir,
    matrices_subdir  = "train_test_all_omics",
    reference_LFs,
    reference_axes,
    id_col       = NULL,
    group        = "group1",
    python_bin,
    out_plot_dir = NULL,
    prefix       = ""
) {

  # ---- input checks --------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  reticulate::use_python(python_bin, required = TRUE)

  if (length(reference_axes) < 2)
    stop("At least two reference_axes are required.")

  sample_dirs <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)
  if (!length(sample_dirs))
    stop("No sample subdirectories found in: ", models_dir)

  message("Found ", length(sample_dirs), " sample dir(s) in: ", models_dir)

  # ---- load reference LFs --------------------------------------------------

  if (is.character(reference_LFs)) {
    ref_df <- read.csv(reference_LFs, check.names = FALSE)
  } else {
    ref_df <- as.data.frame(reference_LFs)
  }

  if (!all(reference_axes %in% colnames(ref_df)))
    stop("Some reference_axes not found in reference_LFs: ",
         paste(setdiff(reference_axes, colnames(ref_df)), collapse = ", "))

  if (is.null(id_col))
    id_col <- if ("Sample" %in% colnames(ref_df)) "Sample" else colnames(ref_df)[1]

  if (!(id_col %in% colnames(ref_df)))
    stop("id_col '", id_col, "' not found in reference_LFs.")

  ref_ids           <- as.character(ref_df[[id_col]])
  ref_mat           <- as.matrix(ref_df[, reference_axes, drop = FALSE])
  rownames(ref_mat) <- ref_ids

  # ---- greedy factor picker ------------------------------------------------

  .greedy_pick <- function(cor_mat, k) {
    pick         <- rep(NA_integer_, k)
    used_factors <- integer(0)
    for (j in seq_len(k)) {
      ord <- order(abs(cor_mat[, j]), decreasing = TRUE)
      for (idx in ord) {
        if (!(idx %in% used_factors)) {
          pick[j]      <- idx
          used_factors <- c(used_factors, idx)
          break
        }
      }
    }
    pick
  }

  # ---- alignment QC PDF helper ---------------------------------------------

  .write_alignment_QC_pdf <- function(
    out_pdf, sample_id, model_id,
    MOFA.LFs, cor_mat, pick,
    ref_sub, Zk, common_samples,
    reference_axes
  ) {
    pdf(out_pdf, width = 6, height = 6)
    on.exit(dev.off())

    # -- page 1: text summary ------------------------------------------------
    mofa_factor_names <- colnames(MOFA.LFs)[pick]
    r_chosen <- sapply(seq_along(pick), function(j) cor_mat[pick[j], j])

    summary_lines <- c(
      paste0("Sample : ", sample_id),
      paste0("Model  : ", model_id),
      paste0("Mode   : all-samples (n = ", length(common_samples), ")"),
      "",
      "Factor assignment (reference axis <- new MOFA factor):",
      strrep("\u2500", 50),
      mapply(function(ref_ax, mofa_ax, r)
        sprintf("  %-12s  <-  %-12s  r = %+.3f", ref_ax, mofa_ax, r),
        reference_axes, mofa_factor_names, r_chosen)
    )

    p_text <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 1,
                        label = paste(summary_lines, collapse = "\n"),
                        hjust = 0, vjust = 1, size = 3.6, family = "mono") +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))
    print(p_text)

    # -- page 2: correlation heatmap -----------------------------------------
    cor_df           <- as.data.frame(as.table(cor_mat))
    colnames(cor_df) <- c("MOFA_Factor", "Ref_Axis", "Correlation")
    cor_df$MOFA_Factor <- factor(cor_df$MOFA_Factor, levels = rev(rownames(cor_mat)))
    cor_df$Ref_Axis    <- factor(cor_df$Ref_Axis, levels = reference_axes)

    chosen_df <- data.frame(
      MOFA_Factor = factor(rownames(cor_mat)[pick], levels = levels(cor_df$MOFA_Factor)),
      Ref_Axis    = factor(reference_axes,          levels = levels(cor_df$Ref_Axis)),
      stringsAsFactors = FALSE
    )

    p_heat <- ggplot2::ggplot(cor_df,
                              ggplot2::aes(x = Ref_Axis, y = MOFA_Factor, fill = Correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::geom_tile(data = chosen_df,
                         ggplot2::aes(x = Ref_Axis, y = MOFA_Factor),
                         inherit.aes = FALSE, fill = NA, color = "black", linewidth = 1.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(Correlation, digits = 2, format = "f")),
        size = 2.8, color = "grey10") +
      ggplot2::geom_text(data = chosen_df,
                         ggplot2::aes(x = Ref_Axis, y = MOFA_Factor, label = "\u2605"),
                         inherit.aes = FALSE, size = 4, vjust = -0.5, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = "#4575B4", mid = "white", high = "#D73027",
        midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
      ggplot2::labs(
        title    = paste0(sample_id, "  \u2014  Factor Matching correlation heatmap"),
        subtitle = "bold border = Assigned factor for that reference axis",
        x = "Reference MOFA axes", y = "Retrained MOFA factors") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(angle = 35, hjust = 1),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40"))
    print(p_heat)

    # -- pages 3+: per-axis scatter plots ------------------------------------
    alignment_metrics <- vector("list", length(reference_axes))

    for (j in seq_along(reference_axes)) {
      axis_name <- reference_axes[j]
      ref_vals  <- ref_sub[, j]
      new_vals  <- Zk[common_samples, j]

      r_val    <- cor(ref_vals, new_vals, use = "pairwise.complete.obs")
      rmse_val <- sqrt(mean((ref_vals - new_vals)^2, na.rm = TRUE))

      alignment_metrics[[j]] <- data.frame(
        Model = model_id, Sample = sample_id, Axis = axis_name,
        Correlation = r_val, RMSE = rmse_val, stringsAsFactors = FALSE)

      df_sc   <- data.frame(Reference = ref_vals, Aligned = new_vals)
      lim_min <- min(c(df_sc$Reference, df_sc$Aligned), na.rm = TRUE)
      lim_max <- max(c(df_sc$Reference, df_sc$Aligned), na.rm = TRUE)

      p_sc <- ggplot2::ggplot(df_sc, ggplot2::aes(x = Reference, y = Aligned)) +
        ggplot2::geom_point(size = 3, shape = 21,
                            fill = "grey40", colour = "black", alpha = 0.8) +
        ggplot2::geom_smooth(method = "lm", se = FALSE,
                             colour = "#B81330", linewidth = 0.9) +
        ggplot2::coord_equal(
          xlim = c(lim_min, lim_max), ylim = c(lim_min, lim_max), expand = TRUE) +
        ggplot2::annotate("label",
                          x = lim_min, y = lim_max, hjust = 0, vjust = 1,
                          label = paste0("r = ", sprintf("%.3f", r_val),
                                         "\nRMSE = ", sprintf("%.3f", rmse_val)),
                          size = 3.3, label.size = 0.4, fill = "white") +
        ggplot2::labs(
          title = paste0(sample_id, " \u2013 ", axis_name),
          x = "Reference LF", y = "Matched retrained LF") +
        ggplot2::theme_classic(base_size = 12)
      print(p_sc)
    }

    do.call(rbind, alignment_metrics)
  }

  # ---- main loop -----------------------------------------------------------

  n_models    <- length(sample_dirs)
  all_metrics <- vector("list", n_models)

  for (i in seq_along(sample_dirs)) {

    sdir       <- sample_dirs[i]
    sample_id  <- basename(sdir)
    inputs_dir <- file.path(sdir, matrices_subdir)

    plot_dir <- if (!is.null(out_plot_dir)) out_plot_dir else
      file.path(inputs_dir, "plots")

    if (!dir.exists(plot_dir))
      dir.create(plot_dir, recursive = TRUE)

    model_path <- file.path(inputs_dir, paste0("MOFA-", sample_id, ".hdf5"))

    message("")
    message(strrep("=", 45))
    message(sprintf("  [%d / %d]  %s", i, n_models, sample_id))
    message(strrep("=", 45))

    if (!dir.exists(inputs_dir)) {
      message("  [SKIP] No '", matrices_subdir, "/' subfolder found.")
      next
    }

    if (!file.exists(model_path)) {
      message("  [SKIP] HDF5 not found: ", model_path)
      next
    }

    model_id <- sub("\\.hdf5$", "", basename(model_path))

    # -- load model ----------------------------------------------------------
    MOFAmodel     <- MOFA2::load_model(model_path)
    MOFA.LFs      <- as.data.frame(MOFAmodel@expectations$Z[[group]])
    model_samples <- rownames(MOFA.LFs)

    # -- common samples ------------------------------------------------------
    common_samples <- intersect(model_samples, rownames(ref_mat))
    message("  Common samples   : ", length(common_samples))

    if (length(common_samples) < 3) {
      message("  [SKIP] Too few common samples (< 3). Skipping.")
      next
    }

    # -- correlation matrix (all common samples) -----------------------------
    ref_sub <- ref_mat[common_samples, , drop = FALSE]

    cor_mat <- cor(
      MOFA.LFs[common_samples, , drop = FALSE],
      ref_sub,
      use = "pairwise.complete.obs"
    )

    # -- greedy factor assignment --------------------------------------------
    k    <- length(reference_axes)
    pick <- .greedy_pick(cor_mat, k)

    if (any(is.na(pick))) {
      message("  [SKIP] Could not uniquely assign all reference axes. Skipping.")
      next
    }

    names(pick) <- reference_axes

    # -- build Zk ------------------------------------------------------------
    Zk           <- as.matrix(MOFA.LFs[, pick, drop = FALSE])
    colnames(Zk) <- reference_axes

    for (j in seq_len(k)) {
      if (!is.na(cor_mat[pick[j], j]) && cor_mat[pick[j], j] < 0)
        Zk[, j] <- -Zk[, j]
    }

    message("  Zk               : ", nrow(Zk), " samples x ", ncol(Zk), " axes")
    message("  Writing outputs to: ", plot_dir)

    # -- alignment QC PDF + metrics ------------------------------------------
    align_pdf    <- file.path(plot_dir, paste0(prefix, sample_id, "_alignment_QC.pdf"))
    alignment_df <- .write_alignment_QC_pdf(
      out_pdf        = align_pdf,
      sample_id      = sample_id,
      model_id       = model_id,
      MOFA.LFs       = MOFA.LFs,
      cor_mat        = cor_mat,
      pick           = pick,
      ref_sub        = ref_sub,
      Zk             = Zk,
      common_samples = common_samples,
      reference_axes = reference_axes
    )

    all_metrics[[i]] <- alignment_df

    # -- CSVs ----------------------------------------------------------------
    write.csv(
      cbind(Sample = rownames(Zk), as.data.frame(Zk)),
      file.path(plot_dir, paste0(prefix, sample_id, "_aligned_LFs_all_samples.csv")),
      row.names = FALSE
    )

    write.csv(
      alignment_df,
      file.path(plot_dir, paste0(prefix, sample_id, "_quality_check_metrics.csv")),
      row.names = FALSE
    )

    # -- 2D projection PDF ---------------------------------------------------
    df_plot        <- as.data.frame(Zk[common_samples, , drop = FALSE])
    df_plot$Sample <- rownames(df_plot)
    combos         <- combn(reference_axes, 2, simplify = FALSE)
    proj_pdf       <- file.path(plot_dir, paste0(prefix, sample_id, "_projection.pdf"))

    local({
      pdf(proj_pdf, width = 10, height = 8)
      on.exit(dev.off())
      for (pair in combos) {
        ax1 <- pair[1]; ax2 <- pair[2]
        p <- ggplot2::ggplot(df_plot,
                             ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]])) +
          ggplot2::geom_point(size = 3, shape = 21,
                              fill = "grey70", colour = "black", alpha = 0.8) +
          ggplot2::coord_equal() +
          ggplot2::labs(x = ax1, y = ax2,
                        title = paste0(sample_id, "  |  ", ax1, " vs ", ax2)) +
          ggplot2::theme_bw(base_size = 12)
        print(p)
      }
    })

    message("  [DONE] ", sample_id)
  }

  # ---- aggregated metrics --------------------------------------------------

  all_metrics_df <- do.call(rbind, Filter(Negate(is.null), all_metrics))
  agg_path <- file.path(models_dir,
                        paste0(prefix, "all_samples_quality_check_metrics.csv"))
  write.csv(all_metrics_df, agg_path, row.names = FALSE)

  message("")
  message(strrep("=", 45))
  message("  All ", n_models, " model(s) processed.")
  message("  Aggregated metrics: ", agg_path)
  message("  Outputs in        : ", models_dir)
  message(strrep("=", 45))

  invisible(all_metrics_df)
}

s3_plot_all_samples_mofa(
  models_dir      = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO",
  matrices_subdir = "inputs",
  reference_LFs   = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/lungNEN_references/lungNEN_LFs.csv",
  reference_axes  = c("Factor1", "Factor2",
                      "Factor5"),
  python_bin      = "/home/lipikal/miniconda3/envs/mofa_env/bin/python"
)
