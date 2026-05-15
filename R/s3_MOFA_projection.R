#' Step 3:  Align retrained MOFA latent factors to reference space and project query sample
#'
#' @description
#' Matches and aligns latent factors from retrained MOFA models to a reference
#' MOFA factor space and projects the query sample into that reference space.
#'
#' @details
#' For each MOFA2 HDF5 model in `models_dir`, this function:
#' \itemize{
#'   \item Loads the retrained MOFA model
#'   \item Computes correlations between retrained latent factors and reference latent factors
#'   \item Assigns each reference axis to the best-matching retrained factor
#'   \item Aligns factor signs to ensure consistency with the reference space
#'   \item Projects the query sample into the aligned reference latent factor space
#' }
#'
#' The function generates multiple outputs including:
#' \itemize{
#'   \item Quality control PDF with factor assignment and correlation heatmaps
#'   \item 2D projection plots of all samples with the query sample highlighted
#'   \item CSV files containing aligned latent factors and projection results
#' }
#'
#' @name s3_plot_query_samples_mofa
#'
#' @param models_dir Root directory folder. Same as `s1_add_sample_to_mofa() outdir`.
#'
#' @param matrices_subdir Folder name where `.hdf5` files are stored (`inputs` by default).
#'
#' @param query_sample Character. Sample ID of the query sample.
#'
#' @param reference_LFs data.frame or path to CSV containing reference latent
#'   factors. Must contain \code{id_col} and all \code{reference_axes} columns.
#'
#' @param reference_axes Character vector of reference axis column names (latent factors) we need to match and align.
#'
#' @param id_col Character. Column name in \code{reference_LFs} that contains
#'   sample identifiers. These IDs must match the sample names used in the
#'   MOFA model. If \code{NULL}, the function uses the \code{"Sample"} column
#'   if present, otherwise the first column.
#'
#' @param group MOFA group name (default \code{"group1"}).
#'
#' @param python_bin Path to the Python binary used by the MOFA
#' environment via the `reticulate` package.
#'
#' @param output_dir Output directory. Defaults to each sample folder's \code{outputs/} directory.
#'
#' @param prefix Optional character prefix for all written files.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{metrics}: data frame of alignment metrics.
#'   \item \code{samples}: nested list of per-sample metrics, ggplot objects,
#'   and output file paths.
#' }
#'
#' @examples
#' mofa_dir <- system.file("extdata/MESOMICS_references", package = "PACMOS")
#' reference_LFs <- system.file("extdata/MESOMICS_references",
#'                             "MESOMICS_latent_factors.csv", package = "PACMOS")
#' query_csv <- system.file("extdata/test_data", "MESOMICS_test_expr.csv", package = "PACMOS")
#' out_dir <- file.path(tempdir(), "pacmos_s3")
#' dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#' s1_add_sample_to_mofa(
#'   query_matrix_path = query_csv,
#'   mofa_dir = mofa_dir,
#'   value_data_types = "D_exprB_MOFA",
#'   outdir = out_dir,
#'   python_bin = Sys.getenv("PACMOS_PYTHON", unset = "")
#' )
#'
#' s2_run_mofa(
#'   models_dir = out_dir,
#'   matrices_subdir = "inputs",
#'   num_factors = 2,
#'   maxiter = 5,
#'   python_bin = Sys.getenv("PACMOS_PYTHON", unset = ""),
#'   views_map = c(RNA = "D_exprB_MOFA")
#' )
#'
#' sample_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
#' sample_dirs <- sample_dirs[dir.exists(file.path(sample_dirs, "inputs"))]
#' sample_id <- basename(sample_dirs[1])
#'
#' s3_plot_query_samples_mofa(
#'   models_dir = out_dir,
#'   matrices_subdir = "inputs",
#'   query_sample = sample_id,
#'   reference_LFs = reference_LFs,
#'   reference_axes = c("Morphology_LF", "Adaptive-response_LF"),
#'   python_bin = Sys.getenv("PACMOS_PYTHON", unset = "")
#' )
#'
#' @export
s3_plot_query_samples_mofa <- function(
    models_dir,
    matrices_subdir = "inputs",
    query_sample,
    reference_LFs,
    reference_axes,
    id_col       = NULL,
    group        = "group1",
    python_bin   = NULL,
    output_dir= NULL,
    prefix       = ""
) {

  # ---- dependency check ------------------------------------------------------
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required.", call. = FALSE)
  }

  if (!is.null(python_bin) &&
      length(python_bin) == 1L &&
      nzchar(python_bin)) {

    if (!file.exists(python_bin)) {
      stop("python_bin does not exist: ", python_bin, call. = FALSE)
    }

    reticulate::use_python(python_bin, required = TRUE)
  }

  # ---- input checks --------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (missing(query_sample) || is.null(query_sample) ||
      !nzchar(query_sample))
    stop("query_sample must be provided (LOO query sample ID).")

  if (length(reference_axes) < 2)
    stop("At least two reference_axes are required.")

  all_dirs    <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- all_dirs[
    dir.exists(file.path(all_dirs, matrices_subdir))
  ]
  if (!length(sample_dirs))
    stop("No sample subdirectories found in: ", models_dir)

  #message("Found ", length(sample_dirs), " sample dir(s) in: ", models_dir)

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

  ref_ids <- as.character(ref_df[[id_col]])
  ref_mat <- as.matrix(ref_df[, reference_axes, drop = FALSE])
  rownames(ref_mat) <- ref_ids

  # ---- alignment QC PDF ---------------------------------------------

  .write_alignment_QC_pdf <- function(
    out_pdf, sample_id, model_id,
    MOFA.LFs, cor_mat, pick,
    ref_sub, Zk, train_samples,
    reference_axes
  ) {

    pdf(out_pdf, width = 6, height = 6)
    on.exit(dev.off())

    # -- page 1: text summary ------------------------------------------------
    mofa_factor_names <- colnames(MOFA.LFs)[pick]
    r_chosen <- sapply(seq_along(pick),
                       function(j) cor_mat[pick[j], j])

    summary_lines <- c(
      paste0("Sample : ", sample_id),
      paste0("Model  : ", model_id),
      "",
      "Factor assignment (reference axis <- new MOFA factor):",
      strrep("-", 50),
      mapply(function(ref_ax, mofa_ax, r)
        sprintf("  %-12s  <-  %-12s  r = %+.3f",
                ref_ax, mofa_ax, r),
        reference_axes, mofa_factor_names, r_chosen)
    )

    p_text <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 1,
        label    = paste(summary_lines, collapse = "\n"),
        hjust = 0, vjust = 1, size = 3.6, family = "mono"
      ) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))

    print(p_text)

    # -- page 2: correlation heatmap -----------------------------------------
    cor_df           <- as.data.frame(as.table(cor_mat))
    colnames(cor_df) <- c("MOFA_Factor", "Ref_Axis", "Correlation")
    cor_df$MOFA_Factor <- factor(cor_df$MOFA_Factor,
                                 levels = rev(rownames(cor_mat)))
    cor_df$Ref_Axis    <- factor(cor_df$Ref_Axis, levels = reference_axes)

    chosen_df <- data.frame(
      MOFA_Factor = factor(rownames(cor_mat)[pick],
                           levels = levels(cor_df$MOFA_Factor)),
      Ref_Axis    = factor(reference_axes,
                           levels = levels(cor_df$Ref_Axis)),
      stringsAsFactors = FALSE
    )

    p_heat <- ggplot2::ggplot(
      cor_df,
      ggplot2::aes(x = Ref_Axis, y = MOFA_Factor, fill = Correlation)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::geom_tile(
        data        = chosen_df,
        ggplot2::aes(x = Ref_Axis, y = MOFA_Factor),
        inherit.aes = FALSE,
        fill = NA, color = "black", linewidth = 1.5
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(Correlation, digits = 2, format = "f")),
        size = 2.8, color = "grey10"
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#4575B4", mid = "white", high = "#D73027",
        midpoint = 0, limits = c(-1, 1), name = "Pearson r"
      ) +
      ggplot2::labs(
        title    = paste0(sample_id),
        subtitle = "bold border = Assigned factor for that reference axis",
        x = "Reference MOFA axes", y = "Retrained MOFA factors"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x   = ggplot2::element_text(angle = 35, hjust = 1),
        plot.subtitle = ggplot2::element_text(size = 8, color = "grey40")
      )

    print(p_heat)

    # -- pages 3+: per-axis scatter plots ------------------------------------
    alignment_metrics <- vector("list", length(reference_axes))
    scatter_plots <- list()

    for (j in seq_along(reference_axes)) {

      axis_name <- reference_axes[j]
      ref_vals <- ref_sub[, j]
      new_vals <- Zk[train_samples, j]

      r_val    <- cor(ref_vals, new_vals, use = "pairwise.complete.obs")
      rmse_val <- sqrt(mean((ref_vals - new_vals)^2, na.rm = TRUE))

      alignment_metrics[[j]] <- data.frame(
        Model       = model_id,
        Sample      = sample_id,
        Axis        = axis_name,
        Correlation = r_val,
        RMSE        = rmse_val,
        stringsAsFactors = FALSE
      )

      df_sc   <- data.frame(Reference = ref_vals, Aligned = new_vals)
      lim_min <- min(c(df_sc$Reference, df_sc$Aligned), na.rm = TRUE)
      lim_max <- max(c(df_sc$Reference, df_sc$Aligned), na.rm = TRUE)

      p_sc <- ggplot2::ggplot(
        df_sc, ggplot2::aes(x = Reference, y = Aligned)
      ) +
        ggplot2::geom_point(size = 3, shape = 21,
                            fill = "grey40", colour = "black", alpha = 0.8) +
        ggplot2::geom_smooth(method = "lm", se = FALSE,
                             colour = "#B81330", linewidth = 0.9) +
        ggplot2::coord_equal(
          xlim = c(lim_min, lim_max),
          ylim = c(lim_min, lim_max),
          expand = TRUE
        ) +
        ggplot2::annotate(
          "label",
          x = lim_min, y = lim_max, hjust = 0, vjust = 1,
          label = paste0("r = ",    sprintf("%.3f", r_val),
                         "\nRMSE = ", sprintf("%.3f", rmse_val)),
          size = 3.3, fill = "white"
        ) +
        ggplot2::labs(
          title = paste0(sample_id, " - ", axis_name),
          x = "Reference LF", y = "Matched retrained LF"
        ) +
        ggplot2::theme_classic(base_size = 12)

      print(p_sc)
      scatter_plots[[axis_name]] <- p_sc

    }

    #do.call(rbind, alignment_metrics)
    list(
      metrics = do.call(rbind, alignment_metrics),
      plots = list(
        summary = p_text,
        heatmap = p_heat,
        scatter_by_axis = scatter_plots
      )
    )
  }

  # ---- main loop -----------------------------------------------------------

  n_models        <- length(sample_dirs)
  all_metrics     <- vector("list", n_models)
  all_sample_results <- list()

  for (i in seq_along(sample_dirs)) {

    sdir        <- sample_dirs[i]
    sample_id   <- basename(sdir)
    if (sample_id != query_sample) next

    inputs_dir  <- file.path(sdir, matrices_subdir)

    output_dir<- if (!is.null(output_dir)) output_dir else
      file.path(sdir, "outputs")

    if (!dir.exists(output_dir))
      dir.create(output_dir, recursive = TRUE)

    # locate HDF5
    model_path <- file.path(inputs_dir,
                            paste0("MOFA-", sample_id, ".hdf5"))

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

    model_id <- sub("\\.hdf5$", "", basename(model_path))   # "MOFA-MESO_001_T"

    # -- load model ----------------------------------------------------------
    MOFAmodel     <- MOFA2::load_model(model_path)
    MOFA.LFs      <- as.data.frame(MOFAmodel@expectations$Z[[group]])
    model_samples <- rownames(MOFA.LFs)

    # -- validate query sample -------------------------------------------
    if (!(query_sample %in% model_samples)) {
      message("  [SKIP] '", query_sample, "' not in model. Skipping.")
      next
    }

    train_samples <- setdiff(model_samples, query_sample)
    message("  Training samples : ", length(train_samples))
    message("  Query sample      : ", query_sample)

    missing_ref <- setdiff(train_samples, rownames(ref_mat))
    if (length(missing_ref)) {
      message("  [SKIP] ", length(missing_ref),
              " training sample(s) missing from reference. Skipping.")
      next
    }

    # -- correlation matrix (training samples only) --------------------------
    ref_sub <- ref_mat[train_samples, , drop = FALSE]

    cor_mat <- cor(
      MOFA.LFs[train_samples, , drop = FALSE],
      ref_sub,
      use = "pairwise.complete.obs"
    )

    # -- factor assignment --------------------------------------------

    k      <- length(reference_axes)
    cost_t <- t(1 - abs(cor_mat)^2)       # [k × n_factors]

    assignment <- clue::solve_LSAP(cost_t)
    pick       <- as.integer(assignment)
    names(pick) <- reference_axes

    r_chosen <- cor_mat[cbind(pick, seq_len(k))] # rows= retraiend, col = ref

    if (any(is.na(pick)) || any(pick < 1) || any(pick > nrow(cor_mat))) {
      message("  [SKIP] Invalid LSAP assignment for: ", sample_id)
      next
    }
    if (any(abs(r_chosen) < 0.3)) {
      warning("Weak factor match (|r| < 0.3): ",
              paste(reference_axes[abs(r_chosen) < 0.3], collapse = ", "))
    }

    # -- build Zk: aligned factor matrix ---------------
    Zk           <- as.matrix(MOFA.LFs[, pick, drop = FALSE])
    colnames(Zk) <- reference_axes

    for (j in seq_len(k)) {
      if (!is.na(cor_mat[pick[j], j]) && cor_mat[pick[j], j] < 0)
        Zk[, j] <- -Zk[, j]
    }

    # -- build stable_input --------------------------------------------------
    # Reference coordinates are fixed for all reference samples except
    # query_sample. The query sample row comes from Zk (freshly projected).
    stable_input <- rbind(
      ref_mat[setdiff(rownames(ref_mat), query_sample), , drop = FALSE],
      Zk[query_sample, , drop = FALSE]
    )

    #message("  Stable input     : ", nrow(stable_input), " samples x ",
    #        ncol(stable_input), " axes")
    message("  Writing outputs to: ", output_dir)

    # -- alignment QC PDF + metrics ------------------------------------------
    align_pdf <- file.path(
      output_dir,
      paste0(prefix, sample_id, "_quality_check_metrics.pdf")
    )

    #alignment_df <- .write_alignment_QC_pdf(
    #  out_pdf       = align_pdf,
    #  sample_id     = sample_id,
    #  model_id      = model_id,
    #  MOFA.LFs      = MOFA.LFs,
    #  cor_mat       = cor_mat,
    #  pick          = pick,
    #  ref_sub       = ref_sub,
    #  Zk            = Zk,
    #  train_samples = train_samples,
    #  reference_axes = reference_axes
    #)

    #all_metrics[[i]] <- alignment_df

    qc_result <- .write_alignment_QC_pdf(
      out_pdf       = align_pdf,
      sample_id     = sample_id,
      model_id      = model_id,
      MOFA.LFs      = MOFA.LFs,
      cor_mat       = cor_mat,
      pick          = pick,
      ref_sub       = ref_sub,
      Zk            = Zk,
      train_samples = train_samples,
      reference_axes = reference_axes
    )

    alignment_df <- qc_result$metrics
    all_metrics[[i]] <- alignment_df

    # -- CSVs ----------------------------------------------------------------
    write.csv(
      cbind(Sample = query_sample,
            as.data.frame(Zk[query_sample, , drop = FALSE])),
      file.path(output_dir,
                paste0(prefix, sample_id, "_query_sample_LFs.csv")),
      row.names = FALSE
    )

    write.csv(
      cbind(Sample = rownames(Zk), as.data.frame(Zk)),
      file.path(output_dir,
                paste0(prefix, sample_id, "_retrained_LFs_all_samples.csv")),
      row.names = FALSE
    )

    write.csv(
      cbind(Sample = rownames(stable_input), as.data.frame(stable_input)),
      file.path(output_dir,
                paste0(prefix, sample_id, "_stable_input.csv")),
      row.names = FALSE
    )

    write.csv(
      alignment_df,
      file.path(output_dir,
                paste0(prefix, sample_id, "_quality_check_metrics.csv")),
      row.names = FALSE
    )

    # -- 2D projection PDF ---------------------------------------------------
    df_plot         <- as.data.frame(stable_input)
    df_plot$Sample  <- rownames(stable_input)
    combos          <- combn(reference_axes, 2, simplify = FALSE)

    proj_pdf <- file.path(
      output_dir,
      paste0(prefix, sample_id, "_projection.pdf")
    )

    pdf(proj_pdf, width = 10, height = 8)
    projection_plots <- list()

    for (pair in combos) {
      ax1 <- pair[1]; ax2 <- pair[2]

      p <- ggplot2::ggplot(
        df_plot,
        ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]])
      ) +
        ggplot2::geom_point(
          size = 3, shape = 21,
          fill = "grey70", colour = "black", alpha = 0.8
        ) +
        ggplot2::geom_point(
          data   = df_plot[df_plot$Sample == query_sample, ],
          fill   = "red", colour = "black",
          shape  = 21, size = 5
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(
          x     = ax1, y = ax2,
          title = paste0(sample_id, "  -  ", ax1, " vs ", ax2)
        ) +
        ggplot2::theme_bw(base_size = 12)

      print(p)

      plot_name <- paste(pair, collapse = "_vs_")
      projection_plots[[plot_name]] <- p
    }

    dev.off()
    all_sample_results[[sample_id]] <- list(
      metrics = alignment_df,
      plots = list(
        alignment_summary = qc_result$plots$summary,
        alignment_heatmap = qc_result$plots$heatmap,
        scatter_by_axis = qc_result$plots$scatter_by_axis,
        projection_pairs = projection_plots
      ),
      files = list(
        alignment_qc_pdf = align_pdf,
        projection_pdf = proj_pdf,
        output_dir= output_dir
      )
    )
    message("  [DONE] ", sample_id)
  }

  # ---- aggregated alignment metrics ----------------------------------------

  all_metrics_df <- do.call(rbind, Filter(Negate(is.null), all_metrics))

  #agg_path <- file.path(models_dir,
  ##                      paste0(prefix, "all_samples_quality_check_metrics.csv"))
  #write.csv(all_metrics_df, agg_path, row.names = FALSE)

  #message("")
  message(strrep("=", 45))
  #message("  All ", n_models, " model(s) processed.")
  #message("  Aggregated metrics: ", agg_path)
  message("  Outputs in        : ", output_dir)
  message(strrep("=", 45))

  #invisible(all_metrics_df)
  invisible(list(
    metrics = all_metrics_df,
    plots = all_sample_results[[query_sample]]$plots,
    files = all_sample_results[[query_sample]]$files
  ))
}
