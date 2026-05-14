#' Infer Archetype Weights
#'
#' @description
#' Computes archetype mixture weights for samples based on MOFA latent factors
#' derived from Step 3 (`s3_plot_query_samples_mofa`).
#'
#' @details
#' For each sample folder in `models_dir`, this function loads the aligned
#' latent factor matrix produced by `s3_plot_query_samples_mofa()` —
#' either `<sample_id>_stable_input.csv` or
#' `<sample_id>_retrained_LFs_all_samples.csv` (controlled via `input_type`) —
#' and infers archetype mixture weights using constrained least squares.
#'
#'
#' For each sample, two CSV files are written:
#' \itemize{
#'   \item `<sample_id>_archetype_weights.csv` — query sample only
#'   \item `<sample_id>_archetype_weights_all_samples.csv` — all samples
#' }
#'
#' Additionally, an aggregated CSV (query samples only) is written to `out_dir`.
#'
#' @name infer_fuzzy_weights
#'
#' @param models_dir      Character. Root directory folder. Same as `s1_add_sample_to_mofa() outdir`.
#'
#' @param coord           data.frame of archetype coordinates. Either:
#'   (A) first column = archetype names, remaining columns = coordinates, or
#'   (B) rownames = archetype names, all columns = coordinates.
#'
#' @param n_archetypes    Integer. Expected number of archetypes. Validated
#'   against \code{nrow(coord)}.
#'
#' @param input_type Character. Which aligned matrix to use as input.
#'   Either \code{"stable"} (reference + query; from \code{_stable_input.csv})
#'   or \code{"retrained"} (all retrained model samples; from
#'   \code{_retrained_LFs_all_samples.csv}).
#'
#' @param reference_axes Character vector specifying which latent factor
#'   columns (e.g. "Factor1", "Factor2") to use for archetype weight inference.
#'   These must match column names in the input matrices and correspond to the
#'   same dimensions used to define the archetype coordinates in \code{coord}.
#'   If \code{NULL}, all numeric latent factor columns are used.
#'
#' @param out_dir         Output directory for aggregated CSV. Defaults to
#'   \code{models_dir/archetype_weights/}.
#'
#' @param prefix          Optional character prefix for output files.
#'
#' @return Invisibly returns a data.frame with columns
#'   \code{Sample} + one column per archetype.
#'
#' @examples
#' mofa_dir <- system.file("extdata/MESOMICS_references", package = "PACMOS")
#' reference_LFs <- system.file(
#'   "extdata/MESOMICS_references",
#'   "MESOMICS_latent_factors.csv",
#'   package = "PACMOS"
#' )
#' query_csv <- system.file(
#'   "extdata/test_data",
#'   "MESOMICS_test_expr.csv",
#'   package = "PACMOS"
#' )
#' out_dir <- file.path(tempdir(), "pacmos_s4")
#' dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#'
#' s1_add_sample_to_mofa(
#'   query_matrix_path = query_csv,
#'   mofa_dir = mofa_dir,
#'   value_data_types = "D_exprB_MOFA",
#'   outdir = out_dir,
#'   python_bin = Sys.which("/home/lipikal/miniconda3/envs/pacmos_env/bin/python")
#' )
#'
#' s2_run_mofa(
#'   models_dir = out_dir,
#'   matrices_subdir = "inputs",
#'   num_factors = 2,
#'   convergence_mode = "fast",
#'   maxiter = 5,
#'   use_basilisk = FALSE,
#'   skip_existing = TRUE,
#'   python_bin = Sys.which("/home/lipikal/miniconda3/envs/pacmos_env/bin/python"),
#'   views_map = c(RNA = "D_exprB_MOFA"),
#'   binary_views = NULL
#' )
#'
#' sample_id <- basename(list.dirs(out_dir, recursive = FALSE, full.names = TRUE)[1])
#'
#' s3_plot_query_samples_mofa(
#'   models_dir = out_dir,
#'   matrices_subdir = "inputs",
#'   query_sample = sample_id,
#'   reference_LFs = reference_LFs,
#'   reference_axes = c("Morphology_LF", "Adaptive-response_LF"),
#'   group = "group1",
#'   python_bin = Sys.which("/home/lipikal/miniconda3/envs/pacmos_env/bin/python")
#' )
#'
#' archetype_coords <- data.frame(
#'   Archetype = c("Cell division", "Tumor-immune-interaction", "Acinar"),
#'   Morphology_LF = c(-3.5973, -1.7981, 3.8586),
#'   `Adaptive-response_LF` = c(-2.3960, 3.5635, -0.8256),
#'   stringsAsFactors = FALSE
#' )
#'
#' infer_fuzzy_weights(
#'   models_dir = out_dir,
#'   coord = archetype_coords,
#'   n_archetypes = 3,
#'   reference_axes = c("Morphology_LF", "Adaptive-response_LF"),
#'   input_type = "stable"
#' )
#'
#' @export
infer_fuzzy_weights <- function(
    models_dir,
    coord,
    n_archetypes,
    reference_axes = NULL,
    input_type = c("stable", "retrained"),
    out_dir  = file.path(models_dir, "archetype_weights"),
    prefix   = ""
) {

  # ---- checks --------------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (!requireNamespace("limSolve", quietly = TRUE))
    stop("Package 'limSolve' is required. Please install it.")

  coord <- as.data.frame(coord)

  if (is.character(coord[[1]]) || is.factor(coord[[1]])) {
    archetype_names <- as.character(coord[[1]])
    Z               <- as.matrix(coord[, -1])
  } else {
    archetype_names <- rownames(coord)
    Z               <- as.matrix(coord)
  }

  if (is.null(archetype_names) || any(is.na(archetype_names)))
    stop("Could not determine archetype names from coord. ",
         "Provide names as first column or rownames.")

  if (!missing(n_archetypes) && nrow(Z) != n_archetypes)
    stop("n_archetypes = ", n_archetypes,
         " but coord has ", nrow(Z), " rows.")

  k <- nrow(Z)
  message("Archetypes (", k, "): ", paste(archetype_names, collapse = ", "))

  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # ---- inline solver -------------------------------------------------------

  .solve_weights <- function(x, Z_t, k, archetype_names) {
    E    <- matrix(1, nrow = 1, ncol = k)
    f_eq <- 1
    G    <- diag(k)
    h    <- rep(0, k)
    res  <- limSolve::lsei(
      A = Z_t, B = x,
      E = E,   F = f_eq,
      G = G,   H = h,
      verbose = FALSE
    )
    alpha        <- res$X
    names(alpha) <- archetype_names
    alpha
  }

  Z_t <- t(Z)   # dim x k

  # ---- discover valid sample dirs ------------------------------------------

  sample_dirs <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)

  input_type <- match.arg(input_type)

  # map input_type to the correct file suffix
  file_suffix <- switch(input_type,
                        stable    = "_stable_input.csv",
                        retrained = "_retrained_LFs_all_samples.csv"
  )

  sample_dirs <- sample_dirs[
    file.exists(
      file.path(sample_dirs, "outputs",
                paste0(basename(sample_dirs), file_suffix))
    )
  ]

  if (!length(sample_dirs))
    stop("No valid sample dirs found. Check models_dir and sample outputs folders.")

  message("Found ", length(sample_dirs), " sample dir(s).")

  all_weights <- vector("list", length(sample_dirs))

  # ---- main loop -----------------------------------------------------------

  for (i in seq_along(sample_dirs)) {

    tryCatch(
      {
        sdir      <- sample_dirs[i]
        sample_id <- basename(sdir)

        message("")
        message(strrep("=", 45))
        message(sprintf("  [%d / %d]  %s", i, length(sample_dirs), sample_id))
        message(strrep("=", 45))

        si_path <- file.path(
          sdir, "outputs",
          paste0(sample_id, file_suffix)
        )

        input_df <- read.csv(si_path, check.names = FALSE)

        if (!("Sample" %in% colnames(input_df))) {
          stop("Column 'Sample' not found in: ", basename(si_path))
        }

        if (is.null(reference_axes)) {
          numeric_cols <- vapply(input_df, is.numeric, logical(1))
          numeric_cols["Sample"] <- FALSE
          input_df <- input_df[, c("Sample", names(input_df)[numeric_cols]), drop = FALSE]
        } else {
          missing_axes <- setdiff(reference_axes, colnames(input_df))
          if (length(missing_axes)) {
            stop("reference_axes not found in ", basename(si_path), ": ",
                 paste(missing_axes, collapse = ", "))
          }
          input_df <- input_df[, c("Sample", reference_axes), drop = FALSE]
        }

        if ((ncol(input_df) - 1) != ncol(Z)) {
          stop(
            "Dimension mismatch in ", basename(si_path), ": input has ",
            ncol(input_df) - 1, " feature(s) but coord expects ", ncol(Z),
            ". Check reference_axes or coord."
          )
        }

        sid_col <- input_df[["Sample"]]

        # ---- solve weights for ALL rows in stable_input ----------------------

        all_rows_weights <- vector("list", nrow(input_df))

        for (r in seq_len(nrow(input_df))) {

          s_id <- sid_col[r]
          x_r  <- as.numeric(input_df[r, -1])

          if (length(x_r) != ncol(Z)) {
            message("    [SKIP row] Dimension mismatch for: ", s_id,
                    " (length=", length(x_r), ", expected=", ncol(Z), ")")
            next
          }

          alpha_r <- tryCatch(
            .solve_weights(x_r, Z_t, k, archetype_names),
            error = function(e) {
              message("    [SKIP row] Solver error for '", s_id, "': ",
                      conditionMessage(e))
              NULL
            }
          )

          if (is.null(alpha_r)) next

          row_df <- as.data.frame(t(c(Sample = s_id, alpha_r)),
                                  stringsAsFactors = FALSE)
          for (col in archetype_names)
            row_df[[col]] <- as.numeric(row_df[[col]])

          all_rows_weights[[r]] <- row_df
        }

        all_samples_df <- do.call(rbind,
                                  Filter(Negate(is.null), all_rows_weights))

        if (is.null(all_samples_df) || nrow(all_samples_df) == 0) {
          message("  [SKIP] No weights computed for any sample in: ", sample_id)
          next
        }

        sample_plot_dir <- file.path(sdir, "outputs")

        # ---- CSV 1: query sample only ----------------------------------------

        query_df <- all_samples_df[all_samples_df$Sample == sample_id, ,
                                   drop = FALSE]

        if (nrow(query_df) == 0) {
          message("  [SKIP] '", sample_id, "' not found in ",
                  basename(si_path), " rows.")
          next
        }

        sample_csv <- file.path(sample_plot_dir,
                                paste0(sample_id,
                                       "_archetype_weights.csv"))
        write.csv(query_df, sample_csv, row.names = FALSE)
        message("  Saved (query only) : ", basename(sample_csv))

        # ---- CSV 2: all samples in stable_input ------------------------------

        all_csv <- file.path(sample_plot_dir,
                             paste0(sample_id,
                                    "_archetype_weights_all_samples.csv"))
        write.csv(all_samples_df, all_csv, row.names = FALSE)
        message("  Saved (all samples): ", basename(all_csv))

        # store query row for aggregation
        all_weights[[i]] <- c(Sample = sample_id,
                              setNames(as.numeric(query_df[, archetype_names]),
                                       archetype_names))
      }, error = function(e) {
        warning("Error in sample ", basename(sample_dirs[i]), ": ",
                conditionMessage(e))
      })

  }

  # ---- aggregated CSV (query rows only) ------------------------------------

  results_df <- do.call(
    rbind,
    lapply(Filter(Negate(is.null), all_weights), function(x) {
      as.data.frame(t(x), stringsAsFactors = FALSE)
    })
  )

  if (is.null(results_df) || nrow(results_df) == 0) {
    stop("No archetype weights were computed for any sample.")
  }

  for (col in archetype_names)
    results_df[[col]] <- as.numeric(results_df[[col]])


  agg_csv <- file.path(out_dir,
                       paste0(prefix,
                              "_archetype_weights.csv"))
  write.csv(results_df, agg_csv, row.names = FALSE)

  message("")
  message(strrep("=", 45))
  message("  Processed  : ", nrow(results_df), " samples")
  message("  Aggregated : ", agg_csv)
  message(strrep("=", 45))

  invisible(results_df)
}
