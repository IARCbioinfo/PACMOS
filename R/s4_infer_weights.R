#' Infer Archetype Weights
#'
#' @description
#' Computes archetype mixture weights for samples based on MOFA latent factors
#' derived from Step 3 (`s3_plot_query_samples_mofa`).
#'
#' @details
#' For each sample folder in `models_dir`, this function loads the
#' `<sample_id>_stable_input.csv` produced by
#' `s3_plot_query_samples_mofa()` and infers archetype mixture weights
#' using constrained least squares.
#'
#' Weights are computed for all samples in the `stable_input` matrix,
#' including both reference and query samples.
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
#' @param models_dir      Character. Root directory
#'
#' @param matrices_subdir Character. Folder name where `.RData` files are stored.
#'
#' @param coord           data.frame of archetype coordinates. Either:
#'   (A) first column = archetype names, remaining columns = coordinates, or
#'   (B) rownames = archetype names, all columns = coordinates.
#'
#' @param n_archetypes    Integer. Expected number of archetypes. Validated
#'   against \code{nrow(coord)}.
#'
#' @param out_dir         Output directory for aggregated CSV. Defaults to
#'   \code{models_dir/archetype_weights/}.
#'
#' @param prefix          Optional character prefix for output files.
#'
#' @return Invisibly returns a data.frame with columns
#'   \code{Sample} + one column per archetype (query samples only).
#'
#'
#' @export
infer_fuzzy_weights <- function(
    models_dir,
    matrices_subdir,
    coord,
    n_archetypes,
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
  sample_dirs <- sample_dirs[
    file.exists(
      file.path(sample_dirs, matrices_subdir, "plots",
                paste0(basename(sample_dirs), "_stable_input.csv"))
    )
  ]

  if (!length(sample_dirs))
    stop("No valid sample dirs found. Check models_dir and matrices_subdir.")

  message("Found ", length(sample_dirs), " sample dir(s).")

  all_weights <- vector("list", length(sample_dirs))

  # ---- main loop -----------------------------------------------------------

  for (i in seq_along(sample_dirs)) {

    sdir      <- sample_dirs[i]
    sample_id <- basename(sdir)

    message("")
    message(strrep("=", 45))
    message(sprintf("  [%d / %d]  %s", i, length(sample_dirs), sample_id))
    message(strrep("=", 45))

    si_path <- file.path(sdir, matrices_subdir, "plots",
                         paste0(sample_id, "_stable_input.csv"))

    stable_input <- read.csv(si_path, check.names = FALSE)
    sid_col      <- stable_input[[1]]

    # ---- solve weights for ALL rows in stable_input ----------------------

    all_rows_weights <- vector("list", nrow(stable_input))

    for (r in seq_len(nrow(stable_input))) {

      s_id <- sid_col[r]
      x_r  <- as.numeric(stable_input[r, -1])

      if (length(x_r) != ncol(Z)) {
        message("    [SKIP row] Dimension mismatch for: ", s_id)
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

    sample_plot_dir <- file.path(sdir, matrices_subdir, "plots")

    # ---- CSV 1: query sample only ----------------------------------------

    query_df <- all_samples_df[all_samples_df$Sample == sample_id, ,
                               drop = FALSE]

    if (nrow(query_df) == 0) {
      message("  [SKIP] '", sample_id, "' not found in stable_input rows.")
      next
    }

    sample_csv <- file.path(sample_plot_dir,
                            paste0(prefix, sample_id,
                                   "_archetype_weights.csv"))
    write.csv(query_df, sample_csv, row.names = FALSE)
    message("  Saved (query only) : ", basename(sample_csv))

    # ---- CSV 2: all samples in stable_input ------------------------------

    all_csv <- file.path(sample_plot_dir,
                         paste0(prefix, sample_id,
                                "_archetype_weights_all_samples.csv"))
    write.csv(all_samples_df, all_csv, row.names = FALSE)
    message("  Saved (all samples): ", basename(all_csv))

    # store query row for aggregation
    all_weights[[i]] <- c(Sample = sample_id,
                          setNames(as.numeric(query_df[, archetype_names]),
                                   archetype_names))
  }

  # ---- aggregated CSV (query rows only) ------------------------------------

  results_df <- do.call(
    rbind,
    lapply(Filter(Negate(is.null), all_weights), function(x) {
      as.data.frame(t(x), stringsAsFactors = FALSE)
    })
  )

  for (col in archetype_names)
    results_df[[col]] <- as.numeric(results_df[[col]])

  agg_csv <- file.path(out_dir,
                       paste0(prefix, matrices_subdir,
                              "_archetype_weights.csv"))
  write.csv(results_df, agg_csv, row.names = FALSE)

  message("")
  message(strrep("=", 45))
  message("  Processed  : ", nrow(results_df), " samples")
  message("  Aggregated : ", agg_csv)
  message(strrep("=", 45))

  invisible(results_df)
}
