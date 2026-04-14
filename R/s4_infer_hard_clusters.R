#' Infer K-means Cluster Labels
#'
#' For each sample folder in \code{models_dir}, loads the
#' \code{<sample_id>_stable_input.csv}, runs k-means clustering
#' on ALL samples, and saves one CSV per sample:
#' \itemize{
#'   \item \code{<sample_id>_sample_clusters.csv} - sample | kmeans_cluster
#' }
#' Cluster-to-label mapping can be done separately afterwards.
#'
#' @name infer_kmeans_clusters
#'
#' @param models_dir      Character. Root directory.
#'
#' @param matrices_subdir Character. Folder name where `.RData` files are stored.
#'
#' @param k               Integer. Number of k-means clusters.
#'
#' @param lf_cols         Character vector. LF columns to use as features for k-means.
#'
#' @param prefix          Optional character prefix for output files.
#'
#' @param seed            Integer. Random seed for reproducibility.
#'   Default \code{42}.
#'
#' @return Invisibly returns NULL.
#'
#' @export
infer_kmeans_clusters <- function(
    models_dir,
    matrices_subdir,
    k,
    lf_cols,
    prefix  = "",
    seed    = 42
) {

  # ---- checks --------------------------------------------------------------

  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)

  if (missing(k) || !is.numeric(k) || k < 2)
    stop("k must be an integer >= 2.")

  if (missing(lf_cols) || !is.character(lf_cols) || !length(lf_cols))
    stop("lf_cols must be a non-empty character vector of column names.")

  message("K-means settings - k: ", k,
          " | features: ", paste(lf_cols, collapse = ", "),
          " | seed: ", seed)

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

    stable_input      <- read.csv(si_path, check.names = FALSE)
    sid_col           <- as.character(stable_input[[1]])

    # ---- validate LF columns -----------------------------------------------

    missing_cols <- setdiff(lf_cols, colnames(stable_input))
    if (length(missing_cols)) {
      message("  [SKIP] Missing LF columns: ",
              paste(missing_cols, collapse = ", "))
      next
    }

    lf_data           <- stable_input[, lf_cols, drop = FALSE]
    rownames(lf_data) <- sid_col

    if (nrow(lf_data) < k) {
      message("  [SKIP] Fewer samples (", nrow(lf_data), ") than k=", k)
      next
    }

    # ---- k-means on ALL samples in stable_input ----------------------------

    set.seed(seed)
    km <- tryCatch(
      kmeans(lf_data, centers = k, nstart = 25, iter.max = 100),
      error = function(e) {
        message("  [SKIP] kmeans error: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(km)) next

    # ---- CSV: sample - kmeans_cluster --------------------------------------

    cluster_df <- data.frame(
      sample         = sid_col,
      kmeans_cluster = as.integer(km$cluster),
      stringsAsFactors = FALSE
    )

    csv_out <- file.path(sdir, matrices_subdir, "plots",
                         paste0(prefix, sample_id, "_sample_clusters.csv"))
    write.csv(cluster_df, csv_out, row.names = FALSE)
    message("  Saved: ", basename(csv_out))
  }

  message("")
  message(strrep("=", 45))
  message("  Done. Cluster CSVs saved per sample.")
  message(strrep("=", 45))

  invisible(NULL)
}
