#' Step 2: Train MOFA model for reference + query sample
#'
#' @description
#' This function trains a MOFA2 model using the matrices generated in
#' \code{s1_add_sample_to_mofa()}.
#'
#' @details
#' The function loads layer-specific `.RData` matrices, constructs a MOFA object,
#' trains the model, and writes the resulting model to an HDF5 file.
#'
#' Each sample is processed independently. For each sample directory, the
#' corresponding matrices are loaded, aligned across views, and used to train
#' a MOFA model.
#'
#' If `skip_existing = TRUE`, samples for which a model file already exists
#' will be skipped.
#'
#'
#' @name s2_run_mofa
#'
#' @param models_dir Root folder containing `.RData` matrices created by
#' `s1_add_sample_to_mofa()`.
#'
#' @param matrices_subdir Folder name where `.RData` files are stored.
#'
#' @param num_factors Integer. Number of latent factors.
#'
#' @param convergence_mode MOFA convergence mode (`"slow"` or `"fast"`).
#'
#' @param maxiter Maximum number of training iterations.
#'
#' @param use_basilisk Logical. Whether to run MOFA inside basilisk environment.
#'
#' @param skip_existing Logical. If TRUE, skip samples for which output
#' HDF5 already exists.
#'
#' @param python_bin Path to the Python executable used by MOFA.
#'
#' @param outfile_prefix Optional prefix for output files.
#'
#' @param views_map Named character vector mapping view names to matrix
#' object names.
#'
#' Example:
#' \preformatted{
#' views_map = c(
#'   RNA = "D_expr_MOFA",
#'   CNV = "D_cnv_MOFA",
#'   ALT = "D_alt_MOFA"
#' )
#' }
#'
#' @param binary_views Character vector specifying which views contain
#' binary data (e.g. mutation or alteration layers).
#'
#' @return Invisibly returns TRUE when MOFA training completes. Trained MOFA models are
#' written to disk as `.hdf5` files.
#'
#'
#' @export
s2_run_mofa <- function(
    models_dir,
    matrices_subdir  = "train_query_all_omics",
    num_factors      = 10,
    convergence_mode = "slow",
    maxiter          = 10000,
    use_basilisk     = FALSE,
    skip_existing    = TRUE,
    python_bin,
    outfile_prefix   = "",
    views_map,
    binary_views     = NULL
) {

  # ---- Check directories ------------------------------------------------------
  if (!dir.exists(models_dir))
    stop("models_dir does not exist: ", models_dir)


  # ---- dependency check ------------------------------------------------------
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required.")

  reticulate::use_python(python_bin, required = TRUE)

  if (is.null(views_map) || is.null(names(views_map)))
    stop("views_map must be a named character vector")


  # ---- Discover sample subdirectories ------------------------------------------------------
  sample_dirs <- list.dirs(models_dir, recursive = FALSE, full.names = TRUE)

  if (length(sample_dirs) == 0)
    stop("No sample subdirectories found in: ", models_dir)

  message("Found ", length(sample_dirs), " sample(s) in: ", models_dir)

  results <- list()


  # ---- loop over samples ------------------------------------------------------
  for (sdir in sample_dirs) {

    sample_name <- basename(sdir)
    inputs_dir  <- file.path(sdir, matrices_subdir)   # sampleA/train_plus_query*/

    if (!dir.exists(inputs_dir)) {
      warning("Skipping '", sample_name, "': no '", matrices_subdir, "/' subfolder found.")
      results[[sample_name]] <- "SKIPPED"
      next
    }

    # ---- outfile saved INSIDE the sample's train_plus_query folder ------------
    outfile <- file.path(
      inputs_dir,                                      # sampleA/train_plus_query*/
      paste0(outfile_prefix, "MOFA-", sample_name, ".hdf5")
    )

    if (skip_existing && file.exists(outfile)) {
      message("Skipping existing model: ", outfile)
      results[[sample_name]] <- "SKIPPED"
      next
    }

    cat("\n==============================\n")
    cat("Processing sample: ", sample_name, "\n")
    cat("==============================\n")

    tryCatch({

      # ---- Locate and load matrices ------------------------------------------------------
      rdata_files <- list.files(inputs_dir, pattern = "\\.RData$", full.names = TRUE)

      if (length(rdata_files) == 0)
        stop("No .RData matrices found in: ", inputs_dir)

      tenv <- new.env(parent = emptyenv())
      for (f in rdata_files)
        load(f, envir = tenv)


      # ---- Build view list ------------------------------------------------------
      view_list <- list()

      for (view_nm in names(views_map)) {

        obj_nm <- views_map[[view_nm]]

        if (!exists(obj_nm, envir = tenv)) {
          message("View not found: ", view_nm)
          next
        }

        m <- get(obj_nm, envir = tenv)
        if (!is.matrix(m)) m <- as.matrix(m)
        storage.mode(m) <- "double"
        view_list[[view_nm]] <- m
      }

      if (length(view_list) == 0)
        stop("No valid views found for sample: ", sample_name)


      # ---- Ensure identical sample columns across views ------------------------------------------------------
      all_cols <- Reduce(union, lapply(view_list, colnames))

      view_list <- lapply(view_list, function(m) {
        missing_cols <- setdiff(all_cols, colnames(m))
        if (length(missing_cols)) {
          na_mat <- matrix(NA, nrow = nrow(m), ncol = length(missing_cols),
                           dimnames = list(rownames(m), missing_cols))
          m <- cbind(m, na_mat)
        }
        m[, all_cols, drop = FALSE]
      })


      # ---- Create and train MOFA model ------------------------------------------------------
      MOFAobject <- MOFA2::create_mofa(view_list)

      data_opts  <- MOFA2::get_default_data_options(MOFAobject)
      model_opts <- MOFA2::get_default_model_options(MOFAobject)
      model_opts$num_factors <- as.integer(num_factors)

      if (!is.null(binary_views)) {
        for (bv in binary_views)
          if (bv %in% names(view_list))
            model_opts$likelihoods[bv] <- "bernoulli"
      }

      train_opts                  <- MOFA2::get_default_training_options(MOFAobject)
      train_opts$convergence_mode <- convergence_mode
      train_opts$maxiter          <- as.integer(maxiter)

      MOFAobject <- MOFA2::prepare_mofa(MOFAobject, data_opts, model_opts, train_opts)

      MOFA2::run_mofa(MOFAobject, outfile = outfile, use_basilisk = use_basilisk)

      message("- Saved: ", outfile)
      results[[sample_name]] <- "OK"

    }, error = function(e) {
      warning("Failed for sample '", sample_name, "': ", e$message)
      results[[sample_name]] <<- paste0("ERROR: ", e$message)
    })
  }

  # ----Summary ------------------------------------------------------
  message("\n- Run summary -")
  for (nm in names(results))
    message("  ", nm, ": ", results[[nm]])

  invisible(results)
}
