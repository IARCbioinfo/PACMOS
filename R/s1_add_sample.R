#' Step 1: Add a query sample to reference MOFA input matrices
#'
#' @description
#' Reads one or more query sample CSV matrices and appends query sample(s)
#' to reference MOFA input matrices.
#'
#' @details
#' Each query matrix must correspond to a specific MOFA input layer. The mapping
#' between query matrices and MOFA layers is defined through the
#' `value_data_types` argument, where each entry corresponds to the same index
#' in `query_matrix_path`.
#'
#' The updated matrices are saved as individual `.RData` files that can be used
#' for retraining MOFA models.
#'
#' For the matched layer, real values from the query sample are inserted.
#' For all other MOFA layers, a column containing `NA` values is added for that
#' sample to maintain consistent sample structure across layers.
#'
#'
#' @name s1_add_sample_to_mofa
#'
#' @param query_matrix_path Character vector of CSV file paths containing query
#' sample matrices. First column is expected as gene_id, if not first column is
#' assumed as gene_id/feature_id (depending on the omic type).
#'
#' @param mofa_dir Character. Directory containing reference MOFA `.RData`
#' matrices (e.g. `D_expr_MOFA.RData`, `D_alt_MOFA.RData`). These objects must
#' have rownames corresponding to gene IDs/feature_ids (depending on the omic type).
#'
#' @param value_data_types Character vector maping each query matrix to its corresponding reference MOFA data layer.
#' Must be the same length as `query_matrix_path`. Each element should exactly match the name of a MOFA input object
#' as it appears when the reference .RData files are loaded into R.
#'
#' @param outdir Character. Root directory where output will be stored.
#'
#' @param python_bin Character. Path to the Python binary used by the MOFA
#' environment via the `reticulate` package.
#'
#' @return
#' Saves updated MOFA matrices as `.RData` files in `outdir`. One file is saved
#' per `(data layer, sample)` combination.
#'
#' @examples
#' mofa_dir <- system.file("extdata/MESOMICS_references", package = "PACMOS")
#' query_csv <- system.file("extdata/test_data", "MESOMICS_test_expr.csv", package = "PACMOS")
#' out_dir <- file.path(tempdir(), "pacmos_s1")
#' dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#' s1_add_sample_to_mofa(
#'   query_matrix_path = query_csv,
#'   mofa_dir = mofa_dir,
#'   value_data_types = "D_exprB_MOFA",
#'   outdir = out_dir,
#'   python_bin = Sys.which("/home/lipikal/miniconda3/envs/pacmos_env/bin/python")
#' )
#'
#' @export
s1_add_sample_to_mofa <- function(query_matrix_path,
                                  mofa_dir = system.file("extdata/", package = "PACMOS"),
                                  value_data_types = NULL,
                                  outdir = "output/",
                                  python_bin) {

  # ---- dependency check ------------------------------------------------------
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required.")

  reticulate::use_python(python_bin, required = TRUE)


  # ---- ensure correct input types -------------------------------------------
  query_matrix_path <- as.character(query_matrix_path)
  value_data_types <- as.character(value_data_types)


  # ---- create output directory if needed ------------------------------------
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)


  # ---- validate layer mapping ------------------------------------------------
  if (is.null(value_data_types) || !length(value_data_types))
    stop("value_data_types must be provided (same length as query_matrix_path).")


  if (length(query_matrix_path) != length(value_data_types)) {
    stop(
      "Input data mismatch: found ", length(query_matrix_path),
      " query_matrix_path item(s) but ",
      length(value_data_types), " value_data_types. ",
      "Each query matrix path must align by index with a data type."
    )
  }


  # ---- Load reference MOFA matrices -----------------------------------------
  mofa_inputs_sample <- list()

  rdata_files <- list.files(
    mofa_dir,
    pattern = "\\.RData$",
    full.names = TRUE
  )

  for (f in rdata_files) {

    tenv <- new.env(parent = emptyenv())
    obj_names <- load(f, envir = tenv)

    for (nm in obj_names) {

      obj <- get(nm, envir = tenv)

      # ensure data.frame structure
      mofa_inputs_sample[[nm]] <- as.data.frame(obj, check.names = FALSE)

    }
  }


  # ---- display pairing information ------------------------------------------
  cat("- Pairing (query_matrix_path -> data_type):\n")

  for (i in seq_along(query_matrix_path)) {

    cat(
      "  - [", i, "] ",
      basename(query_matrix_path[i]),
      " -> ",
      value_data_types[i],
      "\n",
      sep = ""
    )

  }


  # ---- read all query matrices first ----------------------------------------
  query_matrices <- list()

  for (i in seq_along(query_matrix_path)) {

    tm_path <- query_matrix_path[i]
    paired_type <- value_data_types[i]

    if (!file.exists(tm_path))
      stop("query matrix file does not exist: ", tm_path)

    query_matrix <- read.csv(
      tm_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (ncol(query_matrix) < 2) {
      stop(
        "query matrix '", tm_path,
        "' must have 'gene_id' + at least one sample column."
      )
    }

    # ---- detect gene_id column --------------------------------------------

    if ("gene_id" %in% colnames(query_matrix)) {

      # ensure it's first column
      query_matrix <- query_matrix[, c("gene_id",
                                       setdiff(colnames(query_matrix), "gene_id")),
                                   drop = FALSE]

    } else {

      # fallback: assume first column is gene_id
      colnames(query_matrix)[1] <- "gene_id"

      warning(
        "Column 'gene_id' not found in '", basename(tm_path),
        "'. Assuming first column is gene_id."
      )
    }

    query_matrices[[paired_type]] <- query_matrix

    sample_names <- setdiff(names(query_matrix), "gene_id")

    cat(
      "- Processing file #", i, ": ",
      basename(tm_path),
      " (", length(sample_names), " sample(s)) against data_type '",
      paired_type,
      "'\n",
      sep = ""
    )
  }


  # ---- identify all sample columns across all query matrices ----------------
  all_sample_names <- unique(unlist(
    lapply(query_matrices, function(x) setdiff(names(x), "gene_id"))
  ))


  # ---- iterate through each query sample ------------------------------------
  for (sample_name in all_sample_names) {

    mofa_inputs_current <- lapply(mofa_inputs_sample, function(x) x)

    sample_dir <- file.path(outdir, sample_name)
    inputs_dir <- file.path(sample_dir, "inputs")

    dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)

    cat("\n==============================\n")
    cat("Processing sample: ", sample_name, "\n")
    cat("==============================\n")


    # ---- iterate through all MOFA layers ------------------------------------
    for (data_type in names(mofa_inputs_sample)) {

      mat <- as.data.frame(
        mofa_inputs_sample[[data_type]],
        check.names = FALSE
      )


      # ---- matching layer: insert real values if provided -------------------
      if (data_type %in% names(query_matrices) &&
          sample_name %in% colnames(query_matrices[[data_type]])) {

        sample_column <- query_matrices[[data_type]][, c("gene_id", sample_name)]

        if (is.null(rownames(mat))) {
          stop(
            "MOFA object '",
            data_type,
            "' has no rownames to match against 'gene_id'."
          )
        }

        # align gene order
        idx <- match(rownames(mat), sample_column$gene_id)

        matched   <- !is.na(idx)
        n_total   <- nrow(mat)
        n_match   <- sum(matched)
        n_missing <- n_total - n_match

        order_ok <- (n_missing == 0) &&
          isTRUE(all(sample_column$gene_id[idx[!is.na(idx)]] == rownames(mat)[!is.na(idx)]))

        if (order_ok) {

          cat(
            "--[",
            data_type,
            " | ",
            sample_name,
            "] gene order OK (",
            n_match,
            "/",
            n_total,
            " matched)\n",
            sep = ""
          )

        } else {

          cat(
            "-- [",
            data_type,
            " | ",
            sample_name,
            "] reordered to MOFA row order; ",
            n_missing,
            " genes missing in query (filled NA)\n",
            sep = ""
          )

        }

        mat[[sample_name]] <- sample_column[idx, sample_name]

      } else {

        # ---- non-matching layers: add NA if column absent -------------------
        if (!(sample_name %in% colnames(mat))) {
          mat[[sample_name]] <- NA_real_
        }

      }

      mofa_inputs_current[[data_type]] <- mat

      # ---- save matrix -------------------------------------------------------
      obj_to_save <- data_type

      assign(obj_to_save, mat)

      save(
        list = obj_to_save,
        file = file.path(
          inputs_dir,
          paste0(data_type, "_", sample_name, ".RData")
        )
      )

      cat("   |_ Saved: ", data_type, "_", sample_name, ".RData\n", sep = "")
    }

    cat(
      "\n Completed sample: ", sample_name,
      "\n\n",
      sep = ""
    )
  }

}
