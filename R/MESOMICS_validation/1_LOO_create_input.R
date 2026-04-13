load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_met.proB_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_met.enhB_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_met.bodB_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_loh_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_exprB_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_cnv_MOFA.RData")
load("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/D_alt_MOFA.RData")

gene_counts <- read.csv(
  "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/inst/extdata/MESOMICS_references/gene_count_matrix_1pass.csv",
  row.names = 1
)

gene_counts <- as.matrix(gene_counts)
storage.mode(gene_counts) <- "numeric"


library(DESeq2)

############################################################
# Function 1: Gene Expression processing
############################################################

process_expression_layer <- function(expr_counts, train_samples, test_sample, gene_set) {

  # ── TRAIN ───────────────────────────────────────────────────────────────────

  train_expr_samples <- intersect(train_samples, colnames(expr_counts))
  counts_train       <- expr_counts[, train_expr_samples, drop = FALSE]

  col_train <- data.frame(
    dummy     = rep(1, length(train_expr_samples)),
    row.names = train_expr_samples
  )

  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_train),
    colData   = col_train,
    design    = ~1
  )

  # size factors from training samples
  dds_train <- DESeq2::estimateSizeFactors(dds_train)

  # geometric means from training (used to anchor test size factors)
  counts.train   <- DESeq2::counts(dds_train)
  geoMeans_train <- apply(counts.train, 1, function(x) {
    if (all(x == 0)) return(NA_real_)
    exp(mean(log(x[x > 0])))
  })

  # dispersions on training — stored inside dds_train
  dds_train <- DESeq2::estimateDispersionsGeneEst(dds_train, quiet = TRUE)
  dds_train <- DESeq2::estimateDispersionsFit(dds_train,     quiet = TRUE)

  # RUN VST using frozen dispersion — blind=FALSE reads stored dispersionFunction
  vsd_train_full <- SummarizedExperiment::assay(
    DESeq2::varianceStabilizingTransformation(dds_train, blind = FALSE)
  )

  # scaling parameters from TRAIN only
  mu  <- rowMeans(vsd_train_full)
  sdv <- apply(vsd_train_full, 1, sd)
  sdv[sdv == 0 | is.na(sdv)] <- 1

  train_scaled_full <- t(scale(t(vsd_train_full), center = mu, scale = sdv))

  # ── TEST ────────────────────────────────────────────────────────────────────

  if (test_sample %in% colnames(expr_counts)) {

    counts_test <- expr_counts[, test_sample, drop = FALSE]

    col_test <- data.frame(dummy = 1, row.names = test_sample)

    dds_test <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(counts_test),
      colData   = col_test,
      design    = ~1
    )

    counts.test <- DESeq2::counts(dds_test)

    # size factor anchored to training geometric means
    valid <- !is.na(geoMeans_train) & geoMeans_train > 0
    sf    <- DESeq2::estimateSizeFactorsForMatrix(
      counts.test[valid, , drop = FALSE],
      geoMeans = geoMeans_train[valid]
    )
    sizeFactors(dds_test) <- sf

    # transfer frozen dispersion function from training
    dispersionFunction(dds_test) <- dispersionFunction(dds_train)

    # RUN VST using frozen dispersion — blind=FALSE reads stored dispersionFunction
    vsd_test_full <- SummarizedExperiment::assay(
      DESeq2::varianceStabilizingTransformation(dds_test, blind = FALSE)
    )

    # align gene order to training before scaling
    vsd_test_full <- vsd_test_full[rownames(vsd_train_full), , drop = FALSE]

    # scale with frozen TRAIN mu/sdv
    test_scaled_full <- t(scale(t(vsd_test_full), center = mu, scale = sdv))

  } else {

    test_scaled_full <- matrix(
      NA,
      nrow     = nrow(train_scaled_full),
      ncol     = 1,
      dimnames = list(rownames(train_scaled_full), test_sample)
    )
  }

  # ── Restrict to MOFA gene set (after normalization) ─────────────────────────

  gene_set <- intersect(gene_set, rownames(train_scaled_full))

  train_scaled <- train_scaled_full[gene_set, , drop = FALSE]
  test_scaled  <- test_scaled_full[gene_set, , drop = FALSE]

  list(
    train = train_scaled,
    test  = test_scaled
  )
}

############################################################
# Function 2: Other omic layer processing
############################################################

process_continuous_layer <- function(mat, train_samples, test_sample, scale_layer = FALSE){

  if(is.null(mat)){
    return(list(train=NULL, test=NULL))
  }

  # Identify training samples present in this layer
  train_present <- intersect(train_samples, colnames(mat))

  train <- mat[, train_present, drop = FALSE]

  # Compute scaling parameters from TRAIN only
  if(scale_layer){

    mu  <- rowMeans(train, na.rm = TRUE)
    sdv <- apply(train, 1, sd, na.rm = TRUE)

    sdv[sdv == 0 | is.na(sdv)] <- 1

    train <- t(scale(t(train), center = mu, scale = sdv))
  }

  # Add NA columns for missing training samples
  missing_train <- setdiff(train_samples, colnames(train))

  if(length(missing_train) > 0){

    na_mat <- matrix(
      NA,
      nrow = nrow(mat),
      ncol = length(missing_train),
      dimnames = list(rownames(mat), missing_train)
    )

    train <- cbind(train, na_mat)
  }

  # Ensure column order matches train_samples
  train <- train[, train_samples, drop = FALSE]

  # Extract TEST sample
  if(test_sample %in% colnames(mat)){

    test <- mat[, test_sample, drop = FALSE]

    if(scale_layer){
      test <- t(scale(t(test), center = mu, scale = sdv))
    }

  } else {

    test <- matrix(
      NA,
      nrow = nrow(mat),
      ncol = 1,
      dimnames = list(rownames(mat), test_sample)
    )
  }

  # Return processed matrices
  list(train = train, test = test)
}

############################################################
# Function 3: split binary layers (CNV / LOH / ALT)
############################################################

process_binary_layer <- function(mat, train_samples, test_sample){

  if(is.null(mat)){
    return(list(train=NULL, test=NULL))
  }

  # ----------------------------
  # TRAIN
  # ----------------------------

  train <- mat[, intersect(train_samples, colnames(mat)), drop = FALSE]

  # ----------------------------
  # TEST
  # ----------------------------

  if(test_sample %in% colnames(mat)){

    test <- mat[, test_sample, drop = FALSE]

  } else {

    # create NA column if sample missing

    test <- matrix(
      NA,
      nrow = nrow(mat),
      ncol = 1,
      dimnames = list(rownames(mat), test_sample)
    )
  }

  list(train=train, test=test)
}
############################################################
# Main LOO pipeline
############################################################
expr_gene_set <- rownames(D_exprB_MOFA)

generate_LOO_inputs_MESOMICS <- function(
    expr_counts = NULL,
    met_bod = NULL,
    met_enh = NULL,
    met_pro = NULL,
    cnv = NULL,
    loh = NULL,
    alt = NULL,
    output_dir = "Mesomics_LOO"
){

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ############################################################
  # Collect all omic layers
  ############################################################

  omics_list <- list(
    expr = expr_counts,
    met_bod = met_bod,
    met_enh = met_enh,
    met_pro = met_pro,
    cnv = cnv,
    loh = loh,
    alt = alt
  )

  omics_list <- omics_list[!sapply(omics_list, is.null)]

  ############################################################
  # Get UNION of samples across all layers
  ############################################################

  samples_all <- sort(unique(unlist(lapply(omics_list, colnames))))

  presence_matrix <- sapply(omics_list, function(mat) {
    samples_all %in% colnames(mat)
  })

  presence_df <- as.data.frame(presence_matrix)
  rownames(presence_df) <- samples_all

  presence_counts <- rowSums(presence_df)

  # Keep samples present in >= 2 layers
  samples <- samples_all[presence_counts >= 2]

  ############################################################
  # LOO loop
  ############################################################

  for(test_sample in samples){

    cat("Processing sample:", test_sample, "\n")

    train_samples <- setdiff(samples, test_sample)

    ##########################################################
    # Create sample directory
    ##########################################################

    sample_dir <- file.path(output_dir, test_sample)
    input_dir  <- file.path(sample_dir, "inputs")

    dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

    ##########################################################
    # Containers for train/test layers
    ##########################################################

    train_list <- list()
    test_list  <- list()

    ##########################################################
    # Expression processing
    ##########################################################

    if(!is.null(expr_counts)){

      expr_processed <- process_expression_layer(
        expr_counts,
        train_samples,
        test_sample,
        gene_set = expr_gene_set
      )

      train_list$D_exprB_MOFA <- expr_processed$train
      test_list$D_exprB_MOFA  <- expr_processed$test
    }

    ##########################################################
    # Methylation layers (continuous scaling)
    ##########################################################

    if(!is.null(met_bod)){

      res <- process_continuous_layer(met_bod, train_samples, test_sample, scale_layer = TRUE)

      train_list$D_met.bodB_MOFA <- res$train
      test_list$D_met.bodB_MOFA  <- res$test
    }

    if(!is.null(met_enh)){

      res <- process_continuous_layer(met_enh, train_samples, test_sample, scale_layer = TRUE)

      train_list$D_met.enhB_MOFA <- res$train
      test_list$D_met.enhB_MOFA  <- res$test
    }

    if(!is.null(met_pro)){

      res <- process_continuous_layer(met_pro, train_samples, test_sample, scale_layer = TRUE)

      train_list$D_met.proB_MOFA <- res$train
      test_list$D_met.proB_MOFA  <- res$test
    }

    ##########################################################
    # CNV / LOH / ALT (no scaling, just split)
    ##########################################################

    if(!is.null(cnv)){

      res <- process_binary_layer(
        cnv,
        train_samples,
        test_sample
      )

      train_list$D_cnv_MOFA <- res$train
      test_list$D_cnv_MOFA  <- res$test
    }

    if(!is.null(loh)){

      res <- process_binary_layer(
        loh,
        train_samples,
        test_sample
      )

      train_list$D_loh_MOFA <- res$train
      test_list$D_loh_MOFA  <- res$test
    }

    if(!is.null(alt)){

      res <- process_binary_layer(
        alt,
        train_samples,
        test_sample
      )

      train_list$D_alt_MOFA <- res$train
      test_list$D_alt_MOFA  <- res$test
    }

    ##########################################################
    # Save TRAIN matrices as individual RData
    ##########################################################

    for(layer in names(train_list)){

      train_mat <- train_list[[layer]]

      assign(layer, train_mat)

      save(
        list = layer,
        file = file.path(input_dir, paste0(layer, ".RData"))
      )
    }

    ##########################################################
    # Save TEST matrices as CSV
    ##########################################################

    for(layer in names(test_list)){

      test_mat <- test_list[[layer]]

      df <- data.frame(
        gene_id = rownames(test_mat),
        value = as.numeric(test_mat[,1]),
        check.names = FALSE
      )

      colnames(df)[2] <- colnames(test_mat)

      write.csv(
        df,
        file = file.path(
          input_dir,
          paste0(layer, "_", test_sample, "_test.csv")
        ),
        row.names = FALSE
      )
    }

    cat("Saved inputs for:", test_sample, "\n\n")
  }

  cat("All LOO inputs generated.\n")
}

generate_LOO_inputs_MESOMICS(
  expr_counts = gene_counts,
  met_bod = D_met.bodB_MOFA,
  met_enh = D_met.enhB_MOFA,
  met_pro = D_met.proB_MOFA,
  cnv = D_cnv_MOFA,
  loh = D_loh_MOFA,
  alt = D_alt_MOFA,

  output_dir = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/Mesomics_LOO"
)
