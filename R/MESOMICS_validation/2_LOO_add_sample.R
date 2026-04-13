run_add_sample_to_mofa_MESOMICS_validation <- function(
    input_dirs,
    layers = "all",
    output_subdir = "train_test_all_omics",
    python_bin="/home/lipikal/miniconda3/envs/mofa_env/bin/python"
){

  for(sub_input_dir in input_dirs){

    cat("=====================================\n")
    cat("Processing:", sub_input_dir, "\n")


    ############################################################
    # Directory containing train/test matrices
    ############################################################

    dir_needed <- file.path(sub_input_dir, "inputs")

    if(!dir.exists(dir_needed)){
      cat("Directory missing:", dir_needed, "\n")
      next
    }


    ############################################################
    # Create output directory
    ############################################################

    outdir <- file.path(sub_input_dir, output_subdir)

    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


    ############################################################
    # List training matrices (.RData)
    ############################################################

    list_train_files <- list.files(
      dir_needed,
      pattern="^D_.*_MOFA\\.RData$",
      full.names=TRUE
    )


    ############################################################
    # Extract layer names from RData
    ############################################################

    layer_names <- gsub("\\.RData$", "", basename(list_train_files))


    ############################################################
    # Subset layers if user specified
    ############################################################

    if(!identical(layers, "all")){
      layer_names <- intersect(layer_names, layers)
    }


    if(length(layer_names) == 0){
      cat("No valid layers found\n")
      next
    }


    ############################################################
    # Detect sample name
    ############################################################

    sample_name <- basename(normalizePath(sub_input_dir))

    ############################################################
    # Match test CSV files to layers
    ############################################################

    list_test_files <- as.character(sapply(layer_names, function(layer){

      f <- list.files(
        dir_needed,
        pattern = paste0("^", layer, "_", sample_name, "_test\\.csv$"),
        full.names = TRUE
      )

      if(length(f) == 0){
        return(NA)
      }

      return(f[1])

    }))

    ############################################################
    # Remove layers without test files
    ############################################################

    missing_layers <- layer_names[is.na(list_test_files)]

    if(length(missing_layers)){
      cat("Missing test matrices for layers:",
          paste(missing_layers, collapse=", "), "\n")
    }

    valid <- !is.na(list_test_files)

    layer_names <- layer_names[valid]
    list_test_files <- list_test_files[valid]

    if(length(list_test_files) == 0){
      cat("No matching test CSV files found\n")
      next
    }


    cat("Layers used:", paste(layer_names, collapse=", "), "\n")


    ############################################################
    # Run add_sample_to_mofa
    ############################################################

    s1_add_sample_to_mofa(
      test_matrix_path = list_test_files,
      mofa_dir = dir_needed,
      value_data_types = layer_names,
      outdir = outdir,
      python_bin = python_bin
    )


    cat("Finished:", sub_input_dir, "\n\n")

  }

  cat("=====================================\n")
  cat("All folds processed.\n")
}

input_dirs <- list.dirs(
  "Analysis_010426/Mesomics_LOO",
  recursive = FALSE,
  full.names = TRUE
)

#change this for different omics combination
run_add_sample_to_mofa_MESOMICS_validation(
  input_dirs = input_dirs,
  layers =c('D_exprB_MOFA', 'D_cnv_MOFA', 'D_alt_MOFA', 'D_loh_MOFA'),
  output_subdir = "train_test_Expr_Genomic_only"
)
