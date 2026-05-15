## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github(
#    "IARCbioinfo/PACMOS",
#    dependencies = TRUE
#  )

## ----set-python, eval=FALSE---------------------------------------------------
#  ## Linux/macOS
#  Sys.setenv(
#    PACMOS_PYTHON =
#      "/home/lipikal/miniconda3/envs/pacmos_env/bin/python"
#  ) # replace with your path
#  
#  ## Windows
#  Sys.setenv(
#    PACMOS_PYTHON =
#      "C:/Users/YOUR_USERNAME/miniconda3/envs/pacmos_env/python.exe"
#  ) # replace with your path
#  
#  file.exists(Sys.getenv("PACMOS_PYTHON"))
#  # [1] TRUE

## ----build-vignettes, eval=FALSE----------------------------------------------
#  devtools::build_vignettes()

## ----setup--------------------------------------------------------------------
library(PACMOS)

## ----setup-paths--------------------------------------------------------------
mofa_dir <- system.file("extdata/lungNEN_references", package = "PACMOS")
reference_LFs  <- system.file("extdata/lungNEN_references", "lungNEN_LFs.csv", package = "PACMOS")
query_csv <- system.file("extdata/test_data", "lungNEN_test_expr.csv", package = "PACMOS")


out_dir <- file.path('.', "lungNEN_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir


python_path <- Sys.getenv("PACMOS_PYTHON", unset = "")

if (!nzchar(python_path) || !file.exists(python_path)) {
  stop(
    "PACMOS_PYTHON is not set to a valid Python path.\n",
    "Before running this vignette, set PACMOS_PYTHON to the Python executable ",
    "from your pacmos_env Conda environment."
  )
}

## ----run-step1, results='markup', warning=FALSE-------------------------------

## -----------------------------------------------------------------------------
print('Starting Step 1')
PACMOS::s1_add_sample_to_mofa(
  query_matrix_path = query_csv,
  mofa_dir = mofa_dir,
  value_data_types = c('D_expr_MOFA'),
  outdir = out_dir,
  python_bin = python_path
)
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step2, warning=FALSE-------------------------------------------------
print('Starting Step 2')

PACMOS::s2_run_mofa(
  models_dir      = out_dir,
  matrices_subdir = "inputs",
  num_factors = 10,
  convergence_mode = "slow",
  maxiter = 10000,
  use_basilisk = FALSE,
  skip_existing = TRUE,
  python_bin      = python_path,
  views_map = c(
    RNA     = "D_expr_MOFA",
    Meth    = "D_met_MOFA",
    CNV     = "D_cnv_MOFA",
    Alt     = "D_alt_MOFA"
  ),
  binary_views = "Alt"
)
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step3,warning=FALSE--------------------------------------------------
print('Starting Step 3')

sample_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[dir.exists(file.path(sample_dirs, "inputs"))]
sample_ids <- basename(sample_dirs)

s3_results <- list()

for (id in sample_ids) {
  s3_results[[id]] <- PACMOS::s3_plot_query_samples_mofa(
    models_dir      = out_dir,
    matrices_subdir = "inputs",
    query_sample    = id,
    reference_LFs   = reference_LFs,
    reference_axes  = c(
      "Factor1",
      "Factor2",
      "Factor5"
    ),
    group = "group1",
    python_bin = python_path
  )
}
cat("Outputs are saved in:", out_dir, "\n")


## ----step3-qc-heatmap, fig.width=7, fig.height=5------------------------------
sample_to_show <- sample_ids[1]

s3_results[[sample_to_show]]$plots$alignment_heatmap

## ----step3-qc-scatter, fig.width=6, fig.height=5------------------------------
s3_results[[sample_to_show]]$plots$scatter_by_axis$Factor1

## ----step3-projection, fig.width=6, fig.height=5------------------------------
s3_results[[sample_to_show]]$plots$projection_pairs[[1]]

## ----run-step4, warning=FALSE-------------------------------------------------
print('Starting Step 4')

archetype_coords <- data.frame(
  Archetype = c("sc-enriched", "Ca A1", "Ca B", "Ca A2"),
  Factor1   = c( 0.115975, -2.808502,  3.587424, -0.009624),
  Factor2   = c(-0.663702, -1.362303, -0.980481,  2.544914),
  Factor5   = c(-5.090584,  1.194101,  0.919848,  0.810954),
  stringsAsFactors = FALSE
)

PACMOS::infer_fuzzy_weights(
  models_dir      = out_dir,
  coord           = archetype_coords,
  n_archetypes    = 4,
  input_type      = "stable"
)
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step5, warning=FALSE-------------------------------------------------
print('Starting Step 5')

fuzzy_results <- PACMOS::plot_fuzzy_query_sample(
  models_dir      = out_dir
)
cat("Outputs are saved in:", out_dir, "\n")


## ----s5-ternary, fig.width=20, fig.height=80,out.width="700px"----------------
sample_to_show <- sample_ids[1]

knitr::include_graphics(
  fuzzy_results$samples[[sample_to_show]]$files$projection_pdf
)

## ----session-info-------------------------------------------------------------
sessionInfo()

