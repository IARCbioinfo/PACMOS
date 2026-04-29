## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PACMOS)

## ----setup-paths--------------------------------------------------------------
mofa_dir <- system.file("extdata/lungNEN_references", package = "PACMOS")
reference_LFs  <- system.file("extdata/lungNEN_references", "lungNEN_LFs.csv", package = "PACMOS")
query_csv <- system.file("extdata/test_data", "lungNEN_test_expr.csv", package = "PACMOS")


out_dir <- file.path('.', "lungNEN_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir

Sys.setenv(MOFA_PYTHON = "/home/lipikal/miniconda3/envs/mofa_env/bin/python") #replace this with your 'which python' output
python_path <- Sys.getenv("MOFA_PYTHON", unset = NA)

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

sample_ids <- list.files(out_dir)

for (id in sample_ids) {
  PACMOS::s3_plot_query_samples_mofa(
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
  matrices_subdir = "inputs",
  coord           = archetype_coords,
  n_archetypes    = 4,
  input_type      = "stable"
)
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step5, warning=FALSE-------------------------------------------------
print('Starting Step 5')

PACMOS::plot_fuzzy_query_sample(
  models_dir      = out_dir,
  matrices_subdir = "inputs"
)
cat("Outputs are saved in:", out_dir, "\n")


