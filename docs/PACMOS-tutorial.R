## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PACMOS)

## ----setup-paths--------------------------------------------------------------
mofa_dir <- system.file("extdata/MESOMICS_references", package = "PACMOS")
reference_LFs  <- system.file("extdata/MESOMICS_references", "MESOMICS_latent_factors.csv", package = "PACMOS")
query_csv <- system.file("extdata/test_data", "test-normalised-gene_count.csv", package = "PACMOS")


out_dir <- file.path('.', "query_output")
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
      value_data_types = c('D_exprB_MOFA'),
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
    RNA     = "D_exprB_MOFA",
    MethPro = "D_met.proB_MOFA",
    MethBod = "D_met.bodB_MOFA",
    MethEnh = "D_met.enhB_MOFA",
    Total   = "D_cnv_MOFA",
    Minor   = "D_loh_MOFA",
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
      "Morphology_LF",
      "Adaptive-response_LF"
    ),
    group = "group1",
    python_bin = python_path
  )
}
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step4, warning=FALSE-------------------------------------------------
print('Starting Step 4')

archetype_coords <- data.frame(
  Archetype = c(
    "Cell division",
    "Tumor-immune-interaction",
    "Acinar"
  ),
  Morphology_factor.MESOMICS = c(
    -3.59733245242571,
    -1.79807036681759,
    3.85856701073313
  ),
  Adaptive.response_factor.MESOMICS = c(
    -2.39600999565261,
    3.56348237698459,
    -0.825586204203109
  ),
  stringsAsFactors = FALSE
)

PACMOS::infer_fuzzy_weights(
  models_dir      = out_dir,
  matrices_subdir = "inputs",
  coord           = archetype_coords,
  n_archetypes    = 3
)
cat("Outputs are saved in:", out_dir, "\n")


## ----run-step5, warning=FALSE-------------------------------------------------
print('Starting Step 5')

PACMOS::plot_fuzzy_query_sample(
  models_dir      = out_dir,
  matrices_subdir = "inputs"
)
cat("Outputs are saved in:", out_dir, "\n")


