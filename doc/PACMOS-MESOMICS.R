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
#      "/home/user/miniconda3/envs/pacmos_env/bin/python"
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
mofa_dir <- system.file("extdata/MESOMICS_references", package = "PACMOS")
reference_LFs  <- system.file("extdata/MESOMICS_references", "MESOMICS_latent_factors.csv", package = "PACMOS")
query_csv <- system.file("extdata/test_data", "MESOMICS_test_expr.csv", package = "PACMOS")


out_dir <- file.path('.', "MESOMICS_output")
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
    reference_axes  = c("Morphology_LF", "Adaptive-response_LF"),
    group           = "group1",
    python_bin      = python_path
  )
}
cat("Outputs are saved in:", out_dir, "\n")


## ----step3-qc-heatmap, fig.width=7, fig.height=5------------------------------
sample_to_show <- sample_ids[1]

s3_results[[sample_to_show]]$plots$alignment_heatmap

## ----step3-qc-scatter, fig.width=6, fig.height=5------------------------------
s3_results[[sample_to_show]]$plots$scatter_by_axis$Morphology_LF

## ----step3-projection, fig.width=6, fig.height=5------------------------------
s3_results[[sample_to_show]]$plots$projection_pairs[[1]]

## ----run-step4, warning=FALSE-------------------------------------------------
print('Starting Step 4')

archetype_coords <- data.frame(
  Archetype = c("Cell division", "Tumor-immune-interaction", "Acinar"),
  Morphology_LF = c(-3.59733245242571, -1.79807036681759, 3.85856701073313),
  `Adaptive-response_LF` = c(-2.39600999565261, 3.56348237698459, -0.825586204203109),
  stringsAsFactors = FALSE,
  check.names = FALSE)

PACMOS::infer_fuzzy_weights(
  models_dir      = out_dir,
  coord           = archetype_coords,
  n_archetypes    = 3,
  input_type      = "stable",
  reference_axes = c("Morphology_LF", "Adaptive-response_LF")
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

