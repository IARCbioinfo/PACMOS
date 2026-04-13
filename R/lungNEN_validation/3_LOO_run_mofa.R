devtools::load_all("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS")

run_mofa_MESOMICS_validation <- function(
    project_dir,
    train_test_subdir = "train_test_all_omics",
    python_bin        = "/home/lipikal/miniconda3/envs/mofa_env/bin/python"
) {

  views_map <- c(
    RNA     = "D_expr_MOFA",
    Meth    = "D_met_MOFA",
    CNV     = "D_cnv_MOFA",
    Alt     = "D_alt_MOFA"
  )

  s2_run_mofa(
    project_dir     = project_dir,
    matrices_subdir = train_test_subdir,
    python_bin      = python_bin,
    views_map       = views_map,
    binary_views    = "Alt"
  )
}

run_mofa_MESOMICS_validation(
  project_dir       = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO",
  train_test_subdir = "train_test_Genomic_Meth_only",
  python_bin        = "/home/lipikal/miniconda3/envs/mofa_env/bin/python"
)
