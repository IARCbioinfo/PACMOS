devtools::load_all('/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/')

loo_root   <- "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/Mesomics_LOO"
all_sids    <- basename(list.dirs(loo_root, recursive = FALSE))
sample_ids <- all_sids[
  file.exists(file.path(loo_root, all_sids, "inputs", "plots",
                        paste0(all_sids, "_aligned_LFs_all_samples.csv")))
]

for (sid in sample_ids) {
  s3_plot_test_samples_mofa(
    models_dir      = loo_root,
    matrices_subdir = "train_test_Meth_only",
    query_sample    = sid,
    reference_LFs   = file.path(loo_root, sid, "inputs", "plots",
                                paste0(sid, "_aligned_LFs_all_samples.csv")),
    reference_axes  = c("Morphology_LF",
                        "Adaptive-response_LF"),
    python_bin      = "/home/lipikal/miniconda3/envs/mofa_env/bin/python",
  )
}
