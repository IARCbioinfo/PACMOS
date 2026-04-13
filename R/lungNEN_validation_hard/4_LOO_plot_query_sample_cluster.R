devtools::load_all("/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS")

loo_root   <- "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_hard_LOO"
all_sids   <- basename(list.dirs(loo_root, recursive = FALSE))
sample_ids <- all_sids[
  file.exists(file.path(loo_root, all_sids,
                        "train_test_Meth_only", "plots",
                        paste0(all_sids, "_sample_clusters.csv")))
]

for (sid in sample_ids) {

  ref_file <- file.path(loo_root, sid, "inputs", "plots",
                        paste0(sid, "_kmeans_results.csv"))

  if (!file.exists(ref_file)) {
    message("Skipping ", sid, " — ref file not found")
    next
  }

  plot_kmeans_query_sample(
    models_dir      = loo_root,
    matrices_subdir = "train_test_Meth_only",
    lf_cols         = c("Factor1", "Factor2", "Factor5"),
    ref_labels_path = ref_file,
    ref_cols        = c("sample", "bio_label"),
    query_sample    = sid
  )
}
