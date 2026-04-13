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
  `Adaptive-response_factor.MESOMICS` = c(
    -2.39600999565261,
    3.56348237698459,
    -0.825586204203109
  ),
  stringsAsFactors = FALSE
)



weights_df <- infer_fuzzy_weights(
  models_dir      = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/Mesomics_LOO",
  matrices_subdir = "train_test_all_omics",
  coord           = archetype_coords,
  n_archetypes    = 3
)
