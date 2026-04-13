archetype_coords <- data.frame(
  Archetype = c("SC-enriched", "Ca A1", "Ca B", "Ca A2"),
  Factor1   = c( 0.115975, -2.808502,  3.587424, -0.009624),
  Factor2   = c(-0.663702, -1.362303, -0.980481,  2.544914),
  Factor5   = c(-5.090584,  1.194101,  0.919848,  0.810954),
  stringsAsFactors = FALSE
)

weights_df <- infer_fuzzy_weights(
  models_dir      = "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_LOO",
  matrices_subdir = "train_test_Meth_only",
  coord           = archetype_coords,
  n_archetypes    = 4
)
