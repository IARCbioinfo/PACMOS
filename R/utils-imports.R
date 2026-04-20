#' @importFrom grDevices dev.off pdf
#' @importFrom graphics axis layout legend mtext par rect segments title
#' @importFrom stats ave cor kmeans setNames
#' @importFrom utils combn read.csv write.csv
#' @importFrom rlang .data
NULL
utils::globalVariables(c(
  "bio_label", "kmeans_cluster", "pct", "n",
  "Ref_Axis", "MOFA_Factor", "Correlation",
  "Reference", "Aligned"
))
