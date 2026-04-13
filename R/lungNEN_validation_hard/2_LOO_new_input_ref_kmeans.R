library(tidyverse)
library(ggplot2)

# ── 1. LOAD REFERENCE LABELS ─────────────────────────────────
ref_labels <- read.csv(
  "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_hard_ref/sample_cluster_labels.csv"
) %>% select(sample, bio_label)

# ── 2. PATHS ─────────────────────────────────────────────────
loo_dir <- "/data/mesomics/work/mesomics2/lipikal/lipikal-2026/Rpackage/PACMOS/PACMOS/Analysis_010426/lungNEN_hard_LOO"
loo_folders <- list.dirs(loo_dir, recursive = FALSE)

# ── 3. LOOP OVER EACH LOO SAMPLE ─────────────────────────────
for (folder in loo_folders) {

  left_out <- basename(folder)
  lf_file  <- file.path(folder, "inputs", "plots",
                        paste0(left_out, "_aligned_LFs_all_samples.csv"))

  if (!file.exists(lf_file)) {
    message("Skipping ", left_out, " — file not found")
    next
  }

  # ── Read LFs ───────────────────────────────────────────────
  lf_data <- read.csv(lf_file, row.names = 1)
  # columns: Factor1, Factor2, Factor5

  # ── Keep only samples present in reference labels ──────────
  common_samples <- intersect(rownames(lf_data), ref_labels$sample)

  if (length(common_samples) < 4) {
    message("Skipping ", left_out, " — too few common samples")
    next
  }

  lf_data   <- lf_data[common_samples, , drop = FALSE]
  bio_ref   <- ref_labels %>% filter(sample %in% common_samples)

  # ── K-means ────────────────────────────────────────────────
  set.seed(42)
  km <- kmeans(lf_data, centers = 4, nstart = 25, iter.max = 100)

  result_df <- as.data.frame(lf_data) %>%
    mutate(
      sample         = rownames(.),
      left_out_sample = left_out,
      kmeans_cluster  = as.factor(km$cluster)
    ) %>%
    left_join(ref_labels, by = "sample")   # attach bio_label from reference

  # ── Majority: map k-means cluster → bio label ─────────
  cluster_map <- result_df %>%
    count(kmeans_cluster, bio_label) %>%
    group_by(kmeans_cluster) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(kmeans_cluster, mapped_bio_label = bio_label)

  result_df <- result_df %>%
    left_join(cluster_map, by = "kmeans_cluster") %>%
    select(left_out_sample, sample, Factor1, Factor2, Factor5,
           kmeans_cluster, bio_label, mapped_bio_label)

  # ── Save CSV ───────────────────────────────────────────────
  out_plots <- file.path(folder, "inputs", "plots")

  write.csv(result_df,
            file.path(out_plots, paste0(left_out, "_kmeans_results.csv")),
            row.names = FALSE, quote = FALSE)

  # ── Build overlap for heatmap ──────────────────────────────
  overlap <- result_df %>%
    count(kmeans_cluster, bio_label) %>%
    complete(kmeans_cluster, bio_label, fill = list(n = 0)) %>%
    group_by(kmeans_cluster) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()

  # Order for diagonal
  cluster_order <- overlap %>%
    group_by(kmeans_cluster) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    arrange(bio_label) %>%
    pull(kmeans_cluster)

  label_order <- overlap %>%
    group_by(kmeans_cluster) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    arrange(bio_label) %>%
    pull(bio_label) %>%
    as.character() %>%
    unique()

  overlap <- overlap %>%
    mutate(
      kmeans_cluster = factor(kmeans_cluster, levels = cluster_order),
      bio_label      = factor(bio_label,      levels = label_order)
    )

  # ── Heatmap ────────────────────────────────────────────────
  p_heatmap <- ggplot(overlap, aes(x = bio_label, y = kmeans_cluster, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      data = overlap %>% filter(n > 0),
      aes(label = paste0(n, "\n(", pct, "%)")),
      size = 3, fontface = "bold"
    ) +
    scale_fill_gradient(low = "white", high = "#2166AC") +
    labs(
      title = paste0("LOO: ", left_out, " — K-means vs. Biological Labels"),
      x     = "Biological Label",
      y     = "K-means Cluster",
      fill  = "% of cluster"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(
    file.path(out_plots, paste0(left_out, "_kmeans_heatmap.pdf")),
    p_heatmap, width = 6, height = 4.5
  )

  message("Done: ", left_out)
}

# ── 4. SUMMARY: confusion matrix accuracy per LOO run ────────
all_results <- map_dfr(loo_folders, function(folder) {
  left_out <- basename(folder)
  res_file <- file.path(folder, "inputs", "plots",
                        paste0(left_out, "_kmeans_results.csv"))
  if (file.exists(res_file)) read.csv(res_file) else NULL
})

out_summary_dir <- file.path(loo_dir, "..", "lungNEN_hard_ref")

# Accuracy per LOO run = diagonal sum / total
# i.e. samples where mapped_bio_label == bio_label / all samples in that run
per_run_accuracy <- all_results %>%
  filter(!is.na(bio_label), !is.na(mapped_bio_label)) %>%
  group_by(left_out_sample) %>%
  summarise(
    n_total      = n(),
    n_correct    = sum(mapped_bio_label == bio_label),   # diagonal sum
    accuracy_pct = round(100 * n_correct / n_total, 1),  # overall accuracy
    .groups = "drop"
  ) %>%
  arrange(desc(accuracy_pct))

write.csv(per_run_accuracy,
          file.path(out_summary_dir, "LOO_new_ref_inputs_per_run_accuracy.csv"),
          row.names = FALSE, quote = FALSE)

cat("\n── LOO Confusion Matrix Accuracy Summary ──\n")
cat("Mean accuracy across all runs: ", round(mean(per_run_accuracy$accuracy_pct), 1), "%\n")
cat("Min accuracy:                  ", min(per_run_accuracy$accuracy_pct), "%\n")
cat("Max accuracy:                  ", max(per_run_accuracy$accuracy_pct), "%\n")
cat("\nRuns with accuracy < 90%:\n")
print(per_run_accuracy %>% filter(accuracy_pct < 90))

message("\n── All LOO samples processed ──")

message("\n─ All LOO samples processed ──")
