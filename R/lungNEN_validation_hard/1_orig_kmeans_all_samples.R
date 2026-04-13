library(tidyverse)
library(ggplot2)

# ── 1. LOAD DATA ─────────────────────────────────────────────
df <- read.csv("inst/extdata/lungNEN_references/lungNEN_LFs.csv", row.names = 1)
df <- df[, -1]

label <- read.delim("inst/extdata/lungNEN_references/lungNEN_LF_K4_with_label.txt", row.names = 1)

# ── 2. SELECT LFs & ALIGN ────────────────────────────────────
selected_lfs   <- c("Factor1", "Factor2", "Factor5")
lf_data        <- df[, selected_lfs, drop = FALSE]
common_samples <- intersect(rownames(lf_data), rownames(label))
lf_data        <- lf_data[common_samples, , drop = FALSE]
bio_labels     <- label[common_samples, "archetype_k4_LF3_label"]

# ── 3. K-MEANS ───────────────────────────────────────────────
set.seed(42)
k <- 4

km_result <- kmeans(lf_data, centers = k, nstart = 25, iter.max = 100)

result_df <- as.data.frame(lf_data) %>%
  mutate(
    sample         = rownames(.),
    kmeans_cluster = as.factor(km_result$cluster),
    bio_label      = factor(bio_labels)
  )

# ── 4. SCATTER: shape = k-means cluster, color = bio label ───
p_scatter <- ggplot(result_df, aes(x = Factor1, y = Factor2,
                                   color = bio_label,
                                   shape = kmeans_cluster)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = paste0("K-means (k=", k, ") on MOFA Latent Factors"),
    x     = "Factor1",
    y     = "Factor2",
    color = "Biological Label",
    shape = "K-means Cluster"
  ) +
  theme_bw(base_size = 13)

print(p_scatter)
#ggsave("kmeans_lf_scatter.pdf", p_scatter, width = 8, height = 5.5)

# ── 5. BUILD OVERLAP TABLE ────────────────────────────────────
overlap <- result_df %>%
  count(kmeans_cluster, bio_label) %>%
  complete(kmeans_cluster, bio_label, fill = list(n = 0)) %>%
  group_by(kmeans_cluster) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()

# ── 6. ORDER FOR DIAGONAL: map each cluster to its dominant label
cluster_order <- overlap %>%
  group_by(kmeans_cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(bio_label) %>%           # sort by bio label alphabetically
  pull(kmeans_cluster)

label_order <- overlap %>%
  group_by(kmeans_cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(bio_label) %>%
  pull(bio_label) %>%
  as.character() %>%
  unique()

overlap <- overlap %>%
  mutate(
    kmeans_cluster = factor(kmeans_cluster, levels = cluster_order),
    bio_label      = factor(bio_label,      levels = label_order)
  )

# ── 7. HEATMAP ────────────────────────────────────────────────
p_heatmap <- ggplot(overlap, aes(x = bio_label, y = kmeans_cluster, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    data = overlap %>% filter(n > 0),
    aes(label = paste0(n, "\n(", pct, "%)")),
    size = 3.5, fontface = "bold"
  ) +
  scale_fill_gradient(low = "white", high = "#2166AC") +
  labs(
    title = "K-means Clusters vs. Biological Labels",
    x     = "Biological Label",
    y     = "K-means Cluster",
    fill  = "% of cluster"
  ) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p_heatmap)
#ggsave("cluster_label_heatmap.pdf", p_heatmap, width = 6, height = 4.5)

# ── 8. CLUSTER MAPPING SUMMARY ───────────────────────────────
mapping <- overlap %>%
  group_by(kmeans_cluster) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  select(kmeans_cluster, dominant_label = bio_label, n_dominant = n, pct)

cat("\n── K-means Cluster → Biological Label Mapping ──\n")
print(mapping)

#### save the outputs ###########

# ── 9. CREATE OUTPUT FOLDER & SAVE ───────────────────────────
out_dir <- "Analysis_010426/lungNEN_hard_ref"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── SAVE PLOTS ───────────────────────────────────────────────
ggsave(file.path(out_dir, "kmeans_lf_scatter.pdf"),      p_scatter, width = 8,   height = 5.5)
ggsave(file.path(out_dir, "cluster_label_heatmap.pdf"),  p_heatmap, width = 6,   height = 4.5)

# ── CSV 1: sample + biological label + k-means cluster only ──
sample_labels <- result_df %>%
  select(sample, bio_label, kmeans_cluster)

write.csv(sample_labels,
          file.path(out_dir, "sample_cluster_labels.csv"),
          row.names = FALSE, quote = FALSE)

# ── CSV 2: full details — factors, labels, cluster, centroids distance ──
centroids <- as.data.frame(km_result$centers)
colnames(centroids) <- selected_lfs

# Distance of each sample to its assigned centroid
dist_to_centroid <- mapply(function(i, cl) {
  sqrt(sum((lf_data[i, ] - centroids[cl, ])^2))
}, i = seq_len(nrow(lf_data)), cl = km_result$cluster)

full_details <- result_df %>%
  mutate(dist_to_centroid = dist_to_centroid) %>%
  select(sample, bio_label, kmeans_cluster, all_of(selected_lfs), dist_to_centroid)

write.csv(full_details,
          file.path(out_dir, "sample_full_details.csv"),
          row.names = FALSE, quote = FALSE)

# ── CSV 3: cluster centroids ──────────────────────────────────
centroid_df <- centroids %>%
  mutate(kmeans_cluster = as.factor(rownames(.))) %>%
  left_join(
    mapping %>% select(kmeans_cluster, dominant_label),
    by = "kmeans_cluster"
  ) %>%
  relocate(kmeans_cluster, dominant_label)

write.csv(centroid_df,
          file.path(out_dir, "cluster_centroids.csv"),
          row.names = FALSE, quote = FALSE)

cat("\n── Files saved to:", out_dir, "──\n")
cat("  - kmeans_lf_scatter.pdf\n")
cat("  - cluster_label_heatmap.pdf\n")
cat("  - sample_cluster_labels.csv\n")
cat("  - sample_full_details.csv\n")
cat("  - cluster_centroids.csv\n")
