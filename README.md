
# PACMOS

<img src="man/figures/PACMOS_logo.png" align="right" width="130" style="margin-left:50px;"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/lipikakalson/PACMOS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lipikakalson/PACMOS/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License:
MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)
<!-- badges: end -->

This R package provides a streamlined workflow to integrate query samples into reference MOFA (Multi-Omics Factor Analysis) models, perform inference (fuzzy/hard clustering), and visualize biological patterns.

The package is designed for reproducible multi-omics datasets, enabling:
<ul>
  <li>Integration of query samples into existing MOFA inputs.</li>
  <li>MOFA model retraining.</li>
  <li>Projection of query samples into reference latent spaces. </li>
  <li>Fuzzy/Hard clustering.</li>
  <li>Visualization of clustering.</li>
  </ul>
  
## Installation

You can install the development version of PACMOS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lipikakalson/PACMOS")
```

## Tutorial
A tutorial on the usage of PACMOS is available in the vignette, which is available at this link[].

## Functions

### Step 1
This function adds one or more query samples to reference MOFA input matrices in an incremental manner. Each query matrix is matched to a specific MOFA data layer, and the corresponding values are inserted into the reference matrices. For all other layers, NA values are added to maintain consistent structure across views.

The updated matrices (reference + one query sample) are saved as `.RData` files
```
PACMOS::s1_add_sample_to_mofa(
      query_matrix_path = c('query_expr.csv', 'query_met.csv'),
      mofa_dir = mofa_dir,
      value_data_types = c('D_expr_MOFA', 'D_met_MOFA'),
      outdir = out_dir,
      python_bin = python_path
    )
```
where,
<ol>
  <li>query_matrix_path =  Character vector of CSV file paths containing query sample matrices.</li>
  <li>mofa_dir          =  Directory containing reference MOFA `.RData` matrices.</li>
  <li>value_data_types  = Character vector of MOFA object names indicating which layer each query matrix corresponds to. Must be the same length as `query_matrix_path`.</li>
  <li>outdir            = Directory where updated `.RData` matrices will be written.</li>
  <li>python_bin        = Path to the Python binary used by the MOFA environment via the `reticulate` package.</li>
</ol>

### Step 2
This function trains a MOFA2 model using the matrices generated in `s1_add_sample_to_mofa()`.
```
PACMOS::s2_run_mofa(
  models_dir      = out_dir,
  matrices_subdir = "inputs",
  python_bin      = python_path,
  views_map = c(
    RNA     = "D_exprB_MOFA",
    Meth    = "D_met_MOFA",
    Total   = "D_cnv_MOFA",
    Minor   = "D_loh_MOFA",
    Alt     = "D_alt_MOFA"
  ),
  binary_views = "Alt"
)
```
where,
<ol>
  <li>models_dir      =  Root folder containing `.RData` matrices created by `s1_add_sample_to_mofa()`.</li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>python_bin      =  Path to the Python binary used by the MOFA environment via the `reticulate` package.</li>
  <li>views_map       =  Named character vector mapping view names to matrix object names.</li>
  <li>binary_views    =  Character vector specifying which views contain binary data</li>
</ol>

### Step 3
This function projects query samples into the reference latent factor space and generates quality control and projection visualizations.
```
PACMOS::s3_plot_query_samples_mofa(
    models_dir      = out_dir,
    matrices_subdir = "inputs",
    query_sample    = id,
    reference_LFs   = reference_LFs,
    reference_axes  = c(
      "Morphology_LF",
      "Adaptive-response_LF"
    ),
    python_bin = python_path
  )
```
where,
<ol>
  <li>models_dir      =  Root folder containing output of `s2_run_MOFA()`. (out_dir) </li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>query_sample    =  Sample ID of the query sample.</li>
  <li>reference_LFs   =  data.frame or path to CSV containing reference latent factors.</li>
  <li>reference_axes  =  Character vector of reference axis column names we need to match and align.</li>
  <li>binary_views    =  Character vector specifying which views contain binary data</li>
  <li>python_bin      =  Path to the Python binary used by the MOFA environment via the `reticulate` package.</li>
</ol>

## FUZZY CLUSTERING
### Step 4
This function estimates the degree to which each query sample belongs to each biological archetype predefined in reference multiomics data using a fuzzy weighting approach.
```
archetype_coords <- data.frame(
  Archetype = c("a", "b", "c"), # archetypes
  LF1 = c(
    -3.59733245242571,
    -1.79807036681759,
    3.85856701073313
  ),
  LF2 = c(
    -2.39600999565261,
    3.56348237698459,
    -0.825586204203109
  ), # archetype coordinates in LF space
  stringsAsFactors = FALSE
)

PACMOS::infer_fuzzy_weights(
  models_dir      = out_dir,
  matrices_subdir = "inputs",
  coord           = archetype_coords,
  n_archetypes    = 3
```
where,
<ol>
  <li>models_dir      =  Root folder. (out_dir) </li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>coord           =  Data frame of archetype coordinates</li>
  <li>n_archetypes    =  Expected number of archetypes.</li>
</ol>

### Step 5
This function visualizes the archetype composition of query samples in the reference archetypal space based on the inferred fuzzy weights.
```
PACMOS::plot_fuzzy_query_sample(
  models_dir        = out_dir,
  matrices_subdir   = "inputs"
)
```
where,
<ol>
  <li>models_dir      =  Root folder. (out_dir) </li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
</ol>

## HARD CLUSTERING
### Step 4
This function assigns each query sample to a discrete cluster using k-means clustering in the latent factor space.
```
PACMOS::infer_kmeans_labels(
  models_dir      = out_dir",
  matrices_subdir = "inputs",
  k               = 4, # number of clusters
  lf_cols         = c("LF1", "LF2", "LF3")
)

```
where,
<ol>
  <li>models_dir      =  Root folder. (out_dir) </li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>k               =  Number of k-means clusters.</li>
  <li>lf_cols         =  LF columns to use as features for k-means.</li>
</ol>

### Step 5
This function visualizes the hard clustering results obtained from infer_kmeans_clusters().
```
PACMOS::plot_kmeans_query_sample(
    models_dir      = loo_root,
    matrices_subdir = "train_test_Meth_only",
    lf_cols         = c("Factor1", "Factor2", "Factor5"),
    ref_labels_path = ref_file,
    ref_cols        = c("sample", "bio_label"),
    query_sample    = sid
  )
```
where,
<ol>
  <li>models_dir      =  Root folder. (out_dir) </li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>lf_cols         =  LF columns to use as features for k-means. (used for plotting) </li>
  <li>ref_labels_path =  Path to CSV containing reference labels.
  <li>ref_cols        =  Character vector of length 2. Column names in the `ref_labels_path` CSV to use. E.g. c("sample", "bio_label")</li>
  <li>query_sample    =  query_sample </li>
</ol>

## Output directory structure
```
Output directory structure:
  <outdir>/
    <query_sample>/
      <inputs>/
        value_data_type_query_sample.RData (Step 1 output)
        MOFA-query_sample.hdf5 (Step 2 output)
        <plots>/
          #--Step 3 outputs--
          query_sample_projection.pdf
          query_sample__quality_check_metrics.csv
          query_sample__quality_check_metrics.pdf
          query_sample_query_sample_LFs.csv
          query_sample__retrained_LFs_all_samples.csv
          query_sample__stable_input.csv

          #--Step 4 output--
          query_sample_archetype_weights_all_samples.csv
          query_sample__archetype_weights.csv

          #--Step5 output--
          query_sample__archetype_projection.pdf

```
