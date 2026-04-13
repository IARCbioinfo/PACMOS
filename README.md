
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

The package is designed for reproducible, modular analysis of multi-omics datasets, enabling:
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
  <li>project_dir     =  Root folder containing `.RData` matrices created by `s1_add_sample_to_mofa()`.</li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>python_bin      = Path to the Python binary used by the MOFA environment via the `reticulate` package.</li>
  <li>views_map       = Named character vector mapping view names to matrix object names.</li>
  <li>binary_views    = Character vector specifying which views contain binary data</li>
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
  <li>models_dir =  Root folder containing `.RData` matrices created by `s1_add_sample_to_mofa()`.</li>
  <li>matrices_subdir =  Folder name where `.RData` files are stored.</li>
  <li>query_sample    =  Sample ID of the query sample.</li>
  <li>reference_LFs   =  data.frame or path to CSV containing reference latent factors.</li>
  <li>reference_axes  =  Character vector of reference axis column names we need to match and align.</li>
  <li>binary_views    =  Character vector specifying which views contain binary data</li>
  <li>python_bin      =  Path to the Python binary used by the MOFA environment via the `reticulate` package.</li>

</ol>
