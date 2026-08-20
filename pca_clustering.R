#!/usr/bin/env Rscript

# ============================================================
# PCA and unsupervised clustering of transcription factors
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ============================================================
# Command-line arguments
# ============================================================
#
# Usage:
# Rscript pca_clustering.R \
#   SelectedGeneConstant_updated.csv \
#   output_directory
#
# Input format:
#   Rows    = transcription factors
#   Columns = kinetic parameters or other numerical features
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript pca_clustering.R ",
    "<input.csv> <output_directory>"
  )
}

input_file <- args[1]
output_dir <- args[2]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Parameters
# ============================================================

max_pcs <- 4
k_neighbors <- 5
cluster_resolution <- 0.8
random_seed <- 42

# ============================================================
# Load input data
# ============================================================

df <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(df) < 2) {
  stop("At least two transcription factors are required.")
}

if (ncol(df) < 2) {
  stop("At least two numerical features are required.")
}

if (!all(vapply(df, is.numeric, logical(1)))) {
  stop(
    "All feature columns in the input file must contain numerical values."
  )
}

if (anyNA(df)) {
  stop("Input data contain missing values.")
}

# Rows = transcription factors
# Columns = features
#
# Seurat expects features x samples, so the matrix is transposed.
mat <- t(as.matrix(df))

if (any(!is.finite(mat))) {
  stop("Input data contain non-finite values.")
}

cat(
  "Matrix dimensions (features x transcription factors):",
  nrow(mat), "x", ncol(mat), "\n"
)

# ============================================================
# Create Seurat object
# ============================================================

# Kinetic parameters are continuous measurements rather than
# sequencing counts. A zero-valued count matrix is therefore used
# only to initialize the Seurat object, while the original kinetic
# parameter matrix is stored in the data layer.

dummy_counts <- matrix(
  0,
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
)

seu <- CreateSeuratObject(
  counts = dummy_counts,
  min.cells = 0,
  min.features = 0,
  project = "TF_clustering"
)

if (packageVersion("Seurat") >= "5.0.0") {

  LayerData(
    seu,
    assay = "RNA",
    layer = "data"
  ) <- mat

} else {

  seu <- SetAssayData(
    seu,
    assay = "RNA",
    slot = "data",
    new.data = mat
  )
}

# ============================================================
# Feature scaling
# ============================================================

seu <- ScaleData(
  seu,
  features = rownames(seu),
  verbose = FALSE
)

# ============================================================
# PCA
# ============================================================

n_pcs_compute <- min(
  max_pcs,
  nrow(mat),
  ncol(mat) - 1
)

if (n_pcs_compute < 2) {
  stop("Insufficient matrix dimensions for PCA.")
}

seu <- RunPCA(
  seu,
  features = rownames(seu),
  npcs = n_pcs_compute,
  seed.use = random_seed,
  verbose = FALSE
)

actual_pcs <- ncol(
  Embeddings(
    seu,
    reduction = "pca"
  )
)

use_pcs <- min(
  actual_pcs,
  max_pcs
)

cat(
  "PCs computed:", actual_pcs,
  "| PCs used for clustering:", use_pcs,
  "\n"
)

# ============================================================
# SNN graph construction
# ============================================================

k_use <- min(
  k_neighbors,
  ncol(seu) - 1
)

seu <- FindNeighbors(
  seu,
  dims = seq_len(use_pcs),
  k.param = k_use,
  verbose = FALSE
)

# ============================================================
# Louvain clustering
# ============================================================

seu <- FindClusters(
  seu,
  resolution = cluster_resolution,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = random_seed,
  verbose = FALSE
)

cat("\nCluster sizes:\n")
print(table(seu$seurat_clusters))

# ============================================================
# PCA coordinates
# ============================================================

pca_embeddings <- Embeddings(
  seu,
  reduction = "pca"
)

pca_coordinates <- as.data.frame(
  pca_embeddings[, 1:2, drop = FALSE]
)

colnames(pca_coordinates) <- c(
  "PC1",
  "PC2"
)

pca_coordinates$Gene <- rownames(
  pca_coordinates
)

pca_coordinates$Cluster <- factor(
  seu$seurat_clusters
)

write.csv(
  pca_coordinates,
  file = file.path(
    output_dir,
    "PCA_coordinates.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Cluster membership
# ============================================================

cluster_df <- data.frame(
  Gene = colnames(seu),
  Cluster = factor(seu$seurat_clusters),
  stringsAsFactors = FALSE
)

cluster_df <- cluster_df[
  order(
    as.integer(cluster_df$Cluster),
    cluster_df$Gene
  ),
  ,
  drop = FALSE
]

write.csv(
  cluster_df,
  file = file.path(
    output_dir,
    "cluster_membership.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

cat("\nCluster membership:\n")
print(
  cluster_df,
  row.names = FALSE
)

# ============================================================
# PCA visualization
# ============================================================

cluster_levels <- levels(
  pca_coordinates$Cluster
)

color_palette <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628",
  "#F781BF",
  "#999999"
)

if (length(cluster_levels) > length(color_palette)) {
  stop(
    "The number of clusters exceeds the available color palette."
  )
}

cluster_colors <- setNames(
  color_palette[
    seq_along(cluster_levels)
  ],
  cluster_levels
)

p_pca <- ggplot(
  pca_coordinates,
  aes(
    x = PC1,
    y = PC2,
    color = Cluster
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = Gene),
    size = 2.8,
    hjust = -0.15,
    vjust = 0.4,
    show.legend = FALSE,
    check_overlap = FALSE
  ) +
  scale_color_manual(
    values = cluster_colors,
    name = "Cluster"
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.08, 0.22)
    )
  ) +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(
      color = "black"
    ),
    axis.text = element_text(
      color = "black"
    )
  )

print(p_pca)

# ============================================================
# Save PCA figures
# ============================================================

ggsave(
  filename = file.path(
    output_dir,
    "PCA_TF_clustering.pdf"
  ),
  plot = p_pca,
  width = 7,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "PCA_TF_clustering.tiff"
  ),
  plot = p_pca,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    output_dir,
    "PCA_TF_clustering.png"
  ),
  plot = p_pca,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nPCA and clustering analysis completed.\n")
cat("Results saved to:", output_dir, "\n")
