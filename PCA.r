# ============================================================
# PCA and unsupervised clustering of 29 transcription factors
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

input_file <- "/path/to/input_directory"
output_dir <- "/path/to/output_directory"

df <- read.csv(
  input_file,
  row.names = 1,
  check.names = FALSE
)

# rows = genes, columns = features
# transpose to features x genes for Seurat
mat <- t(as.matrix(df))

cat(
  "Matrix dimensions (features x genes):",
  nrow(mat), "x", ncol(mat), "\n"
)

# ------------------------------------------------------------
# 2. Create Seurat object
# ------------------------------------------------------------

seu <- CreateSeuratObject(
  counts = mat,
  min.cells = 0,
  min.features = 0,
  project = "TF_clustering"
)

# Use the input kinetic parameters directly as the data layer
if (packageVersion("Seurat") >= "5.0.0") {

  LayerData(
    seu,
    assay = "RNA",
    layer = "data"
  ) <- mat

} else {

  seu <- SetAssayData(
    seu,
    slot = "data",
    new.data = mat
  )
}

# ------------------------------------------------------------
# 3. Feature scaling
# ------------------------------------------------------------

seu <- ScaleData(
  seu,
  features = rownames(seu),
  verbose = FALSE
)

# ------------------------------------------------------------
# 4. PCA
# ------------------------------------------------------------

n_pcs_compute <- min(
  nrow(mat) - 1,
  ncol(mat) - 1,
  4
)

seu <- RunPCA(
  seu,
  features = rownames(seu),
  npcs = n_pcs_compute,
  seed.use = 42,
  verbose = FALSE
)

actual_pcs <- ncol(
  Embeddings(seu, reduction = "pca")
)

use_pcs <- min(actual_pcs, 4)

cat(
  "PCs computed:", actual_pcs,
  "| PCs used for clustering:", use_pcs,
  "\n"
)

# ------------------------------------------------------------
# 5. SNN graph
# ------------------------------------------------------------

seu <- FindNeighbors(
  seu,
  dims = 1:use_pcs,
  k.param = min(5, ncol(seu) - 1),
  verbose = FALSE
)

# ------------------------------------------------------------
# 6. Louvain clustering
# ------------------------------------------------------------

seu <- FindClusters(
  seu,
  resolution = 0.8,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  verbose = FALSE
)

cat("\nCluster sizes:\n")
print(table(seu$seurat_clusters))

# ------------------------------------------------------------
# 7. Extract PCA coordinates
# ------------------------------------------------------------

pca_coordinates <- as.data.frame(
  Embeddings(seu, reduction = "pca")[, 1:2]
)

colnames(pca_coordinates) <- c("PC1", "PC2")

pca_coordinates$Gene <- rownames(pca_coordinates)

pca_coordinates$Cluster <- factor(
  seu$seurat_clusters
)

# Source data
write.csv(
  pca_coordinates,
  file.path(output_dir, "PCA_coordinates.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------
# 8. Cluster membership
# ------------------------------------------------------------

cluster_df <- data.frame(
  Gene = colnames(seu),
  Cluster = seu$seurat_clusters
)

cluster_df <- cluster_df[
  order(cluster_df$Cluster),
]

write.csv(
  cluster_df,
  file.path(output_dir, "cluster_membership.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\nCluster membership:\n")
print(cluster_df, row.names = FALSE)

# ------------------------------------------------------------
# 9. PCA visualization
# ------------------------------------------------------------

n_clusters <- nlevels(
  factor(seu$seurat_clusters)
)

cluster_colors <- c(
  "#E41A1C",
  "#377EB8",
  "#4DAF4A",
  "#984EA3",
  "#FF7F00",
  "#A65628"
)[1:n_clusters]

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

# ------------------------------------------------------------
# 10. Save PCA figure
# ------------------------------------------------------------

ggsave(
  file.path(output_dir, "PCA_29TFs.pdf"),
  plot = p_pca,
  width = 7,
  height = 6,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "PCA_29TFs.tiff"),
  plot = p_pca,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(output_dir, "PCA_29TFs.png"),
  plot = p_pca,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

cat("\nPCA analysis completed.\n")
cat("Outputs:\n")
cat("  PCA_29TFs.pdf\n")
cat("  PCA_29TFs.tiff\n")
cat("  PCA_29TFs.png\n")
cat("  PCA_coordinates.csv\n")
cat("  cluster_membership.csv\n")
