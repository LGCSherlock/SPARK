#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

# =========================
# Command-line arguments
# =========================
#
# Usage:
# Rscript go_enrichment_by_cluster.R \
#   cluster_membership.csv \
#   output_directory
#
# Input CSV must contain two columns:
#   Gene    Gene symbol
#   Cluster Cluster assignment
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript go_enrichment_by_cluster.R ",
    "<cluster_membership.csv> <output_directory>"
  )
}

input_file <- args[1]
output_dir <- args[2]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Parameters
# =========================

top_n <- 3
bar_color <- "#bdcbe6"

# =========================
# Load input data
# =========================

gene_data <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

required_columns <- c("Gene", "Cluster")

if (!all(required_columns %in% colnames(gene_data))) {
  stop(
    "Input file must contain the following columns: ",
    paste(required_columns, collapse = ", ")
  )
}

gene_data <- gene_data %>%
  select(Gene, Cluster) %>%
  filter(
    !is.na(Gene),
    !is.na(Cluster),
    Gene != ""
  ) %>%
  distinct()

# =========================
# Convert gene symbols to Entrez IDs
# =========================

id_mapping <- bitr(
  unique(gene_data$Gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_cluster_df <- gene_data %>%
  inner_join(
    id_mapping,
    by = c("Gene" = "SYMBOL")
  ) %>%
  distinct(Gene, Cluster, ENTREZID)

# Save gene ID mapping for reproducibility
write.csv(
  gene_cluster_df,
  file = file.path(output_dir, "gene_id_mapping.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Report genes that could not be mapped
unmapped_genes <- gene_data %>%
  anti_join(
    id_mapping,
    by = c("Gene" = "SYMBOL")
  )

if (nrow(unmapped_genes) > 0) {
  write.csv(
    unmapped_genes,
    file = file.path(output_dir, "unmapped_genes.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

# =========================
# GO Biological Process enrichment
# =========================

cluster_levels <- unique(gene_data$Cluster)

result_list <- list()

for (cluster_id in cluster_levels) {

  cluster_genes <- gene_cluster_df %>%
    filter(Cluster == cluster_id) %>%
    pull(ENTREZID) %>%
    unique()

  if (length(cluster_genes) == 0) {
    message(
      "No mapped genes found for Cluster ",
      cluster_id,
      ". Skipping."
    )
    next
  }

  ego <- enrichGO(
    gene = cluster_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  cluster_result <- as.data.frame(ego)

  if (nrow(cluster_result) == 0) {
    message(
      "No significant GO terms found for Cluster ",
      cluster_id,
      "."
    )
    next
  }

  cluster_result <- cluster_result %>%
    mutate(
      log10_adjusted_p = -log10(p.adjust),
      Cluster = as.character(cluster_id),
      ClusterLabel = paste("Cluster", cluster_id)
    )

  result_list[[as.character(cluster_id)]] <- cluster_result
}

if (length(result_list) == 0) {
  stop("No significant GO enrichment results were obtained.")
}

all_results <- bind_rows(result_list)

# =========================
# Select top GO terms for visualization
# =========================

top_results <- all_results %>%
  group_by(Cluster, ClusterLabel) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup()

# Generate facet-specific term identifiers to preserve
# term ordering independently within each cluster
plot_data <- top_results %>%
  group_by(ClusterLabel) %>%
  arrange(log10_adjusted_p, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    term_id = paste(
      Description,
      ClusterLabel,
      seq_len(n()),
      sep = "___"
    ),
    term_id = factor(term_id, levels = term_id)
  )

# =========================
# Plot
# =========================

p <- ggplot(
  plot_data,
  aes(
    x = log10_adjusted_p,
    y = term_id
  )
) +
  geom_col(
    fill = bar_color,
    width = 0.75
  ) +
  facet_wrap(
    ~ClusterLabel,
    scales = "free_y",
    ncol = 1
  ) +
  scale_y_discrete(
    labels = function(x) sub("___.*$", "", x)
  ) +
  labs(
    x = expression(-log[10]("adjusted p-value")),
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_blank(),
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      color = "grey85"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p)

# =========================
# Save results
# =========================

write.csv(
  all_results,
  file = file.path(
    output_dir,
    "GO_BP_enrichment_all_results.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  top_results,
  file = file.path(
    output_dir,
    "GO_BP_enrichment_top_terms.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

figure_height <- max(
  4,
  2.0 * length(unique(plot_data$ClusterLabel))
)

ggsave(
  filename = file.path(
    output_dir,
    "GO_BP_enrichment_plot.png"
  ),
  plot = p,
  width = 7,
  height = figure_height,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "GO_BP_enrichment_plot.pdf"
  ),
  plot = p,
  width = 7,
  height = figure_height,
  units = "in",
  bg = "white"
)

message("GO enrichment analysis completed.")
message("Results saved to: ", output_dir)
