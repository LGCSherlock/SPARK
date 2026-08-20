#!/usr/bin/env Rscript

# ============================================================
# Heatmap visualization of kinetic parameters
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(scales)
})

# ============================================================
# Command-line arguments
# ============================================================
#
# Usage:
# Rscript parameter_heatmaps.R \
#   experimental.csv \
#   control.csv \
#   output_directory
#
# Input files must contain:
#   gene, gamma, lambda, alpha, beta
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript parameter_heatmaps.R ",
    "<experimental.csv> <control.csv> <output_directory>"
  )
}

experimental_file <- args[1]
control_file <- args[2]
output_dir <- args[3]

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Parameters
# ============================================================

parameters <- c(
  "gamma",
  "lambda",
  "alpha",
  "beta"
)

panel_labels <- c(
  gamma = "bold(gamma)",
  lambda = "bold(lambda)",
  alpha = "bold(alpha)",
  beta = "bold(beta)"
)

group_levels <- c(
  "DEX-treated",
  "Control"
)

figure_width <- 7.1
figure_height <- 9.3
figure_resolution <- 1200

# ============================================================
# Input validation
# ============================================================

read_and_validate <- function(path, label) {

  if (!file.exists(path)) {
    stop(label, " file not found: ", path)
  }

  dat <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns <- c(
    "gene",
    parameters
  )

  missing_columns <- setdiff(
    required_columns,
    names(dat)
  )

  if (length(missing_columns) > 0) {
    stop(
      label,
      " file is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat <- dat[
    ,
    required_columns,
    drop = FALSE
  ]

  if (anyNA(dat$gene) || any(dat$gene == "")) {
    stop(
      label,
      " file contains missing or empty gene names."
    )
  }

  if (anyDuplicated(dat$gene)) {
    stop(
      label,
      " file contains duplicated gene names."
    )
  }

  numeric_check <- vapply(
    dat[parameters],
    is.numeric,
    logical(1)
  )

  if (!all(numeric_check)) {
    stop(
      label,
      " file contains non-numeric parameter values."
    )
  }

  if (anyNA(dat[parameters])) {
    stop(
      label,
      " file contains missing parameter values."
    )
  }

  if (any(!is.finite(as.matrix(dat[parameters])))) {
    stop(
      label,
      " file contains non-finite parameter values."
    )
  }

  dat
}

# ============================================================
# Load data
# ============================================================

experimental <- read_and_validate(
  experimental_file,
  "Experimental"
)

control <- read_and_validate(
  control_file,
  "Control"
)

if (!setequal(
  experimental$gene,
  control$gene
)) {
  stop(
    "Experimental and control files must contain the same genes."
  )
}

# Reorder the control dataset to match the experimental dataset
control <- control[
  match(
    experimental$gene,
    control$gene
  ),
  ,
  drop = FALSE
]

# ============================================================
# Prepare plotting data
# ============================================================

# Preserve the input gene order in all panels.
gene_levels <- rev(
  experimental$gene
)

heatmap_data <- do.call(
  rbind,
  lapply(
    parameters,
    function(parameter) {

      rbind(
        data.frame(
          Gene = experimental$gene,
          Parameter = parameter,
          Group = group_levels[1],
          Value = experimental[[parameter]],
          stringsAsFactors = FALSE
        ),
        data.frame(
          Gene = control$gene,
          Parameter = parameter,
          Group = group_levels[2],
          Value = control[[parameter]],
          stringsAsFactors = FALSE
        )
      )
    }
  )
)

heatmap_data$Gene <- factor(
  heatmap_data$Gene,
  levels = gene_levels
)

heatmap_data$Group <- factor(
  heatmap_data$Group,
  levels = group_levels
)

heatmap_data$Parameter <- factor(
  heatmap_data$Parameter,
  levels = parameters
)

# Save source data used for figure generation
write.csv(
  heatmap_data,
  file = file.path(
    output_dir,
    "parameter_heatmap_source_data.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Legend formatting
# ============================================================

format_legend <- function(parameter) {

  if (parameter == "alpha") {

    label_scientific(
      digits = 2
    )

  } else {

    label_number(
      accuracy = 0.001,
      trim = TRUE
    )
  }
}

# ============================================================
# Heatmap function
# ============================================================

make_heatmap <- function(parameter) {

  panel_data <- heatmap_data[
    heatmap_data$Parameter == parameter,
    ,
    drop = FALSE
  ]

  values <- panel_data$Value

  base_plot <- ggplot(
    panel_data,
    aes(
      x = Group,
      y = Gene,
      fill = Value
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.25
    ) +
    scale_x_discrete(
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      expand = c(0, 0),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = parse(
        text = panel_labels[[parameter]]
      )[[1]],
      fill = "Value"
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = unit(
          28,
          "mm"
        ),
        barwidth = unit(
          3.2,
          "mm"
        ),
        ticks = TRUE,
        frame.colour = "black",
        frame.linewidth = 0.35
      )
    ) +
    theme_minimal(
      base_size = 7,
      base_family = "sans"
    ) +
    theme(
      panel.grid = element_blank(),

      axis.text.x = element_text(
        color = "black",
        size = 6.5,
        face = "bold",
        margin = margin(t = 3)
      ),

      axis.text.y = element_text(
        color = "black",
        size = 5.3,
        margin = margin(r = 2)
      ),

      axis.ticks = element_blank(),

      plot.title = element_text(
        color = "black",
        size = 10,
        hjust = 0.5,
        margin = margin(b = 4)
      ),

      legend.position = "right",

      legend.title = element_text(
        color = "black",
        size = 6.5
      ),

      legend.text = element_text(
        color = "black",
        size = 5.5
      ),

      legend.ticks = element_line(
        color = "black",
        linewidth = 0.35
      ),

      legend.margin = margin(l = 2),

      plot.margin = margin(
        5,
        5,
        7,
        5
      )
    )

  # Alpha is displayed using a sequential scale because its values
  # are non-negative in this analysis.
  if (parameter == "alpha") {

    value_limits <- range(
      values,
      finite = TRUE
    )

    if (value_limits[1] == value_limits[2]) {
      padding <- max(
        abs(value_limits[1]) * 0.01,
        .Machine$double.eps
      )

      value_limits <- c(
        value_limits[1] - padding,
        value_limits[2] + padding
      )
    }

    legend_breaks <- c(
      value_limits[1],
      mean(value_limits),
      value_limits[2]
    )

    base_plot +
      scale_fill_gradientn(
        colours = c(
          "#F7FBFF",
          "#6BAED6",
          "#08306B"
        ),
        limits = value_limits,
        breaks = legend_breaks,
        labels = format_legend(parameter),
        oob = squish
      )

  } else {

    max_absolute <- max(
      abs(values),
      na.rm = TRUE
    )

    if (max_absolute == 0) {
      max_absolute <- 1
    }

    base_plot +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        limits = c(
          -max_absolute,
          max_absolute
        ),
        breaks = c(
          -max_absolute,
          0,
          max_absolute
        ),
        labels = format_legend(parameter),
        oob = squish
      )
  }
}

# ============================================================
# Generate individual heatmaps
# ============================================================

plots <- lapply(
  parameters,
  make_heatmap
)

# ============================================================
# Combine heatmaps
# ============================================================

draw_combined_figure <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        widths = unit(
          c(1, 1),
          "null"
        ),
        heights = unit(
          c(1, 1),
          "null"
        )
      )
    )
  )

  for (i in seq_along(plots)) {

    row_position <- (
      (i - 1L) %/% 2L
    ) + 1L

    column_position <- (
      (i - 1L) %% 2L
    ) + 1L

    print(
      plots[[i]],
      vp = viewport(
        layout.pos.row = row_position,
        layout.pos.col = column_position
      ),
      newpage = FALSE
    )
  }

  popViewport()
}

# ============================================================
# Save figures
# ============================================================

cairo_pdf(
  filename = file.path(
    output_dir,
    "parameter_heatmaps.pdf"
  ),
  width = figure_width,
  height = figure_height,
  family = "sans"
)

draw_combined_figure()
dev.off()

tiff(
  filename = file.path(
    output_dir,
    "parameter_heatmaps.tiff"
  ),
  width = figure_width,
  height = figure_height,
  units = "in",
  res = figure_resolution,
  compression = "lzw",
  type = "cairo"
)

draw_combined_figure()
dev.off()

png(
  filename = file.path(
    output_dir,
    "parameter_heatmaps.png"
  ),
  width = figure_width,
  height = figure_height,
  units = "in",
  res = figure_resolution,
  type = "cairo"
)

draw_combined_figure()
dev.off()

cat(
  "Parameter heatmaps generated for",
  nrow(experimental),
  "paired genes.\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)
