# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- eda-preparation ----
# Load visualization packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(ggplot2)
library(pheatmap)
# Verify cleaned expression matrix
stopifnot(
  ncol(expr_gene) == 63,
  nrow(metadata) == 63,
  !anyNA(expr_gene)
)
# Match metadata to expression samples
metadata <- metadata[
  match(
    colnames(expr_gene),
    metadata$Sample_ID
  ),
  ,
  drop = FALSE
]
# Verify sample alignment
stopifnot(
  identical(
    colnames(expr_gene),
    metadata$Sample_ID
  )
)
# Display group distribution
print(table(metadata$Group))
# Define group colours
group_colors <- c(
  "Control" = "#2C7BB6",
  "Stroke" = "#D7191C"
)
# Ensure figures directory exists
dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cat("EDA preparation completed successfully.\n")

## ---- expression-boxplot ----
# Create sample colours based on clinical groups
sample_colors <- unname(
  group_colors[as.character(metadata$Group)]
)
# Define output file
boxplot_file <- file.path(
  figure_dir,
  "fig1_expression_boxplot.png"
)
# Save the boxplot with a white background
png(
  filename = boxplot_file,
  width = 2800,
  height = 1500,
  res = 200,
  bg = "white"
)

par(mar = c(11, 4, 4, 2) + 0.1)

boxplot(
  expr_gene,
  names = metadata$Sample_ID,
  las = 2,
  col = sample_colors,
  outline = FALSE,
  cex.axis = 0.55,
  ylab = "Log2 Gene Expression",
  xlab = "",
  main = "Gene Expression Distribution Across 63 Samples"
)

legend(
  "topright",
  legend = names(group_colors),
  fill = group_colors,
  title = "Clinical Group",
  cex = 0.8
)

dev.off()
# Display the saved figure in the report
knitr::include_graphics(boxplot_file)

## ---- expression-density ----
# Set seed for reproducibility
set.seed(123)
# Randomly select 3,000 genes
selected_genes <- sample(
  seq_len(nrow(expr_gene)),
  size = min(3000, nrow(expr_gene))
)
# Calculate density for each sample
density_data <- do.call(
  rbind,
  lapply(seq_len(ncol(expr_gene)), function(i) {
    sample_density <- density(
      expr_gene[selected_genes, i],
      na.rm = TRUE
    )
    data.frame(
      Expression = sample_density$x,
      Density = sample_density$y,
      Sample_ID = colnames(expr_gene)[i],
      Group = as.character(metadata$Group[i])
    )
  })
)
# Set clinical group order
density_data$Group <- factor(
  density_data$Group,
  levels = c("Control", "Stroke")
)
# Create density plot
density_plot <- ggplot2::ggplot(
  density_data,
  ggplot2::aes(
    x = Expression,
    y = Density,
    group = Sample_ID,
    colour = Group
  )
) +
  ggplot2::geom_line(
    alpha = 0.4,
    linewidth = 0.5
  ) +
  ggplot2::facet_wrap(~Group, ncol = 1) +
  ggplot2::scale_colour_manual(values = group_colors) +
  ggplot2::labs(
    title = "Gene Expression Density Distribution",
    subtitle = "Distribution across 3,000 selected genes per sample",
    x = "Log2 Gene Expression",
    y = "Density",
    colour = "Clinical Group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold"),
    # White background for the complete figure
    plot.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    # White background for the plotting panels
    panel.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    )
  )
# Save the figure with a white background
density_file <- file.path(
  figure_dir,
  "fig2_expression_density.png"
)

ggplot2::ggsave(
  filename = density_file,
  plot = density_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
# Display the saved figure in the report
knitr::include_graphics(density_file)

## ---- pca-analysis ----
# Calculate the variance of each gene
gene_variances <- apply(
  expr_gene,
  1,
  var
)
# Select the 500 most variable genes
top_variable_genes <- order(
  gene_variances,
  decreasing = TRUE
)[seq_len(min(500, length(gene_variances)))]

pca_matrix <- expr_gene[
  top_variable_genes,
  ,
  drop = FALSE
]
# Perform PCA using samples as observations
pca_result <- prcomp(
  t(pca_matrix),
  center = TRUE,
  scale. = FALSE
)
# Calculate percentage of variance explained
variance_explained <- (
  pca_result$sdev^2 /
  sum(pca_result$sdev^2)
) * 100
# Extract PC1 and PC2 scores
pca_scores <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)
# Attach clinical group information
pca_scores$Group <- metadata$Group[
  match(pca_scores$Sample_ID, metadata$Sample_ID)
]

pca_scores$Group <- factor(
  pca_scores$Group,
  levels = c("Control", "Stroke")
)
# Verify that all PCA samples have matching clinical groups
stopifnot(
  nrow(pca_scores) == 63,
  !anyNA(pca_scores$Group)
)
# Display variance explained
cat(
  "PC1 variance explained:",
  round(variance_explained[1], 2), "%\n"
)

cat(
  "PC2 variance explained:",
  round(variance_explained[2], 2), "%\n"
)
# Create PCA plot
pca_plot <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = Group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey70"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey70"
  ) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.85
  ) +
  ggplot2::scale_colour_manual(values = group_colors) +
  ggplot2::labs(
    title = "PCA of Gene Expression Profiles",
    subtitle = "Based on the 500 most variable genes",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 2),
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 2),
      "% variance)"
    ),
    colour = "Clinical Group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    # White background for the complete figure
    plot.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    # White background for the plotting panel
    panel.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    )
  )
# Save PCA plot with a white background
pca_file <- file.path(
  figure_dir,
  "fig3_pca_plot.png"
)

ggplot2::ggsave(
  filename = pca_file,
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)
# Export PCA scores
write.csv(
  pca_scores,
  file.path(
    data_dir,
    "GSE16561_PCA_scores.csv"
  ),
  row.names = FALSE
)
# Display plot in the report
knitr::include_graphics(pca_file)

## ---- sample-correlation ----
# Calculate Pearson correlation between samples
sample_correlation <- cor(
  expr_gene,
  method = "pearson"
)
# Create sample group annotation
sample_annotation <- data.frame(
  Group = metadata$Group
)

rownames(sample_annotation) <- metadata$Sample_ID
# Ensure annotation follows correlation matrix order
sample_annotation <- sample_annotation[
  rownames(sample_correlation),
  ,
  drop = FALSE
]
# Define output file
correlation_file <- file.path(
  figure_dir,
  "fig4_sample_correlation_heatmap.png"
)
# Generate a cleaner correlation heatmap
correlation_heatmap <- pheatmap(
  sample_correlation,
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  breaks = seq(-1, 1, length.out = 101),
  # Show group annotation only once, across the top
  annotation_col = sample_annotation,
  annotation_colors = list(
    Group = group_colors
  ),
  annotation_names_col = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 40,
  treeheight_col = 40,
  main = "Sample-to-Sample Pearson Correlation",
silent = TRUE
)
# Save the heatmap using an explicitly managed device
grDevices::png(
  filename = correlation_file,
  width = 2400,
  height = 2200,
  res = 200,
  bg = "white"
)

grid::grid.newpage()
grid::grid.draw(correlation_heatmap$gtable)

grDevices::dev.off()
# Export correlation matrix
write.csv(
  sample_correlation,
  file.path(
    data_dir,
    "GSE16561_sample_correlation.csv"
  )
)
# Display heatmap in the report
knitr::include_graphics(correlation_file)

## ---- eda-summary ----
# Calculate pairwise sample correlations
correlation_values <- sample_correlation[
  upper.tri(sample_correlation)
]
# Create EDA summary table
eda_summary <- data.frame(
  Metric = c(
    "Total unique genes",
    "Total clinical samples",
    "Control samples",
    "Stroke samples",
    "Mean gene expression",
    "Standard deviation of gene expression",
    "Minimum gene expression",
    "Maximum gene expression",
    "Number of genes used for PCA",
    "PC1 variance explained (%)",
    "PC2 variance explained (%)",
    "Mean sample correlation",
    "Minimum sample correlation",
    "Maximum sample correlation"
  ),
  Value = as.character(c(
    nrow(expr_gene),
    ncol(expr_gene),
    sum(metadata$Group == "Control"),
    sum(metadata$Group == "Stroke"),
    round(mean(expr_gene), 4),
    round(sd(as.vector(expr_gene)), 4),
    round(min(expr_gene), 4),
    round(max(expr_gene), 4),
    length(top_variable_genes),
    round(variance_explained[1], 2),
    round(variance_explained[2], 2),
    round(mean(correlation_values), 4),
    round(min(correlation_values), 4),
    round(max(correlation_values), 4)
  )),
  stringsAsFactors = FALSE
)
# Display summary table
knitr::kable(
  eda_summary,
  caption = "Exploratory Data Analysis Summary"
)
# Export summary table
write.csv(
  eda_summary,
  file.path(
    data_dir,
    "GSE16561_EDA_summary.csv"
  ),
  row.names = FALSE
)
# Verify generated figures and result files
eda_required_files <- c(
  file.path(figure_dir, "fig1_expression_boxplot.png"),
  file.path(figure_dir, "fig2_expression_density.png"),
  file.path(figure_dir, "fig3_pca_plot.png"),
  file.path(figure_dir, "fig4_sample_correlation_heatmap.png"),
  file.path(data_dir, "GSE16561_PCA_scores.csv"),
  file.path(data_dir, "GSE16561_sample_correlation.csv"),
  file.path(data_dir, "GSE16561_EDA_summary.csv")
)
# Check whether each expected output file exists
eda_file_exists <- file.exists(eda_required_files)
# Get file sizes safely, including when a file is missing
eda_file_sizes <- rep(0, length(eda_required_files))

if (any(eda_file_exists)) {
  eda_file_sizes[eda_file_exists] <- file.info(
    eda_required_files[eda_file_exists]
  )$size
}

eda_output_status <- data.frame(
  File = basename(eda_required_files),
  Exists = eda_file_exists,
  Size_KB = round(eda_file_sizes / 1024, 2)
)

print(eda_output_status)
# Stop if any expected output is missing or empty
verify_output_files(eda_required_files)

cat("All EDA output files verified successfully.\n")