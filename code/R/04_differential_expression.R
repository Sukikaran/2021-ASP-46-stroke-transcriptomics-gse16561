# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- limma-design-setup ----
# Set Control as the reference group
metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Stroke")
)
# Verify sample alignment and group information
stopifnot(
  ncol(expr_gene) == nrow(metadata),
  identical(colnames(expr_gene), metadata$Sample_ID),
  !anyNA(metadata$Group),
  !anyNA(expr_gene),
  sum(metadata$Group == "Control") == 24,
  sum(metadata$Group == "Stroke") == 39
)
# Create the experimental design matrix
design_matrix <- model.matrix(
  ~ Group,
  data = metadata
)
# Assign sample IDs to the design matrix rows
rownames(design_matrix) <- metadata$Sample_ID
# Display design matrix information
cat("Design matrix dimensions:\n")
print(dim(design_matrix))

cat("\nDesign matrix columns:\n")
print(colnames(design_matrix))

cat("\nFirst six rows of the design matrix:\n")
print(head(design_matrix))
# Define the Stroke versus Control comparison
contrast_matrix <- limma::makeContrasts(
  Stroke_vs_Control = GroupStroke,
  levels = design_matrix
)
# Display the contrast matrix
cat("\nContrast matrix:\n")
print(contrast_matrix)
# Export the design matrix
design_export <- data.frame(
  Sample_ID = rownames(design_matrix),
  design_matrix,
  check.names = FALSE
)

write.csv(
  design_export,
  file.path(data_dir, "GSE16561_limma_design_matrix.csv"),
  row.names = FALSE
)
# Export the contrast matrix
write.csv(
  contrast_matrix,
  file.path(data_dir, "GSE16561_limma_contrast_matrix.csv"),
  row.names = TRUE
)
# Verify the generated files
limma_setup_files <- c(
  file.path(data_dir, "GSE16561_limma_design_matrix.csv"),
  file.path(data_dir, "GSE16561_limma_contrast_matrix.csv")
)

verify_output_files(limma_setup_files)

cat("\nLimma design and contrast setup completed successfully.\n")

## ---- limma-model-fitting ----
# Verify the experimental design and contrast
stopifnot(
  identical(colnames(expr_gene), rownames(design_matrix)),
  nrow(design_matrix) == 63,
  ncol(design_matrix) == 2,
  nrow(contrast_matrix) == 2,
  ncol(contrast_matrix) == 1,
  colnames(contrast_matrix)[1] == "Stroke_vs_Control"
)
# Fit a linear model for each gene
fit_initial <- limma::lmFit(
  expr_gene,
  design_matrix
)
# Apply the Stroke versus Control contrast
fit_contrasted <- limma::contrasts.fit(
  fit_initial,
  contrast_matrix
)
# Apply robust empirical Bayes moderation
fit_ebayes <- limma::eBayes(
  fit_contrasted,
  robust = TRUE
)
# Verify model results
stopifnot(
  nrow(fit_ebayes$coefficients) == nrow(expr_gene),
  ncol(fit_ebayes$coefficients) == 1,
  identical(
    colnames(fit_ebayes$coefficients),
    "Stroke_vs_Control"
  ),
  all(is.finite(fit_ebayes$coefficients)),
  all(is.finite(fit_ebayes$p.value)),
  all(fit_ebayes$p.value >= 0),
  all(fit_ebayes$p.value <= 1)
)
# Display model information
cat("Limma model fitting summary:\n")

cat(
  "Number of genes analysed:",
  nrow(fit_ebayes$coefficients), "\n"
)

cat(
  "Number of clinical samples:",
  ncol(expr_gene), "\n"
)

cat(
  "Comparison:",
  colnames(fit_ebayes$coefficients), "\n"
)

cat(
  "Median residual degrees of freedom:",
  median(fit_ebayes$df.residual), "\n"
)

cat(
  "Number of genes with valid p-values:",
  sum(is.finite(fit_ebayes$p.value)), "\n"
)

cat("Robust empirical Bayes moderation: Applied\n")
# Save the fitted model for subsequent analysis
model_rds_file <- file.path(
  data_dir,
  "GSE16561_limma_ebayes_fit.rds"
)

saveRDS(
  fit_ebayes,
  model_rds_file
)
# Create a summary of the model fitting stage
limma_model_summary <- data.frame(
  Metric = c(
    "Number of genes analysed",
    "Number of clinical samples",
    "Control samples",
    "Stroke samples",
    "Comparison",
    "Statistical method",
    "Empirical Bayes moderation",
    "Robust moderation",
    "Median residual degrees of freedom"
  ),
  Value = c(
    nrow(fit_ebayes$coefficients),
    ncol(expr_gene),
    sum(metadata$Group == "Control"),
    sum(metadata$Group == "Stroke"),
    "Stroke versus Control",
    "Limma linear model",
    "Applied",
    "Applied",
    round(median(fit_ebayes$df.residual), 2)
  )
)
# Export the model summary
write.csv(
  limma_model_summary,
  file.path(data_dir, "GSE16561_limma_model_summary.csv"),
  row.names = FALSE
)
# Verify generated output files
limma_model_files <- c(
  model_rds_file,
  file.path(data_dir, "GSE16561_limma_model_summary.csv")
)

verify_output_files(limma_model_files)

cat("\nLimma model fitting completed successfully.\n")

## ---- deg-identification ----
# Define the false discovery rate significance threshold
fdr_threshold <- 0.05
# Extract differential expression results for all genes
deg_results <- limma::topTable(
  fit_ebayes,
  coef = "Stroke_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P",
  p.value = 1
)
# Add gene symbols as an explicit column
deg_results$Gene_Symbol <- rownames(deg_results)
# Arrange the result columns
deg_results <- deg_results[
  ,
  c(
    "Gene_Symbol",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]
# Verify the differential expression results
stopifnot(
  nrow(deg_results) == nrow(expr_gene),
  !anyDuplicated(deg_results$Gene_Symbol),
  all(is.finite(deg_results$P.Value)),
  all(is.finite(deg_results$adj.P.Val)),
  all(deg_results$P.Value >= 0 & deg_results$P.Value <= 1),
  all(deg_results$adj.P.Val >= 0 & deg_results$adj.P.Val <= 1)
)
# Identify statistically significant genes
is_significant <- deg_results$adj.P.Val < fdr_threshold
# Classify the direction of differential expression
deg_results$Regulation <- ifelse(
  !is_significant,
  "Not significant",
  ifelse(
    deg_results$logFC > 0,
    "Upregulated in Stroke",
    ifelse(
      deg_results$logFC < 0,
      "Downregulated in Stroke",
      "Significant (zero logFC)"
    )
  )
)
# Extract significant differentially expressed genes
significant_degs <- deg_results[
  is_significant,
  ,
  drop = FALSE
]
# Extract significant upregulated genes
upregulated_degs <- deg_results[
  deg_results$Regulation == "Upregulated in Stroke",
  ,
  drop = FALSE
]
# Extract significant downregulated genes
downregulated_degs <- deg_results[
  deg_results$Regulation == "Downregulated in Stroke",
  ,
  drop = FALSE
]
# Display differential expression summary
cat("Differential Expression Analysis Summary\n")
cat("-----------------------------------------\n")

cat("Total genes tested:", nrow(deg_results), "\n")

cat("FDR significance threshold:", fdr_threshold, "\n")

cat("Significant DEGs:", nrow(significant_degs), "\n")

cat("Upregulated genes:", nrow(upregulated_degs), "\n")

cat("Downregulated genes:", nrow(downregulated_degs), "\n")

cat(
  "Genes not significant:",
  sum(!is_significant),
  "\n"
)
# Display the top 10 genes ranked by p-value
cat("\nTop 10 genes ranked by p-value:\n")

print(
  head(deg_results, 10),
  row.names = FALSE
)
# Export all differential expression results
write.csv(
  deg_results,
  file.path(
    data_dir,
    "GSE16561_differential_expression_all_genes.csv"
  ),
  row.names = FALSE
)
# Export significant differentially expressed genes
write.csv(
  significant_degs,
  file.path(
    data_dir,
    "GSE16561_significant_DEGs.csv"
  ),
  row.names = FALSE
)
# Create a summary table
deg_summary <- data.frame(
  Metric = c(
    "Total genes tested",
    "FDR significance threshold",
    "Significant differentially expressed genes",
    "Upregulated genes in Stroke",
    "Downregulated genes in Stroke",
    "Genes not significant",
    "Minimum adjusted p-value"
  ),
  Value = c(
    as.character(nrow(deg_results)),
    as.character(fdr_threshold),
    as.character(nrow(significant_degs)),
    as.character(nrow(upregulated_degs)),
    as.character(nrow(downregulated_degs)),
    as.character(sum(!is_significant)),
    format(
      min(deg_results$adj.P.Val),
      scientific = TRUE,
      digits = 4
    )
  ),
  stringsAsFactors = FALSE
)
# Export the DEG summary
write.csv(
  deg_summary,
  file.path(
    data_dir,
    "GSE16561_DEG_summary.csv"
  ),
  row.names = FALSE
)
# Verify all output files
deg_output_files <- c(
  file.path(
    data_dir,
    "GSE16561_differential_expression_all_genes.csv"
  ),
  file.path(
    data_dir,
    "GSE16561_significant_DEGs.csv"
  ),
  file.path(
    data_dir,
    "GSE16561_DEG_summary.csv"
  )
)

verify_output_files(deg_output_files)

cat("\nAll differential expression output files verified.\n")

cat("\nDifferential expression analysis completed successfully.\n")

## ---- volcano-plot ----
# Prepare data for volcano plot
volcano_data <- deg_results
# Calculate negative log10 of adjusted p-values
# A small lower limit prevents infinite values when p-values are zero
volcano_data$neg_log10_FDR <- -log10(
  pmax(
    volcano_data$adj.P.Val,
    .Machine$double.xmin
  )
)
# Set the order of regulation categories
volcano_data$Regulation <- factor(
  volcano_data$Regulation,
  levels = c(
    "Upregulated in Stroke",
    "Downregulated in Stroke",
    "Not significant",
    "Significant (zero logFC)"
  )
)
# Create volcano plot
volcano_plot <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(
    x = logFC,
    y = neg_log10_FDR,
    colour = Regulation
  )
) +
  # Display genes as points
  ggplot2::geom_point(
    alpha = 0.65,
    size = 1.3
  ) +
  # Add horizontal FDR significance threshold
  ggplot2::geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed",
    colour = "black"
  ) +
  # Add vertical reference line at zero log fold change
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "black"
  ) +
  # Define colours for regulation categories
  ggplot2::scale_colour_manual(
    values = c(
      "Upregulated in Stroke" = "#D7191C",
      "Downregulated in Stroke" = "#2C7BB6",
      "Not significant" = "grey75",
      "Significant (zero logFC)" = "#984EA3"
    ),
    drop = TRUE
  ) +
  # Add plot titles and labels
  ggplot2::labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    subtitle = paste(
      "Stroke versus Control | FDR threshold:",
      fdr_threshold
    ),
    x = "Log2 Fold Change (Stroke vs Control)",
    y = "-log10(BH-adjusted p-value)",
    colour = "Gene Regulation"
  ) +
  # Apply plot formatting
  ggplot2::theme_minimal(base_size = 12) +

  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    legend.position = "bottom",
    # White background for the complete plot
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
# Define volcano plot output file
volcano_plot_file <- file.path(
  figure_dir,
  "fig5_volcano_plot.png"
)
# Save volcano plot with a white background
ggplot2::ggsave(
  filename = volcano_plot_file,
  plot = volcano_plot,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)
# Verify volcano plot data and output
stopifnot(
  nrow(volcano_data) == nrow(deg_results),
  all(is.finite(volcano_data$neg_log10_FDR))
)
verify_output_files(volcano_plot_file)
# Display plot in the report
print(volcano_plot)

cat("\nVolcano plot generated successfully.\n")

cat(
  "Total genes displayed:",
  nrow(volcano_data),
  "\n"
)

cat(
  "Significant genes:",
  sum(volcano_data$adj.P.Val < fdr_threshold),
  "\n"
)

cat(
  "Upregulated genes:",
  nrow(upregulated_degs),
  "\n"
)

cat(
  "Downregulated genes:",
  nrow(downregulated_degs),
  "\n"
)