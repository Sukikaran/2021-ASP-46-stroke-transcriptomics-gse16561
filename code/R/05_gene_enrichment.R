# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- gene-id-mapping ----
# Install Bioconductor packages if they are missing
annotation_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

for (pkg in annotation_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )

  }

}
# Load human gene annotation database
library(AnnotationDbi)
library(org.Hs.eg.db)
# Retrieve unique gene symbols from all tested genes
gene_symbols <- unique(
  as.character(deg_results$Gene_Symbol)
)
# Verify gene symbols
stopifnot(
  length(gene_symbols) > 0,
  !anyNA(gene_symbols),
  all(nzchar(gene_symbols))
)

cat("Total unique gene symbols:", length(gene_symbols), "\n")
# Retrieve Entrez IDs and official gene names
annotation_table <- AnnotationDbi::select(
  x = org.Hs.eg.db::org.Hs.eg.db,
  keys = gene_symbols,
  columns = c("ENTREZID", "GENENAME"),
  keytype = "SYMBOL"
)
# Remove duplicate annotation records
annotation_table <- unique(annotation_table)
# Return one record per symbol; discard one-to-many Entrez mappings
map_single_gene <- function(symbol) {
  rows <- annotation_table[!is.na(annotation_table$SYMBOL) & annotation_table$SYMBOL == symbol, , drop = FALSE]
  entrez <- unique(na.omit(as.character(rows$ENTREZID)))
  entrez <- entrez[nzchar(entrez)]
  gene_names <- unique(na.omit(as.character(rows$GENENAME)))
  gene_names <- gene_names[nzchar(gene_names)]
  status <- if (!length(entrez)) "Unmapped" else if (length(entrez) > 1) "Ambiguous" else "Mapped"

  data.frame(
    Gene_Symbol = symbol,
    Gene_Name = if (status == "Mapped" && length(gene_names) == 1) gene_names[1] else NA_character_,
    ENTREZID = if (status == "Mapped") entrez[1] else NA_character_,
    Mapping_Status = status,
    stringsAsFactors = FALSE
  )
}
# Map all tested gene symbols
mapping_results <- lapply(
  gene_symbols,
  map_single_gene
)

gene_id_mapping <- do.call(
  rbind,
  mapping_results
)

rownames(gene_id_mapping) <- NULL
# Verify gene identifier mapping
stopifnot(
  nrow(gene_id_mapping) == length(gene_symbols),
  !anyDuplicated(gene_id_mapping$Gene_Symbol),
  all(
    gene_id_mapping$Mapping_Status %in%
      c("Mapped", "Unmapped", "Ambiguous")
  )
)
# Add gene identifiers to the complete DEG results
mapping_index <- match(
  deg_results$Gene_Symbol,
  gene_id_mapping$Gene_Symbol
)

deg_annotated_results <- deg_results

deg_annotated_results$Gene_Name <-
  gene_id_mapping$Gene_Name[mapping_index]

deg_annotated_results$ENTREZID <-
  gene_id_mapping$ENTREZID[mapping_index]

deg_annotated_results$Mapping_Status <-
  gene_id_mapping$Mapping_Status[mapping_index]
# Prepare background gene set from all tested genes
background_entrez <- unique(
  na.omit(gene_id_mapping$ENTREZID)
)

background_entrez <- sort(
  as.character(background_entrez)
)
# Prepare significant Entrez gene identifiers
significant_entrez <- unique(
  na.omit(
    deg_annotated_results$ENTREZID[
      deg_annotated_results$adj.P.Val < fdr_threshold
    ]
  )
)
# Prepare upregulated Entrez gene identifiers
upregulated_entrez <- unique(
  na.omit(
    deg_annotated_results$ENTREZID[
      deg_annotated_results$Regulation ==
        "Upregulated in Stroke"
    ]
  )
)
# Prepare downregulated Entrez gene identifiers
downregulated_entrez <- unique(
  na.omit(
    deg_annotated_results$ENTREZID[
      deg_annotated_results$Regulation ==
        "Downregulated in Stroke"
    ]
  )
)
# Check the enrichment gene sets
if (length(background_entrez) == 0) {
  stop("No genes were successfully mapped to Entrez IDs.")
}

if (length(significant_entrez) == 0) {
  stop("No significant DEGs have valid Entrez IDs.")
}

stopifnot(
  !anyDuplicated(background_entrez),
  !anyDuplicated(significant_entrez),
  all(significant_entrez %in% background_entrez),
  all(upregulated_entrez %in% background_entrez),
  all(downregulated_entrez %in% background_entrez)
)
# Store gene sets for subsequent enrichment analysis
enrichment_gene_sets <- list(
  significant = significant_entrez,
  upregulated = upregulated_entrez,
  downregulated = downregulated_entrez,
  universe = background_entrez
)
# Create an enrichment gene-set table
enrichment_gene_sets_table <- deg_annotated_results[
  ,
  c(
    "Gene_Symbol",
    "Gene_Name",
    "ENTREZID",
    "Regulation",
    "adj.P.Val",
    "Mapping_Status"
  )
]

enrichment_gene_sets_table$In_Enrichment_Background <-
  !is.na(enrichment_gene_sets_table$ENTREZID)

enrichment_gene_sets_table$In_Significant_DEG_Set <-
  enrichment_gene_sets_table$Regulation != "Not significant" &
  enrichment_gene_sets_table$Mapping_Status == "Mapped"
# Display mapping summary
cat("\nGene ID Mapping Summary\n")
cat("-----------------------\n")

cat(
  "Total genes tested:",
  nrow(gene_id_mapping), "\n"
)

cat(
  "Successfully mapped genes:",
  sum(gene_id_mapping$Mapping_Status == "Mapped"), "\n"
)

cat(
  "Unmapped genes:",
  sum(gene_id_mapping$Mapping_Status == "Unmapped"), "\n"
)

cat(
  "Ambiguous mappings:",
  sum(gene_id_mapping$Mapping_Status == "Ambiguous"), "\n"
)

cat(
  "Unique Entrez IDs in background:",
  length(background_entrez), "\n"
)

cat(
  "Significant DEGs with Entrez IDs:",
  length(significant_entrez), "\n"
)

cat(
  "Upregulated genes with Entrez IDs:",
  length(upregulated_entrez), "\n"
)

cat(
  "Downregulated genes with Entrez IDs:",
  length(downregulated_entrez), "\n"
)
# Create mapping summary table
mapping_summary <- data.frame(
  Metric = c(
    "Total genes tested",
    "Successfully mapped genes",
    "Unmapped genes",
    "Ambiguous mappings",
    "Unique Entrez IDs in background",
    "Significant DEGs with Entrez IDs",
    "Upregulated genes with Entrez IDs",
    "Downregulated genes with Entrez IDs"
  ),
  Value = as.character(c(
    nrow(gene_id_mapping),
    sum(gene_id_mapping$Mapping_Status == "Mapped"),
    sum(gene_id_mapping$Mapping_Status == "Unmapped"),
    sum(gene_id_mapping$Mapping_Status == "Ambiguous"),
    length(background_entrez),
    length(significant_entrez),
    length(upregulated_entrez),
    length(downregulated_entrez)
  )),
  stringsAsFactors = FALSE
)
# Export gene identifier mapping
write.csv(
  gene_id_mapping,
  file.path(
    data_dir,
    "GSE16561_gene_id_mapping.csv"
  ),
  row.names = FALSE
)
# Export annotated differential expression results
write.csv(
  deg_annotated_results,
  file.path(
    data_dir,
    "GSE16561_annotated_DEG_results.csv"
  ),
  row.names = FALSE
)
# Export enrichment gene-set information
write.csv(
  enrichment_gene_sets_table,
  file.path(
    data_dir,
    "GSE16561_enrichment_gene_sets.csv"
  ),
  row.names = FALSE
)
# Export mapping summary
write.csv(
  mapping_summary,
  file.path(
    data_dir,
    "GSE16561_enrichment_mapping_summary.csv"
  ),
  row.names = FALSE
)
# Save gene sets for subsequent analysis
saveRDS(
  enrichment_gene_sets,
  file.path(
    data_dir,
    "GSE16561_enrichment_gene_sets.rds"
  )
)
# Verify all output files
mapping_output_files <- c(
  file.path(data_dir, "GSE16561_gene_id_mapping.csv"),
  file.path(data_dir, "GSE16561_annotated_DEG_results.csv"),
  file.path(data_dir, "GSE16561_enrichment_gene_sets.csv"),
  file.path(data_dir, "GSE16561_enrichment_mapping_summary.csv"),
  file.path(data_dir, "GSE16561_enrichment_gene_sets.rds")
)

verify_output_files(mapping_output_files)

cat("\nAll gene annotation output files verified.\n")

cat("\nGene ID annotation completed successfully.\n")

## ---- functional-enrichment-analysis ----
# Validate inputs and prepare the three gene sets
stopifnot(
  requireNamespace("limma", quietly = TRUE),
  requireNamespace("org.Hs.eg.db", quietly = TRUE),
  requireNamespace("ggplot2", quietly = TRUE),
  exists("background_entrez"), exists("significant_entrez"),
  exists("upregulated_entrez"), exists("downregulated_entrez")
)

enrichment_fdr_threshold <- 0.05
gene_sets <- list(
  All = significant_entrez,
  Upregulated = upregulated_entrez,
  Downregulated = downregulated_entrez
)
gene_sets <- lapply(gene_sets, function(x) unique(na.omit(as.character(x))))
background_entrez <- unique(na.omit(as.character(background_entrez)))
significant_entrez <- gene_sets$All
upregulated_entrez <- gene_sets$Upregulated
downregulated_entrez <- gene_sets$Downregulated
stopifnot(
  length(background_entrez) > 0,
  length(gene_sets$All) > 0,
  all(unlist(gene_sets) %in% background_entrez)
)
# GO enrichment; adjust FDR separately within each ontology
run_go_enrichment <- function(gene_ids, analysis_name) {
  result <- as.data.frame(limma::goana(de = gene_ids, universe = background_entrez, species = "Hs"))
  if (!nrow(result)) return(result)
  result$GO_ID <- rownames(result)
  result$Analysis <- analysis_name
  result$Input_Genes <- length(gene_ids)
  result$FDR <- NA_real_
  for (ontology in unique(result$Ont)) {
    idx <- result$Ont == ontology
    result$FDR[idx] <- p.adjust(result$P.DE[idx], method = "BH")
  }
  result$Significant <- !is.na(result$FDR) & result$FDR < enrichment_fdr_threshold
  rownames(result) <- NULL
  result
}

go_results <- setNames(
  lapply(names(gene_sets), function(nm) run_go_enrichment(gene_sets[[nm]], paste(nm, if (nm == "All") "significant DEGs" else "genes"))),
  names(gene_sets)
)
go_bp_all <- subset(go_results$All, Ont == "BP")
go_mf_all <- subset(go_results$All, Ont == "MF")
go_cc_all <- subset(go_results$All, Ont == "CC")
go_bp_upregulated <- subset(go_results$Upregulated, Ont == "BP")
go_bp_downregulated <- subset(go_results$Downregulated, Ont == "BP")
cat("\nGO enrichment analysis completed.\n")
# KEGG enrichment; retain status information when the REST service is unavailable
run_kegg_enrichment <- function(gene_ids, analysis_name) {
  tryCatch({
    result <- as.data.frame(limma::kegga(de = gene_ids, universe = background_entrez, species.KEGG = "hsa"))
    if (nrow(result)) {
      result$KEGG_ID <- rownames(result)
      result$Term <- result$Pathway
      result$Analysis <- analysis_name
      result$Input_Genes <- length(gene_ids)
      result$FDR <- p.adjust(result$P.DE, method = "BH")
      result$Significant <- !is.na(result$FDR) & result$FDR < enrichment_fdr_threshold
      rownames(result) <- NULL
    }
    list(table = result, status = "Completed")
  }, error = function(e) {
    cat("\nKEGG analysis unavailable for", analysis_name, ":", conditionMessage(e), "\n")
    list(table = data.frame(), status = paste("Unavailable:", conditionMessage(e)))
  })
}

kegg_results <- setNames(
  lapply(names(gene_sets), function(nm) run_kegg_enrichment(gene_sets[[nm]], paste(nm, if (nm == "All") "significant DEGs" else "genes"))),
  names(gene_sets)
)
kegg_status <- setNames(
  vapply(kegg_results, `[[`, character(1), "status"),
  paste0("KEGG_", names(kegg_results))
)

enrichment_results <- list(
  GO_BP_All = go_bp_all,
  GO_MF_All = go_mf_all,
  GO_CC_All = go_cc_all,
  GO_BP_Upregulated = go_bp_upregulated,
  GO_BP_Downregulated = go_bp_downregulated,
  KEGG_All = kegg_results$All$table,
  KEGG_Upregulated = kegg_results$Upregulated$table,
  KEGG_Downregulated = kegg_results$Downregulated$table
)
# Export each enrichment table and a single summary table
enrichment_csv_files <- setNames(
  file.path(data_dir, paste0("GSE16561_", tolower(names(enrichment_results)), "_enrichment.csv")),
  names(enrichment_results)
)
for (nm in names(enrichment_results)) {
  tab <- enrichment_results[[nm]]
  if (!ncol(tab)) tab <- data.frame(Note = "No results were returned. Check the analysis status.")
  write.csv(tab, enrichment_csv_files[[nm]], row.names = FALSE)
}

enrichment_summary <- do.call(rbind, lapply(names(enrichment_results), function(nm) {
  tab <- enrichment_results[[nm]]
  fdr <- if ("FDR" %in% names(tab)) tab$FDR else numeric()
  data.frame(
    Analysis = nm,
    Returned_Terms = nrow(tab),
    Significant_Terms_FDR_0.05 = sum(!is.na(fdr) & fdr < enrichment_fdr_threshold),
    Status = if (startsWith(nm, "KEGG_")) unname(kegg_status[nm]) else "Completed",
    stringsAsFactors = FALSE
  )
}))
enrichment_summary_file <- file.path(data_dir, "GSE16561_enrichment_summary.csv")
write.csv(enrichment_summary, enrichment_summary_file, row.names = FALSE)
# Create a dot plot for significant terms, or a clear empty-result panel
create_enrichment_plot <- function(result_table, plot_title, output_file) {
  required <- c("FDR", "Term", "DE", "Input_Genes")
  significant <- if (nrow(result_table) && all(required %in% names(result_table))) {
    result_table[!is.na(result_table$FDR) & result_table$FDR < enrichment_fdr_threshold, , drop = FALSE]
  } else data.frame()

  if (nrow(significant)) {
    significant$GeneRatio <- significant$DE / significant$Input_Genes
    significant <- head(significant[order(significant$FDR), , drop = FALSE], 15)
    plot <- ggplot2::ggplot(significant, ggplot2::aes(GeneRatio, reorder(Term, GeneRatio))) +
      ggplot2::geom_point(ggplot2::aes(size = DE, color = FDR), alpha = 0.85) +
      ggplot2::scale_color_gradient(low = "#D7191C", high = "#2C7BB6") +
      ggplot2::labs(title = plot_title, x = "Gene Ratio", y = NULL, color = "Adjusted p-value", size = "Gene Count") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
  } else {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No significantly enriched terms (FDR < 0.05).") +
      ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = plot_title) + ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      )
  }
  ggplot2::ggsave(output_file, plot, width = 10, height = 7, dpi = 300, bg = "white")
  plot
}

figure_names <- c(
  "fig6_GO_BP_all_DEGs.png", "fig7_GO_MF_all_DEGs.png", "fig8_GO_CC_all_DEGs.png",
  "fig9_GO_BP_upregulated_DEGs.png", "fig10_GO_BP_downregulated_DEGs.png",
  "fig11_KEGG_all_DEGs.png", "fig12_KEGG_upregulated_DEGs.png", "fig13_KEGG_downregulated_DEGs.png"
)
enrichment_plot_files <- setNames(file.path(figure_dir, figure_names), names(enrichment_results))
enrichment_plots <- setNames(lapply(names(enrichment_results), function(nm) {
  p <- create_enrichment_plot(enrichment_results[[nm]], paste("Enrichment Analysis:", gsub("_", " ", nm)), enrichment_plot_files[[nm]])
  print(p)
  p
}), names(enrichment_results))

enrichment_rds_file <- file.path(data_dir, "GSE16561_functional_enrichment_results.rds")
saveRDS(enrichment_results, enrichment_rds_file)
verify_output_files(c(enrichment_csv_files, enrichment_summary_file, enrichment_plot_files, enrichment_rds_file))

cat("\nFunctional Enrichment Summary\n-----------------------------\n")
print(enrichment_summary, row.names = FALSE)
cat("\nMapped background genes:", length(background_entrez),
    "\nSignificant genes analysed:", length(significant_entrez),
    "\nUpregulated genes analysed:", length(upregulated_entrez),
    "\nDownregulated genes analysed:", length(downregulated_entrez),
    "\nAll enrichment output files verified.\n")