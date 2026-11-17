# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- metadata-cleaning ----
# Extract GEO sample identifiers
sample_ids <- as.character(
  pdata_raw$geo_accession
)
# Match metadata to the order of the expression matrix
sample_index <- match(
  colnames(expr_raw),
  sample_ids
)

if (anyNA(sample_index)) {
  stop("Could not match sample identifiers.")
}
# Reorder metadata to match expression columns
pdata_ordered <- pdata_raw[
  sample_index,
  ,
  drop = FALSE
]

sample_ids <- sample_ids[sample_index]
# Combine sample description fields
text_columns <- vapply(
  pdata_ordered,
  function(x) is.character(x) || is.factor(x),
  logical(1)
)

sample_description <- apply(
  as.data.frame(
    lapply(
      pdata_ordered[, text_columns, drop = FALSE],
      as.character
    )
  ),
  1,
  paste,
  collapse = " "
)
# Avoid interpreting non-stroke as stroke
sample_description <- gsub(
  "non[- ]stroke",
  " ",
  sample_description,
  ignore.case = TRUE
)
# Identify stroke and control samples
stroke_flag <- grepl(
  "\\bstroke\\b|ischemic",
  sample_description,
  ignore.case = TRUE
)

control_flag <- grepl(
  "\\bcontrol\\b|healthy",
  sample_description,
  ignore.case = TRUE
)
# Detect ambiguous sample classifications
if (any(stroke_flag == control_flag)) {

  print(
    sample_description[stroke_flag == control_flag]
  )

  stop("Some samples could not be classified correctly.")

}
# Create the clinical metadata table
metadata <- data.frame(
  Sample_ID = sample_ids,
  Title = as.character(pdata_ordered$title),
  Group = ifelse(
    stroke_flag,
    "Stroke",
    "Control"
  ),
  stringsAsFactors = FALSE
)

metadata$Group <- factor(
  metadata$Group,
  levels = c("Control", "Stroke")
)

rownames(metadata) <- metadata$Sample_ID
# Display group distribution
cat("Clinical Sample Distribution:\n")

group_counts <- table(metadata$Group)

print(group_counts)
# Verify expected sample counts
if (
  group_counts["Control"] != 24 ||
  group_counts["Stroke"] != 39
) {

  stop("Expected 24 Control and 39 Stroke samples.")

}

cat("Clinical metadata verification successful.\n")

## ---- missing-values ----
# Count original missing values
na_count <- sum(is.na(expr_raw))

missing_percentage <- (
  na_count / length(expr_raw)
) * 100

cat(
  "Missing values before filtering:",
  na_count, "\n"
)

cat(
  "Missing percentage:",
  round(missing_percentage, 3), "%\n"
)
# Calculate missing-value proportion for each probe
missing_fraction <- rowMeans(
  is.na(expr_raw)
)
# Identify probes containing at least one valid observation
nonempty_probes <- rowSums(
  !is.na(expr_raw)
) > 0
# Retain probes with no more than 20% missing values
keep_probes <- (
  missing_fraction <= 0.20
) & nonempty_probes

expr_clean <- expr_raw[
  keep_probes,
  ,
  drop = FALSE
]

cat(
  "Probes before filtering:",
  nrow(expr_raw), "\n"
)

cat(
  "Probes retained:",
  nrow(expr_clean), "\n"
)

cat(
  "Probes removed:",
  nrow(expr_raw) - nrow(expr_clean), "\n"
)
# Calculate median expression for each probe
row_medians <- apply(
  expr_clean,
  1,
  median,
  na.rm = TRUE
)
# Impute each missing cell using the median of its probe (row)
missing_cells <- which(is.na(expr_clean), arr.ind = TRUE)
if (nrow(missing_cells)) {
  expr_clean[missing_cells] <- row_medians[missing_cells[, 1]]
}
# Verify missing-value handling
cat(
  "Missing values after imputation:",
  sum(is.na(expr_clean)), "\n"
)

if (anyNA(expr_clean)) {
  stop("Missing values remain after preprocessing.")
}

## ---- expression-inspection ----
# Calculate expression value quantiles
expression_quantiles <- quantile(
  expr_clean,
  probs = c(0, 0.25, 0.50, 0.75, 0.99, 1),
  na.rm = TRUE
)

cat("Expression Value Quantiles:\n")

print(expression_quantiles)

cat(
  "Minimum expression:",
  min(expr_clean), "\n"
)

cat(
  "Maximum expression:",
  max(expr_clean), "\n"
)

cat(
  "Mean expression:",
  mean(expr_clean), "\n"
)
# Retain the GEO-processed expression values
expr_processed <- expr_clean

cat(
  "Existing log2 and quantile-normalized values retained.\n"
)

cat(
  "Additional normalization was not applied.\n"
)

## ---- platform-annotation ----
# Platform annotation file inside the project
gpl_file <- file.path(
  data_dir,
  "GPL6883.annot.gz"
)
# NCBI GEO annotation download URL
gpl_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/",
  "GPL6nnn/GPL6883/annot/GPL6883.annot.gz"
)
# Download annotation only if missing
if (!file.exists(gpl_file)) {

  cat("Downloading GPL6883 annotation...\n")

  download.file(
    url = gpl_url,
    destfile = gpl_file,
    mode = "wb",
    method = "libcurl"
  )

}

if (
  !file.exists(gpl_file) ||
  file.info(gpl_file)$size < 1000
) {

  stop("GPL6883 annotation file is missing or incomplete.")

}

cat("GPL6883 annotation file located successfully.\n")
# Read the annotation header
con <- gzfile(
  gpl_file,
  open = "rt"
)

header_lines <- readLines(
  con,
  n = 500
)

close(con)

header_index <- grep(
  "^#?ID\\t",
  header_lines
)[1]

if (is.na(header_index)) {
  stop("GPL6883 annotation header was not found.")
}
# Import the complete annotation table
con <- gzfile(
  gpl_file,
  open = "rt"
)

gpl_table <- read.delim(
  con,
  header = TRUE,
  sep = "\t",
  skip = header_index - 1,
  quote = "",
  comment.char = "",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

close(con)
# Clean column names
colnames(gpl_table) <- sub(
  "^#",
  "",
  colnames(gpl_table)
)

cat(
  "Platform annotation dimensions:",
  nrow(gpl_table), "x", ncol(gpl_table), "\n"
)

cat("Available annotation columns:\n")

print(colnames(gpl_table))

## ---- gene-annotation ----
# Identify annotation columns
id_col <- "ID"

symbol_col <- "Gene symbol"
# Verify columns
stopifnot(
  id_col %in% colnames(gpl_table),
  symbol_col %in% colnames(gpl_table)
)
# Match expression probes with GPL6883 annotation
probe_index <- match(
  rownames(expr_processed),
  as.character(gpl_table[[id_col]])
)
# Extract gene symbols
probe_symbols <- as.character(
  gpl_table[[symbol_col]][probe_index]
)

probe_symbols <- trimws(probe_symbols)

cat(
  "Successfully matched probes:",
  sum(!is.na(probe_symbols)), "\n"
)

cat(
  "Unmatched probes:",
  sum(is.na(probe_symbols)), "\n"
)
# Handle multiple gene symbols assigned to probes
gene_symbols <- sub(
  "\\s*///.*$",
  "",
  probe_symbols
)

gene_symbols <- trimws(gene_symbols)
# Identify valid gene annotations
valid_probes <- !is.na(gene_symbols) &
  gene_symbols != "" &
  !toupper(gene_symbols) %in% c(
    "---", "NA", "NONE", "NULL"
  )
# Retain successfully annotated probes
expr_annotated <- expr_processed[
  valid_probes,
  ,
  drop = FALSE
]

gene_symbols_valid <- gene_symbols[
  valid_probes
]

cat(
  "Valid annotated probes:",
  nrow(expr_annotated), "\n"
)
# Summarise duplicate probes using mean expression
expr_gene <- limma::avereps(
  expr_annotated,
  ID = gene_symbols_valid
)

cat(
  "Final Cleaned Gene Matrix (Unique Genes x Samples):",
  nrow(expr_gene), "x", ncol(expr_gene), "\n"
)

## ---- final-quality-checks ----
# Match metadata to expression matrix columns
metadata <- metadata[
  match(
    colnames(expr_gene),
    metadata$Sample_ID
  ),
  ,
  drop = FALSE
]
# Verify sample matching
if (anyNA(metadata$Sample_ID)) {
  stop("Sample matching failed.")
}
# Check final matrix
cat(
  "Final missing values:",
  sum(is.na(expr_gene)), "\n"
)

cat(
  "Duplicate gene symbols:",
  sum(duplicated(rownames(expr_gene))), "\n"
)

cat(
  "Final gene count:",
  nrow(expr_gene), "\n"
)

cat(
  "Final sample count:",
  ncol(expr_gene), "\n"
)

cat("Final sample distribution:\n")

print(table(metadata$Group))
# Validate final dataset
stopifnot(
  is.matrix(expr_gene),
  is.numeric(expr_gene),
  !anyNA(expr_gene),
  !anyDuplicated(rownames(expr_gene)),
  ncol(expr_gene) == 63,
  sum(metadata$Group == "Control") == 24,
  sum(metadata$Group == "Stroke") == 39
)

cat("All quality checks completed successfully.\n")

## ---- save-cleaned-data ----
# Export clinical metadata
write.csv(
  metadata,
  file.path(
    data_dir,
    "GSE16561_clinical_metadata.csv"
  ),
  row.names = FALSE
)
# Export annotated gene expression matrix
write.csv(
  expr_gene,
  file.path(
    data_dir,
    "GSE16561_gene_expression_matrix.csv"
  ),
  quote = FALSE
)
# Save combined cleaned dataset
saveRDS(
  list(
    expr = expr_gene,
    metadata = metadata
  ),
  file.path(
    data_dir,
    "GSE16561_cleaned_data.rds"
  )
)
# Create quality control summary
qc_summary <- data.frame(
  Metric = c(
    "Original probes",
    "Original samples",
    "Original missing values",
    "Missing value percentage",
    "Probes retained after filtering",
    "Annotated probes",
    "Unique genes",
    "Control samples",
    "Stroke samples",
    "Final missing values"
  ),
  Value = c(
    nrow(expr_raw),
    ncol(expr_raw),
    na_count,
    round(missing_percentage, 3),
    nrow(expr_clean),
    nrow(expr_annotated),
    nrow(expr_gene),
    24,
    39,
    sum(is.na(expr_gene))
  )
)

write.csv(
  qc_summary,
  file.path(
    data_dir,
    "GSE16561_QC_summary.csv"
  ),
  row.names = FALSE
)

cat("Cleaned data exported successfully.\n")

cat("Generated files:\n")

print(
  list.files(
    data_dir,
    full.names = FALSE
  )
)