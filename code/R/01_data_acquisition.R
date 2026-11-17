# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- locate-data ----
options(timeout = 600)
# Possible dataset files inside the project data folder
possible_paths <- file.path(
  data_dir,
  c(
    "GSE16561_series_matrix.txt",
    "GSE16561_series_matrix.txt.gz"
  )
)

available_paths <- possible_paths[
  file.exists(possible_paths)
]
# Download from NCBI GEO only if the file is missing
if (length(available_paths) == 0) {

  cat("Dataset not found locally. Downloading from NCBI GEO...\n")

  geo_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    "GSE16nnn/GSE16561/matrix/",
    "GSE16561_series_matrix.txt.gz"
  )

  downloaded_file <- file.path(
    data_dir,
    "GSE16561_series_matrix.txt.gz"
  )

  download.file(
    url = geo_url,
    destfile = downloaded_file,
    mode = "wb",
    method = "libcurl"
  )

  if (!file.exists(downloaded_file)) {
    stop("Dataset download failed.")
  }

  file_to_use <- downloaded_file

} else {

  file_to_use <- available_paths[1]

}

cat("Dataset found at:", file_to_use, "\n")

## ---- import-data ----
# Cache the imported ExpressionSet inside data/
gse_rds <- file.path(
  data_dir,
  "GSE16561_raw_eSet.rds"
)
# Import the dataset or load the existing cache
if (file.exists(gse_rds)) {

  cat("Loading cached GEO dataset...\n")

  gse <- readRDS(gse_rds)

} else {

  cat("Importing GSE16561 Series Matrix...\n")

  gse_object <- GEOquery::getGEO(
    filename = file_to_use,
    getGPL = FALSE
  )

  if (inherits(gse_object, "ExpressionSet")) {

    gse <- gse_object

  } else if (is.list(gse_object) &&
             length(gse_object) > 0) {

    gse <- gse_object[[1]]

  } else {

    stop("Unexpected GEO dataset format.")

  }

  saveRDS(
    gse,
    file = gse_rds
  )

}
# Extract expression and phenotype data
expr_raw <- Biobase::exprs(gse)

pdata_raw <- Biobase::pData(gse)
# Ensure a numeric expression matrix
expr_raw <- as.matrix(expr_raw)

storage.mode(expr_raw) <- "numeric"

cat(
  "Raw Expression Matrix Dimensions (Probes x Samples):",
  nrow(expr_raw), "x", ncol(expr_raw), "\n"
)

cat(
  "Number of clinical samples:",
  nrow(pdata_raw), "\n"
)