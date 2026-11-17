# Modular code chunks for stroke_biomarker_analysis.Rmd.
# These named chunks are loaded with knitr::read_chunk();
# knit the R Markdown entry point rather than source() this file directly.

## ---- load-packages ----
# Install Bioconductor manager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Required Bioconductor packages
bioc_pkgs <- c(
  "GEOquery",
  "Biobase",
  "limma"
)
# Install missing packages
for (pkg in bioc_pkgs) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )

  }

}
# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
})

cat("Required packages loaded successfully.\n")

cat("GEOquery version:",
    as.character(packageVersion("GEOquery")), "\n")

cat("Biobase version:",
    as.character(packageVersion("Biobase")), "\n")

cat("limma version:",
    as.character(packageVersion("limma")), "\n")