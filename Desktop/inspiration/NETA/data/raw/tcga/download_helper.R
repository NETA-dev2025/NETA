
# TCGA Data Download Helper
# Run this separately if you have TCGAbiolinks installed

if (!require("TCGAbiolinks")) {
  cat("Please install TCGAbiolinks from Bioconductor:\n")
  cat("  BiocManager::install(\"TCGAbiolinks\")\n")
} else {
  library(TCGAbiolinks)
  library(tidyverse)

  # Example: Download PCPG (pure neuroendocrine tumors)
  query <- GDCquery(
    project = "TCGA-PCPG",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # GDCdownload(query)
  # data <- GDCprepare(query)

  cat("Use this template to download real TCGA data\n")
}

