library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query)

datatable(as.data.frame(colData(data)),
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
