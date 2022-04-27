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


query_maf <- GDCquery(project = "TCGA-COAD",
                data.category = "Simple Nucleotide Variation",
                data.type = "Masked Somatic Mutation",
                legacy = F)
GDCdownload(query_maf, directory = "GDCdata/")
maf = GDCprepare(query_maf, directory = "GDCdata/")
# Only first 50 to make render faster
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

