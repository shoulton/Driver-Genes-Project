library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)


query <- GDCquery(
  project = c("TCGA-GBM",'TCGA-COAD'),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)

patients <- GDCprepare(query)

datatable(as.data.frame(colData(patients)),
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

#data.type = "Masked Somatic Mutation",

query_maf <- GDCquery(project = "TCGA-COAD",
                data.category = "Simple Nucleotide Variation",
                legacy = F)
GDCdownload(query_maf, directory = "GDCdata/")
maf = GDCprepare(query_maf, directory = "GDCdata/")
# Only first 50 to make render faster

driver_genes <- c("MGA", "GTF2I", "MSH6", "RRAS2", "RPL22", "PIK3R2",
                  "TAF1", "SOS1", "PDS5B", "ZFHX3", "ZMYM2", "ATF7IP",
                  "INPPL1", "ATR", "ARID5B", "SOX17", "MYCN", "CTNND1",
                  "FOXA2", "DICER1", "SCAF4", "SIN3A", "KMT2B", "ACVR1",
                  "CCND1", "CHD3", "AMER1")

only_gene <- filter(maf, (Hugo_Symbol %in% driver_genes))
only_gene <- transform(only_gene, patient_ID = substr(only_gene$Tumor_Sample_Barcode,1,12))

driver_gene_sample_barcodes = only_gene$patient_ID

patients_df <- as.data.frame(colData(patients))

patients_with_driver_genes = filter(as.data.frame(colData(patients)), patient %in% driver_gene_sample_barcodes)

datatable(only_gene,
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(patients_with_driver_genes,
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

Survival_Function = survfit(Surv(patients_with_driver_genes$days_to_death,
                                 patients_with_driver_genes$vital_status == "Dead")~1)

Survival_Function

plot(Survival_Function)

`%notin%` <- Negate(`%in%`)

patients_without_driver_genes = filter(as.data.frame(colData(patients)), patient %notin% driver_gene_sample_barcodes)

datatable(patients_without_driver_genes,
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

Survival_Function = survfit(Surv(patients_without_driver_genes$days_to_death,
                                 patients_without_driver_genes$vital_status == "Dead")~1)

Survival_Function

plot(Survival_Function)

