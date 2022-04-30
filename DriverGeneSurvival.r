library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)

#Define notin for use later
`%notin%` <- Negate(`%in%`)


# Query Patient database, in this case COAD
query <- GDCquery(
  project = c('TCGA-COAD'),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)

patients_COAD <- GDCprepare(query)

#Query Maf database, specifically COAD
query_maf_COAD <- GDCquery(project = "TCGA-COAD",
                data.category = "Simple Nucleotide Variation",
                data.type = "Masked Somatic Mutation",
                legacy = F)

GDCdownload(query_maf_COAD, directory = "GDCdata/")
maf_COAD = GDCprepare(query_maf_COAD, directory = "GDCdata/")

#driver_genes_COAD <- c("MGA", "GTF2I", "MSH6", "RRAS2", "RPL22", "PIK3R2",
#                  "TAF1", "SOS1", "PDS5B", "ZFHX3", "ZMYM2", "ATF7IP",
#                  "INPPL1", "ATR", "ARID5B", "SOX17", "MYCN", "CTNND1",
#                  "FOXA2", "DICER1", "SCAF4", "SIN3A", "KMT2B", "ACVR1",
#                  "CCND1", "CHD3", "AMER1")

#List the Driver Genes for this specific cancer type
COAD_driver_genes <- c("AMER1", "TCFZL2", "ZFP36L2", "PCBP1", "SOX9", "TGIF1")

#Find maf files with these driver genes
COAD_gene <- filter(maf_COAD, (Hugo_Symbol %in% COAD_driver_genes))

#Generate patient ID based on tumor sample barcode
COAD_gene <- transform(COAD_gene, patient_ID = substr(COAD_gene$Tumor_Sample_Barcode,1,12))

#Find patient IDs where mutations are present
driver_gene_sample_barcodes = COAD_gene$patient_ID

#Select patients that have the driver gene mutation based off of tumor sample barcodes as found above
COAD_patients_with_driver_genes = filter(as.data.frame(colData(patients_COAD)), patient %in% driver_gene_sample_barcodes)

#Fit a survival function for patients with these specific driver gene mutations and plot
Survival_Function = survfit(Surv(COAD_patients_with_driver_genes$days_to_last_follow_up,
                                 COAD_patients_with_driver_genes$vital_status == "Alive")~1)

Survival_Function

svg(filename = "COAD with driver genes survival.svg")
plot(Survival_Function)
dev.off()

#Fit a survival model for patients without these specific driver gene mutations and plot
COAD_patients_without_driver_genes = filter(as.data.frame(colData(patients_COAD)), patient %notin% driver_gene_sample_barcodes)

Survival_Function = survfit(Surv(COAD_patients_without_driver_genes$days_to_last_follow_up,
                                 COAD_patients_without_driver_genes$vital_status == "Alive")~1)
Survival_Function

svg(filename = "COAD without driver genes survival.svg")
plot(Survival_Function)
dev.off()



### Repeat with LUAD

# Query Patient database, in this case COAD
query <- GDCquery(
  project = c('TCGA-LUAD'),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)

patients_LUAD <- GDCprepare(query)

#Query Maf database, specifically COAD
query_maf_LUAD <- GDCquery(project = "TCGA-LUAD",
                           data.category = "Simple Nucleotide Variation",
                           data.type = "Masked Somatic Mutation",
                           legacy = F)

GDCdownload(query_maf_LUAD, directory = "GDCdata/")
maf_LUAD = GDCprepare(query_maf_LUAD, directory = "GDCdata/")

#List the Driver Genes for this specific cancer type
LUAD_driver_genes <- c("MGA", "RIT1", "TP53")

#Find maf files with these driver genes
LUAD_gene <- filter(maf_LUAD, (Hugo_Symbol %in% LUAD_driver_genes))

#Generate patient ID based on tumor sample barcode
LUAD_gene <- transform(LUAD_gene, patient_ID = substr(LUAD_gene$Tumor_Sample_Barcode,1,12))

#Find patient IDs where mutations are present
driver_gene_sample_barcodes = LUAD_gene$patient_ID

LUAD_patients_df <- as.data.frame(colData(patients_LUAD))

#Select patients that have the driver gene mutation based off of tumor sample barcodes as found above
LUAD_patients_with_driver_genes = filter(as.data.frame(colData(patients_LUAD)), patient %in% driver_gene_sample_barcodes)


#Fit a survival function for patients with these specific driver gene mutations and plot
Survival_Function = survfit(Surv(LUAD_patients_with_driver_genes$days_to_last_follow_up,
                                 LUAD_patients_with_driver_genes$vital_status == "Alive")~1)

Survival_Function

svg(filename = "LUAD with driver genes survival.svg")
plot(Survival_Function)
dev.off()

#Fit a survival model for patients without these specific driver gene mutations and plot
LUAD_patients_without_driver_genes = filter(as.data.frame(colData(patients_LUAD)), patient %notin% driver_gene_sample_barcodes)

Survival_Function = survfit(Surv(LUAD_patients_without_driver_genes$days_to_last_follow_up,
                                 LUAD_patients_without_driver_genes$vital_status == "Alive")~1)
Survival_Function

svg(filename = "LUAD without driver genes survival.svg")
plot(Survival_Function)
dev.off()



### Repeat with UCEC

# Query Patient database, in this case COAD
query <- GDCquery(
  project = c('TCGA-UCEC'),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)

patients_UCEC <- GDCprepare(query)

#Query Maf database, specifically COAD
query_maf_UCEC <- GDCquery(project = "TCGA-UCEC",
                           data.category = "Simple Nucleotide Variation",
                           data.type = "Masked Somatic Mutation",
                           legacy = F)

GDCdownload(query_maf_UCEC, directory = "GDCdata/")
maf_UCEC = GDCprepare(query_maf_UCEC, directory = "GDCdata/")

#List the Driver Genes for this specific cancer type
UCEC_driver_genes <- c("MSH6", "RRAS2", "RPL22", "PIK3R2",
                       "TAF1", "SOS1", "PDS5B", "ZFHX3", "ZMYM2", "ATF7IP",
                       "INPPL1", "ATR", "ARID5B", "SOX17", "MYCN", "CTNND1",
                       "FOXA2", "DICER1", "SCAF4", "SIN3A", "KMT2B", "ACVR1",
                       "CCND1", "CHD3")

#Find maf files with these driver genes
UCEC_gene <- filter(maf_UCEC, (Hugo_Symbol %in% UCEC_driver_genes))

#Generate patient ID based on tumor sample barcode
UCEC_gene <- transform(UCEC_gene, patient_ID = substr(UCEC_gene$Tumor_Sample_Barcode,1,12))

#Find patient IDs where mutations are present
driver_gene_sample_barcodes = UCEC_gene$patient_ID

UCEC_patients_df <- as.data.frame(colData(patients_UCEC))

#Select patients that have the driver gene mutation based off of tumor sample barcodes as found above
UCEC_patients_with_driver_genes = filter(as.data.frame(colData(patients_UCEC)), patient %in% driver_gene_sample_barcodes)


#Fit a survival function for patients with these specific driver gene mutations and plot
Survival_Function = survfit(Surv(UCEC_patients_with_driver_genes$days_to_last_follow_up,
                                 UCEC_patients_with_driver_genes$vital_status == "Alive")~1)

Survival_Function

svg(filename = "UCEC with driver genes survival.svg")
plot(Survival_Function)
dev.off()

#Fit a survival model for patients without these specific driver gene mutations and plot
UCEC_patients_without_driver_genes = filter(as.data.frame(colData(patients_UCEC)), patient %notin% driver_gene_sample_barcodes)

Survival_Function = survfit(Surv(UCEC_patients_without_driver_genes$days_to_last_follow_up,
                                 UCEC_patients_without_driver_genes$vital_status == "Alive")~1)
Survival_Function

svg(filename = "UCEC without driver genes survival.svg")
plot(Survival_Function)
dev.off()
