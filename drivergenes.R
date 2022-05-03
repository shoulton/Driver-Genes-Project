library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)

# LUAD
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  ## platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)

GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query)

datatable(as.data.frame(colData(data)),
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# test
clin <- GDCquery_clinic("TCGA-LUAD","clinical")
my_data <- data.frame(clin, stringsAsFactors = FALSE, check.names = FALSE)
colnames(my_data)[colnames(my_data)=="submitter_id"] <-  "Tumor_Sample_Barcode"
colnames(my_data)[colnames(my_data)=="vital_status"] <-  "Overall_Survival_Status"
my_data$Overall_Survival_Status <- ifelse(my_data$Overall_Survival_Status == 'Alive', 1, 0)

laml <- read.maf(maf, clinicalData = my_data, isTCGA = TRUE)

# original
# clinical <- as.data.frame(colData(data))
# clinical$Tumor_Sample_Barcode <- clinical$barcode
# clinical <- subset(clinical, !is.na(days_to_last_follow_up))
# clinical$days_to_last_follow_up <- as.numeric(clinical$days_to_last_follow_up)
# clinical$vital_status <- ifelse(clinical$vital_status == 'Alive', 1, 0)

query_maf <- GDCquery(project = "TCGA-LUAD",
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     legacy = F)
GDCdownload(query_maf, directory = "GDCdata/")
maf <- GDCprepare(query_maf, directory = "GDCdata/")
maf_file <- read.maf(maf = maf, clinicalData = clinical)

# Only first 50 to make render faster
datatable(maf[,10],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


# visualize maf data

plotmafSummary(laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(laml, top = 20)

# survival
#Survival analysis based on grouping of RIT1 mutation status
mafSurvival(maf = laml, genes = c('TP53', 'RIT1', 'MGA', 'PIK3CA', 'KRAS'), time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)


# COAD
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)

GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query)

datatable(as.data.frame(colData(data)),
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# test
clin <- GDCquery_clinic("TCGA-COAD","clinical")
my_data <- data.frame(clin, stringsAsFactors = FALSE, check.names = FALSE)
colnames(my_data)[colnames(my_data)=="submitter_id"] <-  "Tumor_Sample_Barcode"
colnames(my_data)[colnames(my_data)=="vital_status"] <-  "Overall_Survival_Status"
my_data$Overall_Survival_Status <- ifelse(my_data$Overall_Survival_Status == 'Alive', 1, 0)

laml <- read.maf(maf, clinicalData = my_data, isTCGA = TRUE)

# original
# clinical <- as.data.frame(colData(data))
# clinical$Tumor_Sample_Barcode <- clinical$barcode
# clinical <- subset(clinical, !is.na(days_to_last_follow_up))
# clinical$days_to_last_follow_up <- as.numeric(clinical$days_to_last_follow_up)
# clinical$vital_status <- ifelse(clinical$vital_status == 'Alive', 1, 0)

query_maf <- GDCquery(project = "TCGA-COAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf, directory = "GDCdata/")
maf <- GDCprepare(query_maf, directory = "GDCdata/")
# maf_file <- read.maf(maf = maf, clinicalData = clinical)

# Only first 50 to make render faster
datatable(maf[,10],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


# visualize maf data

plotmafSummary(laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(laml, top = 10)

# survival
#Survival analysis based on grouping of RIT1 mutation status
# mafSurvival(maf = laml, genes = c('TP53', 'RIT1', 'MGA', 'PIK3CA', 'KRAS'), time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvival(maf = laml, genes = c('TCFZL2', 'AMER1', 'ZFP36L2', 'PCBP1', 'SOX9', 'TGIF1'), time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)

# ('AMER1', 'TCFZL2', 'ZFP36L2', 'PCBP1', 'SOX9', 'TGIF1')

# MESO
query <- GDCquery(
  project = "TCGA-MESO",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)

GDCdownload(query, method = "api", files.per.chunk = 10)

data <- GDCprepare(query)

datatable(as.data.frame(colData(data)),
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# test
clin <- GDCquery_clinic("TCGA-MESO","clinical")
my_data <- data.frame(clin, stringsAsFactors = FALSE, check.names = FALSE)
colnames(my_data)[colnames(my_data)=="submitter_id"] <-  "Tumor_Sample_Barcode"
colnames(my_data)[colnames(my_data)=="vital_status"] <-  "Overall_Survival_Status"
my_data$Overall_Survival_Status <- ifelse(my_data$Overall_Survival_Status == 'Alive', 1, 0)


# original
# clinical <- as.data.frame(colData(data))
# clinical$Tumor_Sample_Barcode <- clinical$barcode
# clinical <- subset(clinical, !is.na(days_to_last_follow_up))
# clinical$days_to_last_follow_up <- as.numeric(clinical$days_to_last_follow_up)
# clinical$vital_status <- ifelse(clinical$vital_status == 'Alive', 1, 0)

query_maf <- GDCquery(project = "TCGA-MESO",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf, directory = "GDCdata/")
maf <- GDCprepare(query_maf, directory = "GDCdata/")
# maf_file <- read.maf(maf = maf, clinicalData = clinical)

# Only first 50 to make render faster
datatable(maf[,10],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

laml <- read.maf(maf, clinicalData = my_data, isTCGA = TRUE)

# visualize maf data

plotmafSummary(laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(laml, top = 10)

# survival
#Survival analysis based on grouping of RIT1 mutation status
# mafSurvival(maf = laml, genes = c('TP53', 'RIT1', 'MGA', 'PIK3CA', 'KRAS'), time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvival(maf = laml, genes = c( 'TP53','LATS2','BAP1',), time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', isTCGA = TRUE)

# ('AMER1', 'TCFZL2', 'ZFP36L2', 'PCBP1', 'SOX9', 'TGIF1')




