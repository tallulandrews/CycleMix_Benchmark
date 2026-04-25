# Idea: Characterize proliferating Kupffer cells responsible for self-renewal in liver kupffer, lung alveolar, and splenic macrophages 

require(Seurat)
require(CycleMix)

source("~/My_R_Packages/CycleMix/R/Analysis.R")
require("mclust")
require("Matrix")

#counts <- Read10X("Guilliams/rawData_human/countTable_human/")

counts <- readMM("Guilliams/rawData_human/countTable_human/matrix.mtx.gz")
cells <- read.table("Guilliams/rawData_human/countTable_human/barcodes.tsv.gz")
genes <- read.table("Guilliams/rawData_human/countTable_human/features.tsv.gz")
rownames(counts) <- genes[,1]
colnames(counts) <- cells[,1]

norm <- log2(t(t(counts)/colSums(counts)*median(colSums(counts)))+1)
CycleMix <- classifyCells(norm, HGeneSets$Tirosh)
saveRDS(CycleMix, "Guilliams_Cyclemix.rds")


Andrews <- readRDS("Andrews2024_Meyloid.rds")
CycleMix<- classifyCells(norm, HGeneSets$Tirosh, expr_name="RNA", symbol_column=NULL)
saveRDS(CycleMix, "Andrews_Cyclemix.rds")
