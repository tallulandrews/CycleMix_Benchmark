source("~/scripts/My_R_Scripts.R")

obj <- readRDS("/home/tandrew6/scratch/CycleMix_Manuscript/GSE121265_eset-final.rds")
counts <- exprs(obj); rownames(counts) <- fData(obj)$name
GT <- pData(obj)$theta



Hsiao_counts <- counts
Hsiao_meta <- meta

