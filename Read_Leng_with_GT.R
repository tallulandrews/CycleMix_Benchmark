data <- read.delim("/home/tandrew6/scratch/CycleMix_Manuscript/GSE64016_H1andFUCCI_normalized_EC.csv.gz", sep=",")

counts <- data
rownames(counts) <- counts[,1]
counts <- counts[,-1]

meta <- sapply(strsplit(colnames(counts), "_"), function(x){x[[1]]})
meta[meta == "H1"] <- "Unknown"

Leng_counts <- counts
Leng_meta <- meta


