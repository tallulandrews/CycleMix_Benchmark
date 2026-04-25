source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")


mat <- t(read.delim("GSE146773_Counts.csv.gz", sep=",", header=T))
tmp <- mat[1,]
mat <- mat[-1,]
mat <- t(apply(mat, 1, as.numeric))
colnames(mat) <- tmp

gene_symbols <- General_Map(rownames(mat), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
source("~/scripts/My_R_Scripts.R")

mat <- remove_duplicate_rows(mat, gene_symbols, method="max")

# timepoint 2 classification: https://www.nature.com/articles/s41586-021-03232-9#Fig1
# G1 = 0->15
# S&G2m = 15->25

meta <- read.delim("GSE146773_fucci_coords.csv.gz", sep=",")
Madhessian_meta <- rep("G0", ncol(mat))
Madhessian_meta[(meta$fucci_time_hrs > 5)] <- "G1"
Madhessian_meta[(meta$fucci_time_hrs > 15)] <- "S"
Madhessian_meta[(meta$fucci_time_hrs > 21)] <- "G2M"

Madhessian_counts <- Matrix(as.matrix(mat))


