source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")


mat <- readRDS("RAW.UMI.counts.BC.cell.lines.rds")
gene_symbols <- General_Map(rownames(mat), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
source("~/scripts/My_R_Scripts.R")

mat <- remove_duplicate_rows(mat, gene_symbols, method="max")

BC_celllines_counts <- mat
BC_cellines_line <- unlist(sapply(strsplit(colnames(mat), "_"), function(x){x[[1]]}))


