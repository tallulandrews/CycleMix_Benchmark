files <- Sys.glob("/home/tandrew6/scratch/CycleMix_Manuscript/E-MTAB-2805/*_counts.txt")

source("~/scripts/My_R_Scripts.R")

counts <- NULL
meta <- c()
for (f in files) {
	mat <- read.table(f, header=T)
	mat <- mat[!is.na(mat[,3]),]
	mat <- remove_duplicate_rows(mat[,-1*(1:4)], mat[,3], method="max")
	if (is.null(counts)){
		counts <- mat
	} else {
		common_genes <- intersect(rownames(mat), rownames(counts))
		mat <- mat[match(common_genes, rownames(mat)),]
		counts <- counts[match(common_genes, rownames(counts)),]
		counts <- cbind(counts, mat)
	}
	meta <- c(meta, sapply(strsplit(colnames(mat), "_"), function(x){x[[1]]}))
}
Buttner_counts <- counts
Buttner_meta <- meta

