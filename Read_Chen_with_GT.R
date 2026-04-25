files <- Sys.glob("/home/tandrew6/scratch/CycleMix_Manuscript/GSE107644_RAW/*_UMIcount.txt.gz")

source("~/scripts/My_R_Scripts.R")

anno <- read.table("/home/tandrew6/scratch/CycleMix_Manuscript/GSE107644_barcode_information.txt.gz")

counts <- NULL
meta <- c()
for (f in files) {
	mat <- read.table(f, header=T)
	rownames(mat) <- mat[,1]
	mat <- mat[,-1]
	if (is.null(counts)){
		counts <- mat
	} else {
		common_genes <- intersect(rownames(mat), rownames(counts))
		mat <- mat[match(common_genes, rownames(mat)),]
		counts <- counts[match(common_genes, rownames(counts)),]
		counts <- cbind(counts, mat)
	}
	#meta <- c(meta, sapply(strsplit(colnames(mat), "_"), function(x){x[[1]]}))
	meta <- c(meta, colnames(mat))
}

phase <- rep("Unknown", length(meta));
phase[grep("TypeBG2M", meta)] <- "G2M"
phase[grep("_G1_", meta)] <- "G1"
phase[grep("_TypeBS_", meta)] <- "S"


Chen_counts <- counts
Chen_meta <- phase

