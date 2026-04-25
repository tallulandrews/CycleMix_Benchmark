data <- read.delim("/home/tandrew6/scratch/CycleMix_Manuscript/McDavid_DataSetS2.csv.gz", sep="\t")

cellID <- sub("\\..+$", "", data$cellID)
geneID <- unique(data$primerid)

mat <- matrix(0, nrow=length(geneID), ncol=length(unique(cellID)))
matn <- matrix(0, nrow=length(geneID), ncol=length(unique(cellID)))
rownames(mat) <- geneID
rownames(matn) <- geneID
colnames(mat) <- unique(cellID)
colnames(matn) <- unique(cellID)
for (i in 1:nrow(data)) {
	cell <- cellID[i]
	gene <- data$primerid[i]
	expr <- 2^data$lCount[i]
	norm <- data$nCount[i]
	mat[gene, cell] <- expr
	matn[gene, cell] <- norm
}
	
meta <- data$cycle[match(colnames(mat), cellID)]

MacDavid_counts <- mat
MacDavid_lognorm <- matn
MacDavid_meta <- meta


