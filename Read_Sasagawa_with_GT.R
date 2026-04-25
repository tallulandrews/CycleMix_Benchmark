meta <- read.delim("/home/tandrew6/scripts/CycleMix_Benchmark/GSE42268-GPL13112_series_matrix.txt", "\t", header=F)
meta <- meta[c(2,8,10,11,13),]

meta[5,] <- sub("cell cycle phase: ", "", meta[5,])
meta[4,] <- sub("sample type: ", "", meta[4,])

keep <- meta[4,] == "single cell" | meta[4,] == "Averaged single cell lysis from 12 cells"

meta <- meta[,keep]

files <- Sys.glob("/home/tandrew6/scratch/CycleMix_Manuscript/GSE42268_RAW/*.txt.gz")
files_id <- sapply(strsplit(files, "/"), function(x){x[[7]]})
files_id <- sapply(strsplit(files_id, "_"), function(x){x[[1]]})
files <- files[match(meta[1,], files_id)]

fpkms <- NULL;
for (f in files){
	expr <- read.table(f, header=TRUE)
	if (is.null(fpkms)){
		fpkms <- expr[,2:3]
	} else {
		common_genes <- intersect(fpkms[,1], expr[,2])
		fpkms <- fpkms[match(common_genes, fpkms[,1]),]
		expr <- expr[match(common_genes, expr[,2]),]
		fpkms <- cbind(fpkms, expr[,3])
	}
}
rownames(fpkms) <- fpkms[,1]
fpkms <- fpkms[,-1]
sample_ids <- sub("^.+\\/","", files)
sample_ids <- sub("_.+$","", sample_ids)
colnames(fpkms) <- sample_ids

all_meta <- meta
meta <- meta[,meta[1,] %in% colnames(fpkms)]
fpkms <- fpkms[,colnames(fpkms) %in% meta[1,]]

meta <- meta[,match(colnames(fpkms), meta[1,])]


Sasagawa_meta <- meta[5,]
Sasagawa_fpkms <- fpkms
