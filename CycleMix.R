require(CycleMix)
require(SingleCellExperiment)
require(scater)
require(scran)
source("/home/tandrew6/My_R_Packages/CycleMix/R/Analysis.R")

run_CycleMix <- function(counts, gene_table) {
	SCE <- SingleCellExperiment(assays=list(counts=counts), rowData=data.frame(feature_symbol=rownames(counts)))
	logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(counts)/median(colSums(counts)))
	print("Cyclemix is classifying cells")
	#start.time <- Sys.time()
	elapsed.time <- system.time(out <- classifyCells(SCE, gene_table))[3]
	#end.time <- Sys.time()
	#elapsed.time=end.time - start.time
	print(paste("Elapsed time:", elapsed.time))
	return(list(phase=out$phase, time=elapsed.time))
}

run_CycleMix_scores <- function(counts, gene_table) {
        SCE <- SingleCellExperiment(assays=list(counts=counts), rowData=data.frame(feature_symbol=rownames(counts)))
        logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(counts)/median(colSums(counts)))
	start.time <- Sys.time()
        out <- classifyCells(SCE, gene_table)
	end.time <- Sys.time()
	print(paste("Elapsed time:", end.time - start.time))
	rownames(out$scores) <-out$phase
        return(out$scores)
}

