require(Seurat)

run_Seurat <- function(counts, gene_table) {
	counts <- counts[rownames(counts)!="",]
	s.genes <- gene_table[(gene_table[,2] == "S" | gene_table[,2] == "G1S") & gene_table[,3] == 1,1]
	g2m.genes <- gene_table[(gene_table[,2] == "G2M" | gene_table[,2] == "M" | gene_table[,2] == "G2") & gene_table[,3] == 1,1]
	seur <- CreateSeuratObject(counts)
	seur <- NormalizeData(seur)
	print("Seurat is classifying cells")
	#start.time <- Sys.time()
	elapsed.time <- system.time(
	seur <- tryCatch({out = CellCycleScoring(seur, s.features=s.genes, g2m.features=g2m.genes)},
	error = function(cond) {
            message(paste("few genes detected rerunning with 'ctrl=10'"))
            message("Here's the original error message:")
            message(conditionMessage(cond))
	    out <- CellCycleScoring(seur, s.features=s.genes, g2m.features=g2m.genes, ctrl=10)
	    return(out)
        }
        )
	)
	#end.time <- Sys.time()
	#elapsed.time=end.time - start.time
	print(paste("Elapsed time:", elapsed.time))
	return(list(phase=seur@meta.data$Phase, time=elapsed.time[3]))
}

run_Seurat_scores <- function(counts, gene_table) {
        s.genes <- gene_table[(gene_table[,2] == "S" | gene_table[,2] == "G1S") & gene_table[,3] == 1,1]
        g2m.genes <- gene_table[(gene_table[,2] == "G2M" | gene_table[,2] == "M" | gene_table[,2] == "G2") & gene_table[,3] == 1,1]
        seur <- CreateSeuratObject(counts)
        seur <- NormalizeData(seur)
	start.time <- Sys.time()
	tryCatch({out = CellCycleScoring(seur, s.features=s.genes, g2m.features=g2m.genes); return(out)},
	error = function(cond) {
            message(paste("few genes detected rerunning with 'ctrl=10'"))
            message("Here's the original error message:")
            message(conditionMessage(cond))
	    out <- CellCycleScoring(seur, s.features=s.genes, g2m.features=g2m.genes, ctrl=10)
	    return(out)
        }
        )
        seur <- CellCycleScoring(seur, s.features=s.genes, g2m.features=g2m.genes)
	end.time <- Sys.time()
	print(paste("Elapsed time:", end.time - start.time))
	out <- as.matrix(seur@meta.data[,c("S.Score", "G2M.Score")])
	rownames(out) <- seur@meta.data[,"Phase"]
        return(seur@meta.data$Phase)
}

