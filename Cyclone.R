require("scran")

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")

convert_gene_names <- function(pairs_obj, species=c("Mmus", "Hsap")){
	for (i in 1:length(pairs_obj)) {
		ensg_pairs <- pairs_obj[[i]]
		symbol_pairs <- ensg_pairs
		symbol_pairs[,1] <- General_Map(symbol_pairs[,1], in.org=species[1], out.org=species[1], in.name="ensg", out.name="symbol")
		symbol_pairs[,2] <- General_Map(symbol_pairs[,2], in.org=species[1], out.org=species[1], in.name="ensg", out.name="symbol")
		symbol_pairs <- symbol_pairs[symbol_pairs[,1] != symbol_pairs[,2] & symbol_pairs[,1] != "" & symbol_pairs[,2] != "",]
		pairs_obj[[i]] <- symbol_pairs
	}
	return(pairs_obj)
}

mm.pairs <- convert_gene_names(mm.pairs, species="Mmus")
hs.pairs <- convert_gene_names(hs.pairs, species="Hsap")

run_cyclone <- function(counts, species=c("mmus", "hsap")) {
	#start.time <- Sys.time()
	elapsed.time = -1
	print("Cyclone is classifying cells")
	if (species == "mmus") {
		elapsed.time <- system.time(assignments <- cyclone(counts, mm.pairs))[3]
	} else if (species == "hsap") {
		elapsed.time <- system.time(assignments <- cyclone(counts, hs.pairs))[3]
	} else {
		stop(paste(species, "is not recognized, please specify either 'mmus' or 'hsap' for mouse or human respectively."))
	}
	#end.time <- Sys.time()
	#elapsed.time=end.time - start.time
	print(paste("Elapsed time:", elapsed.time))
	out <- assignments$phases
	return(list(phase=out, time=elapsed.time))
}

run_cyclone_scores <- function(counts, species=c("mmus", "hsap")) {
	start.time <- Sys.time()
	if (species == "mmus") {
		assignments <- cyclone(counts, mm.pairs)
	} else if (species == "hsap") {
		 assignments <- cyclone(counts, hs.pairs)
	} else {
		stop(paste(species, "is not recognized, please specify either 'mmus' or 'hsap' for mouse or human respectively."))
	}
	end.time <- Sys.time()
	print(paste("Elapsed time:", end.time - start.time))
	out <- as.matrix(assignments$scores)
	rownames(out) <- assignments$phases
	return(out)
}



