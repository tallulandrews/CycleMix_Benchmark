source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

Calculate_accuracy_covariate <- function(pred, truth, covariate, window=20) {
	# Sync names.
	if ("None" %in% pred) {
		pred <- as.character(pred)
		pred[pred=="None"] <- "G1"
	}
	if ("G1S" %in% pred) {
		pred <- as.character(pred)
		pred[pred=="G1S"] <- "S"
	}
	if ("G2" %in% pred) {
		pred <- as.character(pred)
		pred[pred=="G2"] <- "G2M"
	}
	TP <- pred==truth
	reorder <- order(covariate)
	covariate<- covariate[reorder]
	TP<- TP[reorder]
	cx <- c(0,cumsum(TP))
	n <- window
	rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
	return(cbind(covariate, rsum))
}


Calculate_accuracy_trio <- function(pred, truth) {
	tab <- table(pred, truth)
	if (nrow(tab) < 3) {
		for (i in c("G1", "S", "G2M")) {
			if (!i %in% rownames(tab)) {
				t_ns <- rownames(tab)
				tab <- rbind(tab, rep(0, ncol(tab)));
				rownames(tab) <- c(t_ns, i)
			}
		}
	}
	if ("None" %in% rownames(tab)) {
		if ("G1" %in% rownames(tab)) {
			tab["G1",] <- tab["G1",] + tab["None",]
		} else {
			tmp <- rownames(tab)
			tmp[tmp=="None"] <- "G1"
			rownames(tab) <- tmp
		}
		tab <- tab[rownames(tab) != "None",]
	}
	if ("G1S" %in% rownames(tab)) {
		tmp <- rownames(tab)
		tmp[tmp=="G1S"] <- "S"
		rownames(tab) <- tmp
	}
	if ("G2" %in% colnames(tab)) {
		tmp <- colnames(tab)
		tmp[tmp=="G2"] <- "G2M"
		colnames(tab) <- tmp
	}
	if (nrow(tab) < 3) {
		for (i in c("G1", "S", "G2M")) {
			if (!i %in% rownames(tab)) {
				t_ns <- rownames(tab)
				tab <- rbind(tab, rep(0, ncol(tab)));
				rownames(tab) <- c(t_ns, i)
			}
		}
	}
	tab <- tab[,c("G1", "G2M", "S")]
	TPs <- c(tab["G1", "G1"], tab["G2M", "G2M"], tab["S","S"])
	FPs <- c(tab["G1", "G2M"]+tab["G1","S"], 
		 tab["G2M", "G1"]+tab["G2M","S"], 
		 tab["S","G1"]+tab["S","G2M"])
	FNs <- c(tab["G2M", "G1"]+tab["S","G1"],
                 tab["G1", "G2M"]+tab["S","G2M"],
                 tab["G1","S"]+tab["G2M","S"])
	names(TPs) <- names(FPs) <- names(FNs) <- c("G1", "G2M", "S")
	total <- sum(tab)
	#return(list(TPs=TPs, FPs=FPs, FNs=FNs, total=total))
	return(sum(TPs)/total)
}
set.seed(3891)

#test <- readRDS("/home/tandrew6/projects/def-tandrew6/tandrew6/External_Data/Merixtell/Ma_Malignant.rds")
#ngenes <- colSums(test@assays$RNA@counts > 0)
#lib.sizes <- colSums(test@assays$RNA@counts)
test1 <-  readRDS("Wang_ileum.rds")
test2 <-  readRDS("Wang_colon.rds")
test3 <-  readRDS("Gu_lymph.rds")

ngenes <- c(colSums(test1@assays$RNA@counts > 0), colSums(test2@assays$RNA@counts > 0), colSums(test3@assays$RNA@counts > 0))
lib.sizes <- c(colSums(test1@assays$RNA@counts), colSums(test1@assays$RNA@counts), colSums(test1@assays$RNA@counts))
resample <- function(count_mat, ncells=4000, lib.n=lib.sizes) {
	require(dplyr)
	prop_mat <- t(t(count_mat)/colSums(count_mat))
	#cumsum_mat <- apply(count_mat,2,cumsum)
	#cumprop_mat <- t(t(cumsum_mat)/colSums(count_mat))
	cell_id_expr = sample(colnames(count_mat), size=ncells, replace=TRUE)
	cell_id_lib = sample(1:length(lib.n), size=ncells, replace=TRUE)
	libs <- lib.sizes[cell_id_lib]
	#genes <- ngenes[cell_id_lib]
	generate_cell <- function(x, gene_names) {
		gene_read_id <- sample(rownames(prop_mat), libs[x], prob=prop_mat[,cell_id_expr[x]], replace=TRUE)
		this_counts <- table(gene_read_id)
		reorder_counts <- this_counts[match(gene_names, names(this_counts))]
		reorder_counts[is.na(reorder_counts)] <- 0
		names(reorder_counts) <- gene_names
		return(reorder_counts)
	}
	new_mat <- sapply(1:length(libs), generate_cell, gene_names=rownames(prop_mat))
}

# Scalability

set.seed(101)
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Buttner_with_GT.R")
Buttner_orig = list(counts=Buttner_counts, meta=Buttner_meta)

G1 <- resample(Buttner_orig$counts[,Buttner_orig$meta=="G1"], ncells=1000)
S <- resample(Buttner_orig$counts[,Buttner_orig$meta=="S"], ncells=1000)
G2M <- resample(Buttner_orig$counts[,Buttner_orig$meta=="G2M"], ncells=1000)

time <- c()
output <- c();

# 300
Buttner_counts <- cbind(G1[,1:100], S[,1:100], G2M[,1:100]); Buttner_meta = rep(c("G1", "S", "G2M"), each=100); colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)
Butt_cyclone_acc_vs_depth <- Calculate_accuracy_covariate(Butt_cyclone$phase, Buttner_meta, colSums(Buttner_counts))
Butt_cyclemix_acc_vs_depth <- Calculate_accuracy_covariate(Butt_cyclemix$phase, Buttner_meta, colSums(Buttner_counts))
Butt_seurat_acc_vs_depth <- Calculate_accuracy_covariate(Butt_seurat$phase, Buttner_meta, colSums(Buttner_counts))

# 1500
Buttner_counts <- cbind(G1[,1:500], S[,1:500], G2M[,1:500]); Buttner_meta = rep(c("G1", "S", "G2M"), each=500);colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)

# 3000
Buttner_counts <- cbind(G1, S, G2M); Buttner_meta = rep(c("G1", "S", "G2M"), each=1000);colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)

# 6000
Buttner_counts <- cbind(G1, G1, S, S, G2M, G2M); Buttner_meta = rep(c("G1", "S", "G2M"), each=2000);colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)

# 9000
Buttner_counts <- cbind(G1, G1,G1, S, S,S, G2M, G2M, G2M); Buttner_meta = rep(c("G1", "S", "G2M"), each=3000); colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)

# 12000
Buttner_counts <- cbind(G1,G1, G1,G1, S,S, S,S, G2M,G2M, G2M, G2M); Buttner_meta = rep(c("G1", "S", "G2M"), each=4000); colnames(Buttner_counts) <- paste0("cell",1:ncol(Buttner_counts));

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)

time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)

saveRDS(list(output=output, time=time), file="Resample_Scalability_output.rds")

### Plotting ###
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
stuff <- readRDS("Resample_Scalability_output.rds")

colnames(stuff$time) <- c("Cyclone", "CycleMix", "Seurat")
png("Scalability.png", width=3, height=3, units="in", res=300)
plot(stuff$output[,4],stuff$time[,1], col=colours["Cyclone"], ylab="time (s)", xlab="N cells", log="y", type="l", lwd=2, ylim=c(0.1, max(stuff$time[,1])), main="Runtime")
lines(stuff$output[,4], stuff$time[,2], col=colours["CycleMix"], lwd=2)
lines(stuff$output[,4], stuff$time[,3], col=colours["Seurat"], lwd=2)
legend("bottomright", 
dev.off()


