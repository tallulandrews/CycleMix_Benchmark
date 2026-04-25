source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

Calculate_accuracy_covariate <- function(pred, truth, covariate, window=50) {
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
	if (sum(is.na(pred)) > 0) {
		pred[is.na(pred)] <- "G1"
	}
	TP <- pred==truth
	reorder <- order(covariate)
	covariate<- covariate[reorder]
	TP<- TP[reorder]
	cx <- c(0,cumsum(TP))
	n <- window
	rsum <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
	rsum <- c(rsum,rep(rsum[length(rsum)], times=length(covariate)-length(rsum))) 
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
test1 <-  readRDS("Wang_ileum.rds")
test2 <-  readRDS("Wang_colon.rds")
test3 <-  readRDS("Gu_lymph.rds")

ngenes <- c(colSums(test1@assays$RNA@counts > 0), colSums(test2@assays$RNA@counts > 0), colSums(test3@assays$RNA@counts > 0))
lib.sizes <- c(colSums(test1@assays$RNA@counts), colSums(test2@assays$RNA@counts), colSums(test3@assays$RNA@counts))
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

output <- c();
time <- c();

#source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Sasagawa_with_GT.R")
#Sasagawa_meta <- sub("\\/", "", Sasagawa_meta)
#Sask_cyclone <- run_cyclone(Sasagawa_fpkms, species="mmus")
#Sask_cyclemix <- run_CycleMix(Sasagawa_fpkms, MGeneSets$Cyclone)
#Sask_seurat <- run_Seurat(Sasagawa_fpkms, MGeneSets$Cyclone)
#acc <- c(Calculate_accuracy_trio(Sask_cyclone$phase, Sasagawa_meta), 41/length(Sasagawa_meta), Calculate_accuracy_trio(Sask_seurat$phase, Sasagawa_meta), n=sum(Sasagawa_meta != "Unknown"))
#output <- rbind(output, acc)
#time <- rbind(time, c(Sask_cyclone$time, Sask_cyclemix$time, Sask_seurat$time))

print("Resample Beuttner")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Buttner_with_GT.R")
Buttner_orig = list(counts=Buttner_counts, meta=Buttner_meta)

set.seed(3736)
G1 <- resample(Buttner_counts[,Buttner_meta=="G1"], ncells=1000)
S <- resample(Buttner_counts[,Buttner_meta=="S"], ncells=1000)
G2M <- resample(Buttner_counts[,Buttner_meta=="G2M"], ncells=1000)
Buttner_counts <- cbind(G1, S, G2M); Buttner_meta = rep(c("G1", "S", "G2M"), each=1000)
colnames(Buttner_counts) <- paste0("Cell", 1:ncol(Buttner_counts))

png("Demo_gene_counts.png", width=4, height=4, units="in", res=150)
par(mar=c(4,4,1,1))
my_col1 <- rgb(255, 0, 255, max = 255, alpha = 125, names = "purple50")
my_col2 <- rgb(20, 20, 225, max = 255, alpha = 125, names = "blue50")
hist(ngenes, col=my_col1, xlab="Detected Genes", main="", freq=FALSE)
hist(colSums(Buttner_counts>0), add=TRUE, col=my_col2, xlab="", ylab="", main="", ylim=c(0,600), freq=FALSE)
legend("topright", c("Chromium", "Synthetic"), fill=c(my_col1, my_col2), bty="n")
dev.off()

png("Demo_gene_counts2.png", width=4, height=4, units="in", res=150)
plot(lib.sizes, ngenes, log="xy", xlab="Total Counts", ylab="Detected Genes", col=my_col1, pch=16, cex=0.5)
points(colSums(Buttner_counts), colSums(Buttner_counts>0), col=my_col2, pch=16, cex=0.5)
legend("topleft", c("Chromium", "Synthetic"), fill=c(my_col1, my_col2), bty="n")
dev.off()

detect <- c(rowMeans(test1@assays$RNA@counts > 0), rowMeans(test2@assays$RNA@counts > 0), rowMeans(test3@assays$RNA@counts > 0))
means <- c(rowMeans(test1@assays$RNA@counts), rowMeans(test2@assays$RNA@counts), rowMeans(test3@assays$RNA@counts))

png("Demo_gene_counts3.png", width=4, height=4, units="in", res=150)
plot(means, detect, log="x", xlab="Gene Mean", ylab="Detection rate", col=my_col1, pch=16, cex=0.5)
points(rowMeans(Buttner_counts), rowMeans(Buttner_counts>0), col=my_col2, pch=16, cex=0.5)
legend("topleft", c("Chromium", "Synthetic"), fill=c(my_col1, my_col2), bty="n")
dev.off()

# Distribution of gene counts
Buttner_counts <- Buttner_counts[rowSums(Buttner_counts) > 0,]
test1_counts <- GetAssayData(test1, layer="counts"); test1_counts <- test1_counts[rowSums(test1_counts)>0,]
test2_counts <- GetAssayData(test2, layer="counts"); test2_counts <- test2_counts[rowSums(test2_counts)>0,]
test3_counts <- GetAssayData(test3, layer="counts"); test3_counts <- test3_counts[rowSums(test3_counts)>0,]
Synth_0 <- sum(Buttner_counts==0)/prod(dim(Buttner_counts))
Synth_1 <- sum(Buttner_counts==1)/prod(dim(Buttner_counts))
Synth_2 <- sum(Buttner_counts==2)/prod(dim(Buttner_counts))
Synth_3 <- sum(Buttner_counts==3)/prod(dim(Buttner_counts))
Synth_4 <- sum(Buttner_counts==4)/prod(dim(Buttner_counts))
denom <- sum(prod(dim(test1_counts)), prod(dim(test2_counts)), prod(dim(test3_counts)))
True_0 <- (sum(test1_counts==0)+sum(test2_counts==0)+sum(test3_counts==0))/denom
True_1 <- (sum(test1_counts==1)+sum(test2_counts==1)+sum(test3_counts==1))/denom
True_2 <- (sum(test1_counts==2)+sum(test2_counts==2)+sum(test3_counts==2))/denom
True_3 <- (sum(test1_counts==3)+sum(test2_counts==3)+sum(test3_counts==3))/denom
True_4 <- (sum(test1_counts==4)+sum(test2_counts==4)+sum(test3_counts==4))/denom

png("Demo_gene_counts3.png", width=4, height=4, units="in", res=150)
plot(means, detect, log="x", xlab="Gene Mean", ylab="Detection rate", col=my_col1, pch=16, cex=0.5)
points(rowMeans(Buttner_counts), rowMeans(Buttner_counts>0), col=my_col2, pch=16, cex=0.5)
legend("topleft", c("Chromium", "Synthetic"), fill=c(my_col1, my_col2), bty="n")
dev.off()

Butt_cyclone <- run_cyclone(Buttner_counts, species="mmus")
Butt_cyclemix <- run_CycleMix(Buttner_counts, MGeneSets$Cyclone)
Butt_seurat <- run_Seurat(Buttner_counts, MGeneSets$Cyclone)
acc <- c(Calculate_accuracy_trio(Butt_cyclone$phase, Buttner_meta), Calculate_accuracy_trio(Butt_cyclemix$phase, Buttner_meta), Calculate_accuracy_trio(Butt_seurat$phase, Buttner_meta), n=sum(Buttner_meta != "Unknown"))
output <- rbind(output, acc)
time <- rbind(time, c(Butt_cyclone$time, Butt_cyclemix$time, Butt_seurat$time))

Butt_cyclone_acc_vs_depth <- Calculate_accuracy_covariate(Butt_cyclone$phase, Buttner_meta, colSums(Buttner_counts))
Butt_cyclemix_acc_vs_depth <- Calculate_accuracy_covariate(Butt_cyclemix$phase, Buttner_meta, colSums(Buttner_counts))
Butt_seurat_acc_vs_depth <- Calculate_accuracy_covariate(Butt_seurat$phase, Buttner_meta, colSums(Buttner_counts))

print("Resample Chen")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Chen_with_GT.R")
Chen_orig = list(counts=Chen_counts, meta=Chen_meta)

set.seed(4236)
G1 <- resample(Chen_counts[,Chen_meta=="G1"], ncells=1000)
S <- resample(Chen_counts[,Chen_meta=="S"], ncells=1000)
G2M <- resample(Chen_counts[,Chen_meta=="G2M"], ncells=1000)
Chen_counts <- cbind(G1, S, G2M); Chen_meta = rep(c("G1", "S", "G2M"), each=1000)
colnames(Chen_counts) <- paste0("Cell", 1:ncol(Chen_counts))

Chen_cyclone <- run_cyclone(Chen_counts, species="mmus")
Chen_cyclemix <- run_CycleMix(Chen_counts, MGeneSets$Cyclone)
Chen_seurat <- run_Seurat(Chen_counts, MGeneSets$Cyclone)
acc <- c(Calculate_accuracy_trio(Chen_cyclone$phase, Chen_meta), Calculate_accuracy_trio(Chen_cyclemix$phase, Chen_meta), Calculate_accuracy_trio(Chen_seurat$phase, Chen_meta), n=sum(Chen_meta != "Unknown"))
output <- rbind(output, acc)
time <- rbind(time, c(Chen_cyclone$time, Chen_cyclemix$time, Chen_seurat$time))

Chen_cyclone_acc_vs_depth <- Calculate_accuracy_covariate(Chen_cyclone$phase, Chen_meta, colSums(Chen_counts))
Chen_cyclemix_acc_vs_depth <- Calculate_accuracy_covariate(Chen_cyclemix$phase, Chen_meta, colSums(Chen_counts))
Chen_seurat_acc_vs_depth <- Calculate_accuracy_covariate(Chen_seurat$phase, Chen_meta, colSums(Chen_counts))

print("Resample Leng")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Leng_with_GT.R")
Leng_orig = list(counts=Leng_counts, meta=Leng_meta)

set.seed(8516)
G1 <- resample(Leng_counts[,Leng_meta=="G1"], ncells=1000)
S <- resample(Leng_counts[,Leng_meta=="S"], ncells=1000)
G2M <- resample(Leng_counts[,Leng_meta=="G2"], ncells=1000)
Leng_counts <- cbind(G1, S, G2M); Leng_meta = rep(c("G1", "S", "G2M"), each=1000)
colnames(Leng_counts) <- paste0("Cell", 1:ncol(Leng_counts))

Leng_cyclone <- run_cyclone(Leng_counts, species="hsap")
Leng_cyclemix <- run_CycleMix(Leng_counts, HGeneSets$Seurat)
Leng_seurat <- run_Seurat(Leng_counts, HGeneSets$Seurat)
acc <- c(Calculate_accuracy_trio(Leng_cyclone$phase, Leng_meta), Calculate_accuracy_trio(Leng_cyclemix$phase, Leng_meta), Calculate_accuracy_trio(Leng_seurat$phase, Leng_meta), n=sum(Leng_meta != "Unknown"))
output <- rbind(output, acc)
time <- rbind(time, c(Leng_cyclone$time, Leng_cyclemix$time, Leng_seurat$time))

Leng_cyclone_acc_vs_depth <- Calculate_accuracy_covariate(Leng_cyclone$phase, Leng_meta, colSums(Leng_counts))
Leng_cyclemix_acc_vs_depth <- Calculate_accuracy_covariate(Leng_cyclemix$phase, Leng_meta, colSums(Leng_counts))
Leng_seurat_acc_vs_depth <- Calculate_accuracy_covariate(Leng_seurat$phase, Leng_meta, colSums(Leng_counts))

vsdepth <- list(
	Buttner=list(cyclone=Butt_cyclone_acc_vs_depth, cyclemix=Butt_cyclemix_acc_vs_depth, seurat=Butt_seurat_acc_vs_depth),
	Chen=list(cyclone=Chen_cyclone_acc_vs_depth, cyclemix=Chen_cyclemix_acc_vs_depth, seurat=Chen_seurat_acc_vs_depth),
	Leng=list(cyclone=Leng_cyclone_acc_vs_depth, cyclemix=Leng_cyclemix_acc_vs_depth, seurat=Leng_seurat_acc_vs_depth)
)
saveRDS(vsdepth, file="Resampled_Benchmarking_vsdepth.rds")
#source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_MacDavid_with_GT.R")
#MacDavid_orig = list(counts=MacDavid_counts, meta=MacDavid_meta)

#G1 <- resample(MacDavid_counts[,MacDavid_meta=="G1"], ncells=1000)
#S <- resample(MacDavid_counts[,MacDavid_meta=="S"], ncells=1000)
#G2M <- resample(MacDavid_counts[,MacDavid_meta=="G2M"], ncells=1000)
#MacDavid_counts <- cbind(G1, S, G2M); MacDavid_meta = rep(c("G1", "S", "G2M"), each=1000)
#colnames(MacDavid_counts) <- paste0("Cell", 1:ncol(MacDavid_counts))

#MacDavid_meta[MacDavid_meta == "g0/g1"] <- "G1"
#MacDavid_meta[MacDavid_meta == "g2/m"] <- "G2M"
#MacDavid_meta[MacDavid_meta == "s"] <- "S"
#G1 <- resample(MacDavid_counts[,MacDavid_meta=="G1"], ncells=1000)
#S <- resample(MacDavid_counts[,MacDavid_meta=="S"], ncells=1000)
#G2M <- resample(MacDavid_counts[,MacDavid_meta=="G2M"], ncells=1000)
#MacDavid_cyclone <- run_cyclone(MacDavid_counts, species="hsap")
#MacDavid_cyclemix <- run_CycleMix(MacDavid_counts, HGeneSets$Seurat)
#MacDavid_seurat <- run_Seurat(MacDavid_counts, HGeneSets$Seurat)
#acc <- c(Calculate_accuracy_trio(MacDavid_cyclone$phase, MacDavid_meta), Calculate_accuracy_trio(MacDavid_cyclemix$phase, MacDavid_meta), Calculate_accuracy_trio(MacDavid_seurat$phase, MacDavid_meta), n=sum(MacDavid_meta != "Unknown"))
#output <- rbind(output, acc)
#time <- rbind(time, c(MacDavid_cyclone$time, MacDavid_cyclemix$time, MacDavid_seurat$time))


## Plots! ##
print("Making Plots")
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
#saveRDS(list(output=output, time=time), file="Resampled_Benchmarking_output.rds")

colnames(output) <- c("Cyclone", "CycleMix", "Seurat", "n")
rownames(output) <- c("Buettner", "Chen", "Leng")
colnames(time) <- c("Cyclone", "CycleMix", "Seurat")
rownames(time) <-  c("Buettner", "Chen", "Leng")
output_strerr <- sqrt(output[,1:3]*(1-output[,1:3])/output[,4])

#png("Accuracy_by_method_resampled.png", width=6, height=6, units="in", res=300)
#coords <- barplot(t(output[,1:3]), beside=TRUE, col=colours, ylim=c(0,1), ylab="Accuracy")
#arrows(x0=coords,y0=t(output[,1:3]), x1=coords, y1=t(output[,1:3]+output_strerr*1.96), length=0)
#legend("topleft", c("Cyclone", "CycleMix", "Seurat"), fill=colours, bty="n")
#dev.off()

#labels <- paste0(rownames(time), "\n(",output[,4],")")
#png("Time_vs_N_resampled.png", width=6, height=6, units="in", res=300)
#barplot(t(time), col=colours, ylab="Run Time (s)", beside=TRUE, names=labels)
#dev.off()

# vs depth
vsdepth <- readRDS(file="Resampled_Benchmarking_vsdepth.rds")
Butt_cyclone_acc_vs_depth <- vsdepth$Buttner$cyclone
Butt_cyclemix_acc_vs_depth <- vsdepth$Buttner$cyclemix
Butt_seurat_acc_vs_depth <- vsdepth$Buttner$seurat
Chen_cyclone_acc_vs_depth <- vsdepth$Chen$cyclone
Chen_cyclemix_acc_vs_depth <- vsdepth$Chen$cyclemix
Chen_seurat_acc_vs_depth <- vsdepth$Chen$seurat
Leng_cyclone_acc_vs_depth <- vsdepth$Leng$cyclone
Leng_cyclemix_acc_vs_depth <- vsdepth$Leng$cyclemix
Leng_seurat_acc_vs_depth <- vsdepth$Leng$seurat
png("Resample_vsdepth_Accuracy.png", width=5, height=5, units="in", res=300)
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))
plot(Butt_cyclone_acc_vs_depth[,1], Butt_cyclone_acc_vs_depth[,2], cex=0.2, pch=16, type="p", ylab="Accuracy", xlab="Reads/cell", col=colours["Cyclone"], ylim=c(0,1))
tmp <- smooth.spline(Butt_cyclone_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Cyclone"])

points(Butt_cyclemix_acc_vs_depth[,1], Butt_cyclemix_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["CycleMix"])
tmp <- smooth.spline(Butt_cyclemix_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["CycleMix"])

points(Butt_seurat_acc_vs_depth[,1], Butt_seurat_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["Seurat"])
tmp <- smooth.spline(Butt_seurat_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Seurat"])


plot(Chen_cyclone_acc_vs_depth[,1], Chen_cyclone_acc_vs_depth[,2], cex=0.2, pch=16, type="p", ylab="Accuracy", xlab="Reads/cell", col=colours["Cyclone"], ylim=c(0,1))
tmp <- smooth.spline(Chen_cyclone_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Cyclone"])

points(Chen_cyclemix_acc_vs_depth[,1], Chen_cyclemix_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["CycleMix"])
tmp <- smooth.spline(Chen_cyclemix_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["CycleMix"])

points(Chen_seurat_acc_vs_depth[,1], Chen_seurat_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["Seurat"])
tmp <- smooth.spline(Chen_seurat_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Seurat"])

plot(Leng_cyclone_acc_vs_depth[,1], Leng_cyclone_acc_vs_depth[,2], cex=0.2, pch=16, type="p", ylab="Accuracy", xlab="Reads/cell", col=colours["Cyclone"], ylim=c(0,1))
tmp <- smooth.spline(Leng_cyclone_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Cyclone"])

points(Leng_cyclemix_acc_vs_depth[,1], Leng_cyclemix_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["CycleMix"])
tmp <- smooth.spline(Leng_cyclemix_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["CycleMix"])

points(Leng_seurat_acc_vs_depth[,1], Leng_seurat_acc_vs_depth[,2], cex=0.2, pch=16, col=colours["Seurat"])
tmp <- smooth.spline(Leng_seurat_acc_vs_depth,spar=1.2)
lines(tmp$x, tmp$y, lwd=2, col=colours["Seurat"])

plot.new()
legend("left", bty="n", lwd=2, col=colours, names(colours))
dev.off()

#saveRDS(list(output=output, time=time), file="Resampled_Benchmarking_output.rds")


