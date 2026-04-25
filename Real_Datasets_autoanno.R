source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
source("~/My_R_Packages/CycleMix/R/Analysis.R")
require("mclust")
require("Matrix")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

set.seed(3891)

output <- c();
time <- c();

run_phase_classifications <- function(seur_obj, species=c("hsap", "mmus")) {
	seur_obj <- NormalizeData(seur_obj)
	if (species == "hsap") {
		s.genes <- cc.genes$s.genes
		g2m.genes <- cc.genes$g2m.genes
		gene_table <- HGeneSets$Tirosh
	} else if (species == "mmus") {
		gene_table <- MGeneSets$Cyclone
		s.genes <- gene_table[(gene_table[,2] == "S" | gene_table[,2] == "G1S") & gene_table[,3] == 1,1]
		g2m.genes <- gene_table[(gene_table[,2] == "G2M" | gene_table[,2] == "M" | gene_table[,2] == "G2") & gene_table[,3] == 1,1]
	}
	cyclemix <- run_CycleMix(GetAssayData(seur_obj, layer="counts"), gene_table)
	seur.elapsed.time <- system.time(seur_obj <- CellCycleScoring(seur_obj, s.features=s.genes, g2m.features=g2m.genes))[3]
	seurat <- list(phase=seur_obj@meta.data$Phase, time=seur.elapsed.time)
	if (ncol(seur_obj ) < 80000) {
	cyclone.elapsed.time <- system.time(cyclone <- run_cyclone(GetAssayData(seur_obj, layer="counts"), species=species))
	} else {
	cyclone.elapsed.time <- rep(NA, length(seur.elapsed.time))
	cyclone <- NA
	}
	return(list(metadata=seur_obj@meta.data, 
		cyclemix=cyclemix, seurat=seurat, cyclone=cyclone, 
		time=rbind(cyclemix$time, seur.elapsed.time, cyclone.elapsed.time)))
}

# Calculates the proportion of neighbours of each cell in group1 that are in group2 based on the graph
calc_purity50 <- function(is.group1, is.group2, graph) {
	out <- c()
	for (cell in is.group1) {
		nns <- which(graph[cell,]>0)
		score <- sum(nns %in% seur_g2m)/length(nns)
		out <- c(out, score)
	}
	return(out)
}

calc_purity50_trio <- function(Phases, nn_graph){
	seur_g2m <- which(Phases$seurat$phase == "G2M")
	seur_S <- which(Phases$seurat$phase == "S")
	seur_g2m_purity50 <- calc_purity50(seur_g2m, seur_g2m, nn_graph)
	seur_S_purity50 <- calc_purity50(seur_S, seur_S, nn_graph)

	cyc_g2m <- which(Phases$cyclemix$phase == "G2M")
	cyc_S <- which(Phases$cyclemix$phase == "G1S")
	cyc_g2m_purity50 <- calc_purity50(cyc_g2m, cyc_g2m, nn_graph)
	cyc_S_purity50 <- calc_purity50(cyc_S, cyc_S, nn_graph)

	lone_g2m <- which(Phases$cyclone$phase == "G2M")
	lone_S <- which(Phases$cyclone$phase == "S")
	lone_g2m_purity50 <- calc_purity50(lone_g2m, lone_g2m, nn_graph)
	lone_S_purity50 <- calc_purity50(lone_S, lone_S, nn_graph)

	return(list(seur_g2m_purity50=seur_g2m_purity50,seur_S_purity50=seur_S_purity50, cyc_g2m_purity50=cyc_g2m_purity50, cyc_S_purity50=cyc_S_purity50, lone_g2m_purity50=lone_g2m_purity50, lone_S_purity50=lone_S_purity50))

}


args <- commandArgs(trailingOnly=TRUE)
is.prolif <- list()
is.quies <- list()
type.col <- list()
is.prolif[["Wang"]] <- c("transit amplifying cell of colon")
is.prolif[["Wang1"]] <- c("transit amplifying cell of colon")
is.prolif[["Wang2"]] <- c("transit amplifying cell of colon")
is.quies[["Wang"]] <- c("colon goblet cell", "paneth cell of colon")
type.col[["Wang"]] <- "cell_type"
type.col[["Wang1"]] <- "cell_type"
type.col[["Wang2"]] <- "cell_type"
is.prolif[["Street"]] <- c("PRM1+ Cells")
is.quies[["Street"]] <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
type.col[["Street"]] <- "author_cell_type"
is.prolif[["Gu"]] <- c("Plasmablast-IgM+_cycling","Prolif-Tregs", "GC.BC_DZ.1", "GC.BC_DZ.2", "DC.prolif.")
is.quies[["Gu"]] <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
type.col[["Gu"]] <- "final_annotation"
is.prolif[["Moreno"]] <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
is.quies[["Moreno"]] <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
type.col[["Moreno"]] <- "celltype_original"
is.prolif[["Triana"]] <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")
is.quies[["Triana"]] <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")
type.col[["Triana"]] <- "cell_type"

if (args[1]==1) {
print("Wang1")
# Wang Colon
obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Idents(obj) <- obj@meta.data[,type.col[["Wang"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Wang"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
	this <- markers[markers$cluster==x,];
	this$diff_detect <- this$pct.1-this$pct.2
	this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
	up <- this[this$avg_log2FC > 0,][1:10,];
	dn <- this[this$avg_log2FC < 0,][1:10,];
	return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
	gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
prop_labelled <- sum(out$phase != "None")/ncol(obj)
obj@meta.data$auto_anno <- out$phase
tmp <- obj@meta.data[obj@meta.data$auto_anno != "None",]
tab <- table(tmp[,type.col[["Wang"]]], tmp$auto_anno)
tab <- tab[order(rownames(tab)), order(colnames(tab))]
tab <- tab[rownames(tab) %in% colnames(tab),]
tab <- tab[,colnames(tab) %in% rownames(tab)]
sum(diag(tab))/nrow(tmp)

saveRDS(out, "Wang_colon_autoanno.rds")

}

if (args[1]==2) {
print("Wang2")
# Wang ileum
obj <- readRDS("Wang_ileum.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Idents(obj) <- obj@meta.data[,type.col[["Wang"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Wang"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
        this <- markers[markers$cluster==x,];
        this$diff_detect <- this$pct.1-this$pct.2
        this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
        up <- this[this$avg_log2FC > 0,][1:10,];
        dn <- this[this$avg_log2FC < 0,][1:10,];
        return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
        gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
saveRDS(out, "Wang_ileum_autoanno.rds")

}

if (args[1]==3) {
print("Moreno") ### FAILED ###
# Moreno_glioblastoma
obj <- readRDS("Moreno_core_gliobastoma.rds") 
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
set.seed(998)
counts <- counts[,sample(1:ncol(counts), size=70000)]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data[match(colnames(counts),rownames(obj@meta.data)),])
Idents(obj) <- obj@meta.data[,type.col[["Moreno"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Moreno"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
        this <- markers[markers$cluster==x,];
        this$diff_detect <- this$pct.1-this$pct.2
        this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
        up <- this[this$avg_log2FC > 0,][1:10,];
        dn <- this[this$avg_log2FC < 0,][1:10,];
        return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
        gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
saveRDS(out, "Moreno_glioblastoma_autoanno.rds")


}


if (args[1]==4) {
print("Gu")
# Gu Lymph
obj <- readRDS("Gu_lymph.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Idents(obj) <- obj@meta.data[,type.col[["Gu"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Gu"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
        this <- markers[markers$cluster==x,];
        this$diff_detect <- this$pct.1-this$pct.2
        this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
        up <- this[this$avg_log2FC > 0,][1:10,];
        dn <- this[this$avg_log2FC < 0,][1:10,];
        return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
        gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
saveRDS(out, "Gu_lymph_autoanno.rds")


}

if (args[1]==5) {
print("Triana") 
# Triana bonemarrow
obj <- readRDS("Triana_bonemarrow.rds")
obj@meta.data$UMAP1 <- obj@reductions$bothumap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$bothumap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Idents(obj) <- obj@meta.data[,type.col[["Triana"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Triana"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
        this <- markers[markers$cluster==x,];
        this$diff_detect <- this$pct.1-this$pct.2
        this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
        up <- this[this$avg_log2FC > 0,][1:10,];
        dn <- this[this$avg_log2FC < 0,][1:10,];
        return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
        gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
saveRDS(out, "Triana_bonemarrow_autoanno.rds")

}

if (args[1]==6) {
print("Steets")
# Streats adipose
obj <- readRDS("Streets_adipose.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Idents(obj) <- obj@meta.data[,type.col[["Street"]]]
obj <- NormalizeData(obj)
markers <- FindAllMarkers(obj, min.pct = 0.1, max.cells.per.ident=500, group.by=type.col[["Street"]])
require(dplyr)
gene_tables <- lapply(unique(markers$cluster), function(x){
        this <- markers[markers$cluster==x,];
        this$diff_detect <- this$pct.1-this$pct.2
        this <- this[order(abs(this$diff_detect), decreasing=TRUE),]
        up <- this[this$avg_log2FC > 0,][1:10,];
        dn <- this[this$avg_log2FC < 0,][1:10,];
        return(data.frame(Gene=c(up$gene, dn$gene), Stage=rep(x, 20), Dir=c(rep(1,10),rep(-1,10))))})
gene_table <- c()
for(i in 1:length(gene_tables)){
        gene_table <- rbind(gene_table, gene_tables[[i]])
}

out <- classifyCells(obj, gene_table, expr_name="RNA", symbol_column=NULL)
saveRDS(out, "Streets_glioblastoma_autoanno.rds")

}


# Plots #
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")

# Calculate % of cells annotated (not None), and % correctly annotated.

# This will be supplementary so I can go a bit extravagant

# One plot per dataset with each cell-type with % of cells annotated + % correct for each cell-type as mirrored barplots? - Too complicated & hard to read
# One plot per dataset with Name of cell-type & number of cells per cell-type on one axis, and accuracy on other axis, with number above the bar for (n cells annotated).

# One plot with overall accuracy for each dataset.

dataset2annotations <- list("Gu"="Gu_lymph_autoanno.rds",
			"Triana" = "Triana_bonemarrow_autoanno.rds",
			"Moreno"="Moreno_glioblastoma_autoanno.rds",
			"Wang1" = "Wang_colon_autoanno.rds",
			"Wang2" = "Wang_ileum_autoanno.rds",
			"Street" = "Streets_glioblastoma_autoanno.rds")
dataset2gt_anno <- list("Gu"="Gu_lymph_phases.rds",
			"Moreno"="Moreno_glioblastoma_phases.rds",
			"Street"="Streets_adipose_phases.rds",
			"Triana"="Triana_bonemarrow_phases.rds",
			"Wang1" = "Wang_colon_phases.rds",
			"Wang2" = "Wang_ileum_phases.rds")

create_plot <- function(dataset) {
	# Read Data
	auto_anno <- readRDS(dataset2annotations[[dataset]])$phase
	truth.col <- type.col[[dataset]]
	truth <- readRDS(dataset2gt_anno[[dataset]])$metadata
	if (dataset == "Moreno") {
		set.seed(998)
		truth <- truth[sample(1:nrow(truth),size=70000),]
		keep <- truth[,truth.col] != "unknown"
		truth <- truth[keep,]
		auto_anno <- auto_anno[keep]
	}

	# n cells/celltype
	celltype2Ncells <- table(truth[,truth.col])

	# accuracy/celltype
	celltype_vs_autoanno <- table(truth[,truth.col], auto_anno)
	# exclude the "none" category
	n_not_anno <- celltype_vs_autoanno[,"None"];
	n_not_anno <- n_not_anno[match(names(celltype2Ncells), names(n_not_anno))]
	n_anno <- celltype2Ncells - n_not_anno
	
	# reorder stuff
	celltype_vs_autoanno <- celltype_vs_autoanno[
		match(names(celltype2Ncells), rownames(celltype_vs_autoanno)),
		match(names(celltype2Ncells), colnames(celltype_vs_autoanno))]
	celltype_vs_autoanno[is.na(celltype_vs_autoanno)] <- 0

	# calc scores
	accuracy <- diag(celltype_vs_autoanno)/(celltype2Ncells - n_not_anno)
	TP <- diag(celltype_vs_autoanno)
	FN <- rowSums(celltype_vs_autoanno) - TP
	FP <- colSums(celltype_vs_autoanno) - TP
	TN <- sum(celltype2Ncells)-n_not_anno-TP-FN-FP

	MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

	#Trim to remove rare cell-types
	keep <- celltype2Ncells >= 20
	png(paste(dataset, "autoanno_SFig.png", sep=""), width=8, height=sum(keep)/8+5, units="in", res=100)
	par(mar=c(4,20,2,1))
	coords <- barplot(accuracy[keep], horiz=TRUE, xlim=c(0,1.2), las=1, xlab="Recall", main=dataset)
	abline(v=1, lty=2, col="grey65")
	#text(x=accuracy[keep], y=coords, n_anno[keep], pos=4)
	text(x=rep(1, length(coords)), y=coords, n_anno[keep], pos=4)
	dev.off()

	require(igraph)
	# Exclude None
	keep <- auto_anno != "None"
	return(compare(auto_anno[keep],truth[keep,truth.col], method="adjusted.rand"))
}
	

ARI_Wang1 = create_plot("Wang1")
ARI_Wang2 = create_plot("Wang2")
ARI_Street = create_plot("Street")
ARI_Triana = create_plot("Triana")
ARI_Gu = create_plot("Gu")
ARI_Moreno = create_plot("Moreno")

png("Dataset_ARI.png", width=6, height=6, units="in", res=100)
barplot(c(ARI_Wang1, ARI_Wang2, ARI_Street, ARI_Triana, ARI_Gu, ARI_Moreno), col=dataset_col[c("Wang1", "Wang2", "Streets", "Tiana", "Gu", "Moreno")], ylab="ARI", names=c("Wang1", "Wang2", "Streets", "Triana", "Gu", "Moreno"), ylim=c(0,1))
dev.off()

