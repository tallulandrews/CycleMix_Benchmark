source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
source("~/My_R_Packages/CycleMix/R/Analysis.R")
require("mclust")
require("Matrix")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

HGenesQuiesc <- HGeneSets$Quiesc[,1]
MGenesQuiesc <- General_Map(HGenesQuiesc, in.org="Hsap", out.org="Mmus", in.name="symbol", out.name="symbol")
MGenesQuiesc[MGenesQuiesc == ""] <- HGenesQuiesc[MGenesQuiesc == ""]

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

set.seed(3891)

output <- c();
time <- c();

# Calculates the proportion of neighbours of each cell in group1 that are in group2 based on the graph
calc_purity <- function(is.group1, is.group2, graph) {
	out <- c()
	for (cell in is.group1) {
		nns <- which(graph[cell,]>0)
		score <- sum(nns %in% seur_g2m)/length(nns)
		out <- c(out, score)
	}
	return(out)
}
args <- commandArgs(trailingOnly=TRUE)

print("Moreno")
obj <- readRDS("Moreno_core_gliobastoma.rds") 
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
require(ggplot2)
png("Moreno_celltypes.png", width=7.2, height=6, units="in", res=300)
print(DimPlot(obj, group.by="cell_type")+labs(title=""))
dev.off()
set.seed(101)
obj <- FindVariableFeatures(obj)
obj <- SCTransform(obj, vars.to.regress="donor_id", conserve.memory=TRUE, return.only.var.genes=TRUE, do.correct.umi=FALSE)
obj <- FindNeighbors(obj)
obj <- FindClusters(obj, res=0.8)
png("Moreno_clusters.png", width=7.2, height=6, units="in", res=300)
print(DimPlot(obj, group.by="seurat_clusters", label=TRUE)+labs(title=""))
dev.off()
obj <- FindClusters(obj, res=1.2)
png("Moreno_clusters2.png", width=7.2, height=6, units="in", res=300)
print(DimPlot(obj, group.by="seurat_clusters", label=TRUE)+labs(title=""))
dev.off()
saveRDS("tmp_Moreno_cluster.rds")


is_prolif <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
not_prolif <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
type_column = "celltype_original"

counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Moreno_Glio_CycleMixQ_Phases <- readRDS(file="Moreno_glioblastoma_phases_Quiescence.rds")
Moreno_Glio_Phases <- readRDS(file="Moreno_glioblastoma_phases.rds")
CycleMixQ_stats <- table(Moreno_Glio_CycleMixQ_Phases$cyclemix$phase, Moreno_Glio_CycleMixQ_Phases$metadata[,type_column])

# CycleMix, Seurat, Cyclone accuracy by technology

# Which cell-type do we have enough cells to compare across all datasets:
tab <- table(Moreno_Glio_CycleMixQ_Phases$metadata$assay, Moreno_Glio_CycleMixQ_Phases$metadata$cell_type)
tab <- tab[rownames(tab) != "CEL-seq2",]
cell_types <- tab[,apply(tab,2,min)>50]
assays <- rownames(tab)


source("~/scripts/CycleMix_Benchmark/ColourScheme.R")
for (t in colnames(cell_types)) {
	these_cells <- Moreno_Glio_Phases$metadata$cell_type == t & Moreno_Glio_Phases$metadata$assay %in% assays
	cyclemix_out <- table(Moreno_Glio_Phases$cyclemix$phase[these_cells], Moreno_Glio_Phases$metadata$assay[these_cells])
	seurat_out <- table(Moreno_Glio_Phases$seurat$phase[these_cells], Moreno_Glio_Phases$metadata$assay[these_cells])
	
	dat <- (t(seurat_out)/colSums(seurat_out))
	rownames(dat) <- c("Smartseq2", "STRTseq", "10x 3' v2", "10x 3' v3","CELseq2", "10x 5' v1", "microwell", "BD Rhapsody")
	png(paste("2Assay_Consistency_Seurat", t, "barplot.png", sep="_"), width=6*0.75, height=5*0.75, units="in", res=150)
	par(mar=c(7, 4, 3, 1))
	barplot(t(dat[c(3,4,6,7,8,1,2),])*100, col=phases_col[match(colnames(dat), names(phases_col))], las=2, main=t, ylab="Cells (%)")
	dev.off()
	
	dat <- (t(cyclemix_out)/colSums(cyclemix_out))
	dat <- dat[,c(3,1,2)]
	rownames(dat) <- c("Smartseq2", "STRTseq", "10x 3' v2", "10x 3' v3", "CELseq2","10x 5' v1", "microwell", "BD Rhapsody")
	png(paste("2Assay_Consistency_CycleMix", t, "barplot.png", sep="_"), width=6*0.75, height=5*0.75, units="in", res=150)
	par(mar=c(7, 4, 3, 1))
	barplot(t(dat[c(3,4,6,7,8,1,2),])*100, col=phases_col[match(colnames(dat), names(phases_col))], las=2, main=t, ylab="Cells (%)")
	dev.off()
}

# Plot
require(ggplot2)
require(ggpattern)

#png("Real_CM_Quiescent.png", width=4, height=4, units="in", res=300)
p1 <- ggplot(df, aes(x=tmp, y=G0, pattern=pop))+geom_bar_pattern(stat="identity", pattern_fill="black", 
	fill="white", color="black", pattern_spacing=0.05, 
	pattern_frequency=5, pattern_angle=45)+facet_wrap(~dataset, scales="free_x")+
scale_x_discrete(labels=c("prolif", "quiesc"))+xlab("")+ylab("Quiescent Cells (%)")+ scale_pattern_discrete(choices=c('none', 'stripe'))+ theme_classic() + guides(fill="none", pattern="none")
print(p1)
#dev.off()



