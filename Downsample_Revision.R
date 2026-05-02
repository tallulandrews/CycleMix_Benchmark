source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/Seurat.R")
source("~/My_R_Packages/CycleMix/R/Analysis.R")
 source("~/My_R_Packages/CycleMix/R/Plotting.R")
require("mclust")
require("Matrix")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

set.seed(3891)

output <- c();
time <- c();

source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")
remap_genesets <- function(GeneSet, to.spp = c("hsap", "mmus")){
	species <- "";
	if(sum(GeneSet$Gene %in% ensg_name_map$hgnc_symbol)/nrow(GeneSet) > 0.5) {
		species <- "human"
	} else if (sum(GeneSet$Gene %in% musg_name_map$mgi_symbol)/nrow(GeneSet) > 0.5) {
		species <- "mouse"
	} else {
		stop("Unrecognized species.")
	}
	GeneSet$Orig.Gene <- GeneSet$Gene
	if (species == "human" & to.spp == "mmus") {
		GeneSet$Gene <- General_Map(GeneSet$Orig.Gene, in.name="symbol", out.name="symbol", in.org="Hsap", out.org="Mmus")
	}
	if (species == "mouse" & to.spp == "hsap") {
		GeneSet$Gene <- General_Map(GeneSet$Orig.Gene, in.name="symbol", out.name="symbol", out.org="Hsap", in.org="Mmus")

	}
	GeneSet <- GeneSet[GeneSet$Gene != "",]
	GeneSet <- GeneSet[!is.na(GeneSet$Gene),]
	return(GeneSet)
}

downsample_fit <- function(seur_obj, species=c("hsap","mmus"), plot_prefix, downprop=0.5) {
        if (species == "hsap") {
                s.genes <- cc.genes$s.genes
                g2m.genes <- cc.genes$g2m.genes
                gene_table <- HGeneSets$Tirosh
        } else if (species == "mmus") {
                gene_table <- MGeneSets$Cyclone
                s.genes <- gene_table[(gene_table[,2] == "S" | gene_table[,2] == "G1S") & gene_table[,3] == 1,1]
                g2m.genes <- gene_table[(gene_table[,2] == "G2M" | gene_table[,2] == "M" | gene_table[,2] == "G2") & gene_table[,3] == 1,1]
        }
	counts <- GetAssayData(seur_obj, slot="counts")
	SCE <- SingleCellExperiment(assays=list(counts=counts), rowData=data.frame(feature_symbol=rownames(counts)))
	require("scuttle")
	downsampled <- downsampleMatrix(counts, prop=downprop)
	SCE <- SingleCellExperiment(assays=list(counts=downsampled), rowData=data.frame(feature_symbol=rownames(counts)))
        logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(counts)/median(colSums(counts)))
        print("Cyclemix is classifying cells")
        my_output = classifyCells(SCE, gene_table)
	#print(my_output)
#	for (phase in names(my_output$fits)) {
#		pdf(paste0(plot_prefix,"_",phase,"_", downprop,"_mixturePlot.pdf"), width=4*1.5,height=3*1.5)
#		plotMixture(my_output$fits[[phase]], BIC=TRUE)
#		dev.off()
#	}
	ks_check_fit <- checkFit(my_output, nbootstrap=100, summarize=TRUE)
	saveRDS(list("cyclemix"=my_output,"checkfit"= ks_check_fit), paste0(plot_prefix,"_", downprop,"_CycleMix_fits.rds"))
#	pdf(paste0(plot_prefix,"_G1S_mixturePlot.pdf"), width=4*1.5,height=3*1.5)
#	plotMixture(my_output$fits$G1S, BIC=TRUE)
#	dev.off()
}

args <- commandArgs(trailingOnly=TRUE) # 1 = dataset (1-6), 2 = make plots (TRUE/FALSE)

if (args[1]==1) {
print("Wang1")
Wang_column <- "cell_type"
# Wang Colon
obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, slot="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
dim(obj)
sum(GetAssayData(obj, slot="counts") == 0)/prod(dim(obj))
#Wang_Colon_Phases <- run_phase_classifications(obj, species="hsap")
#saveRDS(Wang_Colon_Phases, file="Wang_colon_phases.rds")

for(prop in seq(from=0.05, to=0.95, by=0.05)) {
	downsample_fit(obj, species="hsap", plot_prefix="Wang_colon", downprop=prop)
}

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="hsap")
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="hsap", type_col=Wang_column)

#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Wang_colon_advphases_allSeurGenes.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="hsap")
#saveRDS(phases_genesets, "Wang_colon_cyclemix_allgenesets.rds")



#Wang_Colon_Phases<- readRDS("Wang_colon_phases.rds")

## Neighbourhood Purity ##
#obj <- ScaleData(obj)
#obj <- FindVariableFeatures(obj)
#obj <- RunPCA(obj, npcs=20)
#obj <- FindNeighbors(obj, k=20)
#nn_graph <- obj@graphs$RNA_nn

#seur_g2m <- which(Wang_Colon_Phases$seurat$phase == "G2M")
#seur_S <- which(Wang_Colon_Phases$seurat$phase == "S")
#seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
#seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)

#cyc_g2m <- which(Wang_Colon_Phases$cyclemix$phase == "G2M")
#cyc_S <- which(Wang_Colon_Phases$cyclemix$phase == "G1S")
#cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
#cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)

#lone_g2m <- which(Wang_Colon_Phases$cyclone$phase == "G2M")
#lone_S <- which(Wang_Colon_Phases$cyclone$phase == "S")
#lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
#lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

#saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Wang_colon_purity20.rds")
}

if (args[1]==2) {
print("Wang2")
Wang_column <- "cell_type"
# Wang ileum
obj <- readRDS("Wang_ileum.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, slot="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
#Wang_Ileum_Phases <- run_phase_classifications(obj, species="hsap")
#saveRDS(Wang_Ileum_Phases, file="Wang_ileum_phases.rds")

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="hsap")
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="hsap", type_col=Wang_column)

#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Wang_ileum_advphases.rds")
#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Wang_ileum_advphases_allSeurGenes.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="hsap")
#saveRDS(phases_genesets, "Wang_ileum_cyclemix_allgenesets.rds")

good.fit <- goodness_of_fit(obj, plot_prefix="Wang_ileum", species="hsap")
saveRDS(good.fit, "Wang_ileum_goodnessoffit.rds")

#Wang_Ileum_Phases <- readRDS(file="Wang_ileum_phases.rds")


## Neighbourhood Purity ##
#obj <- ScaleData(obj)
#obj <- FindVariableFeatures(obj)
#obj <- RunPCA(obj, npcs=20)
#obj <- FindNeighbors(obj, k=20)
#nn_graph <- obj@graphs$RNA_nn

#seur_g2m <- which(Wang_Ileum_Phases$seurat$phase == "G2M")
#seur_S <- which(Wang_Ileum_Phases$seurat$phase == "S")
#seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
#seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)

#cyc_g2m <- which(Wang_Ileum_Phases$cyclemix$phase == "G2M")
#cyc_S <- which(Wang_Ileum_Phases$cyclemix$phase == "G1S")
#cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
#cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)

#lone_g2m <- which(Wang_Ileum_Phases$cyclone$phase == "G2M")
#lone_S <- which(Wang_Ileum_Phases$cyclone$phase == "S")
#lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
#lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

#saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Wang_ileum_purity20.rds")


}

if (args[1]==3) {
print("Moreno") ### FAILED ###

set.seed(3671)
Moreno_column <- "celltype_original"
# Moreno_glioblastoma
#obj <- readRDS("Moreno_core_gliobastoma.rds") 
obj <- readRDS("Moreno_SUBSAMPLE_OBJ.rds")
#obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
#obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
#counts <- GetAssayData(obj, slot="counts")
dim(obj)
sum(GetAssayData(obj, slot="counts") == 0)/prod(dim(obj))

#sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
#keep <- !duplicated(sym) & sym != ""
#counts <- counts[keep,]; rownames(counts) <- sym[keep]
#set.seed(998)
#counts <- counts[,sample(1:ncol(counts), size=70000)]
#obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data[match(colnames(counts),rownames(obj@meta.data)),])
##Moreno_Glio_Phases <- run_phase_classifications(obj, species="hsap")
##saveRDS(Moreno_Glio_Phases, file="Moreno_glioblastoma_phases_SUBSAMPLE.rds")

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="hsap")# FAILS ???
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="hsap", type_col=Moreno_column)

#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Moreno_SUBSAMPLE_advphases.rds")
#saveRDS(list(sct=NULL, smoothed=phases_smoothed), "Moreno_SUBSAMPLE_advphases.rds")
#saveRDS(obj, "Moreno_SUBSAMPLE_OBJ.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="hsap")
#saveRDS(phases_genesets, "Moreno_SUBSAMPLE_cyclemix_allgenesets.rds")

good.fit <- goodness_of_fit(obj, plot_prefix="Moreno_SUBSAMPLE", species="hsap")
saveRDS(good.fit, "Moreno_SUBSAMPLE_goodnessoffit.rds")

#Moreno_Glio_Phases <- readRDS("Moreno_glioblastoma_phases_SUBSAMPLE.rds")

## UMAP with labels
#df <- data.frame(UMAP1=obj@meta.data$UMAP1,
#                        UMAP2=obj@meta.data$UMAP2,
#                        Cyclone=Moreno_Glio_Phases$cyclone$phase)
#p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Cyclone))+geom_point(size=0.2)+scale_color_manual(values=phases_col)+theme_classic()+labs(title="Moreno") + guides(colour = guide_legend(override.aes = list(size=1.5)))
#png(paste("Moreno_Cyclone_UMAP_CC.png", sep="_"), width=8/2,height=4, units="in",res=300)
#        print(p1)
#dev.off()

## Neighbourhood Purity ##
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj, features=VariableFeatures(obj))
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=20)
nn_graph <- obj@graphs$RNA_nn

seur_g2m <- which(Moreno_Glio_Phases$seurat$phase == "G2M")
seur_S <- which(Moreno_Glio_Phases$seurat$phase == "S")
seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)

cyc_g2m <- which(Moreno_Glio_Phases$cyclemix$phase == "G2M")
cyc_S <- which(Moreno_Glio_Phases$cyclemix$phase == "G1S")
cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)

lone_g2m <- which(Moreno_Glio_Phases$cyclone$phase == "G2M")
lone_S <- which(Moreno_Glio_Phases$cyclone$phase == "S")
lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Moreno_glioblastoma_purity20_SUBSAMPLE.rds")

}


if (args[1]==4) {
print("Gu")
Gu_column <- "final_annotation"
# Gu Lymph
obj <- readRDS("Gu_lymph.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, slot="counts")
sym <- General_Map(rownames(counts), in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
dim(obj)
sum(GetAssayData(obj, slot="counts") == 0)/prod(dim(obj))

#Gu_Lymph_Phases <- run_phase_classifications(obj, species="mmus")
#saveRDS(Gu_Lymph_Phases, file="Gu_lymph_phases.rds")

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="mmus")
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="mmus", type_col=Gu_column)

##saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Gu_lymph_advphases.rds")
#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Gu_lymph_advphases_allSeurGenes.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="mmus")
#saveRDS(phases_genesets, "Gu_lymph_cyclemix_allgenesets.rds")

good.fit <- goodness_of_fit(obj, plot_prefix="Gu_Lymph", species="mmus")
saveRDS(good.fit, "Gu_lymph_goodnessoffit.rds")

#Gu_Lymph_Phases <- readRDS("Gu_lymph_phases.rds")

## Neighbourhood Purity ##
#obj <- ScaleData(obj)
#obj <- FindVariableFeatures(obj)
#obj <- RunPCA(obj, npcs=20)
#obj <- FindNeighbors(obj, k=20)
#nn_graph <- obj@graphs$RNA_nn

#seur_g2m <- which(Gu_Lymph_Phases$seurat$phase == "G2M")
#seur_S <- which(Gu_Lymph_Phases$seurat$phase == "S")
#seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
#seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)

#cyc_g2m <- which(Gu_Lymph_Phases$cyclemix$phase == "G2M")
#cyc_S <- which(Gu_Lymph_Phases$cyclemix$phase == "G1S")
#cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
#cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)

#lone_g2m <- which(Gu_Lymph_Phases$cyclone$phase == "G2M")
#lone_S <- which(Gu_Lymph_Phases$cyclone$phase == "S")
#lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
#lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

#saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Gu_lymph_purity20.rds")

}

if (args[1]==5) {
print("Triana") 
Triana_column <- "cell_type"
# Triana bonemarrow
obj <- readRDS("Triana_bonemarrow.rds")
obj@meta.data$UMAP1 <- obj@reductions$bothumap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$bothumap@cell.embeddings[,2]
counts <- GetAssayData(obj, slot="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
dim(obj)
sum(GetAssayData(obj, slot="counts") == 0)/prod(dim(obj))

#Triana_Bone_Phases <- run_phase_classifications(obj, species="hsap")
#saveRDS(Triana_Bone_Phases, file="Triana_bonemarrow_phases.rds")

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="hsap")
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="hsap", type_col=Triana_column)

#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Triana_bonemarrow_advphases.rds")
#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Triana_bonemarrow_advphases_allSeurGenes.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="hsap")

#saveRDS(phases_genesets, "Triana_bonemarrow_cyclemix_allgenesets.rds")

good.fit <- goodness_of_fit(obj, plot_prefix="Triana_bonemarrow", species="hsap")
saveRDS(good.fit, "Triana_bonemarrow_goodnessoffit.rds")

#Triana_Bone_Phases=readRDS("Triana_bonemarrow_phases.rds")

## Neighbourhood Purity ##
#obj <- ScaleData(obj)
#obj <- FindVariableFeatures(obj)
#obj <- RunPCA(obj, npcs=20)
#obj <- FindNeighbors(obj, k=20)
#nn_graph <- obj@graphs$RNA_nn

#seur_g2m <- which(Triana_Bone_Phases$seurat$phase == "G2M")
#seur_S <- which(Triana_Bone_Phases$seurat$phase == "S")
#seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
#seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)

#cyc_g2m <- which(Triana_Bone_Phases$cyclemix$phase == "G2M")
#cyc_S <- which(Triana_Bone_Phases$cyclemix$phase == "G1S")
#cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
#cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)

#lone_g2m <- which(Triana_Bone_Phases$cyclone$phase == "G2M")
#lone_S <- which(Triana_Bone_Phases$cyclone$phase == "S")
#lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
#lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

#saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Triana_bonemarrow_purity20.rds")

}

if (args[1]==6) {
print("Steets")
Street_column <- "author_cell_type"
# Streats adipose
obj <- readRDS("Streets_adipose.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, slot="counts") # should be layer?
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
#Streets_adipose_Phases <- run_phase_classifications(obj, species="hsap")
#saveRDS(Streets_adipose_Phases, file="Streets_adipose_phases.rds")

# Revision CycleMix Variations
# scTransform
#phases_sct <- run_phase_classification_sctransform(obj, species="hsap")
#phases_smoothed <- run_phase_classification_cyclemixsmooth(obj, species="hsap", type_col=Street_column)

#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Streets_adipose_advphases.rds")
#saveRDS(list(sct=phases_sct, smoothed=phases_smoothed), "Streets_adipose_advphases_allSeurGenes.rds")

#phases_genesets <- run_phase_classifications_all_gene_sets(obj, species="hsap")
#saveRDS(phases_genesets, "Streets_adipose_cyclemix_allgenesets.rds")

good.fit <- goodness_of_fit(obj, plot_prefix="Streets_adipose", species="hsap")
saveRDS(good.fit, "Streets_adipose_goodnessoffit.rds")


#Streets_adipose_Phases = readRDS("Streets_adipose_phases.rds")
## Neighbourhood Purity ##
#obj <- FindVariableFeatures(obj)
#obj <- ScaleData(obj)
#obj <- RunPCA(obj, npcs=20)
#obj <- FindNeighbors(obj, k=20)
#nn_graph <- obj@graphs$RNA_nn
#
#seur_g2m <- which(Streets_adipose_Phases$seurat$phase == "G2M")
#seur_S <- which(Streets_adipose_Phases$seurat$phase == "S")
#seur_g2m_purity20 <- calc_purity20(seur_g2m, seur_g2m, nn_graph)
#seur_S_purity20 <- calc_purity20(seur_S, seur_S, nn_graph)
#
#cyc_g2m <- which(Streets_adipose_Phases$cyclemix$phase == "G2M")
#cyc_S <- which(Streets_adipose_Phases$cyclemix$phase == "G1S")
#cyc_g2m_purity20 <- calc_purity20(cyc_g2m, cyc_g2m, nn_graph)
#cyc_S_purity20 <- calc_purity20(cyc_S, cyc_S, nn_graph)
#
#lone_g2m <- which(Streets_adipose_Phases$cyclone$phase == "G2M")
#lone_S <- which(Streets_adipose_Phases$cyclone$phase == "S")
#lone_g2m_purity20 <- calc_purity20(lone_g2m, lone_g2m, nn_graph)
#lone_S_purity20 <- calc_purity20(lone_S, lone_S, nn_graph)

#saveRDS(list(seur_g2m_purity20=seur_g2m_purity20,seur_S_purity20=seur_S_purity20, cyc_g2m_purity20=cyc_g2m_purity20, cyc_S_purity20=cyc_S_purity20, lone_g2m_purity20=lone_g2m_purity20, lone_S_purity20=lone_S_purity20), file="Streets_adipose_purity20.rds")

}


# Plots #
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")
# Calculate the proportion of this cell-type that is cycling (S/G2M) for each method.
get_stats <- function(x, cell_type, type_col="cell_type") {
	cells <- x$metadata[,type_col] %in% cell_type;
	if (length(x$cyclemix[[1]]) != length(cells)) {stop("Cyclemix length no match metadata.")}
	if (length(x$seurat[[1]]) != length(cells)) {stop("Seurat length no match metadata.")}

	cm_freq <- table(factor(x$cyclemix[[1]][cells], levels=c("G1S", "G2M", "None")))/sum(cells)
	names(cm_freq) <- c("S", "G2M", "None"); 
	cm_freq <- cm_freq[c("None", "S", "G2M")]

	if (length(x$cyclone[[1]]) != length(cells)) {
		print("Cyclone length no match metadata.")
		cl_freq <- rep(NA, 3)
		names(cl_freq) <- c("None", "S", "G2M")
	} else {
		cl_freq <- table(factor(x$cyclone[[1]][cells], levels=c("G1", "G2M", "S")))/sum(cells)
		names(cl_freq) <- c("None", "G2M", "S");
		cl_freq <- cl_freq[c("None", "S", "G2M")]
	}

	se_freq <- table(factor(x$seurat[[1]][cells], levels=c("G1","S", "G2M")))/sum(cells)
	names(se_freq) <- c("None", "S", "G2M")
	se_freq <- se_freq[c("None", "S", "G2M")]

	tab <- rbind(cm_freq, cl_freq, se_freq)
	rownames(tab) <- c("CycleMix", "Cyclone", "Seurat")
	return(tab[,2]+tab[,3])
}

calc_simpson <- function(vec) {
	freqs <- table(vec)/length(vec)
	return(sum(freqs^2))
}
rel_simpson <- function(vec_subset, vec_pop) {
	exp_simpson <- calc_simpson(vec_pop)
	if (exp_simpson == 1) {stop("Error: all labels identical at population level.")}
	obs_simpson <- calc_simpson(vec_subset)
	return((obs_simpson-exp_simpson)/(1-exp_simpson))

}
get_simpson <- function(x, cell_type, type_col="cell_type") {
	cells <- x$metadata[,type_col] %in% cell_type;
	if (length(x$cyclemix[[1]]) != length(cells)) {stop("Cyclemix length no match metadata.")}
	if (length(x$seurat[[1]]) != length(cells)) {stop("Seurat length no match metadata.")}

	cm_simp <- rel_simpson(x$cyclemix[[1]][cells], x$cyclemix[[1]])
	cl_simp <- rel_simpson(x$cyclone[[1]][cells], x$cyclone[[1]])
	se_simp <- rel_simpson(x$seurat[[1]][cells], x$seurat[[1]])

	tab <- c(cm_simp, cl_simp, se_simp)
	names(tab) <- c("CycleMix", "Cyclone", "Seurat")
	return(tab)
}

if (args[2] == TRUE) {

Wang_is_prolif <- c("transit amplifying cell of colon")
Wang_not_prolif <- c("colon goblet cell", "paneth cell of colon")
Wang_column <- "cell_type"
Wang1 <- readRDS("Wang_colon_phases.rds")
Wang1_stats_prolif <- get_stats(Wang1, Wang_is_prolif)
Wang1_stats_notprolif <- get_stats(Wang1, Wang_not_prolif)
Wang1_simp_prolif <- get_simpson(Wang1, Wang_is_prolif)
Wang1_simp_notprolif <- get_simpson(Wang1, Wang_not_prolif)
Wang2 <- readRDS("Wang_ileum_phases.rds")
Wang2_stats_prolif <- get_stats(Wang2, Wang_is_prolif)
Wang2_stats_notprolif <- get_stats(Wang2, Wang_not_prolif)
Wang2_simp_prolif <- get_simpson(Wang2, Wang_is_prolif)
Wang2_simp_notprolif <- get_simpson(Wang2, Wang_not_prolif)

Street_is_prolif <- c("PRM1+ Cells")
Street_not_prolif <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
Street_column <- "author_cell_type"
Street <- readRDS("Streets_adipose_phases.rds")
Street_stats_prolif <- get_stats(Street, Street_is_prolif, "author_cell_type")
Street_stats_notprolif <- get_stats(Street, Street_not_prolif, "author_cell_type")
Street_simp_prolif <- get_simpson(Street, Street_is_prolif)
Street_simp_notprolif <- get_simpson(Street, Street_not_prolif)

Gu_is_prolif <- c("Plasmablast-IgM+_cycling","Prolif-Tregs", "GC.BC_DZ.1", "GC.BC_DZ.2", "DC.prolif.")
Gu_not_prolif <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
Gu_column <- "final_annotation"
Gu <- readRDS("Gu_lymph_phases.rds")
Gu_stats_prolif <- get_stats(Gu, Gu_is_prolif, "final_annotation")
Gu_stats_notprolif <- get_stats(Gu, Gu_not_prolif, "final_annotation")
Gu_simp_prolif <- get_simpson(Gu, Gu_is_prolif)
Gu_simp_notprolif <- get_simpson(Gu, Gu_not_prolif)

Moreno_is_prolif <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
Moreno_not_prolif <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
Moreno_column <- "celltype_original"
Moreno <- readRDS("Moreno_glioblastoma_phases.rds")
Moreno_stats_prolif <- get_stats(Moreno, Moreno_is_prolif, "celltype_original")
Moreno_stats_notprolif <- get_stats(Moreno, Moreno_not_prolif, "celltype_original")
Moreno_simp_prolif <- get_simpson(Moreno, Moreno_is_prolif)
Moreno_simp_notprolif <- get_simpson(Moreno, Moreno_not_prolif)

Triana_is_prolif <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")
Triana_not_prolif <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")
Triana_column <- "cell_type"
Triana <- readRDS("Triana_bonemarrow_phases.rds")
Triana_stats_prolif <- get_stats(Triana, Triana_is_prolif, "cell_type")
Triana_stats_notprolif <- get_stats(Triana, Triana_not_prolif, "cell_type")
Triana_simp_prolif <- get_simpson(Triana, Triana_is_prolif)
Triana_simp_notprolif <- get_simpson(Triana, Triana_not_prolif)

## Plot - box plots with lines showing diff between prolif & quiesc


df <- data.frame(method=rep(c("CycleMix","Cyclone","Seurat"),times=12),
		dataset=rep(c("Triana","Moreno","Gu","Street","Wang1","Wang2"),each=6),
		type=rep(c("prolif","prolif","prolif","quiesc","quiesc","quiesc"), times=6),
		prolif_prop= c(Triana_stats_prolif, Triana_stats_notprolif, 
				Moreno_stats_prolif, Moreno_stats_notprolif, 
				Gu_stats_prolif, Gu_stats_notprolif, 
				Street_stats_prolif, Street_stats_notprolif,
				Wang1_stats_prolif, Wang1_stats_notprolif,
				Wang2_stats_prolif, Wang2_stats_notprolif)
		)

df$new <- paste(df$method, df$type, sep="-")
df$prolif_prop <- df$prolif_prop*100

require(ggplot2)
p <- ggplot(df, aes(x = new, y = prolif_prop, color = type)) + 
	geom_boxplot(outlier.shape = NA, fill="grey85") + 
	geom_dotplot(aes(fill = type), binaxis = "y", binwidth = 0.05, 
		stackdir = "center", position = position_dodge(0.75), dotsize=50) + 
	geom_line(aes(group = dataset), color = "grey35") + 
	facet_wrap(~method, scales = "free_x") + 
	scale_x_discrete(labels = c("prolif", "quiesc")) + 
	xlab("") + ylab("Cycling Cells (%)")+theme_classic() + guides(color="none", fill="none")
png("Real_Datasets_prolif_vs_not.png", width=8*0.65, height=6*0.65, units="in", res=300)
print(p)
dev.off()

## Repeat the above for the sctransform & smoothed results ##

## Plot - Simpson Index

## Plots - UMAPs of Cyclemix & Seurat phases

source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
require(ggplot2)
require(patchwork)
Wang1 <- readRDS("Wang_colon_phases.rds")
Wang2 <- readRDS("Wang_ileum_phases.rds")
Gu <- readRDS("Gu_lymph_phases.rds")
Street <- readRDS("Streets_adipose_phases.rds")
Moreno <- readRDS("Moreno_glioblastoma_phases.rds")
Triana <- readRDS("Triana_bonemarrow_phases.rds")

all <- list(Wang_colon=Wang1, Wang_ileum=Wang2, Gu=Gu, Streets=Street, Moreno=Moreno, Triana=Triana)
for (dataset in names(all)) {
	df <- data.frame(UMAP1=all[[dataset]]$metadata$UMAP1, 
			UMAP2=all[[dataset]]$metadata$UMAP2,
			CycleMix=all[[dataset]]$cyclemix$phase,
			Seurat=all[[dataset]]$seurat$phase)
	p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=CycleMix))+geom_point(size=0.2)+scale_color_manual(values=phases_col)+theme_classic()+labs(title=dataset) + guides(colour = guide_legend(override.aes = list(size=1.5)))
	p2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Seurat))+geom_point(size=0.2)+scale_color_manual(values=phases_col)+theme_classic()+labs(title=dataset)+ guides(colour = guide_legend(override.aes = list(size=1.5)))
	png(paste(dataset, "UMAP_CC.png", sep="_"), width=8,height=4, units="in",res=300)
	print(p1+p2)
	dev.off()
}

is_prolif <- list(Wang_colon=Wang_is_prolif, Wang_ileum=Wang_is_prolif, Gu=Gu_is_prolif, Streets=Street_is_prolif, Moreno=Moreno_is_prolif, Triana=Triana_is_prolif)
not_prolif <- list(Wang_colon=Wang_not_prolif, Wang_ileum=Wang_not_prolif, Gu=Gu_not_prolif, Streets=Street_not_prolif, Moreno=Moreno_not_prolif, Triana=Triana_not_prolif)
type_column <- list(Wang_colon=Wang_column, Wang_ileum=Wang_column, Gu=Gu_column, Streets=Street_column, Moreno=Moreno_column, Triana=Triana_column)
for (dataset in names(all)) {
	GT_classification = rep("Unknown", nrow(all[[dataset]]$metadata))
	GT_classification[ (all[[dataset]]$metadata[,type_column[[dataset]]]) %in% is_prolif[[dataset]] ] <- "Proliferative"
	GT_classification[ (all[[dataset]]$metadata[,type_column[[dataset]]]) %in% not_prolif[[dataset]] ] <- "Quiescent"
	df <- data.frame(UMAP1=all[[dataset]]$metadata$UMAP1, 
			UMAP2=all[[dataset]]$metadata$UMAP2,
			GT=GT_classification)
	GT_colours <- c("Proliferative"="purple", "Quiescent"="forestgreen", "Unknown"="grey50")
	p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=GT))+geom_point(size=0.2)+scale_color_manual(values=GT_colours)+theme_classic()+labs(title=dataset) + guides(colour = guide_legend(override.aes = list(size=1.5)))
	png(paste(dataset, "UMAP_GT.png", sep="_"), width=8/2,height=4, units="in",res=300)
	print(p1)
	dev.off()
}

### proliferation neighbourhood purity ###
files <- Sys.glob("*_purity20.rds")
out <- c()
for (f in files) {
	samp_name <- unlist(strsplit(f, "_"))
	samp_name <- paste(samp_name[1], samp_name[2], sep="_")
	stuff <- readRDS(f)
	seur <- cbind(c("Seurat", "Seurat"),c("G1S", "G2M"), c(median(stuff$seur_S_purity20), median(stuff$seur_g2m_purity20)), c(samp_name, samp_name))
	cyc <- cbind(c("CycleMix", "CycleMix"), c("G1S", "G2M"), c(median(stuff$cyc_S_purity20), median(stuff$cyc_g2m_purity20)), c(samp_name, samp_name))
	lone <- cbind(c("Cyclone", "Cyclone"),c("G1S", "G2M"), c(median(stuff$lone_S_purity20), median(stuff$lone_g2m_purity20)), c(samp_name, samp_name))
	out <- rbind(out, seur, cyc, lone)
	seur_test = wilcox.test(stuff$seur_g2m_purity20, stuff$seur_S_purity20)
	cyc_test = wilcox.test(stuff$cyc_g2m_purity20, stuff$cyc_S_purity20)
	lone_test = wilcox.test(stuff$lone_g2m_purity20, stuff$lone_S_purity20)
	print(seur_test$p.value)
	print(cyc_test$p.value)
	print(lone_test$p.value)
}
colnames(out) <- c("method", "phase", "purity", "dataset")
df <- data.frame(out)
df$new <- paste(df$method, df$phase, sep="-")
df$purity <- as.numeric(df$purity)*100


require(ggplot2)
p <- ggplot(df, aes(x = new, y = purity, color = phase)) + 
	geom_boxplot(outlier.shape = NA, fill="grey85") + 
	geom_dotplot(aes(fill = dataset), binaxis = "y", binwidth = 0.05, 
		stackdir = "center", position = position_dodge(0.75), dotsize=50) + 
	facet_wrap(~method, scales = "free_x") + 
	scale_x_discrete(labels = c("S", "G2M")) + 
	scale_color_manual(values=phases_col) +
	xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none", fill="none")
png("Real_Datasets_purity20.png", width=8*0.65, height=6*0.65, units="in", res=300)
print(p)
dev.off()


}


