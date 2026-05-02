source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
require("Seurat")
require("CycleMix")
require("mclust")
require("Matrix")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")

source("~/My_R_Packages/CycleMix/R/Analysis.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

set.seed(3891)

run_phase_classifications <- function(seur_obj, species=c("hsap", "mmus")) {
        if (species == "hsap") {
                s.genes <- cc.genes$s.genes
                g2m.genes <- cc.genes$g2m.genes
                gene_table <- HGeneSets$Tirosh
        } else if (species == "mmus") {
                gene_table <- MGeneSets$Cyclone
                s.genes <- gene_table[(gene_table[,2] == "S" | gene_table[,2] == "G1S") & gene_table[,3] == 1,1]
                g2m.genes <- gene_table[(gene_table[,2] == "G2M" | gene_table[,2] == "M" | gene_table[,2] == "G2") & gene_table[,3] == 1,1]
        }
        seur_obj <- CellCycleScoring(seur_obj, s.features=s.genes, g2m.features=g2m.genes)
}


# Calculates the proportion of neighbours of each cell in group1 that are in group2 based on the graph
calc_purity <- function(is.group1, is.group2, graph) {
	out <- c()
	for (cell in is.group1) {
		nns <- which(graph[cell,]>0)
		score <- sum(nns %in% is.group2)/length(nns)
		out <- c(out, score)
	}
	return(out)
}
calc_purity2 <- function(is.group1, is.group2, graph) {
        out <- c()
	if (length(is.group1) > 1000) {
		is.group1 <- sample(is.group1, 1000)
	}
		
        for (cell in is.group1) {
                nns <- which(graph[cell,]>0)
                score <- sum(nns %in% is.group2)/length(nns)
                out <- c(out, score)
        }
        return(out)
}

calc_purity_cc <- function(Phases, nn_graph, method){
	seur_g2m <- which(Phases[[method]]$phase == "G2M")
	seur_S <- which(Phases[[method]]$phase == "S" | Phases[[method]]$phase == "G1S")
	seur_g2m_purity <- calc_purity2(seur_g2m, seur_g2m, nn_graph)
#	seur_g2m_prolif_purity <- calc_purity(seur_g2m, c(seur_S,seur_g2m), nn_graph)
	seur_prolif_purity <- calc_purity2(c(seur_S,seur_g2m), c(seur_S,seur_g2m), nn_graph)

	return(list(purity=seur_g2m_purity, prolif_purity=seur_prolif_purity))
}

run_CycleMixRegression <- function(obj, phases_column) {
	norm_expr <- GetAssayData(obj, layer="data", assay="RNA")
	out <- regressCyclePartial(norm_expr, obj@meta.data[,phases_column], allow_negative=FALSE, phases=c("G2M", "G1S"), type="norm", method="phase",subsample_cells=5000)
	return(out)
}

run_regressions <- function(obj) {
	# Subset to just the variable features
	obj <- obj[VariableFeatures(obj),]
	# Run SCTransform
	seur.elapsed.time <- system.time(obj <- SCTransform(obj, residual.features=VariableFeatures(obj), do.correct.umi=FALSE, vars.to.regress=c("S.Score", "G2M.Score")))
	obj <- RunPCA(obj, npcs=20)
	obj <- FindNeighbors(obj, k=50, dims=1:20, assay="SCT")
	seur_nn_graph <- obj@graphs$SCT_nn
	obj <- RunUMAP(obj, dims=1:20)
#	print(DimPlot(obj, reduction="umap", group.by="seurat"))
#	print(DimPlot(obj, reduction="umap", group.by="cyclemix"))
	seur_umap <- obj@reductions$umap@cell.embeddings

	# Run CycleMix
	cyclemix.elapsed.time <- system.time(obj@assays$RNA@layers$data<- run_CycleMixRegression(obj, "cyclemix"))
	obj <- ScaleData(obj, assay="RNA")
	obj <- RunPCA(obj, npcs=20, assay="RNA")
	obj <- FindNeighbors(obj, k=50, dims=1:20, assay="RNA")
	cyclemix_nn_graph <- obj@graphs$RNA_nn
	obj <- RunUMAP(obj, dims=1:20)
#	print(DimPlot(obj, reduction="umap", group.by="seurat"))
#	print(DimPlot(obj, reduction="umap", group.by="cyclemix"))
	cyclemix_umap <- obj@reductions$umap@cell.embeddings
	return(list(seur_nn_graph=seur_nn_graph, seur_umap=seur_umap, cyclemix_nn_graph=cyclemix_nn_graph, cyclemix_umap=cyclemix_umap, seur.elapsed.time=seur.elapsed.time, cyclemix.elapsed.time=cyclemix.elapsed.time))
}

args <- commandArgs(trailingOnly=TRUE)

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
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Wang_Colon_Phases<- readRDS("Wang_colon_phases.rds")
obj@meta.data$seurat <- Wang_Colon_Phases$seurat$phase
obj@meta.data$cyclemix <- Wang_Colon_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Wang_uncorrected_graph <- obj@graphs$RNA_nn
Wang_uncorrected_UMAP <- obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Wang_uncorrected_graph
corrected[["orig_umap"]] <- Wang_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$cell_type
saveRDS(corrected, "Wang_colon_corrected.rds")
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
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Wang_Ileum_Phases <- readRDS(file="Wang_ileum_phases.rds")
obj@meta.data$seurat <- Wang_Ileum_Phases$seurat$phase
obj@meta.data$cyclemix <- Wang_Ileum_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Wang_uncorrected_graph <- obj@graphs$RNA_nn
Wang_uncorrected_UMAP <-  obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Wang_uncorrected_graph
corrected[["orig_umap"]] <- Wang_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$cell_type
saveRDS(corrected, "Wang_ileum_corrected.rds")

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
saveRDS(counts, "test.rds")
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")
saveRDS(obj, "test2.rds")

Moreno_Glio_Phases<- readRDS("Moreno_glioblastoma_phases.rds")
obj@meta.data$seurat <- Moreno_Glio_Phases$seurat$phase
obj@meta.data$cyclemix <- Moreno_Glio_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- obj[VariableFeatures(obj),]
obj <- ScaleData(obj, features=VariableFeatures(obj))
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Moreno_uncorrected_graph <- obj@graphs$RNA_nn
Moreno_uncorrected_UMAP <-  obj@meta.data[,c("UMAP1","UMAP2")]
saveRDS(obj, "test3.rds")

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Moreno_uncorrected_graph
corrected[["orig_umap"]] <- Moreno_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$celltype_original
saveRDS(corrected, "Moreno_glioblastoma_corrected.rds")
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
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="mmus")

Gu_Lymph_Phases<- readRDS("Gu_lymph_phases.rds")
obj@meta.data$seurat <- Gu_Lymph_Phases$seurat$phase
obj@meta.data$cyclemix <- as.character(Gu_Lymph_Phases$cyclemix$phase)
obj@meta.data$cyclemix[obj@meta.data$cyclemix=="S"] <- "G1S"
obj@meta.data$cyclemix[obj@meta.data$cyclemix=="G1"] <- "None"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Gu_uncorrected_graph <- obj@graphs$RNA_nn
Gu_uncorrected_UMAP <-  obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Gu_uncorrected_graph
corrected[["orig_umap"]] <- Gu_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$final_annotation
saveRDS(corrected, "Gu_lymph_corrected.rds")
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
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Triana_Marrow_Phases<- readRDS("Triana_bonemarrow_phases.rds")
obj@meta.data$seurat <- Triana_Marrow_Phases$seurat$phase
obj@meta.data$cyclemix <- Triana_Marrow_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Triana_uncorrected_graph <- obj@graphs$RNA_nn
Triana_uncorrected_UMAP <-  obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Triana_uncorrected_graph
corrected[["orig_umap"]] <- Triana_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$cell_type
saveRDS(corrected, "Triana_bonemarrow_corrected.rds")
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
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Streets_Adipose_Phases<- readRDS("Streets_adipose_phases.rds")
obj@meta.data$seurat <- Streets_Adipose_Phases$seurat$phase
obj@meta.data$cyclemix <- Streets_Adipose_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Streets_uncorrected_graph <- obj@graphs$RNA_nn
Streets_uncorrected_UMAP <-  obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- run_regressions(obj)
corrected[["orig_nn_graph"]] <- Streets_uncorrected_graph
corrected[["orig_umap"]] <- Streets_uncorrected_UMAP
corrected[["cell_type"]] <- obj@meta.data$author_cell_type
saveRDS(corrected, "Streets_adipose_corrected.rds")
}


# Plots #
source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
# Calculate the proportion of this cell-type that is cycling (S/G2M) for each method.


if (args[2] == 1) {
	print("Purity of cycling cells")

	Wang1 <- readRDS("Wang_colon_phases.rds")
	Wang2 <- readRDS("Wang_ileum_phases.rds")
	Gu <- readRDS("Gu_lymph_phases.rds")
	Street <- readRDS("Streets_adipose_phases.rds")
	Moreno <- readRDS("Moreno_glioblastoma_phases.rds")
	Triana <- readRDS("Triana_bonemarrow_phases.rds")

	Wang1_corrected <- readRDS("Wang_colon_corrected.rds")
	Wang2_corrected <- readRDS("Wang_ileum_corrected.rds")
	Gu_corrected <- readRDS("Gu_lymph_corrected.rds")
	Moreno_corrected <- readRDS("Moreno_glioblastoma_corrected.rds")
	Triana_corrected <- readRDS("Triana_bonemarrow_corrected.rds")
	Street_corrected <- readRDS("Streets_adipose_corrected.rds")

	# Purity of cycling cells
	calc_df <- function(Phases, Corrected, dataset) {
		purity_orig_mix <- calc_purity_cc(Phases, Corrected$orig_nn_graph, "cyclemix")
		purity_orig_mix2 <- calc_purity_cc(Phases, Corrected$orig_nn_graph, "seurat")
		purity_corr_mix <- calc_purity_cc(Phases, Corrected$cyclemix_nn_graph, "cyclemix")
		purity_corr_mix2 <- calc_purity_cc(Phases, Corrected$cyclemix_nn_graph, "seurat")
		purity_corr_seur2 <- calc_purity_cc(Phases, Corrected$seur_nn_graph, "seurat")
		purity_corr_seur <- calc_purity_cc(Phases, Corrected$seur_nn_graph, "cyclemix")
		df <- data.frame(labels=rep(c("CycleMix", "Seurat"), time=6), 
			correction=rep(c("None", "CycleMix", "Seurat"), each=4), 
			type=rep(c("G2M","G2M", "Prolif", "Prolif"), times=3), 
			dataset=rep(dataset, times=12),
			purity=c(mean(purity_orig_mix[[1]]), mean(purity_orig_mix2[[1]]),
				mean(purity_orig_mix[[2]]), mean(purity_orig_mix2[[2]]),
				mean(purity_corr_mix[[1]]), mean(purity_corr_mix2[[1]]),
				mean(purity_corr_mix[[2]]), mean(purity_corr_mix2[[2]]),
				mean(purity_corr_seur[[1]]), mean(purity_corr_seur2[[1]]),
				mean(purity_corr_seur[[2]]), mean(purity_corr_seur2[[2]])),
			purity_stderr=c(sd(purity_orig_mix[[1]])/sqrt(length(purity_orig_mix[[1]])),
					sd(purity_orig_mix2[[1]])/sqrt(length(purity_orig_mix2[[1]])),
					sd(purity_orig_mix[[2]])/sqrt(length(purity_orig_mix[[2]])),
					sd(purity_orig_mix2[[2]])/sqrt(length(purity_orig_mix2[[2]])),
					sd(purity_corr_mix[[1]])/sqrt(length(purity_corr_mix[[1]])),
					sd(purity_corr_mix2[[1]])/sqrt(length(purity_corr_mix2[[1]])),
					sd(purity_corr_mix[[2]])/sqrt(length(purity_corr_mix[[2]])),
					sd(purity_corr_mix2[[2]])/sqrt(length(purity_corr_mix2[[2]])),
					sd(purity_corr_seur[[1]])/sqrt(length(purity_corr_seur[[1]])),
					sd(purity_corr_seur2[[1]])/sqrt(length(purity_corr_seur2[[1]])),
					sd(purity_corr_seur[[2]])/sqrt(length(purity_corr_seur[[2]])),
					sd(purity_corr_seur2[[2]])/sqrt(length(purity_corr_seur2[[2]]))
			))
		return(df)
	}
	all_df <- rbind(calc_df(Wang1, Wang1_corrected, "Wang_colon"),
			calc_df(Wang2, Wang2_corrected, "Wang_ileum"),
			calc_df(Gu, Gu_corrected, "Gu_lymph"),
			calc_df(Street, Street_corrected, "Streets_adipose"),
			calc_df(Moreno, Moreno_corrected, "Moreno_glio"),
			calc_df(Triana, Triana_corrected, "Triana_bone"))

	saveRDS(all_df, "all_df_regression.rds")

	paired_ttest <- function(correction, labels, type) {
		x <- all_df[all_df$correction == "None" & all_df$labels==labels & all_df$type==type,"purity"]
		y <- all_df[all_df$correction == correction & all_df$labels==labels & all_df$type==type,"purity"]
		return(t.test(x, y, alternative="two.sided", paired=TRUE))
	}
	# y = all_df[all_df$correction == correction & all_df$type==type,"purity"]
	# x = all_df[all_df$correction == "None" & all_df$type==type,"purity"]
	# t.test(x,y, alternative="two.sided", paired=TRUE)

	all_df$new <- paste(all_df$type, all_df$correction)
	all_df$purity <- as.numeric(all_df$purity)*100
	all_df$correction <- factor(all_df$correction, levels=c("None", "Seurat", "CycleMix"))
	p <- ggplot(all_df, aes(x = new, y = purity, color = labels)) +
		geom_boxplot(outlier.shape = NA, fill="grey85") +
		geom_dotplot(aes(fill = labels), binaxis = "y", binwidth = 0.05,
		stackdir = "center", position = position_dodge(0.75), dotsize=50) + 
	facet_wrap(~correction, scales = "free_x") + 
	scale_x_discrete(labels = c("G2M", "Prolif")) + 
	xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none")

	p1 <- ggplot(all_df[all_df$type == "G2M",], aes(x = correction, y = purity, color = labels)) +
		geom_boxplot(outlier.shape = NA, fill="grey85") +
		geom_dotplot(aes(fill = labels), binaxis = "y", binwidth = 0.05,
		stackdir = "center", position = position_dodge(0.75), dotsize=50) +
		facet_wrap(~correction, scales = "free_x") +
		scale_x_discrete(labels = c("G2M")) +
	        xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none")

	p2 <- ggplot(all_df[all_df$type == "Prolif",], aes(x = correction, y = purity, color = labels)) +
		geom_boxplot(outlier.shape = NA, fill="grey85") +
		geom_dotplot(aes(fill = labels), binaxis = "y", binwidth = 0.05,
		stackdir = "center", position = position_dodge(0.75), dotsize=50) +
		facet_wrap(~correction, scales = "free_x") +
		scale_x_discrete(labels = c("Prolif")) +
	        xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none")

	refined_df1 <- all_df[all_df$labels=="CycleMix" & all_df$correction %in% c("None", "CycleMix"),]
	refined_df2 <- all_df[all_df$labels=="Seurat" & all_df$correction %in% c("None", "Seurat"),]
	p3 <- ggplot(refined_df1, aes(x = correction, y = purity, color = type)) + geom_boxplot(outlier.shape=NA, fill="grey85") +
		geom_dotplot(aes(fill = type), binaxis = "y", binwidth = 0.05,
		stackdir = "center", position = position_dodge(0.75), dotsize=20) +
		facet_wrap(~correction, scales = "free_x") +
	        xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none")
		
	p4 <- ggplot(refined_df2, aes(x = correction, y = purity, color = type)) + geom_boxplot(outlier.shape=NA, fill="grey85") +
		geom_dotplot(aes(fill = type), binaxis = "y", binwidth = 0.05,
		stackdir = "center", position = position_dodge(0.75), dotsize=20) +
		facet_wrap(~correction, scales = "free_x") +
	        xlab("") + ylab("Purity (%)")+theme_classic() + guides(color="none")


	png("Real_Datasets_regression_purity_G2M.png", width=6*0.65, height=6*0.65, units="in", res=300)
	print(p1)
	dev.off()
			
	png("Real_Datasets_regression_purity_Prolif.png", width=6*0.65, height=6*0.65, units="in", res=300)
	print(p2)
	dev.off()

	png("Real_Datasets_regression_purity_CycleMix.png", width=6*0.65, height=6*0.65, units="in", res=300)
	print(p3)
	dev.off()
			
	png("Real_Datasets_regression_purity_Seurat.png", width=6*0.65, height=6*0.65, units="in", res=300)
	print(p4)
	dev.off()


}

if (args[2]==2){
	print("Purity of proliferative cell-types")


# Purity of proliferative cell-type
	source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")

	Wang_is_prolif <- c("transit amplifying cell of colon")
	Street_is_prolif <- c("PRM1+ Cells")
	Gu_is_prolif <- c("GC.BC_DZ.1", "GC.BC_DZ.2")
	Moreno_is_prolif <- c("Tumour", "Radial Glia","Neoplastic")
	Triana_is_prolif <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")

	Wang1_corrected <- readRDS("Wang_colon_corrected.rds")
	Wang2_corrected <- readRDS("Wang_ileum_corrected.rds")
	Gu_corrected <- readRDS("Gu_lymph_corrected.rds")
	Moreno_corrected <- readRDS("Moreno_glioblastoma_corrected.rds")
	Triana_corrected <- readRDS("Triana_bonemarrow_corrected.rds")
	Street_corrected <- readRDS("Streets_adipose_corrected.rds")


	type_purity <- function (corrected, is_prolif) {
		orig <- c()
		seurat <- c()
		cyclemix <- c()
		for (t in is_prolif) {
			orig <- c(orig, calc_purity2(which(corrected$cell_type == t), which(corrected$cell_type == t),corrected$orig_nn_graph))
			seurat <- c(seurat,  calc_purity2(which(corrected$cell_type == t), which(corrected$cell_type == t),corrected$seur_nn_graph))
			cyclemix <- c(cyclemix,  calc_purity2(which(corrected$cell_type == t), which(corrected$cell_type == t),corrected$cyclemix_nn_graph))
		}
		scores <- c(mean(orig), mean(seurat), mean(cyclemix))
		stderr <- c(sd(orig)/sqrt(length(orig)), sd(seurat)/sqrt(length(seurat)), sd(cyclemix)/sqrt(length(cyclemix)))
		return(list(scores=scores, stderr=stderr))
	}

	##### I need Error Bars!!! #####
	#Wang1_scores <- type_purity(Wang1_corrected, Wang_is_prolif) # Missing Seur nn
	#Wang2_scores <- type_purity(Wang2_corrected, Wang_is_prolif) # Missing Seur nn
	#Gu_scores <- type_purity(Gu_corrected, Gu_is_prolif) # wrong cell-type labels
	#Triana_scores <- type_purity(Triana_corrected, Triana_is_prolif) # Missing Seur nn
	#Street_scores <- type_purity(Street_corrected, Street_is_prolif) # wrong cell-type labels
	#Moreno_scores <- type_purity(Moreno_corrected, Moreno_is_prolif) # wrong cell-type lavels
	#saveRDS(list(Wang1_scores=Wang1_scores, Wang2_scores=Wang2_scores, Gu_scores=Gu_scores, Triana_scores=Triana_scores, Street_scores=Street_scores, Moreno_scores=Moreno_scores), "regression_type_purity_scores.rds")
	stuff <- readRDS("regression_type_purity_scores.rds")
	Wang1_scores = stuff[["Wang1_scores"]]
	Wang2_scores = stuff[["Wang2_scores"]]
	Gu_scores = stuff[["Gu_scores"]]
	Triana_scores = stuff[["Triana_scores"]]
	Street_scores = stuff[["Street_scores"]]
	Moreno_scores = stuff[["Moreno_scores"]]
	df <- data.frame(dataset=rep(c("Wang1", "Wang2", "Gu", "Triana", "Street", "Moreno"), each=3), 
			correction=rep(c("None", "Seurat", "CycleMix"), times=6), 
			purity = c(Wang1_scores$scores, Wang2_scores$scores, Gu_scores$scores, Triana_scores$scores, Street_scores$scores, Moreno_scores$scores), 
			purity_stderr=c(Wang1_scores$stderr, Wang2_scores$stderr, Gu_scores$stderr, Triana_scores$stderr, Street_scores$stderr, Moreno_scores$stderr))

	colours<-c(colours, c("None"="grey50"))
	df$correction <- factor(df$correction, levels=c("None", "Seurat", "CycleMix"))
	png("Regression_type_purity.png", width=8*0.75, height=4*0.75, units="in", res=300)
	print(ggplot(df, aes(x = dataset, y = purity, fill=correction))+geom_bar(stat="identity", position="dodge") +scale_fill_manual(values=colours) +  xlab("") + ylab("Purity (%)")+theme_classic()) + geom_errorbar(aes(ymin=purity-2*purity_stderr, ymax=purity+2*purity_stderr), width=0.2, position=position_dodge(0.9))
	dev.off()


# Run Time
	source("/home/tandrew6/scripts/CycleMix_Benchmark/ColourScheme.R")
	df <- data.frame(dataset=rep(c("Wang1", "Wang2", "Gu", "Triana", "Street", "Moreno"), each=2), 
		correction=rep(c("Seurat", "CycleMix"), times=6),
		run_time = c(Wang1_corrected$seur.elapsed.time[3], Wang1_corrected$cyclemix.elapsed.time[3], 
				Wang2_corrected$seur.elapsed.time[3], Wang2_corrected$cyclemix.elapsed.time[3], 
				Gu_corrected$seur.elapsed.time[3], Gu_corrected$cyclemix.elapsed.time[3], 
				Triana_corrected$seur.elapsed.time[3], Triana_corrected$cyclemix.elapsed.time[3], 
				Street_corrected$seur.elapsed.time[3], Street_corrected$cyclemix.elapsed.time[3],
				Moreno_corrected$seur.elapsed.time[3], Moreno_corrected$cyclemix.elapsed.time[3]))
        png("Regression_runtime.png", width=4, height=4, units="in", res=300)
	print(ggplot(df, aes(x = dataset, y = run_time, fill=correction))+geom_bar(stat="identity", position="dodge") +scale_fill_manual(values=colours) +  xlab("") + ylab("Run Time (sec)")+theme_classic()) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        dev.off()
}

if (args[2]==3){
# UMAPs

	phases <- list()
        phases[["Wang1"]] <- readRDS("Wang_colon_phases.rds")
        phases[["Wang2"]] <- readRDS("Wang_ileum_phases.rds")
        phases[["Gu"]] <- readRDS("Gu_lymph_phases.rds")
        phases[["Streets"]] <- readRDS("Streets_adipose_phases.rds")
        phases[["Moreno"]] <- readRDS("Moreno_glioblastoma_phases.rds")
        phases[["Triana"]] <- readRDS("Triana_bonemarrow_phases.rds")
	
	corrected <- list()
        corrected[["Wang1"]] <- readRDS("Wang_colon_corrected.rds")
        corrected[["Wang2"]] <- readRDS("Wang_ileum_corrected.rds")
        corrected[["Gu"]] <- readRDS("Gu_lymph_corrected.rds")
        corrected[["Moreno"]] <- readRDS("Moreno_glioblastoma_corrected.rds")
        corrected[["Triana"]] <- readRDS("Triana_bonemarrow_corrected.rds")
        corrected[["Streets"]] <- readRDS("Streets_adipose_corrected.rds")

	prolif <- list()
	prolif[["Wang1"]] <-  c("transit amplifying cell of colon")
	prolif[["Wang2"]] <-  c("transit amplifying cell of colon")
	prolif[["Street"]] <-  c("PRM1+ Cells")
	prolif[["Gu"]] <- c("GC.BC_DZ.1", "GC.BC_DZ.2")
	prolif[["Moreno"]] <- c("Tumour", "Radial Glia","Neoplastic")
	prolif[["Triana"]] <-  c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")

	not_prolif <- list()
	not_prolif[["Wang1"]] <- c("colon goblet cell", "paneth cell of colon")
	not_prolif[["Wang2"]] <- c("colon goblet cell", "paneth cell of colon")
	not_prolif[["Streets"]] <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
	not_prolif[["Gu"]] <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
	not_prolif[["Moreno"]] <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
	not_prolif[["Triana"]] <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")



for (dataset in names(corrected)) {
	# Add whether it is /isn't a prolif cell-type
	is.prolif <- corrected[[dataset]]$cell_type %in% prolif[[dataset]]
	is.notprolif <- corrected[[dataset]]$cell_type %in% not_prolif[[dataset]]
	gt <- rep("Unknown", length(is.prolif))
	gt[is.prolif] <- "Proliferative"	
	gt[is.notprolif] <- "Quiescent"	

        df <- data.frame(UMAP1=corrected[[dataset]]$cyclemix_umap[,1],
                        UMAP2=corrected[[dataset]]$cyclemix_umap[,2],
                        CycleMix=phases[[dataset]]$cyclemix$phase,
                        Seurat=phases[[dataset]]$seurat$phase, GroundTruth=gt)
        p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=CycleMix))+geom_point(size=0.2)+scale_color_manual(values=phases_col)+theme_classic()+labs(title=dataset) + guides(colour = guide_legend(override.aes = list(size=1.5)))
        p2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=Seurat))+geom_point(size=0.2)+scale_color_manual(values=phases_col)+theme_classic()+labs(title=dataset)+ guides(colour = guide_legend(override.aes = list(size=1.5)))
        p3 <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=GroundTruth))+geom_point(size=0.2)+scale_color_manual(values=gt_col)+theme_classic()+labs(title=dataset)+ guides(colour = guide_legend(override.aes = list(size=1.5)))
        png(paste(dataset, "cyclemix_corrected_UMAP.png", sep="_"), width=12,height=4, units="in",res=300)
        print(p1+p2+p3)
        dev.off()
}


}

if (args[3]) {
	# Biological Analysis of corrected data for novel cell subtypes
	# Triana
# Wang Colon
obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Wang_Colon_Phases<- readRDS("Wang_colon_phases.rds")
obj@meta.data$seurat <- Wang_Colon_Phases$seurat$phase
obj@meta.data$cyclemix <- Wang_Colon_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Wang_uncorrected_graph <- obj@graphs$RNA_nn
Wang_uncorrected_UMAP <- obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- readRDS("Wang_colon_corrected.rds")
obj@graphs$CycleMix <- corrected$cyclemix_nn_graph
obj@graphs$SeurReg <- corrected$seur_nn_graph
obj@graphs$UnCor <- corrected$orig_nn_graph
obj@graphs$SeurReg@assay.used <- "RNA"

require(igraph)
obj <- FindClusters(obj, graph.name="CycleMix")
obj@meta.data$clusters_post_cyclemix <- Idents(obj)
obj <- FindClusters(obj, graph.name="SeurReg")
obj@meta.data$clusters_post_sct <- Idents(obj)
obj <- FindClusters(obj, graph.name="UnCor")
obj@meta.data$clusters_uncor <- Idents(obj)

source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")

png("Wang_colon_regression_corrected_clusters.png", width=8, height=6, units="in", res=300)
par(mfcol=c(2,3))

# Uncorrected
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("top", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=TRUE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

# CycleMix
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("top", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=TRUE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

# Seurat
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("top", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=TRUE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

dev.off()


## Wang Ileum
# Wang Colon
obj <- readRDS("Wang_ileum.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Wang_Colon_Phases<- readRDS("Wang_ileum_phases.rds")
obj@meta.data$seurat <- Wang_Colon_Phases$seurat$phase
obj@meta.data$cyclemix <- Wang_Colon_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Wang_uncorrected_graph <- obj@graphs$RNA_nn
Wang_uncorrected_UMAP <- obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- readRDS("Wang_ileum_corrected.rds")
obj@graphs$CycleMix <- corrected$cyclemix_nn_graph
obj@graphs$SeurReg <- corrected$seur_nn_graph
obj@graphs$UnCor <- corrected$orig_nn_graph
obj@graphs$SeurReg@assay.used <- "RNA"

# cell-type lisi
require(lisi)
tmp2 <- compute_lisi(corrected$seur_umap, obj@meta.data, "cell_type")
tmp3 <- compute_lisi(corrected$orig_umap, obj@meta.data, "cell_type")
tmp1 <- compute_lisi(corrected$cyclemix_umap, obj@meta.data, "cell_type")
res <- wilcox.test(tmp1, tmp2)

df <- data.frame(correction=c("None", "CycleMix", "Seurat"), 
		meanLISI = c(mean(tmp3[,1]), mean(tmp1[,1]), mean(tmp2[,1])), 
		medianLISI=c(median(tmp3[,1]), median(tmp1[,1]), median(tmp2[,1]))
p.vals <- c("CycleMix_vs_Seurat"=wilcox.test(tmp1[,1], tmp2[,1])$p.value, 
	    "CycleMix_vs_None"=wilcox.test(tmp1[,1], tmp3[,1])$p.value,
	    "Seurat_vs_None"=wilcox.test(tmp2[,1], tmp3[,1])$p.value)
saveRDS(list(lisi=df, wilcox.pval=p.vals), "Wang_colon_correction_LISI.rds")




require(igraph)
obj <- FindClusters(obj, graph.name="CycleMix")
obj@meta.data$clusters_post_cyclemix <- Idents(obj)
obj <- FindClusters(obj, graph.name="SeurReg")
obj@meta.data$clusters_post_sct <- Idents(obj)
obj <- FindClusters(obj, graph.name="UnCor")
obj@meta.data$clusters_uncor <- Idents(obj)

source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")

png("Wang_ileum_regression_corrected_clusters.png", width=8, height=6, units="in", res=300)
par(mfcol=c(2,3))

# Uncorrected
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topright", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

# CycleMix
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topright", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

# Seurat
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topright", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("transit amplifying cell of colon")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="Transit Amplifying Cells", ylim=c(0,1))

dev.off()


## Streets
obj <- readRDS("Streets_adipose.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
obj <- NormalizeData(obj)
obj <- run_phase_classifications(obj, species="hsap")

Streets_Colon_Phases<- readRDS("Streets_adipose_phases.rds")
obj@meta.data$seurat <- Streets_Colon_Phases$seurat$phase
obj@meta.data$cyclemix <- Streets_Colon_Phases$cyclemix$phase
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
obj <- FindNeighbors(obj, k=50)
Streets_uncorrected_graph <- obj@graphs$RNA_nn
Streets_uncorrected_UMAP <- obj@meta.data[,c("UMAP1","UMAP2")]

corrected <- readRDS("Streets_adipose_corrected.rds")
obj@graphs$CycleMix <- corrected$cyclemix_nn_graph
obj@graphs$SeurReg <- corrected$seur_nn_graph
obj@graphs$UnCor <- corrected$orig_nn_graph
obj@graphs$SeurReg@assay.used <- "RNA"

require(igraph)
obj <- FindClusters(obj, graph.name="CycleMix")
obj@meta.data$clusters_post_cyclemix <- Idents(obj)
obj <- FindClusters(obj, graph.name="SeurReg")
obj@meta.data$clusters_post_sct <- Idents(obj)
obj <- FindClusters(obj, graph.name="UnCor")
obj@meta.data$clusters_uncor <- Idents(obj)

source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")

png("Streets_adipose_regression_corrected_clusters.png", width=8, height=6, units="in", res=300)
par(mfcol=c(2,3))

# Uncorrected
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topleft", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_uncor, obj@meta.data$author_cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("PRM1+ Cells")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="PRM1+ Cells", ylim=c(0,1))

# CycleMix
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topleft", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_cyclemix, obj@meta.data$author_cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("PRM1+ Cells")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="PRM1+ Cells", ylim=c(0,1))

# Seurat
# Plot Clusters vs CellCyclePhase
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$cyclemix)
tab <- tab/rowSums(tab)
tab <- tab[,c("G1S", "G2M")]
barplot(t(tab), col=c(phases_col["G1S"], phases_col["G2M"]), xlab="Cluster", ylab="Cell Proportion", main="Cell Cycle Phase", ylim=c(0,1))
legend("topleft", c("G1S", "G2M"), fill=c(phases_col["G1S"], phases_col["G2M"]), bty="n", horiz=FALSE)

# Plot Clusters vs Proliferating Celltype
tab <- table(obj@meta.data$clusters_post_sct, obj@meta.data$author_cell_type)
tab <- tab/rowSums(tab)
barplot(t(tab[,c("PRM1+ Cells")]), col=c("black"), xlab="Cluster", ylab="Cell Proportion", main="PRM1+ Cells", ylim=c(0,1))

dev.off()

}
