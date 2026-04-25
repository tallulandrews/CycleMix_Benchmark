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

run_phase_classifications <- function(seur_obj, species=c("hsap", "mmus")) {
	seur_obj <- NormalizeData(seur_obj)
	if (species == "hsap") {
		gene_table <- HGeneSets$Tirosh
		gene_table <- rbind(gene_table, HGeneSets$Quiesc)
	} else if (species == "mmus") {
		gene_table <- MGeneSets$Cyclone
		quesc <- HGeneSets$Quiesc
		mouse_g <- General_Map(quesc[,1], in.org="Hsap", out.org="Mmus", in.name="symbol", out.name="symbol")
		mouse_g[mouse_g == ""] <- quesc[mouse_g == "",1]
		quesc[,1] <- mouse_g
		gene_table <- rbind(gene_table,quesc)
	}
	cyclemix <- run_CycleMix(GetAssayData(seur_obj, layer="counts"), gene_table)
	return(list(metadata=seur_obj@meta.data, 
		cyclemix=cyclemix
		))

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

if (args[1]==1) {
print("Wang1")
obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Wang_Colon_Phases <- run_phase_classifications(obj, species="hsap")
saveRDS(Wang_Colon_Phases, file="Wang_colon_phases_Quiescence.rds")

}

if (args[1]==2) {
print("Wang2")
obj <- readRDS("Wang_ileum.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Wang_Ileum_Phases <- run_phase_classifications(obj, species="hsap")
saveRDS(Wang_Ileum_Phases, file="Wang_ileum_phases_Quiescence.rds")

}

if (args[1]==3) {
print("Moreno")
obj <- readRDS("Moreno_core_gliobastoma.rds") 
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Moreno_Glio_Phases <- run_phase_classifications(obj, species="hsap")
saveRDS(Moreno_Glio_Phases, file="Moreno_glioblastoma_phases_Quiescence.rds")

}


if (args[1]==4) {
print("Gu")
obj <- readRDS("Gu_lymph.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Gu_Lymph_Phases <- run_phase_classifications(obj, species="mmus")
saveRDS(Gu_Lymph_Phases, file="Gu_lymph_phases_Quiescence.rds")

}

if (args[1]==5) {
print("Triana") 
obj <- readRDS("Triana_bonemarrow.rds")
obj@meta.data$UMAP1 <- obj@reductions$bothumap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$bothumap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Triana_Bone_Phases <- run_phase_classifications(obj, species="hsap")
saveRDS(Triana_Bone_Phases, file="Triana_bonemarrow_phases_Quiescence.rds")

}

if (args[1]==6) {
print("Steets")
obj <- readRDS("Streets_adipose.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)
Streets_adipose_Phases <- run_phase_classifications(obj, species="hsap")
saveRDS(Streets_adipose_Phases, file="Streets_adipose_phases_Quiescence.rds")
}

if (args[1] == 6) {
Phases <- list()
Phases[["Streets"]]<- readRDS(file="Streets_adipose_phases_Quiescence.rds")
Phases[["Triana"]]<- readRDS(file="Triana_bonemarrow_phases_Quiescence.rds")
Phases[["Wang_colon"]]<- readRDS(file="Wang_colon_phases_Quiescence.rds")
Phases[["Wang_ileum"]]<- readRDS(file="Wang_ileum_phases_Quiescence.rds")
Phases[["Moreno"]]<- readRDS(file="Moreno_glioblastoma_phases_Quiescence.rds")
Phases[["Gu"]]<- readRDS(file="Gu_lymph_phases_Quiescence.rds")

is_prolif <- list()
not_prolif <- list()
is_prolif[["Wang"]] <- c("transit amplifying cell of colon")
not_prolif[["Wang"]] <- c("colon goblet cell", "paneth cell of colon")
is_prolif[["Streets"]] <- c("PRM1+ Cells")
not_prolif[["Streets"]] <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
is_prolif[["Gu"]] <- c("Plasmablast-IgM+_cycling","Prolif-Tregs", "GC.BC_DZ.1", "GC.BC_DZ.2", "DC.prolif.")
not_prolif[["Gu"]] <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
is_prolif[["Moreno"]] <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
not_prolif[["Moreno"]] <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
is_prolif[["Triana"]] <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")
not_prolif[["Triana"]] <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")

type_column = list("Wang"="cell_type","Streets"="author_cell_type","Gu"="final_annotation","Moreno"="celltype_original","Triana"="cell_type")

quesc_stat <- function(x, dataset_name){
	cell_type_prolif <- is_prolif[[dataset_name]]
	cell_type_quiesc <- not_prolif[[dataset_name]]
	type_col <- type_column[[dataset_name]]
	cells1 <- x$metadata[,type_col] %in% cell_type_prolif;
	cells2 <- x$metadata[,type_col] %in% cell_type_quiesc;

	freq1 <- table(factor(x$cyclemix[[1]][cells1], levels=c("G0","G1S", "G2M", "None")))/sum(cells1)
	freq2 <- table(factor(x$cyclemix[[1]][cells2], levels=c("G0", "G1S", "G2M", "None")))/sum(cells2)
	return(rbind(freq1, freq2))

}

quesc_types <- function(x, dataset_name) {
	type_col <- type_column[[dataset_name]]
	tmp <- table(x$cyclemix$phase, x$metadata[,type_col])
	tmp <- t(tmp)/colSums(tmp)
	n_cells <- table(x$metadata[,type_col])
	tmp <- tmp[n_cells > 50,]

	tmp <- tmp[order(tmp[,1], decreasing=TRUE),]
	df <- data.frame(type=rownames(tmp), G0=tmp[,1], G1S=tmp[,2], G2M=tmp[,3], None=tmp[,4])
	return(df[df[,"G0"]>0.7,])
}


df <- c()
for (dataset in names(Phases)) {
	d <- sub("_.*$","", dataset)
	out <- quesc_stat(Phases[[dataset]], d)
	df <- rbind(out,df)
}
df <- df*100
df <- data.frame(df)
df$dataset <- rep(names(Phases), each=2)
df$pop <- rep(c("Prolif", "Quiesc"), times=length(Phases))
df$tmp <- paste(df$dataset, df$pop)

# Plot
require(ggplot2)
require(ggpattern)

png("Real_CM_Quiescent.png", width=4, height=4, units="in", res=300)
p1 <- ggplot(df, aes(x=tmp, y=G0, pattern=pop))+geom_bar_pattern(stat="identity", pattern_fill="black", 
	fill="white", color="black", pattern_spacing=0.05, 
	pattern_frequency=5, pattern_angle=45)+facet_wrap(~dataset, scales="free_x")+
scale_x_discrete(labels=c("prolif", "quiesc"))+xlab("")+ylab("Quiescent Cells (%)")+ scale_pattern_discrete(choices=c('none', 'stripe'))+ theme_classic() + guides(fill="none", pattern="none")
print(p1)
dev.off()




# Which cell-types have highest Quiesc?
# Quiesc vs GT
quesc_types <- function(x, dataset_name) {
        type_col <- type_column[[dataset_name]]
        tmp <- table(x$cyclemix$phase, x$metadata[,type_col])
        tmp <- t(tmp)/colSums(tmp)
        n_cells <- table(x$metadata[,type_col])
        tmp <- tmp[n_cells > 50,]

        tmp <- tmp[order(tmp[,1], decreasing=TRUE),]
        df <- data.frame(type=rownames(tmp), G0=tmp[,1], G1S=tmp[,2], G2M=tmp[,3], None=tmp[,4])
        #return(df[df[,"G0"]>0.6,])
        return(df[1:5,])
}

Qt_Wang1 <- quesc_types(Phases$Wang_colon, "Wang")
Qt_Wang2 <- quesc_types(Phases$Wang_ileum, "Wang")
Qt_Gu <- quesc_types(Phases$Gu, "Gu")
Qt_Triana <- quesc_types(Phases$Triana, "Triana")
Qt_Moreno <- quesc_types(Phases$Moreno, "Moreno")
Qt_Streets <- quesc_types(Phases$Streets, "Streets")

source("~/scripts/CycleMix_Benchmark/ColourScheme.R")
bar_dat <- c(Qt_Wang1[,2],0,Qt_Wang2[,2],0,Qt_Gu[,2],0,Qt_Triana[,2],0,Qt_Moreno[,2],0,Qt_Streets[,2])
bar_names <- c(Qt_Wang1[,1],"",Qt_Wang2[,1],"",Qt_Gu[,1],"",Qt_Triana[,1],"",Qt_Moreno[,1],"",Qt_Streets[,1])
bar_col=c(rep(dataset_col["Wang1"],5),"white", rep(dataset_col["Wang2"],5), "white", rep(dataset_col["Gu"],5), "white", rep(dataset_col["Tiana"],5), "white",rep(dataset_col["Moreno"],5),"white",rep(dataset_col["Streets"],5))

rep(c(dataset_col["Wang1"], dataset_col["Wang2"], dataset_col["Gu"], dataset_col["Triana"], dataset_col["Moreno"], dataset_col["Streets"]), each=5)

bar_names[13] <- "Cxcl12+ Fibroblast"
bar_names[14] <- "Pdgfra-hi Fibroblast"
bar_names[15] <- "Pdgfra-lo Fibroblast"
bar_names[16] <- "Endothelium (Lymphatic)"
bar_names[17] <- "Endosialin+ Fibroblast"
bar_names[19] <- "CD4+ ab memory T"
bar_names[20] <- "CD8+ ab effector memory T"
bar_names[21] <- "CD8+ ab central memory T"
bar_names[22] <- "gd T cell"
bar_names[23] <- "CD4+ ab naive T cell"

png("Quiescent_types_barplot.png", width=4, height=8, units="in", res=150)
par(mar=c(4,12,1,1))
barplot(bar_dat*100, names=bar_names, horiz=T, xlab="% Quiescent Cells", ylab="", las=1, col=bar_col)
abline(v=50, lty=2)
dev.off()
}
# Which genes other than our marker genes are upregulated in the G0 cells in the various datasets?
# Model with G0 vs G1S/G2M with cell-type as a predictor and dataset as batch using nebula


type_column = list("Wang"="cell_type","Streets"="author_cell_type","Gu"="final_annotation","Moreno"="celltype_original","Triana"="cell_type")
source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

run_nebula <- function(counts, meta, type_col) {
	require(nebula)
	design <- model.matrix(~meta$PHASE + meta[,type_col])
	resort <- order(meta$donor_id)
	counts <- counts[,resort]
	meta <- meta[resort,]
	out <- nebula(counts, meta$donor_id, pred=design, ncore=2)
	tmp <- out$summary[,c(grep("G1S", colnames(out$summary)), grep("G2M", colnames(out$summary)), ncol(out$summary))]
	for (p_col in grep("^p_", colnames(tmp))) {
		tmp[,p_col] <- p.adjust(tmp[,p_col], method="fdr")
	}
	return(tmp)
}

run_GLM <- function(counts, meta, type_col) {
	require(MAST)
	counts <- counts[rowSums(counts >0) > 0.01*ncol(counts),]
	lib.size <- colSums(counts)
	tpms <- t(log2(t(counts)/lib.size*median(lib.size)+1))
	mu_G0 = rowMeans(tpms[,meta$PHASE=="G0"])	
	mu_Cyc = rowMeans(tpms[,meta$PHASE!="G0"])	
	filter_genes <- (mu_G0 > 0.25 | mu_Cyc > 0.25) & abs(mu_G0-mu_Cyc) > 0.2
	tpms <- tpms[filter_genes,]
	design <- c()
	if (length(unique(meta[,"donor_id"])) > 1) {
		design <- model.matrix(~meta$PHASE + meta[,type_col] + meta[,"donor_id"])
	} else {
		design <- model.matrix(~meta$PHASE + meta[,type_col])
	}
	design <- design[,colSums(design) > 0]

	output <- apply(tpms, 1, function(x){
			reg <- glm(x~design)	
			return(summary(reg)$coeff)
		})
	rownames(output) <- paste(rep(colnames(design),4), rep(c("Estimate", "se", "tval", "pval"), each=ncol(design)), sep="_")

	logfc_g1s_row <- grep("_Estimate", rownames(output))[2]
	logfc_g2m_row <- grep("_Estimate", rownames(output))[3]
	se_g1s_row <- grep("_se", rownames(output))[2]
	se_g2m_row <- grep("_se", rownames(output))[3]
	pval_g1s_row <- grep("_pval", rownames(output))[2]
	pval_g2m_row <- grep("_pval", rownames(output))[3]
	
	reformated <- data.frame(Estimate_G1S=output[logfc_g1s_row,], se_G1S=output[se_g1s_row,],pval_G1S=output[pval_g1s_row,]*nrow(counts),
				Estimate_G2M=output[logfc_g2m_row,], se_G2M=output[se_g2m_row,],pval_G2M=output[pval_g2m_row,]*nrow(counts),
				gene=colnames(output))

	#sig_G1S <- output[grep("_pval", rownames(output))[2],] < 0.05/nrow(counts)
	#sig_G2M <- output[grep("_pval", rownames(output))[3],] < 0.05/nrow(counts)
	#effect_G1S <- output[grep("_Estimate", rownames(output))[2],]
	#effect_G2M <- output[grep("_Estimate", rownames(output))[3],]
	#dn_genes <- colnames(output)[sig_G1S & sig_G2M & sign(effect_G1S) == sign(effect_G2M) & sign(effect_G1S)>0]
	#up_genes <- colnames(output)[sig_G1S & sig_G2M & sign(effect_G1S) == sign(effect_G2M) & sign(effect_G1S)<0]

	return(reformated)
}	



if (args[2]==1) {
print("Wang1")
obj <- readRDS("Wang_colon.rds")
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Wang_Colon_Phases <- readRDS(file="Wang_colon_phases_Quiescence.rds")
exclude = Wang_Colon_Phases$cyclemix$phase == "None"
meta <- Wang_Colon_Phases$metadata; meta$PHASE = Wang_Colon_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Wang"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Wang_colon_Qnebula.rds")
out <- run_GLM(counts, meta, type_col)
saveRDS(out, "Wang_colon_Qglm.rds")
}

if (args[2]==2) {
print("Wang2")
obj <- readRDS("Wang_ileum.rds")
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Wang_Ileum_Phases <- readRDS(file="Wang_ileum_phases_Quiescence.rds")
exclude = Wang_Ileum_Phases$cyclemix$phase == "None"
meta <- Wang_Ileum_Phases$metadata; meta$PHASE = Wang_Ileum_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Wang"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Wang_ileum_Qnebula.rds")
out <- run_GLM(counts, meta, type_col)
saveRDS(out, "Wang_ileum_Qglm.rds")

}

if (args[2]==3) {
print("Moreno")
obj <- readRDS("Moreno_core_gliobastoma.rds") 
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Moreno_Glio_Phases <- readRDS(file="Moreno_glioblastoma_phases_Quiescence.rds")
exclude = Moreno_Glio_Phases$cyclemix$phase == "None"
meta <- Moreno_Glio_Phases$metadata; meta$PHASE = Moreno_Glio_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Moreno"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Moreno_glioblastoma_Qnebula.rds")
out <- run_GLM(counts, meta, type_col)
saveRDS(out, "Moreno_glioblastoma_Qglm.rds")

}


if (args[2]==4) {
print("Gu")
obj <- readRDS("Gu_lymph.rds")
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Mmus", out.org="Mmus", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Gu_Lymph_Phases <- readRDS(file="Gu_lymph_phases_Quiescence.rds")
exclude = Gu_Lymph_Phases$cyclemix$phase == "None"
meta <- Gu_Lymph_Phases$metadata; meta$PHASE = Gu_Lymph_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Gu"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Gu_lymph_Qnebula.rds")
out <- run_GLM(counts, meta, "cell_type")
saveRDS(out, "Gu_lymph_Qglm.rds")

}

if (args[2]==5) {
print("Triana") 
obj <- readRDS("Triana_bonemarrow.rds")
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Triana_Bone_Phases <- readRDS(file="Triana_bonemarrow_phases_Quiescence.rds")
exclude = Triana_Bone_Phases$cyclemix$phase == "None"
meta <- Triana_Bone_Phases$metadata; meta$PHASE = Triana_Bone_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Triana"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Triana_bonemarrow_Qnebula.rds")
out <- run_GLM(counts, meta, type_col)
saveRDS(out, "Triana_bonemarrow_Qglm.rds")

}

if (args[2]==6) {
print("Steets")
obj <- readRDS("Streets_adipose.rds")
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
Streets_adipose_Phases <- readRDS(file="Streets_adipose_phases_Quiescence.rds")
exclude = Streets_adipose_Phases$cyclemix$phase == "None"
meta <- Streets_adipose_Phases$metadata; meta$PHASE = Streets_adipose_Phases$cyclemix$phase;
counts <- counts[,!exclude]
meta <- meta[!exclude,]
type_col <- type_column[["Streets"]]
#out <- run_nebula(counts, meta, type_col)
#saveRDS(out, "Streets_adipose_Qnebula.rds")
out <- run_GLM(counts, meta, type_col)
saveRDS(out, "Streets_adipose_Qglm.rds")
}


if (args[2]==7) {

DE <- list()
DE[["Streets"]] <- readRDS("Streets_adipose_Qglm.rds")
DE[["Triana"]] <- readRDS("Triana_bonemarrow_Qglm.rds")
DE[["Gu"]] <- readRDS("Gu_lymph_Qglm.rds")
DE[["Moreno"]] <-readRDS("Moreno_glioblastoma_Qglm.rds")
DE[["Wang1"]] <-readRDS("Wang_colon_Qglm.rds")
DE[["Wang2"]] <-readRDS("Wang_ileum_Qglm.rds")

all_genes_up <- c()
all_genes_dn <- c()
for (dataset in names(DE)){
	print(dataset)
	# Criteria - significantly in same direction for G0 vs G1S and G0 vs G2M
	tmp <- DE[[dataset]]; tmp <- tmp[!is.na(tmp[,1]) & !is.na(tmp[,4]) & !is.na(tmp[,3]) & !is.na(tmp[,6]),]
	same_logFC <- sign(tmp[,1]) == sign(tmp[,4]) # consistent direction
	sig <- tmp[,3] < 0.05 | tmp[,6] < 0.05 # FDR < 5%
	dir <- -1*sign(tmp[,1])
	up_genes <- tmp[same_logFC & sig & dir > 0,7]
	dn_genes <- tmp[same_logFC & sig & dir < 0,7]
	all_genes_up <- c(all_genes_up, up_genes)
	print(length(up_genes))
	all_genes_dn <- c(all_genes_dn, dn_genes)
	print(length(dn_genes))
} # Doing this for Wang and Streets found no common genes - not even those I used to define the G0 cells are reproducibly associated across the samples.

reps <- table(all_genes_dn)
consensus_dn_genes <- names(reps[reps >= 4])

reps <- table(all_genes_up)
consensus_up_genes <- names(reps[reps >= 4])

paths <- readRDS("~/projects/def-tandrew6/tandrew6/Human_Pathways.rds")
rich <- MultiPath::do_ora(consensus_dn_genes, pathways=paths, background = rownames(counts))
rich2 <- MultiPath::condense_terms(rich)

tab <- rich2$results[order(rich2$result$FDR),]

png("G0_genes_pathways.png", width=6, height=4, units="in", res=300)
par(mar=c(4,12,1,1))
#barplot(-log10(tab$FDR)[1:9], col="navy", names=rownames(tab)[1:9], horiz=TRUE, xlab="-log(FDR)", las=1)
barplot(-log10(tab$FDR)[1:9], col="navy", names=c("MYC targets", "E2F targets", "Cell Cycle", "Chromosome organization", "Infectious Disease", "DNA Replication", "DNA Geometric Change", "NucEnvelope Reassembly", "mTOR signaling"), horiz=TRUE, xlab="-log(FDR)", las=1)
dev.off()

tab <- rich$results[order(rich$result$FDR),]
png("G0_genes_pathways_uncondensed.png", width=6, height=4, units="in", res=300)
par(mar=c(4,12,1,1))
#barplot(-log10(tab$FDR)[1:9], col="navy", names=rownames(tab)[1:9], horiz=TRUE, xlab="-log(FDR)", las=1)
barplot(-log10(tab$FDR)[1:9], col="navy", names=c("MYC targets", "E2F targets", "G2M Checkpoint","Cell Cycle", "Mitosis", "Chromosome organization", "M Phase", "Metaphase+Anaphase", "mRNA Splicing"), horiz=TRUE, xlab="-log(FDR)", las=1)
dev.off()



Fishers_method_vec <- function(x, y){
	chisq_stat <- -2*(log(x)+log(y))
	pval <- pchisq(chisq_stat, 2*2, lower.tail=FALSE)
	pval[pval < 10^-200] <- 0
	return(signif(pval, digits=1))
}

Gene_table <- matrix(NA, nrow=length(unique(c(all_genes_dn, all_genes_up))), ncol=6) # because moreno is missing
Gene_table_plot <- matrix(NA, nrow=length(unique(c(all_genes_dn, all_genes_up))), ncol=6) # because moreno is missing
colnames(Gene_table) <- names(DE)
rownames(Gene_table) <- unique(c(all_genes_dn, all_genes_up))
colnames(Gene_table_plot) <- names(DE)
rownames(Gene_table_plot) <- unique(c(all_genes_dn, all_genes_up))
for (dataset in names(DE)){
        print(dataset)
        tmp <- DE[[dataset]]; 
	tmp <- tmp[match(rownames(Gene_table), tmp$gene),]
	Eff <- -1*(tmp[,1]+tmp[,4])/2
	Eff[sign(tmp[,1]) != sign(tmp[,4])] <- 0
	Gene_table[,dataset] <- paste(signif(Eff, digits=1), Fishers_method_vec(tmp[,3], tmp[,6]), sep=":")
	Gene_table_plot[,dataset] <- Eff
} 
write.table(Gene_table, "Suppl_Quiescence_Genetable.csv", sep=",")
# Heatmap?
Gene_table_plot[is.na(Gene_table_plot)] <- 0
score <- rowMeans(Gene_table_plot)
head(sort(score)); tail(sort(score))
genes_to_plot <- c(consensus_up_genes, names(head(sort(score[consensus_dn_genes]),20)))

tmp <-  PurpleAndYellow(80)
heat_col <- c(tmp[1:20], "black",tmp[61:80])
dat_range <-  range(Gene_table_plot[genes_to_plot,])
symm_breaks <- seq(from=-1*max(abs(dat_range)), to=max(abs(dat_range)), length=length(heat_col)+1)

is.ref.gene <- genes_to_plot %in% HGeneSets$Quiesc[,1]
df <- data.frame(ref.gene = factor(is.ref.gene))
rownames(df) <- genes_to_plot
tmp <- c("white", "black"); names(tmp) = c(FALSE, TRUE)

require(pheatmap)
png("G0_genes_heatmap.png", width=4, height=6, units="in", res=300)
pheatmap(Gene_table_plot[genes_to_plot,], breaks=symm_breaks, color=heat_col, annotation_row = df, annotation_colors = list(ref.gene=tmp))
dev.off()


}


