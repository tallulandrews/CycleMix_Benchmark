
Wang_is_prolif <- c("transit amplifying cell of colon")
Wang_not_prolif <- c("colon goblet cell", "paneth cell of colon")
Wang_column <- "cell_type"
Wang1 <- readRDS("Wang_colon_corrected.rds")
Wang2 <- readRDS("Wang_ileum_corrected.rds")

Street_is_prolif <- c("PRM1+ Cells")
Street_not_prolif <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
Street_column <- "author_cell_type"
Street <- readRDS("Streets_adipose_corrected.rds")

Gu_is_prolif <- c("Plasmablast-IgM+_cycling","Prolif-Tregs", "GC.BC_DZ.1", "GC.BC_DZ.2", "DC.prolif.")
Gu_not_prolif <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
Gu_column <- "final_annotation"
Gu <- readRDS("Gu_lymph_corrected.rds")

Moreno_is_prolif <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
Moreno_not_prolif <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
Moreno_column <- "celltype_original"
Moreno <- readRDS("Moreno_glioblastoma_corrected.rds")

Triana_is_prolif <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")
Triana_not_prolif <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")
Triana_column <- "cell_type"
Triana <- readRDS("Triana_bonemarrow_corrected.rds")

all_corrected <- list("Wang_colon"=Wang1, "Wang_ileum"=Wang2,
			"Gu_lymph"=Gu, "Streets_adipose"=Street,
			"Moreno_glioblastoma"=Moreno, 
			"Triana_bonemarrow"=Triana)

# cell-type lisi
for (i in 1:length(all_corrected)) {
	corrected <- all_corrected[[i]]
	this_name <- names(all_corrected[i])
	require(lisi)
	df <- data.frame("cell_type" = corrected$cell_type)
	tmp2 <- compute_lisi(corrected$seur_umap, df, "cell_type")
	tmp3 <- compute_lisi(corrected$orig_umap, df, "cell_type")
	tmp1 <- compute_lisi(corrected$cyclemix_umap, df, "cell_type")

	df <- data.frame(correction=c("None", "CycleMix", "Seurat"), 
		meanLISI = c(mean(tmp3[,1]), mean(tmp1[,1]), mean(tmp2[,1])), 
		medianLISI=c(median(tmp3[,1]), median(tmp1[,1]), median(tmp2[,1])))
	p.vals <- c("CycleMix_vs_Seurat"=wilcox.test(tmp1[,1], tmp2[,1])$p.value, 
	    "CycleMix_vs_None"=wilcox.test(tmp1[,1], tmp3[,1])$p.value,
	    "Seurat_vs_None"=wilcox.test(tmp2[,1], tmp3[,1])$p.value)
	saveRDS(list(lisi=df, wilcox.pval=p.vals), paste0(this_name,"_correction_LISI.rds"))
}

