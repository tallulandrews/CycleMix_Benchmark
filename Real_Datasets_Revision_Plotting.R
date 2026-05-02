source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/Seurat.R")
source("~/My_R_Packages/CycleMix/R/Analysis.R")
require("mclust")
require("Matrix")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

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

        se_freq <- table(factor(x$seurat[[1]][cells], levels=c("G1","S", "G2M")))/sum(cells)
        names(se_freq) <- c("None", "S", "G2M")
        se_freq <- se_freq[c("None", "S", "G2M")]

        tab <- rbind(cm_freq, se_freq)
        rownames(tab) <- c("CycleMix", "Seurat")
        return(tab[,2]+tab[,3])
}

get_stats_smooth <- function(x, cell_type, type_col="cell_type") {
	cells <- x$sct$metadata[,type_col] %in% cell_type
	if(length(x$smoothed$knn$phase) != length(cells)){stop("knn length no match")}
	if(length(x$smoothed$cluster$phase) != length(cells)){stop("knn length no match")}
	knn_freq <- table(factor(x$smoothed$knn$phase[cells], levels=c("G1S", "G2M", "None")))/sum(cells)
	names(knn_freq) <- c("S", "G2M", "None")
	knn_freq <- knn_freq[c("None", "S", "G2M")]

	cluster_freq <- table(factor(x$smoothed$cluster$phase[cells], levels=c("G1S", "G2M", "None")))/sum(cells)
	names(cluster_freq) <- c("S", "G2M", "None")
	cluster_freq <- cluster_freq[c("None", "S", "G2M")]

	tab <- rbind(knn_freq, cluster_freq)
	rownames(tab) <- c("KNN", "Cluster")
	return(tab[,2]+tab[,3])
}




set.seed(3891)


Wang_is_prolif <- c("transit amplifying cell of colon")
Wang_not_prolif <- c("colon goblet cell", "paneth cell of colon")
Wang_column <- "cell_type"
Wang1 <- readRDS("Wang_colon_phases.rds")
Wang1_adv <- readRDS("Wang_colon_advphases_allSeurGenes.rds")
Wang1_baseline_prolif <- get_stats(Wang1, Wang_is_prolif)
Wang1_baseline_notprolif <- get_stats(Wang1, Wang_not_prolif)
Wang1_sct_prolif <- get_stats(Wang1_adv$sct, Wang_is_prolif)
Wang1_sct_notprolif <- get_stats(Wang1_adv$sct, Wang_not_prolif)
Wang1_smooth_prolif <- get_stats_smooth(Wang1_adv, Wang_is_prolif)
Wang1_smooth_notprolif <- get_stats_smooth(Wang1_adv, Wang_not_prolif)

Wang2 <- readRDS("Wang_ileum_phases.rds")
Wang2_adv <- readRDS("Wang_ileum_advphases_allSeurGenes.rds")
Wang2_baseline_prolif <- get_stats(Wang2, Wang_is_prolif)
Wang2_baseline_notprolif <- get_stats(Wang2, Wang_not_prolif)
Wang2_sct_prolif <- get_stats(Wang2_adv$sct, Wang_is_prolif)
Wang2_sct_notprolif <- get_stats(Wang2_adv$sct, Wang_not_prolif)
Wang2_smooth_prolif <- get_stats_smooth(Wang2_adv, Wang_is_prolif)
Wang2_smooth_notprolif <- get_stats_smooth(Wang2_adv, Wang_not_prolif)

Street_is_prolif <- c("PRM1+ Cells")
Street_not_prolif <- c("Mature Adipocytes", "Smooth Muscle Cells", "Endothelial Cells")
Street_column <- "author_cell_type"
Street <- readRDS("Streets_adipose_phases.rds")
Street_adv <- readRDS("Streets_adipose_advphases_allSeurGenes.rds")
Street_baseline_prolif <- get_stats(Street, Street_is_prolif, Street_column)
Street_baseline_notprolif <- get_stats(Street, Street_not_prolif, Street_column)
Street_sct_prolif <- get_stats(Street_adv$sct, Street_is_prolif, Street_column)
Street_sct_notprolif <- get_stats(Street_adv$sct, Street_not_prolif, Street_column)
Street_smooth_prolif <- get_stats_smooth(Street_adv, Street_is_prolif, Street_column)
Street_smooth_notprolif <- get_stats_smooth(Street_adv, Street_not_prolif, Street_column)

Gu_is_prolif <- c("Plasmablast-IgM+_cycling","Prolif-Tregs", "GC.BC_DZ.1", "GC.BC_DZ.2", "DC.prolif.")
Gu_not_prolif <- c("PCs.CD19+IgM+", "PCs.CD19-IgM+.1", "PCs.CD19-IgM+.2", "PCs.IgA+","CD4.memory", "Naive.Bcells", "Th17", "Endothelial.cells", "eMBCs", "Naive.CD8+_MLN", "Naive.CD4+_MLN")
Gu_column <- "final_annotation"
Gu <- readRDS("Gu_lymph_phases.rds")
Gu_adv <- readRDS("Gu_lymph_advphases_allSeurGenes.rds")
Gu_baseline_prolif <- get_stats(Gu, Gu_is_prolif, Gu_column)
Gu_baseline_notprolif <- get_stats(Gu, Gu_not_prolif, Gu_column)
Gu_sct_prolif <- get_stats(Gu_adv$sct, Gu_is_prolif, Gu_column)
Gu_sct_notprolif <- get_stats(Gu_adv$sct, Gu_not_prolif, Gu_column)
Gu_smooth_prolif <- get_stats_smooth(Gu_adv, Gu_is_prolif, Gu_column)
Gu_smooth_notprolif <- get_stats_smooth(Gu_adv, Gu_not_prolif, Gu_column)

Moreno_is_prolif <- c("Tumour", "Radial Glia", "Prolif stemcell tumor", "Neoplastic", "Dividing Progenitor")
Moreno_not_prolif <- c("Oligodendrocyte", "Endothelial", "Astrocyte", "Normal brain", "Plasma B")
Moreno_column <- "celltype_original"
Moreno <- readRDS("Moreno_glioblastoma_phases_SUBSAMPLE.rds")
Moreno_adv <- readRDS("Moreno_SUBSAMPLE_advphases_allSeurGenes.rds")
Moreno_adv$sct$metadata <- Moreno$metadata
Moreno_baseline_prolif <- get_stats(Moreno, Moreno_is_prolif, Moreno_column)
Moreno_baseline_notprolif <- get_stats(Moreno, Moreno_not_prolif, Moreno_column)
Moreno_sct_prolif <- get_stats(Moreno_adv$sct, Moreno_is_prolif, Moreno_column)
Moreno_sct_notprolif <- get_stats(Moreno_adv$sct, Moreno_not_prolif, Moreno_column)
Moreno_smooth_prolif <- get_stats_smooth(Moreno_adv, Moreno_is_prolif, Moreno_column) ### ERROR!
Moreno_smooth_notprolif <- get_stats_smooth(Moreno_adv, Moreno_not_prolif, Moreno_column) ### ERROR!

Triana_is_prolif <- c("hematopoietic multipotent progenitor cell", "lymphoid lineage restricted progenitor cell", "common dendritic progenitor", "megakaryocyte-erythroid progenitor cell","erythroid progenitor cell")
Triana_not_prolif <- c("classical monocyte", "central memory CD8-positive, alpha-beta T cell", "unswitched memory B cell", "class switched memory B cell", "CD16-positive, CD56-dim natural killer cell, human", "CD4-positive, alpha-beta memory T cell", "naive B cell", "CD8-positive, alpha-beta memory T cell")
Triana_column <- "cell_type"
Triana <- readRDS("Triana_bonemarrow_phases.rds")
Triana_adv <- readRDS("Triana_bonemarrow_advphases_allSeurGenes.rds")
Triana_baseline_prolif <- get_stats(Triana, Triana_is_prolif, Triana_column)
Triana_baseline_notprolif <- get_stats(Triana, Triana_not_prolif, Triana_column)
Triana_sct_prolif <- get_stats(Triana_adv$sct, Triana_is_prolif, Triana_column)
Triana_sct_notprolif <- get_stats(Triana_adv$sct, Triana_not_prolif, Triana_column)
Triana_smooth_prolif <- get_stats_smooth(Triana_adv, Triana_is_prolif, Triana_column)
Triana_smooth_notprolif <- get_stats_smooth(Triana_adv, Triana_not_prolif, Triana_column)


## Plot 1 : LogNorm vs SCT
df <- data.frame(method=rep(c("CycleMix", "Seurat"),times=20),
		 norm=rep(c("LogNorm", "LogNorm", "SCT", "SCT"), times=10),
		 type=rep(c("prolif","quiesc"), each=20),
		 dataset=rep(rep(c("Triana", "Gu", "Street", "Wang1", "Wang2"),each=4), times=2),
		prolif_prop=c(Triana_baseline_prolif, Triana_sct_prolif,
				Gu_baseline_prolif, Gu_sct_prolif,
				Street_baseline_prolif, Street_sct_prolif,
				Wang1_baseline_prolif, Wang1_sct_prolif,
				Wang2_baseline_prolif, Wang2_sct_prolif,
				Triana_baseline_notprolif, Triana_sct_notprolif,
				Gu_baseline_notprolif, Gu_sct_notprolif,
				Street_baseline_notprolif, Street_sct_notprolif,
				Wang1_baseline_notprolif, Wang1_sct_notprolif,
				Wang2_baseline_notprolif, Wang2_sct_notprolif)*100
		)
df$dataset_type <- paste0(df$dataset, df$type)
require(ggplot2)
				
p <- ggplot(df, aes(x=norm, y=prolif_prop, color=type)) +
	geom_dotplot(aes(fill=type), binaxis="y", binwidth=0.05,
		stackdir="center", dotsize=50) +
	geom_line(aes(group=dataset_type), color="grey35") +
	facet_wrap(~method, scales="free_x") +
	scale_x_discrete(labels=c("LogNorm", "SCT")) +
	xlab("") + ylab("Cycling Cells (%)") + theme_classic() + guides(color="none")

#png("Real_Datasets_LogNorm_vs_SCT.png", width=8*0.65, height=6*0.65, units="in", res=300)
pdf("Real_Datasets_LogNorm_vs_SCT_allseurgenes.pdf", width=8*0.65, height=6*0.65)
print(p)
dev.off()


## Plot 2: CycleMix OG vs Smoothing
df <- data.frame(smoothing=rep(c("None", "KNN", "Cluster"), times=12),
		type=rep(c("prolif", "quiesc"), each=6*3),
		dataset=rep(rep(c("Triana", "Moreno", "Gu", "Street", "Wang1", "Wang2"), each=3), times=2),
		prolif_prop=c(Triana_baseline_prolif[1], Triana_smooth_prolif,
				Moreno_baseline_prolif[1], Moreno_smooth_prolif,
				Gu_baseline_prolif[1], Gu_smooth_prolif,
				Street_baseline_prolif[1], Street_smooth_prolif,
				Wang1_baseline_prolif[1], Wang1_smooth_prolif,
				Wang2_baseline_prolif[1], Wang2_smooth_prolif,
				Triana_baseline_notprolif[1], Triana_smooth_notprolif,
				Moreno_baseline_notprolif[1], Moreno_smooth_notprolif,
				Gu_baseline_notprolif[1], Gu_smooth_notprolif,
				Street_baseline_notprolif[1], Street_smooth_notprolif,
				Wang1_baseline_notprolif[1], Wang1_smooth_notprolif,
				Wang2_baseline_notprolif[1], Wang2_smooth_notprolif)*100
		)
df$smoothing <- factor(df$smoothing, levels=c("None", "KNN", "Cluster"))
df$dataset_type <- paste0(df$dataset, df$type)
require(ggplot2)

p <- ggplot(df, aes(x=smoothing, y=prolif_prop, color=type)) +
	geom_dotplot(aes(fill=type), binaxis="y", binwidth=0.05,
		stackdir="center", dotsize=50) +
	geom_line(aes(group=dataset_type), color="grey35") +
	xlab("") + ylab("Cycling Cells (%)") + theme_classic() + guides(color="none")

#png("Real_Datasets_Smoothing.png", width=8*0.65/4*3, height=6*0.65, units="in", res=300)
pdf("Real_Datasets_Smoothing_allseurgenes.pdf", width=8*0.65/4*3, height=6*0.65)
print(p)
dev.off()


##### Plot 4.1: 
# CycleMix across gene sets.

files <- Sys.glob("*_cyclemix_allgenesets.rds")
# Plot % proliferating for each dataset across all gene sets.
out <- matrix(0, nrow=length(files), ncol=5)
colnames(out) <- c("Seurat", "Tirosh", "Cyclone", "Macosko", "Whit")
for (i in 1:length(files)) {
	f <- files[i]
	dat <- readRDS(f)
	for (set in names(dat)) {
		prolif_prop <- 1-sum(dat[[set]] %in% c("G1", "None"))/length(dat[[set]])
		out[i,set] <- prolif_prop
	}
}
rownames(out) <- c("Gu", "Moreno", "Streets", "Triana", "Wang1", "Wang2")
colnames(out) <- c("Seurat", "Tirosh", "Cyclone", "Macosko", "Whitfield")

out <- round(out, digits=4)
df <- data.frame(dataset=rep(rownames(out), times=ncol(out)),
		 geneset=rep(colnames(out), each=length(files)),
		 prolif=as.numeric(out)
		)
df$geneset <- factor(df$geneset, levels=c("Seurat", "Tirosh", "Cyclone", "Macosko", "Whitfield"))
require(ggplot2)
require(RColorBrewer)
#png("Prolif_prop_vs_geneset.png", width=6, height=4, units="in", res=300)
pdf("Prolif_prop_vs_geneset.pdf", width=6, height=4)
ggplot(df, aes(x=dataset, y=prolif, fill=geneset))+geom_bar(stat="identity", position="dodge")+theme_classic()+scale_fill_manual(values=brewer.pal(5, "Set2")) + ylab("Proportion Proliferative")
dev.off()

# Plot 4.2  GeneSet overlap (genes)
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
require(CycleMix)
cyclone_h <- remap_genesets(MGeneSets$Cyclone, to.spp="hsap")
all_sets <- list("Seurat"=HGeneSets$Seurat, "Tirosh"=HGeneSets$Tirosh, "Macosko"=HGeneSets$Macosko, "Whitfield"=HGeneSets$Whitfield, "Cyclone"=cyclone_h)

out <- matrix(0, nrow=length(all_sets), ncol=length(all_sets))
for (i in 1:length(all_sets)){
	for (j in 1:length(all_sets)){
		olap <- length(intersect(all_sets[[i]]$Gene, all_sets[[j]]$Gene))/max(c(length(unique(all_sets[[i]]$Gene)), length(unique(all_sets[[j]]$Gene))))
		out[i,j] <- olap
	}
}
rownames(out) <- names(all_sets)
colnames(out) <- names(all_sets)
require(pheatmap)
#png("Geneset_vs_geneset.png", width=4, height=4, units="in", res=300)
pdf("Geneset_vs_geneset.pdf", width=4, height=4)
pheatmap(out)
dev.off()

# Plot 4.3 For each pair of gene sets calculate average cells per datasets in each intersections of phases

all_res <- list()
for (f in files) {
	dat <- readRDS(f)
	dat <- dat[c(4,2,5,3,1)]
	for (i in 1:length(dat)){
	        for (j in (i):length(dat)){
			tab <- table(dat[[i]], dat[[j]])/length(dat[[i]])
			if (f == files[1]) {
				all_res[[paste(i,j,sep="_")]] <- tab
			} else {
				tab <- tab[,match(colnames(all_res[[paste(i,j,sep="_")]]), colnames(tab))]
				tab <- tab[match(rownames(all_res[[paste(i,j,sep="_")]]), rownames(tab)),]
				tab[is.na(tab)] <- 0
				all_res[[paste(i,j,sep="_")]] <-  all_res[[paste(i,j,sep="_")]]+tab
			}
		}
	}
}


out <- list()
for (i in 1:(length(dat)-1)){
        for (j in (i+1):length(dat)){
		v1 <- diag(all_res[[paste(i,i,sep="_")]])/6
		v2 <- diag(all_res[[paste(j,j,sep="_")]])/6

		tmp <- all_res[[paste(i,j,sep="_")]]/6
		exp <- v1 %*% t(v2)
		out[[paste(i,j,sep="_")]]<-log2(tmp/exp)
	}
}

# Make Fancy Plot
row_names <- c(rownames(out[["1_2"]]), rownames(out[["2_3"]]), rownames(out[["3_4"]]), rownames(out[["4_5"]]))
col_names <- c(colnames(out[["1_2"]]), colnames(out[["2_3"]]), colnames(out[["3_4"]]), colnames(out[["4_5"]]))
big_mat <- matrix(NA, nrow=length(row_names), ncol=length(col_names)) 

st_col=1
st_row=1
col_step <- 0;
for (i in 1:(length(dat)-1)){
	this_dims <- c(0,0)
        for (j in (i+1):length(dat)){
		this_dims <- dim(out[[paste(i,j,sep="_")]])
		big_mat[st_row:(st_row-1+this_dims[1]), st_col:(st_col-1+this_dims[2])] <- out[[paste(i,j,sep="_")]]
		st_col <- st_col+this_dims[2]
		if (j == i+1) {
			col_step <- col_step + this_dims[2]
		}
	}
	st_col <- 1+col_step
	st_row <- st_row + this_dims[1]
}

colnames(big_mat) <- col_names
rownames(big_mat) <- row_names

reorder_columns <- c(1,2,3,4,7,5,6,8,13,9,10,11,12,15,20,16,17,18,14,19)
reorder_rows <- c(1,2,3,4,5,6,7,10,8,9,11,16,12,13,14,15)

to_plot <- big_mat
to_plot[!is.finite(to_plot)] <- 0
to_plot[to_plot < -10] <- 0
require(pheatmap)
#png("fancy_geneset_vs_geneset_heatmap.png", width=8, height=8, units="in", res=300)
pdf("fancy_geneset_vs_geneset_heatmap.pdf", width=8, height=8)
pheatmap(to_plot[reorder_rows,reorder_columns], cluster_rows=FALSE, cluster_cols=FALSE,  color = colorRampPalette(c("blue", "white", "red"))(21), breaks = c(seq(from=-8, to=-0.01, len=11), seq(from=0.01, to=6, len=11)), border_color=NA, gaps_row=c(3,6,10), gaps_col=c(3,7,13))
dev.off()
	

# Plot - Downsampling
Wang_column <- "cell_type"
# Wang Colon
require(Seurat)
obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=20)
Wang_is_prolif <- c("transit amplifying cell of colon")
Wang_not_prolif <- c("colon goblet cell", "paneth cell of colon")

Wang_get_stats <- function(obj, phases) {
	type_col <- "cell_type"
	Wang_is_prolif <- c("transit amplifying cell of colon")
	Wang_not_prolif <- c("colon goblet cell", "paneth cell of colon")

        cells <- obj@meta.data[,type_col] %in% Wang_is_prolif;
        if (length(phases) != length(cells)) {stop("Cyclemix length no match metadata.")}

        prolif_freq <- table(factor(phases[cells], levels=c("G1S", "G2M", "None")))/sum(cells)
        names(prolif_freq) <- c("S", "G2M", "None");
        prolif_freq <- prolif_freq[c("None", "S", "G2M")]

        cells <- obj@meta.data[,type_col] %in% Wang_not_prolif;
        if (length(phases) != length(cells)) {stop("Cyclemix length no match metadata.")}

        not_freq <- table(factor(phases[cells], levels=c("G1S", "G2M", "None")))/sum(cells)
        names(not_freq) <- c("S", "G2M", "None");
        not_freq <- not_freq[c("None", "S", "G2M")]

        tab <- rbind(prolif_freq, not_freq)
        rownames(tab) <- c("Prolif", "Quiesc")
        return(tab[,2]+tab[,3])
}

downsample_files <- Sys.glob("Wang_colon*fits.rds")
out <- c()
# Plot prolif & quiec vs downsample
for (prop in seq(to=0.95, from=0.05, by=0.05)) {
	f <- readRDS(paste0("Wang_colon_", prop, "_CycleMix_fits.rds"))
	stuff <- c(prop, Wang_get_stats(obj, f$cyclemix$phase))
	out <- rbind(out,stuff)
}


set.seed(3892)
classification <- readRDS("Wang_colon_0.05_CycleMix_fits.rds")$cyclemix
smooth_out <- c()
for (threshold in c(0.5,rev(seq(to=0.25, from=0.05, by=0.05)))) {
	smoothed <- knnSmooth(classification, dim=obj@reductions$pca@cell.embeddings, k=20, min.threshold=threshold)
	stuff <- c(threshold, Wang_get_stats(obj, smoothed$phase))
	smooth_out <- rbind(smooth_out, stuff)
}



# Plot solid lines for downsampling then dotted line for smoothing
#png("Wang_colon_downsample_plot.png", width=8, height=4, units="in", res=300)
pdf("Wang_colon_downsample_plot.pdf", width=8, height=4)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(out[,1], out[,2], type="l", lwd=2, col=gt_col["Proliferative"], ylim=c(0,1), xlim=c(1,0), bty="n", xlab="Downsample (%)", ylab="% of Cells")
lines(out[,1], out[,3], lwd=2, col=gt_col["Quiescent"])

par(mar=c(4,0,1,1))
plot(smooth_out[,1], smooth_out[,2], type="l", lwd=2, col=gt_col["Proliferative"], ylim=c(0,1), xlim=c(0.5,0), yaxt="n", lty=2, xlab="Smoothing threshold", ylab="", bty="n")
lines(smooth_out[,1], smooth_out[,3], lwd=2, col=gt_col["Quiescent"], lty=2)
legend("right", col=c(gt_col["Proliferative"], gt_col["Quiescent"], "black", "black"), lty=c(1,1,1,2), c("Proliferative", "Quiescent", "Original", "knn-smoothed"), bty="n", lwd=2)

dev.off()



####################### Defunct ######################
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




