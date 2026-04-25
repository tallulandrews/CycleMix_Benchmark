source("/home/tandrew6/scripts/CycleMix_Benchmark/Cyclone.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/CycleMix.R")
source("/home/tandrew6/scripts/CycleMix_Benchmark/Seurat.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

obj <- readRDS("developing_immune_system_meyloidcells.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1] 
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
obj <- CreateSeuratObject(counts=counts, meta.data=obj@meta.data)

obj <- NormalizeData(obj)
mat <- GetAssayData(obj, layer="data")

#Fixing Time
# Carnegie stage 13 = 5th week / 28-32 day
# Carnegie stage 19 = 7th week /day 48-51
# Carnegie stage 20 = 7.5th week 51-53 day
# Carnegie stage 22 = 8 week 54-56 day
# Carnegie stage 23 = 8.5th week56-60

age <- rep(0, nrow(obj@meta.data))
age[obj@meta.data$development_stage == "Carnegie stage 13"] <- 32
age[obj@meta.data$development_stage == "Carnegie stage 19"] <- 51
age[obj@meta.data$development_stage == "Carnegie stage 20"] <- 53
age[obj@meta.data$development_stage == "Carnegie stage 22"] <- 56
age[obj@meta.data$development_stage == "Carnegie stage 23"] <- 60
age[obj@meta.data$development_stage == "9th week post-fertilization stage"] <- 7*9
age[obj@meta.data$development_stage == "10th week post-fertilization stage"] <- 7*10
age[obj@meta.data$development_stage == "11th week post-fertilization stage"] <- 7*11
age[obj@meta.data$development_stage == "12th week post-fertilization stage"] <- 7*12
age[obj@meta.data$development_stage == "13th week post-fertilization stage"] <- 7*13
age[obj@meta.data$development_stage == "14th week post-fertilization stage"] <- 7*14
age[obj@meta.data$development_stage == "15th week post-fertilization stage"] <- 7*15
age[obj@meta.data$development_stage == "16th week post-fertilization stage"] <- 7*16
age[obj@meta.data$development_stage == "17th week post-fertilization stage"] <- 7*17
obj@meta.data$age_days <- age



# CycleMix
CycleMix <- classifyCells(mat, HGeneSets$Tirosh, symbol_column=NULL)
obj@meta.data$CycleMix <- CycleMix$phase

# Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


cell_type = "Kupffer cell"
dat <- obj@meta.data[obj@meta.data$cell_type == cell_type,]
dat$Phase <- factor(dat$Phase, c("S", "G2M", "G1"))
tmp <- table(dat$age_days, dat$Phase) #### CycleMix or Phase
toplot1 <- tmp/rowSums(tmp)*100
toplot1 <- cbind(toplot1, as.numeric(rownames(toplot1)))
toplot1 <- cbind(toplot1, as.vector(table(dat$age_days)))
toplot1 <- toplot1[toplot1[,5] > 50,]

cell_type = "macrophage"
dat <- obj@meta.data[obj@meta.data$cell_type == cell_type,]
dat$Phase <- factor(dat$Phase, c("S", "G2M", "G1"))
tmp <- table(dat$age_days, dat$Phase) ### CycleMix or Phase
toplot2 <- tmp/rowSums(tmp)*100
toplot2 <- cbind(toplot2, as.numeric(rownames(toplot2)))
toplot2 <- cbind(toplot2, as.vector(table(dat$age_days)))

png("FetalAtlas_Seurat.png", width=6, height=6, units="in", res=150)
freq <- toplot1[,2]+toplot1[,1]
stderr <- sqrt(freq/100*(1-freq/100)/toplot1[,5])*100
plot(toplot1[,4], freq, ylab="Proliferating (%)", xlab="Age (days)", type="l", lwd=2, col="blue", ylim=c(0,60), main="Human fetal meyloid cells") 
arrows(toplot1[,4], freq, toplot1[,4], freq+stderr, angle=90, length=0.1, lwd=2, col="blue")
arrows(toplot1[,4], freq, toplot1[,4], freq-stderr, angle=90, length=0.1, lwd=2, col="blue")
lines(toplot2[,4], toplot2[,2]+toplot2[,1], col="black", lwd=2)
freq <- toplot2[,2]+toplot2[,1]
stderr <- sqrt(freq/100*(1-freq/100)/toplot2[,5])*100
arrows(toplot2[,4], freq, toplot2[,4], freq+stderr, angle=90, length=0.1, lwd=2, col="black")
arrows(toplot2[,4], freq, toplot2[,4], freq-stderr, angle=90, length=0.1, lwd=2, col="black")
abline(h=56, lty=2, col="grey50", lwd=2)
legend("topright", c("Kupffer cells", "Macrophages"), col=c("blue", "black"), lwd=2)
dev.off()



table(obj@meta.data$CycleMix, obj@meta.data$cell_type)
table(obj@meta.data$CycleMix, obj@meta.data$tissue)
obj@meta.data$tissue_celltype <- paste(obj@meta.data$tissue, obj@meta.data$cell_type, sep="_")
table(obj@meta.data$CycleMix, obj@meta.data$tissue_celltype)








