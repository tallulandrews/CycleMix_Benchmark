require("mclust")
require("Matrix")
require("Seurat")
require("CycleMix")
require(SingleCellExperiment)
require(scater)
require(scran)

source("~/My_R_Packages/CycleMix/R/Plotting.R")
source("~/My_R_Packages/CycleMix/R/Analysis.R")
source("~/projects/def-tandrew6/tandrew6/External_Data/Ensembl/Ensembl_Stuff.R")
load("/home/tandrew6/My_R_Packages/CycleMix/data/HGenes.rda") # get latest gene sets

dir = "/home/tandrew6/scratch/CycleMix_Manuscript"

set.seed(3891)

output <- c();
time <- c();

obj <- readRDS("Wang_colon.rds")
obj@meta.data$UMAP1 <- obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP2 <- obj@reductions$umap@cell.embeddings[,2]
counts <- GetAssayData(obj, layer="counts")
sym <- General_Map(rownames(counts), in.org="Hsap", out.org="Hsap", in.name="ensg", out.name="symbol")
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep,]; rownames(counts) <- sym[keep]
SCE <- SingleCellExperiment(assays=list(counts=counts), rowData=data.frame(feature_symbol=rownames(counts)))
logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(counts)/median(colSums(counts)))
out <- classifyCells(SCE, HGeneSets$Tirosh)


png("Figure1_mixtureG1S_2.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G1S, BIC=TRUE)
dev.off()

png("Figure1_mixtureG2M_2.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G2M, BIC=TRUE)
dev.off()


source("/home/tandrew6/scripts/CycleMix_Benchmark/Read_Buttner_with_GT.R")
SCE <- SingleCellExperiment(assays=list(counts=Buttner_counts), rowData=data.frame(feature_symbol=rownames(Buttner_counts)))
logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(Buttner_counts)/median(colSums(Buttner_counts)))

out <- classifyCells(SCE,MGeneSets$Cyclone)

png("Figure1_mixtureS_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$S, BIC=TRUE)
dev.off()

png("Figure1_mixtureG2M_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G2M, BIC=TRUE)
dev.off()

png("Figure1_mixtureG2M_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
par(mfrow=c(1,2))
h <- hist(out$fits$G2M$data, prob = TRUE, col = "grey85", xlab = "G2M Score",
        main = "", breaks = 20, ylim=c(0,1.2), xlim=c(0,6))
is.G2M <- Buttner_meta == "G2M"
is.S <- Buttner_meta == "S"
points(out$fits$G2M$data[is.G2M], rep(1.1, sum(is.G2M)), pch=18, cex=0.5)
points(out$fits$G2M$data[is.S], rep(1, sum(is.S)), pch=18, cex=0.5)
text(4.5,1.1, "G2M", pos=4)
text(4.5,1, "S", pos=4)

fit = out$fits$G2M
cols=c("blue", "black", "orange")
for (i in 1:fit$G) {
        if (is.na(fit$parameters$variance$sigmasq[i])) {
            fit$parameters$variance$sigmasq[i] <- fit$parameters$variance$sigmasq[1]
        }
        distr_fxn <- function(x) {
            dnorm(x, mean = fit$parameters$mean[i], sd = sqrt(fit$parameters$variance$sigmasq[i])) *
                fit$parameters$pro[i]
        }
        curve(distr_fxn(x), col = cols[i], add = TRUE, lwd = 1.5)
}

h <- hist(out$fits$S$data, prob = TRUE, col = "grey85", xlab = "S Score",
        main = "", breaks = 20, ylim=c(0,1.2), xlim=c(0,6))
is.G2M <- Buttner_meta == "G2M"
is.S <- Buttner_meta == "S"
points(out$fits$S$data[is.G2M], rep(1.1, sum(is.G2M)), pch=18, cex=0.5)
points(out$fits$S$data[is.S], rep(1, sum(is.S)), pch=18, cex=0.5)
text(4.5,1.1, "G2M", pos=4)
text(4.5,1, "S", pos=4)

fit = out$fits$S
cols=c("blue", "black", "orange")
for (i in 1:fit$G) {
        if (is.na(fit$parameters$variance$sigmasq[i])) {
            fit$parameters$variance$sigmasq[i] <- fit$parameters$variance$sigmasq[1]
        }
        distr_fxn <- function(x) {
            dnorm(x, mean = fit$parameters$mean[i], sd = sqrt(fit$parameters$variance$sigmasq[i])) *
                fit$parameters$pro[i]
        }
        curve(distr_fxn(x), col = cols[i], add = TRUE, lwd = 1.5)
}
dev.off()





