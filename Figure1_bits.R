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


png("Figure1_mixtureG1S_Wang.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G1S, BIC=TRUE)
dev.off()

png("Figure1_mixtureG2M_Wang.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G2M, BIC=TRUE)
dev.off()

png("SupFigure1_mixtureG2M_Wang.png", width=4*1.5,height=3*1.5, units="in", res=300)
par(mfrow=c(1,2))
# G1S
h <- hist(out$fits$G1S$data, prob = TRUE, col = "grey85", xlab = "S Score",
        main = "", breaks = 20, ylim=c(0,9), xlim=c(0,0.7))
is.G2M <- out$phase == "G2M"
is.S <- out$phase == "G1S"
points(out$fits$G1S$data[is.G2M], rep(9, sum(is.G2M)), pch=18, cex=0.5)
points(out$fits$G1S$data[is.S], rep(8.5, sum(is.S)), pch=18, cex=0.5)
text(9,0.6, "G2M", pos=2)
text(8.5,0.55, "G1S", pos=2)

fit = out$fits$G1S
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

# G2M
h <- hist(out$fits$G2M$data, prob = TRUE, col = "grey85", xlab = "G2M Score",
        main = "", breaks = 20, ylim=c(0,13), xlim=c(0,0.9))
is.G2M <- out$phase == "G2M"
is.S <- out$phase == "G1S"
points(out$fits$G2M$data[is.G2M], rep(13, sum(is.G2M)), pch=18, cex=0.5)
points(out$fits$G2M$data[is.S], rep(12, sum(is.S)), pch=18, cex=0.5)
text(13,0.8, "G2M", pos=2)
text(12,0.75, "G1S", pos=2)

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

dev.off()


# Alternate Visualizations for Revision
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")

is.G2M <- out$phase == "G2M"
is.S <- out$phase == "G1S"

t1 = table(out$phase, out$fit$G2M$classification)
t2 = table(out$phase, out$fit$G1S$classification)
write.table(t1, "SupFigure1_G2MGMM_vs_phase_Wang.csv", sep=",")
write.table(t2, "SupFigure1_SGMM_vs_phase_Wang.csv", sep=",")

require(ggplot2)
df <- data.frame(scores=c(out$scores[,"G1S"], out$scores[,"G2M"]), phase=rep(c("G1S", "G2M"), each=nrow(out$scores)), GT=c(out$phase, out$phase))

png("SupFigure1_Wang_Scores_vs_Phase.png", width=4*0.75, height=3*0.75, units="in", res=300)
color_scheme <- phases_col[names(phases_col) %in% df$phase]

p <- ggplot(df, aes(x=GT, y=scores, fill=phase)) + theme_classic() + scale_fill_manual(values=color_scheme) + geom_boxplot() + labs(x = "Assigned Phase")
print(p)

dev.off()



############## Buttner ###########

source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/Read_Buttner_with_GT.R")
SCE <- SingleCellExperiment(assays=list(counts=Buttner_counts), rowData=data.frame(feature_symbol=rownames(Buttner_counts)))
logcounts(SCE) <- normalizeCounts(SCE, log=TRUE, pseudo_count=1, size.factors=colSums(Buttner_counts)/median(colSums(Buttner_counts)))

out <- classifyCells(SCE,MGeneSets$Cyclone)

png("Figure1_mixtureS_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$S, BIC=TRUE)
dev.off()

png("Figure1_mixtureG2M_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
plotMixture(out$fits$G2M, BIC=TRUE)
dev.off()

png("SupFigure1_mixtureG2M_Buttner.png", width=4*1.5,height=3*1.5, units="in", res=300)
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


# Alternate Visualizations for Revision
source("/home/tandrew6/scripts/CycleMix_Benchmark_git/CycleMix_Benchmark/ColourScheme.R")

is.G2M <- Buttner_meta == "G2M"
is.S <- Buttner_meta == "S"

t1 = table(Buttner_meta, out$fit$G2M$classification)
t2 = table(Buttner_meta, out$fit$S$classification)
write.table(t1, "SupFigure1_G2MGMM_vs_GT.csv", sep=",")
write.table(t2, "SupFigure1_SGMM_vs_GT.csv", sep=",")

require(ggplot2)
df <- data.frame(scores=c(out$scores[,"S"], out$scores[,"G2M"]), phase=rep(c("S", "G2M"), each=nrow(out$scores)), GT=c(Buttner_meta, Buttner_meta))

png("SupFigure1_Buttner_Scores_vs_Phase.png", width=4*0.75, height=3*0.75, units="in", res=300)
color_scheme <- phases_col[names(phases_col) %in% df$phase]

p <- ggplot(df, aes(x=GT, y=scores, fill=phase)) + theme_classic() + scale_fill_manual(values=color_scheme) + geom_boxplot() + labs(x = "Ground Truth")
print(p)

dev.off()


