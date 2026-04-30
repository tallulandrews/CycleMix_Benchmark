Contamination_genes <- c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")
Hepato_Portal <- c("AGT", "ALB", "ALDOB", "APOA1", "CYP1A2", "CYP2A6", "CYP2A7", "CYP2D6", "FGA", "FGB", "GLS2", "HAL", "HAMP", "SDS")
Hepato_Central <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "CES1", "CYP2E1", "CYP3A4", "CYP3A5", "DCXR", "GSTA1", "OAT")
RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")
Endo_genes_SN_SC <- c("CTGF", "FCGR2B", "S100A13", "FCN2", "FCN3", "LYVE1", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
"F8", "CALCRL", "SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF",
"ENG", "PECAM1", "RAMP3", "INMT", "DNASE1L3", "LIFR", "TIMP3", "C7"
)
Endo_genes <- c("PECAM1","CLEC4M", "CD34", "COL1A2", "ENG", "STAB2", "CLDN5")
Endo_cvLSEC_vs_ppLSEC <- c("CTSD", "CTSL", "CLEC1B", "MS4A6A", 
"STAB1", "CLEC4G", "CRHBP", "DNASE1L3", "FCN2", "FCN3")
Endo_ppLSEC_vs_cvLSEC_portEndo <- c("MGP", "VIM", "ADIRF", "SPARCL1", "CLU", "S100A6", 
"CD9", "CLEC14A", "AQP1", "TM4DF1")
Endo_genes_dot <- unique(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
"SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF", 
"ENG", "PECAM1", "RAMP3", "DNASE1L3", "LIFR", "C7", "TIMP3",
"ACKR1", "WNT2", "RSPO3", "INMT", "PLAC8",
"PODXL", "RBP7", "PLVAP", "GSN", "CD34",
"CCL21", "S100A6", "CST3", "ALB", "APOA1"))
Macrophage_genes <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", 
"CD163", "HLA-DRA", "HLA-DPA1", "CD74", "FABP5",
"PLAUR", "LYZ", "S100A4", "S100A8", "VCAN", "FCN1")
Macrophage_genes_dot <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", "FCGR3A", "FCER1G",
"CD163", "HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "VSIG4", "CD74", 
"GPNMB", "ACTP5", "FABP5", "SPP1", "FABP4", "TREM2", "LGALS3", "CTSB", "PSAP", "APOE", "APOC1", "NPC2",
 "LYZ",  "S100A4", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1", "FTH1", "CD68", 
"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3", "LILRA5", "ILRB2", "PLAC8", "CD52", 
 "LST1", Contamination_genes)
Macrophage_genes_dot <- c(
"MARCO", "CD5L", "LYVE1", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1", 
"HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74", "VSIG4", 
"LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1",
"FABP5", "ACP5", "PLD3", "FTH1", "CD68", "APOE", "PSAP", "CSTB", "LGMN",
"RBP7", "FOLR2", "FCER1G", "MS4A7", "TIMP1",
"JUND", "FOS", "NFKBIA", "ACTG1", "CD14", "CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B",
"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3",
"PLAC8", "CD54", "LST1", "IFITM3", "AIF1", "COTL1",
"DNASE1L3", "FCN2", "CCL14", "FCN3", "SPARC", "CLEC1B", "ENG",
"ALB", "SERPINA1", "APOA1", "HP", "FGA")
label_genes <- function(genes, label) {names(genes) <- rep(label, length(genes)); return(genes)}
Final_Macrophage_genes <- c(
label_genes(c("AREG", "CCL3", "CD83", "CXCL2", "CXCL3", "IL1B", "NAMPT", "PLAUR", "SRGN", "THBS1"), "Activated"),
label_genes(c("C1QA", "C1QB", "C1QC", "CD163", "CD5L", "CTSB", "FTL", "HMOX1", "MARCO", "SEPP1", "SLC40A1", "VCAM1"), "Kupffer"),
label_genes(c("ACP5", "CSTB", "FABP5", "LGMN", "PLD3", "PSAP"), "LAM-like"),
label_genes(c("CD74", "HLA-DPA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1"), "MHCII"),
label_genes(c("FCN1", "LYZ", "MNDA", "S100A12", "S100A4", "S100A6", "S100A8", "S100A9", "VCAN"), "Monocyte"),
label_genes(c("FCER1G", "FOLR2", "MS4A7", "TIMD4", "TIMP1"), "Resident"))
NKT_genes <- c("CD8A", "CD3D", "TRAC", "TRBC2", "TRDC", "GNLY", "GZMB", "GZMA", "CCL5", 
"NKG7", "FCGR3A", "FGFBP2", "CD8B", "IL7R", "CD74", "HLA-DRB1")
NKT_exhaustion = c("PDCD1") 
NKT_Treg = c("ITGAE", "FOXP3")
NKT_genes_dot_old <- c("CD52","NCR1", "CXCR6", "CXCR3", "CD69", "PTPRC", "CD8A", "CD3D", "CD3E", "CD8B", "CD4", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
"IL32", "IL7R", "LTB", "IL17A", "IL18R1", "CD44", "KLRB1", "KLRC1", "KLRF1", "KLRK1", "CCL4", "CCL5", "NKG7",
"FCGR3A", "FGFBP2", "IL2RB", "GZMA", "GZMB", "CSF2", "GNLY", "KLRD1",
"CD74", "CD79A", "CD79B", "HLA-DRB1", "HLA-DRA", "AIF1", "PDCD1",
Prolif_genes, Contamination_genes, RBC_genes)
NKT_genes_dot <- c("CD3D", "CD3E", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
"CD8A", "CD8B", "CCL3", "CCL4", "CCL5", "IL7R", "LTB", "KLRB1", "TPT1",
"KLRC1", "KLRF1", "GZMK", "CMC1", "XCL1", "XCL2", "GZMB", "FCGR3A", "GNLY", "CXCR6", "CD69", "EOMES", "TBX21",
"CD79A", "CD79B", "CD74", "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DQB1",
"AIF1", "VSIG4", "LYZ", "VCAN", "CD163", "C1QC", "ITGAM", "ITGAE", "CST3",
"MKI67", "BRIC5", "TOP2A", "CDK1", "HMGB2", "HBB", "HBA1", "HBA2", "HBD",
"ALB", "SERPINA1", "APOA1", "FGA" 
)
require(Seurat)
require(ggplot2)
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/SubColour_Scheme.R")
flow_cyto_genes <- c("CD68", "MRC1", "CD14", "CD274", "PTPRC", "CD86", "CD163", "IL2RA", "FOXP3",
"NCAM1", "NCR1", "CD3E", "PTPRC", "CD8A", "IL7R", "ICOS", "KLRD1", "IL2RA",
"CD4", "PDCD1", "SELL", "HAVCR2", "CCR3", "CD3E", "IL7R", "CD8A", "PTPRC", 
"LAG3", "CTLA4", "CD27", "CD4")
outname="Redo_Blair"
blair <- readRDS(paste(outname, "Clustered_redo.rds", sep="_"))
DotPlot(blair, feature=c(Final_Macrophage_genes, label_genes(Endo_genes, "Endo"), "APOE", "IPMK", "COX5A", "H3F3B"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
tmp <- RenameIdents(blair, c("0" = "Kupffer", "1" = "Kupffer", "2" = "Activated", "3" = "Activated", "4"="LSEC", "5"="Unk"))
colour <-c("Degraded Cells"="#ef1a1c","Neutrophils"="#377eb8","Kupffer Cells"="#4daf4a","Erythrocytes"="#984ea3","Endothelial Cells"="#ff7f00","LSECs"="#ffff33","Hepatocytes"="#f781bf","NA"="#969696")
mousemap <- readRDS("new_mouse_map_obj.rds")
table(Idents(mousemap))
50/1628
require(CycleMix)
tmp <- as.SingleCellExperiment(blair)
require(scater)
assays(tmp)$logcounts <- scater::normalizeCounts(tmp) <- scater::normalizeCounts(tmp)
rowData(tmp)$feature_symbol <- rownames(tmp)
assay(tmp, "logcounts") <- scater::normalizeCounts(tmp)
rowData(tmp)$feature_symbol <- rownames(tmp)
SCE <- tmp
CC_table <- HGeneSets$Tirosh
expr_name = "logcounts"; do.scale = FALSE
symbol_column = "feature_symbol"
allow.multi = FALSE
assignPhase2 <- function (SCE, CC_table, phase = "G2M", expr_name = "logcounts", 
    do.scale = FALSE, symbol_column = "feature_symbol") 
{
    if (class(SCE)[1] != "SingleCellExperiment") {
        stop("Error: Requires SingleCellExperiment object as input")
    }
    signature <- CC_table[, "Stage"] == phase
    if (sum(signature) <= 1) {
        stop(paste("Error: Insufficient genes associated with phase:", 
            phase))
    }
    signature <- CC_table[signature, ]
    exprmat <- SummarizedExperiment::assays(SCE)[[expr_name]]
    if (do.scale) {
        exprmat <- t(apply(exprmat, 1, scale))
    }
    if (is.null(symbol_column)) {
        rowData(SCE)$CycleMixSym <- rownames(SCE)
        symbol_column <- "CycleMixSym"
    }
    gene_names <- SummarizedExperiment::rowData(SCE)[, symbol_column]
    signature <- signature[signature[, "Gene"] %in% gene_names, 
        ]
    matches <- base::match(signature[, "Gene"], gene_names)
    exprmat <- exprmat[matches, ]
    gene_names <- gene_names[matches]
    keep <- !is.na(gene_names)
    exprmat <- exprmat[keep, ]
    gene_names <- gene_names[keep]
    signature <- signature[keep, ]
    score <- colSums(exprmat * signature[, "Dir"])/sum(abs(signature[, 
        "Dir"]))
    fit <- mclust::Mclust(score, G = 1:3)
    fit$phase <- as.character(fit$classification)
    if (fit$G > 1) {
        tag <- which(fit$parameters$mean == max(fit$parameters$mean))
        fit$phase[fit$phase == tag] <- phase
    }
    fit$phase[fit$phase != phase] <- ""
    return(fit)
}
    if (class(SCE)[1] != "SingleCellExperiment") {
        stop("Error: Requires SingleCellExperiment object as input")
    }
    out_list <- list()
    stages <- as.character(unique(CC_table[, "Stage"]))
    phases <- matrix(nrow = ncol(SCE), ncol = length(stages))
    scores <- matrix(nrow = ncol(SCE), ncol = length(stages))
    for (i in 1:length(stages)) {
        assignment <- assignPhase2(SCE, CC_table, phase = stages[i], 
            do.scale = do.scale, expr_name = expr_name, symbol_column = symbol_column)
        out_list[[stages[i]]] <- assignment
        phases[, i] <- assignment$phase
        scores[, i] <- assignment$data
    }
    if (allow.multi) {
        best <- apply(phases, 1, paste, collapse = "")
    }
    else {
        if (!do.scale) {
            scores <- apply(scores, 2, scale)
        }
        else {
            scores <- apply(scores, 2, scale, center = FALSE)
        }
        best <- sapply(1:nrow(scores), function(i) {
            phases[i, which(scores[i, ] == max(scores[i, ]))]
        })
    }
    best[best == ""] <- "None"
    scores <- as.matrix(scores)
    colnames(scores) <- stages
 if (allow.multi) {
        best <- apply(phases, 1, paste, collapse = "")
    } else {
        if (!do.scale) {
            scores <- apply(scores, 2, scale)
        }
        else {
            scores <- apply(scores, 2, scale, center = FALSE)
        }
        best <- sapply(1:nrow(scores), function(i) {
            phases[i, which(scores[i, ] == max(scores[i, ]))]
        })
    }
    best[best == ""] <- "None"
    scores <- as.matrix(scores)
    colnames(scores) <- stages

CM <- best
blair@meta.data$cyclemix <- CM
blair <- CellCycleScoring(blair, cc.genes$s.genes, cc.genes$g2m.genes)
DimPlot(blair, group.by="cyclemix")
DimPlot(blair, group.by="Phase")
phases_col <- c("G1"="grey50", "G1S"="#7570b3", "None"="grey50", "S"="#7570b3", "G2M"="#d95f02", "G2"="#d95f02")
require(ggplot2)
png("ForTalk_Blair_Macs_Cyclemix.png", width=6*2.1*4/10, height=2*(2+2.4)*4/10, units="in", res=300)
DimPlot(blair, group.by="cyclemix") + scale_color_manual(values=phases_col)
dev.off()
png("ForTalk_Blair_Macs_Seurat.png", width=6*2.1*4/10, height=2*(2+2.4)*4/10, units="in", res=300)
DimPlot(blair, group.by="Phase") + scale_color_manual(values=phases_col)
dev.off()
table(blair@meta.data$cyclemix)
62/ncol(blair)
table(blair@meta.data$Phase)
(128+614)/ncol(blair)
blair <- SCTransform(blair)
FeaturePlot(blair, c("STMN1", "H2AZ1", "HMGN2", "TUBA1B", "HMGB2","NEAT1", "MALAT1", "B2M", "AHNAK"), ncol=4, cols=c("lightgrey", "blue", "black"))


##### Prolif Mac Analysis #####
blair <- readRDS(paste(outname, "Clustered_redo.rds", sep="_"))
blair@meta.data$orig.clusters <- Idents(blair)
p_cells <- blair@reductions$umap@cell.embeddings[,1] > -5.5 & (blair$orig.clusters == 4 |blair$orig.clusters == 5)
blair@meta.data$fixed_clusters <- as.character(blair@meta.data$orig.clusters)
blair@meta.data$fixed_clusters[p_cells] <- "6"
Idents(blair) <- factor(blair@meta.data$fixed_clusters)
blair <- RenameIdents(blair, c("0" = "Kupffer", "1" = "Kupffer", "2" = "Activated", "3" = "Activated", "4"="LSEC", "5"="Unk", "6"="Kupffer2"))
png("TalkFigure_types.png", width=6*2.1/2, height=2*(2+2.4)/2, units="in", res=300)
DimPlot(blair)
dev.off()

png("Revisions_ViolinPlot.png", width=6, height=6, units="in", res=300)
p <- VlnPlot(blair, features=c("STMN1","HMGN2", "TUBA1B", "HMGB2","NEAT1", "MALAT1", "B2M", "AHNAK"), cols=c("grey35", "grey35", "grey35", "grey35", "#d95f02", "#d95f02", "#d95f02", "#d95f02"))
print(p)
dev.off()



colour <-c("Degraded Cells"="#ef1a1c","Neutrophils"="#377eb8","Kupffer Cells"="#4daf4a","Erythrocytes"="#984ea3","Endothelial Cells"="#ff7f00","LSECs"="#ffff33","Hepatocytes"="#f781bf","NA"="#969696")

Idents(blair) <- factor(blair@meta.data$clusters)
#out <- FindAllMarkers(blair, group.by="clusters", logfc.threshold=0)
out <- FindMarkers(blair, group.by="clusters", ident.1="7", logfc.threshold=0)
saveRDS(out, "CC_Mac_markers.rds")

genes <- rownames(out[out$avg_log2FC > 0.5,])
 require(gprofiler2)
rich <- gost(genes, custom_bg=rownames(blair))

plot(blair@reductions$umap@cell.embeddings); abline(v=-5.5)
