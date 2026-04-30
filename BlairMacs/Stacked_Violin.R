library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


##### Prolif Mac Analysis #####
outname="Redo_Blair"
blair <- readRDS(paste(outname, "Clustered_redo.rds", sep="_"))
blair@meta.data$orig.clusters <- Idents(blair)
p_cells <- blair@reductions$umap@cell.embeddings[,1] > -5.5 & (blair$orig.clusters == 4 |blair$orig.clusters == 5)
blair@meta.data$fixed_clusters <- as.character(blair@meta.data$orig.clusters)
blair@meta.data$fixed_clusters[p_cells] <- "6"
Idents(blair) <- factor(blair@meta.data$fixed_clusters)
blair <- RenameIdents(blair, c("0" = "Kupffer", "1" = "Kupffer", "2" = "Activated", "3" = "Activated", "4"="LSEC", "5"="Unk", "6"="Kupffer2"))
png("Revisions_Figure_types.png", width=6*2.1/2, height=2*(2+2.4)/2, units="in", res=300)
DimPlot(blair)
dev.off()

png("Revisions_ViolinPlot.png", width=6, height=6*2, units="in", res=300)
p <- StackedVlnPlot(blair, features=c("STMN1","HMGN2", "TUBA1B", "HMGB2","NEAT1", "MALAT1", "B2M", "AHNAK"))
print(p)
dev.off()

png("Revisions_Dotplot.png", width=6, height=4, units="in", res=300)
p <- DotPlot(blair, features=c("STMN1","HMGN2", "TUBA1B", "HMGB2","NEAT1", "MALAT1", "B2M", "AHNAK"))
print(p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()


