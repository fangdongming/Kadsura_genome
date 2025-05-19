library(getopt)
arg <- matrix(c("input", "i", "1", "character", "input file", 
                "outdir", "o", "1", "character", "outdir", 
                "marker", "m", "1", "character", "marker"
),byrow = T, ncol = 5)

opt = getopt(arg)
if (is.null(opt$input)){
  opt$input <- "input"
}
if (is.null(opt$outdir)){
  opt$outdir <- "outdir"
}
library(Seurat)
library(ggplot2)
library(dplyr)
seurat_object <- readRDS(opt$input)
marker <- read.csv(opt$marker)

pdf("marker_S2.pdf", width = 20, height=10)
dp <- DotPlot(seurat_object, features = unique(marker$Gene_ID)) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
print(dp)
dev.off()

pdf("marker_dotplot_S2_leaf_720.pdf", width = 20, height=10)
match_index <- match(levels(dp$data$features.plot) , marker$Gene_ID)
levels(dp$data$features.plot) <- marker$Gene_Name[match_index]
dp$data$id <- factor(dp$data$id, levels = unique(marker$Cluster))
print(dp)
dev.off()

#pdf("sc_marker_FeaturePlot.pdf")
#FP <- FeaturePlot(seurat_object, features = unique(marker$Gene_ID), slot = "scale.data",label = T, combine = F, order = T)
#print(FP)
#dev.off()




