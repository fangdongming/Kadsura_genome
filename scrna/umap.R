library(SeuratObject)
library(Seurat)

setwd("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/05.hotspot/mingene40_2")

seurat_object <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/04.ann/05.ann/04.ann/harmony__gene2000_8000_lambda1/data.cc.RDS")

module_score <- read.csv("module_scores.csv", header = T)
module <- read.csv("modules.csv",header = T)

n <- max(module$Module)
head(module_score)
colnames(module_score) <- c("cell", paste("module",c(1:n), sep = ""))
rownames(module_score) <- module_score[,"cell"]

cells <- colnames(seurat_object)
module_score <- module_score[cells,]

for (i in 1:n+1) {

  col_name <- colnames(module_score)[i]
  
 
  seurat_object[[col_name]] <- module_score[[col_name]]
}

pdf("umap_module.pdf")
for (i in 1:n) {

  feature_name <- paste0("module", i)
  

  p <- FeaturePlot(seurat_object, features = feature_name)
  print(p)
  }

dev.off()

