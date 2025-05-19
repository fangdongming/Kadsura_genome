library(Seurat)

setwd("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/05.hotspot/mingene40_2")

rds <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/04.ann/05.ann/04.ann/harmony__gene2000_8000_lambda1/data.cc.RDS")


matrix_ <- GetAssayData(rds,assay= "RNA",slot = "data")
transposed_matrix <- t(as.matrix(matrix_))
write.csv(transposed_matrix,"transposed_matrix.csv")
      
