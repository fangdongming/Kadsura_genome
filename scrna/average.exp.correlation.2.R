#library(Seurat)
#library(ggplot2)
#rds <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC_flower_sigleCell/03.cc_soupX_cluster/freeze_harmony_1/harmony__gene2000_8000_lambda1/data.cc.RDS")
#cell.num <- data.frame(table(rds@meta.data[["seurat_clusters"]]))
#write.table(cell.num, paste0("seurat_clusters",".cell.num.txt"))
#Idents(rds) <- rds@meta.data[["seurat_clusters"]]
#ave.exp <- AverageExpression(rds, group.by = "seurat_clusters", slot = "data")[["RNA"]]
#write.table(ave.exp, paste0("seurat_clusters", ".average.expression.txt"))
#cor.exp <- as.data.frame(cor(ave.exp))
#cor.exp[["x"]] <- rownames(cor.exp)
#cor.df <- tidyr::gather(data = cor.exp, y, correlation, -x)
#p1 <- ggplot(cor.df, aes(x, y, fill = correlation)) +
#  geom_tile() + scale_fill_gradient(high = "red", low = "white")
#p1$"data"$"x" <- forcats::fct_inorder(p1$"data"$"x")
#p1$"data"$"y" <- forcats::fct_inorder(p1$"data"$"y")
#pdf(paste0("seurat_clusters", ".correlation.heatmap.pdf"), width = 9.55, height = 8)
#p1
#dev.off()
library(Seurat)  
library(ggplot2)  
library(tidyr)  
library(forcats)  

# 读取两个 RDS 文件  
rds1 <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC_flower_sigleCell/03.cc_soupX_cluster/freeze_harmony_1/harmony__gene2000_8000_lambda1/data.cc.RDS")  
rds2 <- readRDS("/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/04.ann/05.ann/04.ann/harmony__gene2000_8000_lambda1/data.cc.RDS")  

# 计算每个聚类的平均表达量  
ave.exp1 <- AverageExpression(rds1, group.by = "seurat_clusters", slot = "data")[["RNA"]]  
ave.exp2 <- AverageExpression(rds2, group.by = "seurat_clusters", slot = "data")[["RNA"]]  

# 合并两个数据框  
combined.exp <- merge(ave.exp1, ave.exp2, by = "row.names", suffixes = c("_data1", "_data2"))  
rownames(combined.exp) <- combined.exp$Row.names  
combined.exp$Row.names <- NULL  

# 计算相关性  
cor.exp <- as.data.frame(cor(combined.exp))  
cor.exp[["x"]] <- rownames(cor.exp)  
cor.df <- gather(data = cor.exp, y, correlation, -x)  

# 创建热图  
p1 <- ggplot(cor.df, aes(x, y, fill = correlation)) +  
  geom_tile() + scale_fill_gradient(high = "red", low = "white")  

# 重新排序因子  
p1$"data"$"x" <- fct_inorder(p1$"data"$"x")  
p1$"data"$"y" <- fct_inorder(p1$"data"$"y")  

# 保存热图为 PDF  
pdf("combined_correlation_heatmap.pdf", width = 181818181818181818181818181818181818, height = 8)  
print(p1)  
dev.off()
