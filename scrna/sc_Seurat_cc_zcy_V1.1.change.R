library(getopt)
arg <- matrix(c("input", "i","1","character","input file1",
                "outdir","o","1","character","outdir",
                "sample","s","1","character","sample,default=Maize",
                "tissue","t","1","character","tissue,default=Embro",
                "cellcycle","c","1","character","either remove the effect of cell cycle genes or not",
		"regressall","a","1","character","regress all the effects or half",
		"cc.genes","e","1","character","input cell cycle genes list",
		"cc.colname","b","1","character","input cell cycle genes colname",
		"dims","d","1","integer","dims option for FindNeighbors,default=15",
                "resolution","r","1","numeric","resolution option for FindClusters,[0.4-1.2],default=1",
                "help","h","0","logical", "Usage: Rscript runiDrop.R -i <input> -o <outdir> [-s SAMPLE -t TISSUE]",
                "minCG","m","1","integer","minimum gene number for each cell, i.e. nFeature_RNA, default=500",
                "nC_RNA","n","1","integer"," i.e. nCount_RNA, default=1500",
                "rds","f","0","logical", "Save the RDS file",
                "lambda","l","1","numeric"," RUNharmony lambda, default=1",
                "featurePlot","p","1","character","T or F",
                "marker.colname","k","1","character","marker genes colnames",
		"features","g","1","character","plot genes expression"
		),byrow=T,ncol=5)
opt = getopt(arg)
if (is.null(opt$sample)){
        opt$sample <- "Maize"
}
if (is.null(opt$tissue)){
        opt$tissue <- "stem"
}
if (is.null(opt$outdir)){
        opt$outdir <- "output"
}
if (is.null(opt$cellcycle)){
        opt$cellcycle <- F
}
if (is.null(opt$regressall)){
        opt$regressall <- T
}

if (is.null(opt$dims)){
        opt$dims <- 15
}
if (is.null(opt$resolution)){
        opt$resolution <- 1
}
if (is.null(opt$minCG)){
    opt$minCG <- 200
}
if (is.null(opt$nC_RNA)){
    opt$nC_RNA <- 0
}
if (is.null(opt$lambda)){
    opt$lambda <- 1
}
if (is.null(opt$featurePlot)){
 opt$featurePlot <- F
}

library(dplyr)
library(tidyr)
library(harmony)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)


data.harmony <- readRDS(opt$input)
       #input cell cycle genes
cc.col <- opt$cc.colname
cc.genes <- read.csv(opt$cc.genes, header = T)
s.genes <- cc.genes[which(cc.genes[,1]=="S"), cc.col]
s.genes <- unique(unlist(strsplit(s.genes, ",")))
#for(i in 1:length(s.genes)){s.genes[i] <- paste(s.genes,".v3.2")}
g2m.genes <- cc.genes[which(cc.genes[,1] =="G2M"), cc.col]
g2m.genes <- unique(unlist(strsplit(g2m.genes, ",")))
#for(i in 1:length(g2m.genes)){g2m.genes[i] <- paste(g2m.genes,".v3.2")}
##Assign Cell-Cycle Scores
data.harmony <- CellCycleScoring(data.harmony, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
      
pdf("03_cc_umap.cell.stage.pdf", , height = 8, width = 16)
DimPlot(data.harmony, reduction = "umap", group.by = c("Phase","seurat_clusters"),pt.size = 0.5, label = TRUE, repel = TRUE)
dev.off()

#view cell cycle scores and phase assignments
head(data.harmony[[]])

 # Visualize the distribution of cell cycle markers across
#RidgePlot(data.harmony, features = c("AT1G07370", "AT2G29570", "AT3G23890", "AT5G10400"), ncol = 2)


if(opt$regressall){
#Regress out cell cycle scores during data scaling
data.harmony <- ScaleData(data.harmony, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data.harmony))
} else {
data.harmony$CC.Difference <- data.harmony$S.Score - data.harmony$G2M.Score
data.harmony <- ScaleData(data.harmony, vars.to.regress = "CC.Difference", features = rownames(data.harmony))
}

data.harmony <- RunPCA(data.harmony, features = VariableFeatures(data.harmony), nfeatures.print = 10)
        # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
       #data.harmony <- RunPCA(data.harmony, features = c(s.genes, g2m.genes))
        #DimPlot(data.harmony)
		
        ###reharmony
data.harmony <- RunHarmony(data.harmony, 'batch', plot_convergenc=TRUE,lambda = opt$lambda)
data.harmony <- RunUMAP(data.harmony, reduction = "harmony", dims = 1:opt$dims)
data.harmony <- FindNeighbors(data.harmony, reduction = "harmony", dims = 1:opt$dims)
data.harmony <- FindClusters(data.harmony, resolution = opt$resolution)

             
pdf("03_cc_umap.harmony.cc.pdf", height = 8, width = 24)
p3 <- DimPlot(data.harmony, reduction = "umap", group.by = c("Phase","seurat_clusters","old.ident"), label = TRUE, repel = TRUE)
p3
dev.off()
pdf("03_cc_umap.harmony.cc.nfeature.pdf", height = 8, width = 16)
plot1 <- DimPlot(object = data.harmony, reduction = "umap", label = T,pt.size=0.5,group.by = "ident")
plot1.1 <- FeaturePlot(data.harmony, features = "nFeature_RNA", reduction = "umap", label = T, combine = F, order = T)
plot1.2 <- DimPlot(object = data.harmony, reduction = "umap", label = T,pt.size=0.5,group.by = "doublet_info")
dev.off()
	
plot <- DimPlot(object = data.harmony, reduction = "umap", label = T,pt.size=0.5,group.by = c("batch","batch"))
plot.s <- DimPlot(object = data.harmony, reduction = "umap", pt.size=0.5,group.by = "ident", split.by = "batch")
pdf("04.1_cc_umap.metadata.pdf", width = 16, height=8)
plot
dev.off()

pdf("04.2_cc_umap.metadata.split.pdf", width = 16, height=8)
plot.s
dev.off()


data.harmony
saveRDS(data.harmony,"data.cc.RDS")


pdf("05_cc_cell_number.pdf", width = 12, height = 8)
dat<-as.data.frame(table(data.harmony@active.ident))
p<-ggplot(dat,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position = "dodge") + geom_text(aes(label=Freq,vjust=-0.5) )
p
dev.off()

write.table(dat,file = "05_cc_cell_number.txt",sep="\t",quote=F,row.names = F)

###all marker genes
DefaultAssay(data.harmony) <- "RNA"
adata.markers <- FindAllMarkers(data.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allmakers <- adata.markers %>% group_by(cluster)
write.csv(allmakers,file = "06_cc.markergenes_list.csv")

###top10 marker genes heatmap
pdf('06_cc.top10_markergenes.heatmap.integreted.pdf',width = 8 , height = 8 )
if("avg_log2FC" %in% colnames(adata.markers)){
top10 <- adata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
}else{top10 <- adata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)}
DoHeatmap(data.harmony, features = top10$gene)
dev.off()


#?????÷???ùòòμ?±í′?μè??í?
#markergeneheatmap <- read.table('13_markergeneheatmap.txt',header = T)
#p <- VlnPlot(data.harmony, features = markergeneheatmap$gene[1:20])
#ggsave("marker_gene_20.png",p,dpi=600,width=40,height=28)

fplot1 <- c()
heatmap <- c()
DefaultAssay(data.harmony) <- "RNA"
if(opt$featurePlot){
features <- read.csv(opt$features,header = T)
col <- opt$marker.colname
genes <- features[,col]
genes <- unique(unlist(strsplit(genes, ",")))
genes <- unique(unlist(strsplit(genes, " ")))
#for(i in 1:length(genes)){genes[i] <- paste(genes,".v3.2")}
genes <- genes[which(genes%in%rownames(data.harmony))]

h <- length(genes)/5
pdf("07.cc.input.feature.heatmap.pdf",width = 20, height = h)
heatmap <- DoHeatmap(data.harmony, features = genes, slot = "scale.data",)
print(heatmap)
dev.off()

pdf("07.cc.input.featurePlot.pdf")
fplot1 <- FeaturePlot(data.harmony, features = genes, slot = "scale.data",label = T, combine = F, order = T)
print(fplot1)
dev.off()

}

