
##  Metabolic related gene of subcluster(1,2,6,11)
library(Seurat)
library(openxlsx)
library(tidyverse)
library(ComplexHeatmap)
pbmc<-readRDS("pbmc2.rds")
pbmc<-subset(pbmc,seurat_clusters%in%c(1,2,6,11))
DefaultAssay(pbmc)<-"RNA"
expdata<-GetAssayData(pbmc,slot="counts")
selm<-read.xlsx("dre.concern.gl.GO.detail.xlsx")
selm<-subset(selm,sel==1)
marker<-read_tsv("./marker/PMN/marker.all.tsv")
selm<-selm[which(selm$Gene.symbol%in%marker$GeneSymbol),]
cluster_info <- sort(pbmc$seurat_clusters)

for(i in 1:length(unique(selm$Type))){
temp<-subset(selm,Type%in%unique(selm$Type)[i])
temp_data<-as.matrix(expdata[which(rownames(expdata)%in%temp$Gene.symbol),names(cluster_info)])
pdf(str_c("./plots/new/PMN.",unique(selm$Type)[i],".heatmap.pdf"),width=4+0.5*4,height=8)
p<-Heatmap(temp_data,cluster_rows = FALSE,cluster_columns = FALSE,
		show_column_names = FALSE,show_row_names = FALSE)
plot(p)
dev.off()
}


selm<-selm[which(selm$Gene.symbol%in%marker$GeneSymbol),]
oligos.integrated<- RunPCA(object=pbmc,assay="integrated", fearures=selm$Gene.symbol)
oligos.integrated <- RunUMAP(oligos.integrated, dims = 1:30)
oligos.integrated <- RunTSNE(oligos.integrated, dims = 1:30)

spbmc <- FindNeighbors(oligos.integrated, dims = 1:30,reduction = "pca",nn.eps=0,graph.name="new.PMN")
spbmc <- FindClusters(object = spbmc, resolution =1,n.start = 100,algorithm=2,graph.name="new.PMN") # resolution 值越大，cluster越多 1 = original Louvain algorithm; 2 = Louvain algorithm

spbmc$seurat_clusters_pri<-oligos.integrated$seurat_clusters
write_tsv(spbmc@meta.data,file="./result/metadata.tsv")
write.table(spbmc@reductions$umap@cell.embeddings,file="./result/umap.txt",sep="\t")

pdf(str_c("./plots/new/PMN",".umap.cluster3.pdf"),width=4,height=4)
DimPlot(spbmc, reduction = "umap",label = TRUE)
Idents(spbmc)<-"seurat_clusters_pri"
DimPlot(spbmc, reduction = "umap",label = TRUE)
dev.off()



flabel="PMN_sub2"

meta<-read_tsv(file="./result/metadata.addnew.tsv") ### result of find.priOPC.r

dent<-meta$Clu2
names(dent)=meta$samples
spbmc@meta.data$Clu2<-dent[as.character(spbmc@meta.data$samples)]

Idents(spbmc)<-"Clu2"
p2<- DimPlot(spbmc, reduction = "umap",label = TRUE)
p3 <- DimPlot(spbmc, group.by = c("groups"), combine = TRUE)
p2.1<- DimPlot(spbmc, reduction = "tsne",label = TRUE)
p3.2 <- DimPlot(spbmc, split.by = c("groups"), combine = TRUE)

inte_type="int"

ggsave(p2,device="pdf",filename=str_c("plots/new/",flabel,".",inte_type,".UMAP.plot.pdf"),width=6,height=6)
ggsave(p2.1,device="pdf",filename=str_c("plots/new/",flabel,".",inte_type,".TSNE.plot.pdf"),width=6,height=6)
ggsave(p3,device="pdf",filename=str_c("plots/new/",flabel,".",inte_type,".UMAP.Sample.plot.pdf"),width=6,height=6)
ggsave(p3.2,device="pdf",filename=str_c("plots/new/",flabel,".",inte_type,".UMAP.Sample.2.plot.pdf"),width=14,height=6)  ### non.int 中这个图如果没有问题就直接用non.int 方法
