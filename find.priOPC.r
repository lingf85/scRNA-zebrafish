
### find pri-OPC subcluster
library(tidyverse)
library(ggplot2)

meta<-read_tsv("metadata.tsv")  ## result of metabolic.pmn.R
data <-read.table("umap.txt",header=TRUE,sep="\t",quote = "\"",fill=TRUE,stringsAsFactors =FALSE)## result of metabolic.pmn.R
km     <- kmeans(data, centers=4, nstart=3)
ggdata <- as.data.frame(data)
ggdata$Cluster=str_c("c",km$cluster)

ggplot(data = ggdata,aes(x =UMAP_1, y = UMAP_2)) +geom_point(aes(color=Cluster))
sggdata<-ggdata[ggdata$Cluster=="c3",]
data<-sggdata[,1:2]
km     <- kmeans(data, centers=8, nstart=6)
temp<- as.data.frame(data)
temp$Cluster=str_c("c",km$cluster)
temp$Clu2<-rep("other")
temp$Clu2[temp$Cluster %in% c("c1")]<-"sel"

ggdata$Clu2<-ggdata$Cluster
ggdata$Clu2[row.names(ggdata) %in% row.names(temp[temp$Cluster %in% c("c1"),])]<-"sub.c1"

p= ggplot(data = ggdata,aes(x =UMAP_1, y = UMAP_2)) +geom_point(aes(color=Clu2))

#+ geom_text_repel(aes(label =label),size=5,color="deepskyblue4")+mytheme
meta<-meta[match(row.names(ggdata),meta$samples),]
table(row.names(ggdata)==meta$samples)

meta$Clu2<-meta$seurat_clusters_pri
meta$Clu2[ggdata$Clu2=="sub.c1"]<-"sub.c1"

ggdata$Clu3<-meta$Clu2
p= ggplot(data = ggdata,aes(x =UMAP_1, y = UMAP_2)) +geom_point(aes(color=Clu3))
write_tsv(meta,"metadata.addnew.tsv")