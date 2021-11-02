# monocle analysis for pMN:1,2,6,11
library(monocle)
library(tidyverse)
library(colorRamps)
library(pheatmap)

pbmc<-readRDS(file="pbmc2.rds")
counts<-pbmc@assays$RNA@counts  # raw data counts 更适合于monocle 的分析
pd <-  pbmc@meta.data

subclass<-flabel<-"PMN"
selcol<-which(pd$seurat_clusters %in% c(1,2,6,11))
counts<-counts[,selcol]
pd<-pd[selcol,]

fData <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
colnames(pd)

pd <- new("AnnotatedDataFrame", data = pd) #
fd <- new("AnnotatedDataFrame", data =fData)

Obj <- newCellDataSet(
as(counts, "sparseMatrix"),
phenoData = pd,
featureData = fd,lowerDetectionLimit = 0.1,
expressionFamily = negbinomial.size() #gaussianff()
)

Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)

CDE <- differentialGeneTest(Obj[expressed_genes], fullModelFormulaStr = '~seurat_clusters', cores = detectCores()/ 2)
write_tsv(CDE,str_c("monocle/",subclass,"/",subclass,".cde.tsv"))

CDE<-arrange(CDE,pval)
ordering_genes <- row.names(CDE)[order(CDE$qval)][1:500]
Obj1 <- setOrderingFilter(Obj, ordering_genes = c(ordering_genes))
Obj1 <- reduceDimension(Obj1, method = 'DDRTREE') # ICA or DDRTREE # method 很重要。
Obj1<- orderCells(Obj1)

#num<-table(pData(Obj)$seurat_clusters,pData(Obj)$State)

pdf(str_c("monocle/",subclass,"/",subclass,".plot.pdf"),width=4.5,height=4.8)
plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',cell_size = 0.6)
plot_cell_trajectory(Obj1, color_by = 'groups',cell_size = 0.6)
plot_cell_trajectory(Obj1, color_by = 'Cluster',cell_size = 0.6)
dev.off()

pdf(str_c("monocle/",subclass,"/",subclass,".plot.2.pdf"),width=8,height=6)
plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',cell_size = 0.6)+facet_wrap(~seurat_clusters, ncol=3,scales="free_y")
dev.off()


save(Obj1,file=str_c("monocle/",subclass,"/",subclass,".dt.2.Rdata"))

source("/home/dell/project/general/theme.R")
pdf(str_c("monocle/",subclass,"/",subclass,".traj.pdf"),width=6.5,height=4)
	plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',show_tree=FALSE,cell_size = 0.8)
	plot_cell_trajectory(Obj1, color_by = 'Pseudotime',show_tree=FALSE,cell_size = 0.8)
	plot_cell_trajectory(Obj1, color_by = 'groups',show_tree=FALSE,cell_size = 0.3)
dev.off()


### plot ##
TF<-read_tsv("dre_TF.txt")
sCDE<-CDE[1:500,]
sCDE$Gene[sCDE$Gene %in% TF$Symbol]
seltf<-sCDE$Gene[sCDE$Gene %in% TF$Symbol]
seltfde<-sCDE[sCDE$Gene %in% TF$Symbol,]
write_tsv(seltfde,"sel.tf.de.info.2.tsv")
sel<-read_tsv("sel.tf.de.info.2.tsv")
sel<-sel[c(1:50),]

a<-plot_genes_branched_heatmap(Obj1.1[sel$gene_short_name,],
                                          branch_point = 1,
                                          num_clusters = 4,
                                          cores = 1,fontsize=12,
                                          use_gene_short_name = T,
                                         show_rownames = T,branch_labels = c("OPC/OL", "pMN>MN"),return_heatmap=T)
pdf(str_c("monocle/",subclass,"/",subclass,"sel.exp.hp.pdf"),family="Arial",width=6,height=8.0)
a$ph_res
dev.off()

clu<-a$annotation_row
clu$gene<-row.names(clu)
write_tsv(clu,str_c("monocle/PMN/PMN.tf.clu.tsv"))


a<-plot_genes_branched_heatmap(Obj1.1[sCDE$Gene,],
                                          branch_point = 1,
                                          num_clusters = 7,
                                          cores = 1,fontsize=12,
                                          use_gene_short_name = T,
                                         show_rownames = F,branch_labels = c("OPC/OL", "pMN>MN"),return_heatmap=T)
library(extrafont)
pdf(str_c("monocle/PMN/PMN.top500.exp.hp.pdf"),family="Arial",width=6,height=8.0)
a$ph_res
dev.off()

clu<-a$annotation_row
clu$gene<-row.names(clu)
clu$tf<-ifelse(clu$gene %in% sel$gene_short_name,"TF","non-TF")
table(clu$tf)
write_tsv(clu,str_c("monocle/PMN/PMN.top500.clu.tsv"))


blast_genes<-c("neurod1", "olig1", "mnx1", "nkx2.2a","isl1", "sox10")

p1<-plot_genes_in_pseudotime(Obj1[blast_genes,], color_by = "seurat_clusters",ncol=2)+theme(legend.position="right",legend.title= element_blank(),strip.background=element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),strip.text=element_text(size=12,colour="black"))+xlab("")

pdf(str_c("monocle/",subclass,"/",subclass,".pro.m.pseudo.pdf"),width=6,height=7)
p1
dev.off()


### vlnplot ###
cds_exprs <- exprs(cds_subset)
blast_genes<-c("myt1b","myt1a","NPAS3")
blast_genes[blast_genes %in% row.names(cds_exprs)]

plot_genes_violin(cds_subset[blast_genes,],
                  grouping = "seurat_clusters",
                  min_expr = 1.0,color_by = "seurat_clusters",log_scale=FALSE)+theme(legend.position="none",legend.title= element_blank(),strip.background=element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),strip.text=element_text(size=12,colour="black"))+xlab("")

spbmc<-subset(pbmc, seurat_clusters %in% c(6,1,2,11))
spbmc@meta.data$seurat_clusters<-factor(spbmc@meta.data$seurat_clusters,levels=c(6,1,2,11))

pdf(str_c("monocle/",subclass,"/","3genes.violin.pdf"),family="Arial",width=4,height=6)
VlnPlot(spbmc, features = blast_genes,pt.size=0,ncol=1,fill.by="feature",group.by="seurat_clusters")
dev.off()


