#### monocle analysis for 1,4,6； 1,4,6,7， 10
library(monocle)
library(tidyverse)

pbmc<-readRDS(file="pbmc2.rds")
counts<-pbmc@assays$RNA@counts  # raw data counts 更适合于monocle 的分析
pd <-  pbmc@meta.data

subclass<-flabel<-"sub.146"
selcol<-which(pd$seurat_clusters %in% c(1,4,6))
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


pdf(str_c("monocle/",subclass,"/",subclass,".plot.pdf"),width=4.5,height=4.8)
plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',cell_size = 0.6)
plot_cell_trajectory(Obj1, color_by = 'groups',cell_size = 0.6)
plot_cell_trajectory(Obj1, color_by = 'Cluster',cell_size = 0.6)
dev.off()

pdf(str_c("monocle/",subclass,"/",subclass,".plot.2.pdf"),width=8,height=6)
plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',cell_size = 0.6)+facet_wrap(~seurat_clusters, ncol=3,scales="free_y")
dev.off()

save(Obj1,file=str_c("monocle/",subclass,"/",subclass,".dt.2.Rdata"))

pdf(str_c("monocle/",subclass,"/",subclass,".traj.pdf"),width=6.5,height=4)
	plot_cell_trajectory(Obj1, color_by = 'seurat_clusters',show_tree=FALSE,cell_size = 0.8)
	plot_cell_trajectory(Obj1, color_by = 'Pseudotime',show_tree=FALSE,cell_size = 0.8)
	plot_cell_trajectory(Obj1, color_by = 'groups',show_tree=FALSE,cell_size = 0.3)
dev.off()
