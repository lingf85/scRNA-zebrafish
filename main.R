library(Seurat)
library(R.utils)
library(stringr)
library(stringr)
library(scales)
library(RColorBrewer)

### reading cellranger data then normalized data and find cell sub types.
paths<-list.files("raw_data/")
setwd("raw_data/")
SC <- Read10X(data.dir = paths,  unique.features = TRUE)

groups <- strsplit(unlist(SC@Dimnames[2]),"_")
samples <- unlist(lapply(groups,function(x){
            if(is.na(x[2]))
            {
              t=x[1]
              x[1]="28h"
              x[2]=t
            }else if(x[1]=="2"){
              x[1]="42h"              
            }else if(x[1]=="3"){
              x[1]="60h"
            }
      paste0(x[1],"_",x[2])}))
groups<-strsplit(samples,"_")
groups <- unlist(lapply(groups,function(x) x[1]))


emat_10x <- SC

samples <- colnames(SC)
anno_10x <- samples
anno_10x <- as.data.frame(anno_10x,stringsAsFactors = FALSE)
colnames(anno_10x) <- "samples"
rownames(anno_10x)<-samples
anno_10x$groups<-groups

dt<- CreateSeuratObject(emat_10x,meta.data=anno_10x,min.cells = 10, min.features = 200)

dt@meta.data$batchs<-dt@meta.data$groups
dt@meta.data$batchs[which(dt@meta.data$batchs=="42h")]<-"batch1"
dt@meta.data$batchs[which(dt@meta.data$batchs=="28h")]<-"batch2"
dt@meta.data$batchs[which(dt@meta.data$batchs=="60h")]<-"batch2"
dt <- PercentageFeatureSet(dt, pattern = "(^mt-|^MT-)", col.name = "percent.mt")###线粒体
dt <- PercentageFeatureSet(dt, pattern = "(^Rp[sl]|^RP[SL])", col.name = "percent.ribo")###核糖体
summary(dt@meta.data)
dim(dt)


selpheno<-dt@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")]
sta<-apply(selpheno,2,quantile)
write.table(as.data.frame(sta),file=str_c("result/dt",".qc.txt"),col.names=TRUE,sep="\t",quote=F,append=FALSE,row.names=TRUE) 

num<-as.data.frame(table(dt@meta.data$groups))
colnames(num)<- c("groups","numbers")
write.table(num,file=str_c("result/dt",".cellnumber.txt"),col.names=TRUE,sep="\t",quote=F,append=FALSE,row.names=TRUE) 

pdf(str_c("plots/","dt",".qc.plot.group.pdf"))
VlnPlot(dt, group.by = "groups",features = c("percent.ribo"), pt.size = 0)
VlnPlot(dt, group.by = "groups",features = c("percent.mt"), pt.size = 0)
VlnPlot(dt, group.by = "groups",features = c("nFeature_RNA"), pt.size =0)
VlnPlot(dt, group.by = "groups",features = c("nCount_RNA"), pt.size =0)
dev.off()

dt2<- subset(dt, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)##细胞过滤

selpheno2<-dt2@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo")]
sta2<-apply(selpheno2,2,quantile)
write.table(as.data.frame(sta2),file=str_c("result/dt2",".qc.txt"),col.names=TRUE,sep="\t",quote=F,append=FALSE,row.names=TRUE) 

num<-as.data.frame(table(dt@meta.data$groups))
colnames(num)<- c("groups","numbers")
num2<-as.data.frame(table(dt2@meta.data$groups))
num2<-cbind(num2,percent(1-num2[,2]/num[,2]))
colnames(num2)<- c("groups","num_dt2","delete")
write.table(num2,file=str_c("result/dt2",".cellnumber.txt"),col.names=TRUE,sep="\t",quote=F,append=FALSE,row.names=TRUE) 

pdf(str_c("plots/","dt2",".qc.plot.group.pdf"))
VlnPlot(dt2, group.by = "groups",features = c("percent.ribo"), pt.size = 0)
VlnPlot(dt2, group.by = "groups",features = c("percent.mt"), pt.size = 0)
VlnPlot(dt2, group.by = "groups",features = c("nFeature_RNA"), pt.size =0)
VlnPlot(dt2, group.by = "groups",features = c("nCount_RNA"), pt.size =0)
dev.off()


inte_type="int"
oligos<-dt2
oligos.list <- SplitObject(oligos, split.by = "batchs2")
for (i in 1:length(oligos.list)) {
  oligos.list[[i]] <- SCTransform(oligos.list[[i]], verbose = FALSE,vars.to.regress=c("percent.mt","percent.ribo"))
}
oligos.features <- SelectIntegrationFeatures(object.list = oligos.list, nfeatures = 3000)
oligos.list <- PrepSCTIntegration(object.list = oligos.list, anchor.features = oligos.features,verbose = FALSE)
oligos.anchors <- FindIntegrationAnchors(object.list = oligos.list, normalization.method = "SCT", anchor.features = oligos.features, verbose = FALSE)
oligos.integrated <- IntegrateData(anchorset = oligos.anchors, normalization.method = "SCT", verbose = FALSE)

oligos.integrated <- RunPCA(oligos.integrated, verbose = FALSE)
oligos.integrated <- RunUMAP(oligos.integrated, dims = 1:30)
oligos.integrated <- RunTSNE(oligos.integrated, dims = 1:30)
ElbowPlot(oligos.integrated)
pbmc <- FindNeighbors(oligos.integrated, dims = 1:30,reduction = "pca",nn.eps=0)
pbmc <- FindClusters(object = pbmc, resolution =0.3,n.start = 100,algorithm=2)

p2<- DimPlot(pbmc, reduction = "umap",label = TRUE)
p3 <- DimPlot(pbmc, group.by = c("groups"), combine = TRUE)
p2.1<- DimPlot(pbmc, reduction = "tsne",label = TRUE)
p3.2 <- DimPlot(pbmc, split.by = c("groups"), combine = TRUE)

library(ggplot2)
ggsave(p2,device="pdf",filename=str_c("plots/","p2",".UMAP.plot2.pdf"),width=6,height=6)
ggsave(p2.1,device="pdf",filename=str_c("plots/","p2.1",".","TSNE.plot2.pdf"),width=6,height=6)
ggsave(p3,device="pdf",filename=str_c("plots/","p3",".","UMAP.Sample.plot2.pdf"),width=6,height=6)
ggsave(p3.2,device="pdf",filename=str_c("plots/","p3.1",".","UMAP.Sample.2.plot2.pdf"),width=14,height=6)  ### non.int 中这个图如果没有问题就直接用non.int 方法

saveRDS(pbmc,file=str_c("pbmc2",".rds"))


### plot ratio between groups
plot_ratio<-function(num2,label,col) {
  
  if(row.names(num2)[1]!="num_all"){
    num2<-rbind(apply(num2,2,sum),num2)
  }
  ratio<-apply(num2,2,function(x) x[-1]/x[1])
  ratio2<-as.data.frame(ratio);ratio2$celltype<-row.names(ratio2)
  ratio<-melt(ratio)
  colnames(ratio)[c(1,2)]<-c("Var1","Var2")
  p2<-ggplot(ratio, aes(x=Var1,y=value))+geom_bar(aes(fill=Var2),position="fill",stat="identity")+
    mytheme+theme(legend.position="right",legend.title= element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),
                  strip.background=element_blank(),axis.text.x=element_text(angle = 45,vjust =1, hjust = 1))+
    ylab("")+xlab("")+scale_fill_manual(values=col)+
    geom_hline(yintercept=0.5,linetype = "dotted")
  ggsave(p2,filename=str_c("plots/",label,".ratio.pdf"),device="pdf",width=7,height=7)
  return(ratio2)
}
meta2<-pbmc@meta.data
num2<-table(meta2$seurat_cluster,as.character(meta2$groups))
row.names(num2)<-str_c("c",row.names(num2))
flabel="seurat_cluster_groups"
ratio2<-plot_ratio(num2,label=str_c(flabel,".pricell"),col=brewer.pal(6,"Set3"))


pbmc@meta.data$samples2<-str_replace(pbmc@meta.data$samples,"2_|3_","")
spbmc1<-subset(pbmc,groups=="42h")
library(dplyr)
meta2<-left_join(spbmc1@meta.data,pbmc_42h@meta.data,by=c("samples2"="si"))

num2<-table(meta2$seurat_clusters.x,as.character(meta2$leiden))
row.names(num2)<-str_c("c",row.names(num2))
flabel="seurat_cluster_leiden"
ratio2<-plot_ratio(num2,label=str_c(flabel,".pricell"),col=color_all[seq(0,length(color_all),2)])             

## quality control
spbmc<-subset(pbmc,seurat_clusters!="13")

pdf(str_c("plots/","spbmc",".qc.plot.group.pdf"))
VlnPlot(spbmc, group.by = "seurat_clusters",features = c("percent.ribo"), pt.size = 0)
VlnPlot(spbmc, group.by = "seurat_clusters",features = c("percent.mt"), pt.size = 0)
VlnPlot(spbmc, group.by = "seurat_clusters",features = c("nFeature_RNA"), pt.size =0)
VlnPlot(spbmc, group.by = "seurat_clusters",features = c("nCount_RNA"), pt.size =0)
dev.off()

### cell score 
org="dre"
if(org=="hsa" | org=="dre"){
  require(Seurat)
  data(cc.genes)
} else {
  cc.genes<-readRDS(str_c("/home/dell/database/Cycle/",org,"_cell_cycle_genes.rds"))
}

s.genes <-cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

homo<-read_tsv("hs.dre.homo.gene.txt")
s.genes<-unique(homo$gene2[homo$gene1 %in% s.genes])
g2m.genes<-unique(homo$gene2[homo$gene1 %in% g2m.genes])
oligos2<-spbmc
s.genes<-s.genes[s.genes %in% row.names(oligos2)]
g2m.genes<-g2m.genes[g2m.genes %in% row.names(oligos2)]

oligos2 <- CellCycleScoring(oligos2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf(str_c("plots/","spbmc",".qc.plot.cellcycle.pdf"),width=14,height=5)
VlnPlot(oligos2, group.by="seurat_clusters",features = c("S.Score"), pt.size =0)
VlnPlot(oligos2, group.by="seurat_clusters",features = c("G2M.Score"), pt.size =0)
VlnPlot(oligos2, group.by="seurat_clusters",features = c("S.Score"), pt.size =0,split.by = "groups")
VlnPlot(oligos2, group.by="seurat_clusters",features = c("G2M.Score"), pt.size =0,split.by = "groups")
dev.off()

############## plot ratio##
library(readr)
table<-read_tsv("clu.3b.pri.ratio.tsv")
table$Var2[which(table$Var2=="Unknown_9")]="pre_V3IN"
table$Var2[which(table$Var2=="Unknown_6")]="V3IN"
table$Var2[which(table$Var2=="Unknown_12")]="MES/FP"
table<-table[which(table$Var1!=13),]

ratio<-apply(table,1,function(x){
  x<-as.data.frame(t(as.data.frame(x)))
  #print(as.numeric(as.character(x$Freq.x)))
  y=as.numeric(as.character(x$Freq.x))/min(as.numeric(as.character(x$Freq.y)),as.numeric(as.character(x$Freq)))
  #print(y)
  return(y)
  })

data_raw<-cbind(table,ratio)
library(pheatmap)
data<-xtabs(ratio~Var2+Var1,data=data_raw)
pdf(paste0("max_heatmap.pdf"),width=10)
g<-pheatmap(data,cluster_rows=F,cluster_col=F,border_color="grey", color = colorRampPalette(c("white","red"))(256),
            fontsize_row=10,fontsize_col=10,fontsize_number=8,cellwidth=17, 
            display_numbers = matrix(ifelse(data>0.3,round(data,2),""), nrow(data)))
print(g)
dev.off()

#################### find marker ###
spbmc<-readRDS(file=str_c("pbmc2",".rds"))
spbmc<-subset(spbmc,seurat_clusters!="13")  
spbmc <- SCTransform(spbmc, method = "glmGamPoi", vars.to.regress = c("percent.ribo"), verbose = FALSE)

if(! dir.exists("marker/") ) { dir.create("marker/")}
if(! dir.exists(str_c("marker/",flabel) )) { dir.create(str_c("marker/",flabel))}
DefaultAssay(spbmc)<-"SCT"
Idents(spbmc)<-spbmc$seurat_clusters
pri_den<-unique(Idents(spbmc))

go_enrich<-function(DEG,org,overlap_cut,flabel,ont,pcut,ncut,nterm=20,width=6,all_bool=0,red_o_bool=1){ 
 # 如果不指定ont ,则BP,CC,MF都运行
if(missing(pcut)) {pcut=0.05}
if(missing(ncut)) {ncut=5}

OrgDb<-"org.Dr.eg.db"

if(!is.null(ncol(DEG))){ colnames(DEG)[1:2]<-c("GeneID","logFC");genel<-DEG$GeneID;} else  {genel<-DEG}
enrigo<-enrichGO(gene =genel, OrgDb = OrgDb, ont = ont, keyType = "ENTREZID",
                     pvalueCutoff = 1, pAdjustMethod = "BH",
                     qvalueCutoff =1, readable =FALSE,minGSSize = 5,maxGSSize = 3000)
    #enrigo_l[[ont]]<-enrigo 
    ens<-enrigo@result
    if(all_bool==1) { write_tsv(ens,str_c("Func_enri/",flabel,".",ont,".all.tsv"))}
    sens<-ens[ens$pvalue<pcut & ens$Count>=ncut,]
    if(!is.null(ncol(DEG))){ sens<-tack_res(DEG,sens); }
    write_tsv(sens,str_c("Func_enri/",flabel,".",ont,".p005.pri.tsv"))
  
    sens<-ens[ens$p.adjust<pcut & ens$Count>=ncut,]
    if(!is.null(ncol(DEG))){ sens<-tack_res(DEG,sens); }

    print(nrow(sens))
    if(nrow(sens)<100) {sens<-ens[ens$pvalue<pcut & ens$Count>=ncut,];if(!is.null(ncol(DEG))){ sens<-tack_res(DEG,sens); }}
    if(nrow(sens)==0) {next;}
    sens<-sens[order(sens$pvalue),]
    if(nrow(sens)>=300){sens<-sens[1:300,]}
    print(nrow(sens))
    #write_tsv(sens,str_c("Func_enri/",flabel,".",ont,".fdr005.pri.tsv"))

    if(red_o_bool==1){
    result<-reduce_overlap(sens,overlap_cut)
    if(!is.null(ncol(DEG))){ 
        plot_GO_statu(result,str_c("Func_enri/",flabel,".",ont),ncut,nterm=nterm,width=width)   
    } else {
        plot_GO(result,str_c("Func_enri/",flabel,".",ont),ncut,nterm=nterm,width=width)   

    }

   }
}



marker_function(spbmc,pri_den,"MAST",geneinfo,res_ctype,flabel)

marker_function<-function(spbmc,pri_den,test_marker,geneinfo,res_ctype,flabel){
  for(i in 1:length(pri_den)){
    marker2 <- FindMarkers(spbmc, ident.1 = pri_den[i], min.pct = 0.25, logfc.threshold = 0.25,only.pos = TRUE,test.use=test_marker) # 用SCT归一化后的数据做marker筛选
    print(dim(marker2))
    marker2$GeneSymbol=row.names(marker2) 
    marker2<-merge(marker2,geneinfo,by="GeneSymbol")
    
    marker2$sig_bool<-marker2$pct.1>=0.6 & marker2$pct.2<=0.4
    marker2<-merge(marker2,res_ctype,by="GeneSymbol",all.x=TRUE) ## res_ctype 来自cellanno.R
    marker2<-arrange(marker2,desc(avg_logFC))
    
    write_tsv(marker2,str_c("marker/",flabel,"/",pri_den[i],".cm.tsv"))
    
    top<-ifelse(nrow(marker2)>=12,12,nrow(marker2))
    if(nrow(marker2)>=12){
      topm1<-arrange(marker2,desc(avg_logFC))
      topm<-subset(topm1,pct.1>=0.6 & pct.2<=0.4)
      if(nrow(topm)>=12){
        topm<-topm[1:12,]} else{topm<-topm1[1:12,]}
    } else{
      topm<-marker2}
    
    
    
    p1<-FeaturePlot(spbmc, features = topm$GeneSymbol,coord.fixed=FALSE,  cols = c("grey", "red"),label=T)
    ggsave(p1,filename=str_c("marker/",flabel,"/",pri_den[i],".withlabel.top12.png"),device="png",width=12,height=7)
    p2<-FeaturePlot(spbmc, features = topm$GeneSymbol,coord.fixed=FALSE,  cols = c("grey", "red"))
    ggsave(p2,filename=str_c("marker/",flabel,"/",pri_den[i],".top12.png"),device="png",width=12,height=7)
  }
  
  
  
  
  selm<-c();
  markers<-c();
  
  for(i in 1:length(pri_den)){
    marker2<-read_tsv(str_c("marker/",flabel,"/",pri_den[i],".cm.tsv"))
    go_enrich(marker2[,c("GeneID")],"dre",0.85,pri_den[i],"BP",nterm=20,width=6) 
    marker2$celltype<-pri_den[i]
    marker2<-arrange(marker2,p_val)
    marker2$pct.re<-marker2$pct.1/ marker2$pct.2
    marker2<-marker2[marker2$pct.re>2 & marker2$avg_logFC>=0.58,]
    selrow<-ifelse(nrow(marker2)>=3,3,nrow(marker2))
    selm<-unique(c(selm,marker2$GeneSymbol[1:3]))
    markers<-rbind(markers,marker2)
  }
  selm<-na.omit(selm)
  p1<-DotPlot(spbmc, features = selm, cols=c("white","red"), dot.scale = 6) + RotatedAxis()
  ggsave(p1,filename=str_c("marker/",flabel,"/all.top3.pdf"),device="pdf",width=14,height=7)
  write_tsv(markers,str_c("marker/",flabel,"/marker.all.tsv"))  ###保存所有marker的信息
  
}



#####
#plot  known markers
selm<-read_tsv("known.marker.txt")

DefaultAssay(pbmc)<-"SCT" ### 是否统一做了归一化处理,代码需要确认下。
selm<-selm[selm$marker%in%row.names(pbmc),]

flabel<-"all"
pdf(file=str_c("plots/",flabel,".known.marker.hp.pdf"),width=7,height=7)
DoHeatmap(pbmc, features = selm$marker,slot = "data",group.colors = NULL) 
dev.off()
### 3) feature plot
cells<-unique(selm$celltype)

for(j in 1:length(cells)){
  
  sm<-selm[selm$celltype==cells[j],]
  p1<-FeaturePlot(pbmc, features = sm$marker,coord.fixed=FALSE,  cols = c("grey", "red"),ncol=3,label=T)
  cells_temp<-gsub("\\/","|",cells[j])
  ggsave(p1,filename=str_c("marker/",flabel,".",cells_temp,".km.scatter.png"),device="png",width=12,height=4*ceiling(length(sm$marker)/3))
 
}

sm<-selm[which(selm$marker %in% c("NPAS3","myt1a","myt1b")),]
p1<-VlnPlot(pbmc, features = sm$marker,group.by="seurat_clusters",split.by="groups",split.plot =FALSE,pt.size=0,combine=TRUE,ncol=1)
ggsave(p1,filename=str_c("marker/",flabel,".km.vln.png"),device="png",width=30,height=10)


############
### cell type ratio of cell type (1,2,6,11)
pbmc<-readRDS(file="pbmc2.rds")
spbmc<-readRDS(file=str_c("PMN_sub.norm.int.rds"))
plot_ratio<-function(num2,label,col) {
     
  if(row.names(num2)[1]!="num_all"){
    num2<-rbind(apply(num2,2,sum),num2)
  }
  ratio<-apply(num2,2,function(x) x[-1]/x[1])
  ratio2<-as.data.frame(ratio);ratio2$celltype<-row.names(ratio2)
  ratio<-melt(ratio)
  colnames(ratio)[c(1,2)]<-c("Var1","Var2")
  ratio$Var2<-as.factor(ratio$Var2)
  #ratio$Var1<-factor(ratio$Var1,levels=str_c("c",c(0,4,7,2,1,3,5,6,8)))
  p2<-ggplot(ratio, aes(x=Var1,y=value))+geom_bar(aes(fill=Var2),position="fill",stat="identity")+
    mytheme+theme(legend.position="right",legend.title= element_blank(),panel.border=element_rect(fill = NA,colour='black',size=0.8),
                  strip.background=element_blank(),axis.text.x=element_text(angle = 45,vjust =1, hjust = 1))+
    ylab("")+xlab("")+scale_fill_manual(values=col)+
    geom_hline(yintercept=0.5,linetype = "dotted")
  ggsave(p2,filename=str_c("plots/",label,".ratio.pdf"),device="pdf",width=7,height=7)
  return(ratio2)
}

meta<-pbmc@meta.data
num_all<-table(meta$groups)

#spbmc<-Integrated
meta2<-spbmc@meta.data
meta2$seurat_cluster<-as.character(meta2$seurat_cluster)
num2<-table(meta2$seurat_cluster,meta2$groups)
row.names(num2)<-str_c("c",row.names(num2))
flabel="PMN_sub"
num2<-rbind(num_all,num2)
ratio2<-plot_ratio(num2,label=str_c(flabel,".norm"),col=c("lightsalmon","thistle3","skyblue"))
write_tsv(ratio2,str_c(flabel,".celltype.ratio.tsv"))

library(RColorBrewer)
DefaultAssay(pbmc)<-"RNA"
flabel<-"PMN_sub"
clus<-c(2,1,6,11)
spbmc1<-subset(pbmc, seurat_clusters %in% clus)
num2<-table(as.character(meta2$seurat_clusters),as.character(as.factor(spbmc1@meta.data$seurat_clusters)))
row.names(num2)<-str_c("c",row.names(num2))
ratio2<-plot_ratio(num2,label=str_c(flabel,".pricell"),col=brewer.pal(6,"Set3"))  
