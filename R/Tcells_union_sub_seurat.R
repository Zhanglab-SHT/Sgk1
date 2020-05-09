#! /usr/bin/env Rscript

library(tidyverse)
library(Matrix)
library(Seurat)
library(openxlsx)
library(RColorBrewer)

library(scales)
library(pheatmap)
library(scran)
library(SingleCellExperiment)
library(scater)


args=commandArgs(T)
id=args[1]
dim=as.numeric(args[2])
resolution=as.numeric(args[3])
load.data1=args[4]
nor.method=args[5]
source('~/functions.R')
source('/public/home/linli/pbs_template/R_fun/ISnorm_functions.R')


load(load.data1)
sub.clust=c(0,1,8,10,15,18)
seurat.obj=GetSubObj(seurat.obj,sub.clust)


#--------------seurat pipeline--------------
if (nor.method == 'log') {
        print('LogNormalize')
        seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
} else if (nor.method == 'Isnorm') {
	    print('Isnorm')
        seurat.mat=as.matrix(seurat.obj[['RNA']]@counts)
        mt_genes=rownames(seurat.mat)[grepl('mt-',rownames(seurat.mat))]
        gene_dis<-calculate.dis(mat=seurat.mat,detection_rate=0.70,ncore=1,exclude=mt_genes)
        spike_candidate<-dbscan.pick(dis=gene_dis)
        candidate_res<-candidate.norm(mat=seurat.mat,spike_candidate=spike_candidate,ncore=1)
        ISnorm_res<-opt.candidate(mat=seurat.mat,candidate_res=candidate_res)
        normalized=ISnorm_res$normalized
        ncol(normalized)
        normalized=normalized[,colSums(is.na(normalized))==0]
        ncol(normalized)
        sel.cells=colnames(normalized)
        seurat.obj=subset(seurat.obj,cells=sel.cells)
        seurat.obj[['RNA']]@data=log(normalized+1)
        str(gene_dis)
        str(spike_candidate)
        png('candidate_instability_Immune.png')
        boxplot(candidate_res$inst)  ##draw a boxplot to see the instability score of cells for each candidate set
        dev.off()
        identical(colnames(seurat.mat),names(ISnorm_res$size_factor))
        png('counts_vs_df.png')
        plot(colSums(as.matrix(seurat.mat)),ISnorm_res$size_factor)
        dev.off()
        save(spike_candidate,candidate_res,ISnorm_res,file='ISnorm_res.RData')
} else {
        print('scran')
        count.mat=seurat.obj[['RNA']]@counts
        sce=SingleCellExperiment(assays=list(counts=count.mat))
        clusters <- quickCluster(sce)
        sce <- computeSumFactors(sce, clusters=clusters)
        print(ncol(sce))
        sf=sizeFactors(sce)
        filename=paste(id,'_sizefactor.RData',sep='')
        save(sf,file=filename)
        names(sf)=colnames(sce)
        sf=sf[sf>0]
        sce=sce[,names(sf)]
        print(ncol(sce))
        sce=normalize(sce,return_log=F)
	filename=paste(id,'_sce.RData',sep='')
	save(sce,file=filename)
        nor.expr=log(as.matrix(sce@assays$data@listData$normcounts)+1)
        colnames(nor.expr)=colnames(sce)
        rownames(nor.expr)=rownames(sce)
        seurat.obj=subset(seurat.obj,cells=colnames(sce))
        seurat.obj[['RNA']]@data=nor.expr
        print(ncol(seurat.obj))
        png('scran.png')
        plot(colSums(as.matrix(count.mat)),sizeFactors(sce))
        dev.off()


}
filename=paste(id,'_step1.RData',sep='')
save(seurat.obj,file=filename)
#-------------------------seurat analysis------------
seurat.obj <- FindVariableFeatures(seurat.obj,selection.method = "vst") # changeable
seurat.obj <- ScaleData(seurat.obj,features=rownames(seurat.obj))
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj),npcs = 60)
seurat.obj <- FindNeighbors(seurat.obj, dims =1:dim)
seurat.obj <- FindClusters(seurat.obj,resolution = resolution,random.seed=111) # change resultion
seurat.obj=RunUMAP(seurat.obj,dims=1:dim)
seurat.obj=RunTSNE(seurat.obj,perplexity=100,dims=1:dim)

seurat.obj_markers <- FindAllMarkers(seurat.obj, only.pos = T)
group_marker <- subset(seurat.obj_markers, p_val_adj < 0.05)
all_list=list()
for (i in unique(group_marker$cluster)) {
         temseurat.obj=subset(group_marker,cluster==i)
         temseurat.obj=list(temseurat.obj)
         all_list=c(all_list,temseurat.obj)
}

filename=paste(id,'_',resolution,'_',dim,'_cluster.xlsx',sep='')
write.xlsx(all_list,filename)

#------kepp data-------------
filename=paste(id,'_',resolution,'_cluster.RData',sep='')
clusters=Idents(seurat.obj)
save(clusters,file=filename)

filename=paste(id,'_',resolution,'_nor.RData',sep='')
nor.expr=as.matrix(seurat.obj[['RNA']]@data)
save(nor.expr,file=filename)

nor.expr1=nor.expr[rowSums(nor.expr>0)>0,]
nor.expr2=nor.expr[rowSums(nor.expr>0)==0,]

filename=paste(id,'_',resolution,'_nor1.RData',sep='')
save(nor.expr1,file=filename)
filename=paste(id,'_',resolution,'_nor2.RData',sep='')
save(nor.expr2,file=filename)



filename=paste(id,'_',resolution,'_umap.RData',sep='')
s.genes=cc.genes$s.genes
s.genes=paste(toupper(substr(s.genes, 1, 1)), tolower(substr(s.genes, 2, nchar(s.genes))), sep="")
g2m.genes=cc.genes$g2m.genes
g2m.genes=paste(toupper(substr(g2m.genes, 1, 1)), tolower(substr(g2m.genes, 2, nchar(g2m.genes))), sep="")
s.genes=intersect(s.genes,rownames(seurat.obj))
g2m.genes=intersect(g2m.genes,rownames(seurat.obj))
seurat.obj<- CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes)

umap.df <- seurat.obj@reductions$umap@cell.embeddings
tsne.df <- as.data.frame(seurat.obj@reductions$tsne@cell.embeddings)
umap.df=as.data.frame(umap.df)
pati=str_split_fixed(colnames(seurat.obj),'_',2)[,1]
umap.df$pati <- pati
umap.df$clust=Idents(seurat.obj)
umap.df$Gene=seurat.obj$nFeature_RNA
umap.df$Count=seurat.obj$nCount_RNA
umap.df$norma=colSums(nor.expr)
umap.df$mt=seurat.obj$percent.mt
umap.df$s.phase=seurat.obj[[]][,"S.Score"]
umap.df$g2m.phase=seurat.obj[[]][,"G2M.Score"]
umap.df$cellcycle=seurat.obj[[]][,"Phase"]
identical(rownames(tsne.df),rownames(umap.df))
umap.df$TSNE_1=tsne.df[[1]]
umap.df$TSNE_2=tsne.df[[2]]
umap.df$type=seurat.obj@meta.data$type
umap.df$expr=seurat.obj@meta.data$expr
umap.df$batch=seurat.obj@meta.data$batch
rownames(umap.df)=colnames(seurat.obj)
save(umap.df,file=filename)

filename=paste(id,'_',resolution,'_count.RData',sep='')
count.mat=as.matrix(seurat.obj[['RNA']]@counts)
save(count.mat,file=filename)

filename=paste(id,'_',resolution,'_all.RData',sep='')
save(seurat.obj,file=filename)

filename=paste(id,'_',resolution,'_umap','.pdf',sep='')
pdf(filename)
print(DimPlot(seurat.obj,reduction='pca',label=T,pt.size=0.2))

print(DimPlot(seurat.obj,reduction='umap',label=T,pt.size=0.2))
print(DimPlot(seurat.obj,reduction='umap',label=T,group.by = 'type',pt.size=0.2))
print(DimPlot(seurat.obj,reduction='umap',label=T,group.by = 'expr',pt.size=0.2))

print(DimPlot(seurat.obj,reduction='tsne',label=T,pt.size=0.2))
print(DimPlot(seurat.obj,reduction='tsne',label=T,group.by = 'type',pt.size=0.2))
print(DimPlot(seurat.obj,reduction='tsne',label=T,group.by = 'expr',pt.size=0.2))
dev.off()

col_name1='Gene'
col_name2='Count'
col_name3='mt'
col_name4='s.phase'
col_name5='g2m.phase'
col_name6='cellcycle'
col_name7='norma'



plot.df=umap.df %>% group_by(pati,clust) %>% summarise(n=n())
filename=paste(id,'_',resolution,'_static.pdf',sep='')
pdf(filename)
ggplot(plot.df,aes(x=pati,y=n,fill=clust)) + geom_bar(stat='identity',position='fill') +
  theme_bw() + theme(text = element_text(size=20))
ggplot(plot.df,aes(x='clust',y=n,fill=pati)) + geom_bar(stat='identity',position='fill') +
  facet_wrap( ~ clust) + theme_bw() + theme(text = element_text(size=20))
GetFeaturePlot(col_name1,umap.df,'Gene')
GetFeaturePlot(col_name2,umap.df,'Count')
GetFeaturePlot(col_name7,umap.df,'norma')
GetFeaturePlot(col_name3,umap.df,'mt')
GetFeaturePlot(col_name4,umap.df,'s.phase')
GetFeaturePlot(col_name5,umap.df,'g2m.phase')
GetFeaturePlot1(col_name6,umap.df,'cellcycle')
dev.off()

