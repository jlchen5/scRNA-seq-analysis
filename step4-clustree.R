rm(list = ls())
library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
load(file = 'basic.sce.pbmc.Rdata')
DimPlot(pbmc, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
sce=pbmc

# 参考：https://mp.weixin.qq.com/s/WRhMC3Ojy1GWYfLS_4vSeA

#先执行不同resolution 下的分群
library(Seurat)
library(clustree)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1.6,.2))
)
clustree(sce@meta.data, prefix = "RNA_snn_res.")
colnames(sce@meta.data)

p1=DimPlot(sce, reduction = 'umap', group.by = 'RNA_snn_res.0.5',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DimPlot(pbmc, reduction = 'umap',# group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
library(patchwork)
p1+p2

p1=DimPlot(sce, reduction = 'umap', group.by = 'RNA_snn_res.0.5',
           label = TRUE, pt.size = 0.5) + NoLegend()
p2=DimPlot(sce, reduction = 'umap',group.by = 'RNA_snn_res.1.5',
           label = TRUE, pt.size = 0.5) + NoLegend()
library(patchwork)
p1+p2



