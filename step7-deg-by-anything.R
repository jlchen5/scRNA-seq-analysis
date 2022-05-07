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

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE,
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers) 
save(sce.markers,file = 'sce.markers.all_10_celltype.Rdata')

load(file = 'sce.markers.all_10_celltype.Rdata')

head(sce.markers)
table(sce.markers$cluster)
kp=grepl('Mono',sce.markers$cluster)
table(kp)
cg_sce.markers = sce.markers [ kp ,]
DoHeatmap(sce,unique(cg_sce.markers$gene),size=3) 

levels(Idents(sce))
markers_df <- FindMarkers(object = sce, 
                          ident.1 = 'FCGR3A+ Mono',
                          ident.2 = 'CD14+ Mono',
                          #logfc.threshold = 0,
                          min.pct = 0.25)
DoHeatmap(sce,unique(rownames(markers_df)),size=3) 

intersect( rownames(markers_df) ,
          cg_sce.markers$gene)


# drop-out
highCells= colnames(subset(x = sce, subset = FCGR3A > 1,
                           slot = 'counts')) 
highORlow=ifelse(colnames(sce) %in% highCells,'high','low')
table(highORlow)
sce@meta.data$highORlow=highORlow
markers <- FindMarkers(sce, ident.1 = "high", 
                       group.by = 'highORlow' )
head(x = markers)


intersect( rownames(markers_df) ,
           rownames(markers) )


