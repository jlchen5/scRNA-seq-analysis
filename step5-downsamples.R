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

features= c('IL7R', 'CCR7','CD14', 'LYZ',  'IL7R', 'S100A4',"MS4A1", "CD8A",'FOXP3',
            'FCGR3A', 'MS4A7', 'GNLY', 'NKG7',
            'FCER1A', 'CST3','PPBP')

DoHeatmap(subset(sce ), 
          features = features, 
          size = 3
          )
table(Idents(sce))
DoHeatmap(subset(sce, downsample = 15), 
          features = features, size = 3)
# gsva, 细胞通讯，转录因子，拟时序, inferCNV
# 真实项目，10万+ 

# 每个细胞亚群抽10 
allCells=names(Idents(sce))
allType = levels(Idents(sce))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce)== x ]
  cg=sample(cgCells,10)
  cg
}))

cg_sce = sce[, allCells %in% choose_Cells]
cg_sce
table(Idents(cg_sce))


# 如何看原始表达量 
DoHeatmap(subset(sce ), 
          features = features, 
          size = 3,
          slot = 'data'
)


