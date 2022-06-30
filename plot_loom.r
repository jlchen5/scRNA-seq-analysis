library(hdf5r)
library(SeuratDisk)
library(scater)
library(tidyverse)
library(Seurat)

library(ggplot2)
library(cowplot)

# load loom
allcell <- Connect(filename = "./Thienpont_Tumors_52k_v4_R_fixed.loom", mode = "r")
allcell

allcell[["col_attrs/"]]

allcell[["row_attrs/"]]

metadata <- data.frame(
  CellID = allcell[["col_attrs/CellID"]][],
  ClusterName = allcell[["col_attrs/ClusterName"]][],
  CellFromTumor = allcell[["col_attrs/CellFromTumor"]][],
  PatientNumber = allcell[["col_attrs/PatientNumber"]][],
  ClusterID = allcell[["col_attrs/ClusterID"]][],
  nUMI = allcell[["col_attrs/nUMI"]][],
  Embedding = allcell[["col_attrs/Embedding"]][]
)


colnames(metadata)


# plot
# CellFromTumor
p1 <- ggplot(metadata,aes(x = Embedding._X,y = Embedding._Y,
                    color = CellFromTumor)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.1),
        legend.background = element_blank(),
        aspect.ratio = 1) +
  xlab('tSNE 1') + ylab('tSNE 2') +
  scale_color_manual(values = c('0' = '#66CC99','1' = '#336699'),
                    name = '',
                    label = c('0' = 'Non-malignant','1' = 'Tumour')) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# PatientNumber
p2 <- ggplot(metadata,aes(x = Embedding._X,y = Embedding._Y,
                  color = PatientNumber)) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.2),
        legend.background = element_blank(),
        aspect.ratio = 1) +
  xlab('tSNE 1') + ylab('tSNE 2') +
  scale_color_manual(values = c('1' = '#339933','2' = '#FF6600','3' = '#CCFF99',
                                '4' = '#FFCC33','5' = '#336699'),
                     name = '') +
  guides(color = guide_legend(override.aes = list(size = 3)))

# ClusterName
p3 <- ggplot(metadata,aes(x = Embedding._X,y = Embedding._Y,
                    color = ClusterName)) +
  geom_point(size = 0.2,show.legend = F) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.1),
        legend.background = element_blank(),
        aspect.ratio = 1) +
  xlab('tSNE 1') + ylab('tSNE 2')

# nUMI
p4 <- ggplot(metadata,aes(x = Embedding._X,y = Embedding._Y,
                    color = log10(nUMI))) +
  geom_point(size = 0.2) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.1),
        legend.background = element_blank(),
        aspect.ratio = 1) +
  xlab('tSNE 1') + ylab('tSNE 2') +
  scale_color_gradient2(low = 'white',mid = 'grey80',high = 'blue',midpoint = 2.5,
                      name = '',
                      limits = c(2,5),
                      guide = guide_colorbar(direction = 'horizontal'))

# combine
plot_grid(plotlist = list(p1,p2,p3,p4),ncol = 2,align = 'hv')




