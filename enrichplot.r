# 数据预处理

rm(list=ls())
options(stringsAsFactors = F)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(qs)
library(Matrix)
getwd()

# 创建目录
getwd()
gse <- "GSE128033"
dir.create(gse)

# 下载 raw 文件夹
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163558/suppl/GSE163558_RAW.tar
# 使用正则表达式替换：\w{3}$ 匹配每个单词最后的三个字符替换为空字符串，即去掉它们
s <- gsub("(\\w{3}$)", "", gse, perl = TRUE)
s
url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",s,"nnn/",gse,"/suppl/",gse,"_RAW.tar")
url
file <- paste0(gse,"_RAW.tar")
file
# downloader::download(url, destfile = file)

###### step1: 导入数据 ######  
# GSE128033为整理后的每个文件里都是标准的三个文件
samples <- list.dirs("GSE128033/", recursive = F, full.names = F)
samples
scRNAlist <- lapply(samples, function(pro){
  #pro <- samples[1]
  print(pro)
  folder <- file.path("GSE128033/", pro)
  folder
  counts <- Read10X(folder, gene.column = 2)
  sce <- CreateSeuratObject(counts, project=pro, min.features=3)
  return(sce)
})
names(scRNAlist) <-  samples
scRNAlist

# merge
sce.all <- merge(scRNAlist[[1]], y=scRNAlist[-1], add.cell.ids=samples)
sce.all <- JoinLayers(sce.all) # seurat v5
sce.all


# 查看特征
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
sce.all$Sample <- sce.all$orig.ident
sce.all$orig.ident <- str_split(sce.all$orig.ident, pattern = "_", n=2,simplify = T)[,2]
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 

library(qs)
qsave(sce.all, file="GSE128033/sce.all.qs")

## 对每个亚群做差异分析得到差异基因
# 加载包
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

table(Idents(sce.all.int))
sce.all.int <- subset(sce.all.int, ident="Double cells", invert=T)

# 差异分析
sce.markers <- FindAllMarkers(sce.all.int, only.pos = TRUE, min.pct = 0.2, return.thresh = 0.01)
head(sce.markers)

# 查看top10
top10 <- sce.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

# 基因转成ENTREZID
Symbol <- mapIds(get("org.Hs.eg.db"), keys = sce.markers$gene, keytype = "SYMBOL", column="ENTREZID")
head(Symbol)
ids <- bitr(sce.markers$gene,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
head(ids)

# 合并ENTREZID到obj.markers中
data <- merge(sce.markers, ids, by.x="gene", by.y="SYMBOL")
head(data)

gcSample <- split(data$ENTREZID, data$cluster)
gcSample

## 富集分析
## 富集分析 
# KEGG
# xx_kegg <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
# GO
xx_go <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=1, qvalueCutoff=1)
res <- xx_go@compareClusterResult
head(res)

## 将富集结果中的 ENTREZID 重新转为 SYMBOL
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}
head(res)

## 通路筛选Top5
enrich <- res %>% 
  group_by(Cluster) %>% 
  top_n(n = 5, wt = -pvalue) 

dt <- enrich
dt <- dt[order(dt$pvalue,dt$Cluster, decreasing = F), ]
dt$Description <- factor(dt$Description, levels = unique(dt$Description))
colnames(dt)

## 绘图
library(ggforce)
head(dt)
table(dt$Cluster)

colors <- c(
  "Macrophages" = "#B0C4DE",
  "Mast cells" = "#FF69B4", 
  "Endothelial cells" = "#FFB6C1", 
  "ILC1/NK cells" = "#FFFFE0", 
  "B cells" = "#ADD8E6",  
  "T cells" = "#90EE90",   
  "Fibroblasts" = "#FFA07A",  
  "Smooth muscle cells" = "#D3D3D3",
  "Epithelial cells" = "#D8BFD8",   
  "Cycling macrophages" = "#E6E6FA", 
  "Lymphatic endothelial cells" = "#FFD700" 
)

p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
               xend = -log10(pvalue), yend = Description,
               alpha = after_stat(index),
               color = Cluster,
               size = after_stat(index)), 
            n = 500, show.legend = T) 
p

p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        #axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black", size = 0.6), # 加粗x轴和y轴的线条
        axis.text = element_text(face = "bold"), # 加粗x轴和y轴的标签
        axis.title = element_text( size = 13)    # 加粗x轴和y轴的标题
        ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors)
p1
ggsave(filename = "Enrich_cometplot.pdf", width = 15, height = 8, plot = p1)


# 加分面
p2 <- p1 + 
  facet_wrap(~Cluster,scales = "free",ncol = 2) 
p2

ggsave(filename = "Enrich_cometplot_facet.pdf", width = 16, height = 10, plot = p2)

