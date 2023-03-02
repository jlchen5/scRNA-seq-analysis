# 导入需要的库
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# 读入单细胞数据
data = sc.read("data.h5ad")

# 数据预处理
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)
sc.pp.normalize_total(data)
sc.pp.log1p(data)

# 执行PCA降维
sc.pp.highly_variable_genes(data, n_top_genes=2000)
data = data[:, data.var['highly_variable']]
sc.pp.pca(data)

# 聚类
sc.pp.neighbors(data)
sc.tl.umap(data)
sc.tl.leiden(data)

# 可视化
sc.pl.umap(data, color=['leiden'], cmap='viridis')
