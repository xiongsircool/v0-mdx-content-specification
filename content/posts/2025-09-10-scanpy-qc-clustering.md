---
title: "用 Scanpy 进行单细胞质控与聚类"
publishedAt: "2025-09-10"
updatedAt: "2025-09-15"
excerpt: "本文将从 AnnData 对象的创建开始，一步步介绍如何使用 Scanpy 对单细胞数据进行质量控制、标准化、降维和聚类，并探讨流程中的常见问题。"
tags: ["单细胞", "Python", "可视化", "Scanpy"]
authors: ["张三", "李四"]
coverImage: "https://images.unsplash.com/photo-1581093196277-9f6c0cf5a0a5?w=1200&h=630&fit=crop&crop=center&auto=format&q=80"
---

单细胞RNA测序（scRNA-seq）技术的快速发展为我们理解细胞异质性提供了强大的工具。**Scanpy**作为Python生态系统中最受欢迎的单细胞分析包之一，提供了完整的分析流程。本文将详细介绍如何使用Scanpy进行单细胞数据的质量控制和聚类分析。

## 环境准备

首先，我们需要安装必要的Python包：

**使用 pip 安装:**
```bash
pip install scanpy pandas numpy matplotlib seaborn
```

**使用 conda 安装:**
```bash
conda install -c conda-forge scanpy pandas numpy matplotlib seaborn
```

## 数据加载与预处理

### 1. 导入必要的库

\`\`\`python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 设置scanpy参数
sc.settings.verbosity = 3  # 详细输出
sc.settings.set_figure_params(dpi=80, facecolor='white')
\`\`\`

### 2. 加载数据

> **ℹ️ 提示**: Scanpy支持多种数据格式，包括10X Genomics的输出、CSV文件、H5文件等。这里我们以10X数据为例。

\`\`\`python
# 加载10X数据
adata = sc.read_10x_mtx(
    'data/filtered_feature_bc_matrix/',  # 数据路径
    var_names='gene_symbols',            # 使用基因符号作为变量名
    cache=True                           # 缓存数据以加快后续加载
)

# 转置矩阵，使得行为细胞，列为基因
adata = adata.transpose()

# 设置基因名为索引
adata.var_names_unique()
\`\`\`

## 质量控制

质量控制是单细胞分析的关键步骤，主要包括以下几个方面：

### 1. 计算质控指标

\`\`\`python
# 计算线粒体基因比例
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# 计算核糖体基因比例
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
\`\`\`

### 2. 可视化质控指标

![质控指标的小提琴图展示](/placeholder.svg?height=300&width=600)

\`\`\`python
# 绘制质控指标的分布
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
\`\`\`

### 3. 过滤低质量细胞和基因

#### 过滤步骤

1. **过滤基因**: 去除在少于3个细胞中表达的基因
2. **过滤细胞**: 去除表达基因数过少或线粒体基因比例过高的细胞

\`\`\`python
# 过滤基因：至少在3个细胞中表达
sc.pp.filter_genes(adata, min_cells=3)

# 过滤细胞：基于质控指标
sc.pp.filter_cells(adata, min_genes=200)  # 至少表达200个基因
adata = adata[adata.obs.n_genes_by_counts < 5000, :]  # 基因数不超过5000
adata = adata[adata.obs.pct_counts_mt < 20, :]  # 线粒体基因比例小于20%
\`\`\`

> **⚠️ 注意**: 过滤阈值需要根据具体的数据集和实验设计进行调整。过于严格的过滤可能会丢失重要的细胞类型。

## 数据标准化与降维

### 1. 数据标准化

\`\`\`python
# 保存原始数据
adata.raw = adata

# 标准化：每个细胞的总表达量标准化为10,000
sc.pp.normalize_total(adata, target_sum=1e4)

# 对数转换
sc.pp.log1p(adata)
\`\`\`

### 2. 高变基因识别

\`\`\`python
# 识别高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 可视化高变基因
sc.pl.highly_variable_genes(adata)

# 只保留高变基因用于后续分析
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
\`\`\`

### 3. 主成分分析（PCA）

\`\`\`python
# 数据缩放
sc.pp.scale(adata, max_value=10)

# 主成分分析
sc.tl.pca(adata, svd_solver='arpack')

# 可视化PCA结果
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
\`\`\`

## 聚类分析

### 1. 构建邻接图

\`\`\`python
# 计算邻接图
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
\`\`\`

### 2. UMAP降维

\`\`\`python
# UMAP降维
sc.tl.umap(adata)

# 可视化UMAP结果
sc.pl.umap(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'])
\`\`\`

### 3. Leiden聚类

\`\`\`python
# Leiden聚类
sc.tl.leiden(adata, resolution=0.5)

# 可视化聚类结果
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', 
           title='Leiden clustering', frameon=False, save='.pdf')
\`\`\`

![UMAP降维和Leiden聚类结果](/placeholder.svg?height=400&width=600)

## 结果解释与下游分析

### 1. 聚类质量评估

\`\`\`python
# 计算聚类的轮廓系数
from sklearn.metrics import silhouette_score

# 使用PCA空间计算轮廓系数
silhouette_avg = silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'])
print('平均轮廓系数:', round(silhouette_avg, 3))
\`\`\`

### 2. 差异表达分析

\`\`\`python
# 寻找每个聚类的标记基因
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# 可视化标记基因
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
\`\`\`

## 常见问题与解决方案

> **🚨 内存不足**: 对于大型数据集，可以考虑使用`sc.pp.subsample()`进行下采样，或者使用更高效的数据格式如HDF5。

### 1. 聚类分辨率选择

聚类分辨率的选择是一个重要问题：

- **分辨率过低**: 可能会将不同的细胞类型合并
- **分辨率过高**: 可能会过度分割同一细胞类型

建议尝试多个分辨率值（如0.1, 0.3, 0.5, 0.8, 1.0），并结合生物学知识选择最合适的结果。

### 2. 批次效应处理

如果数据存在批次效应，可以使用以下方法：

\`\`\`python
# 使用Harmony进行批次校正
# import scanpy.external as sce
# sce.pp.harmony_integrate(adata, key='batch')

# 或使用scanorama
# import scanorama
# sce.pp.scanorama_integrate(adata, key='batch')
\`\`\`

## 总结

本文介绍了使用Scanpy进行单细胞数据分析的完整流程，包括：

1. **数据加载与预处理**
2. **质量控制与过滤**
3. **数据标准化与高变基因识别**
4. **降维与聚类分析**
5. **结果可视化与解释**

**单细胞分析**是一个迭代的过程，需要根据数据的特点和研究目标不断调整参数。希望本文能为大家的单细胞分析工作提供有用的参考。

---

**参考资料**:
- [Scanpy官方文档](https://scanpy.readthedocs.io/)
- [单细胞分析最佳实践](https://www.sc-best-practices.org/)
\`\`\`
