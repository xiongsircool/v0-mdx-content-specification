---
title: "单细胞RNA测序数据分析入门"
publishedAt: "2025-10-05"
updatedAt: "2025-10-05"
excerpt: "介绍单细胞RNA测序数据的预处理、质量控制、聚类分析和差异表达基因识别等核心流程"
tags: ["单细胞", "RNA-seq", "生物信息学", "Python", "数据分析"]
authors: ["张三", "李四"]
coverImage: "https://images.unsplash.com/photo-1532187863489-5174146d3515?w=1200&h=630&fit=crop&crop=center&auto=format&q=80"
---

# 单细胞RNA测序数据分析入门

单细胞RNA测序（scRNA-seq）技术已经成为现代生物学研究的重要工具，能够揭示细胞异质性和发育过程中的基因表达变化。本文将介绍scRNA-seq数据分析的完整流程。

## 1. 数据质量控制

### 1.1 数据加载

```python
import scanpy as sc
import numpy as np
import pandas as pd

# 读取数据
adata = sc.read_10x_mtx("path/to/filtered_feature_bc_matrix")
adata.var_names_make_unique()
```

### 1.2 质量控制指标

```python
# 计算质量控制指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```

## 2. 数据预处理

### 2.1 过滤低质量细胞

```python
# 过滤条件
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 基于质量控制指标过滤
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
```

### 2.2 数据标准化

```python
# 总计数标准化
sc.pp.normalize_total(adata, target_sum=1e4)

# 对数转换
sc.pp.log1p(adata)

# 识别高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
```

## 3. 降维和聚类

### 3.1 主成分分析

```python
# 缩放数据
sc.pp.scale(adata, max_value=10)

# PCA降维
sc.tl.pca(adata, svd_solver='arpack')
```

### 3.2 聚类分析

```python
# 构建邻域图
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP降维
sc.tl.umap(adata)

# 聚类
sc.tl.leiden(adata, resolution=0.5)
```

## 4. 差异表达分析

### 4.1 寻找标记基因

```python
# 寻找每个聚类的标记基因
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# 查看结果
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
```

## 5. 可视化

### 5.1 UMAP可视化

```python
import matplotlib.pyplot as plt

# UMAP聚类图
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', title='聚类结果')

# 基因表达可视化
sc.pl.umap(adata, color=['CD3D', 'CD8A'], title=['T细胞标记基因'])
```

## 6. 细胞类型注释

```python
# 基于标记基因注释细胞类型
marker_genes = {
    'T细胞': ['CD3D', 'CD8A', 'CD4'],
    'B细胞': ['CD79A', 'MS4A1'],
    'NK细胞': ['GNLY', 'NKG7'],
    '单核细胞': ['CD14', 'LYZ'],
    '树突细胞': ['FCER1A', 'CST3']
}

# 自动注释
for cell_type, genes in marker_genes.items():
    adata.obs[cell_type] = adata[:, genes].X.mean(axis=1)
```

## 7. 功能富集分析

```python
import gseapy as gp

# 对差异表达基因进行GO富集分析
de_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)['0']
enrichment_results = gp.enrichr(gene_list=de_genes.tolist(),
                                gene_sets=['GO_Biological_Process_2023'],
                                organism='human')
```

## 总结

本文介绍了scRNA-seq数据分析的核心流程，包括数据质量控制、预处理、降维聚类、差异表达分析和细胞类型注释。这些分析步骤是单细胞研究的基础，后续还可以进行轨迹分析、细胞通讯分析等更深入的研究。

## 参考资料

1. Scanpy官方文档: https://scanpy.readthedocs.io/
2. Seurat教程: https://satijalab.org/seurat/
3. 单细胞分析最佳实践: https://www.sc-best-practices.org/