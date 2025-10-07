---
title: "单细胞RNA测序数据分析实战"
updatedAt: "2024-09-20"
type: "LearnResource"
level: "进阶"
track: "单细胞分析"
duration: "2周"
prerequisites: ["Python基础", "统计学基础", "R语言基础"]
---

# 单细胞RNA测序数据分析实战

## 项目概述

本项目将带你完成一个完整的单细胞RNA测序数据分析流程，从原始数据处理到生物学解释。我们将使用真实的研究数据进行实战演练。

## 项目目标

完成本项目后，你将掌握：
- 🔬 单细胞数据处理的标准流程
- 📊 质量控制和数据预处理
- 🎯 细胞聚类和类型识别
- 📈 差异表达分析
- 🧬 生物学解释和可视化

## 技术栈

- **Python**: Scanpy, Seurat (via reticulate)
- **R**: Seurat, ggplot2
- **工具**: Jupyter Notebook, RStudio

## 数据集

我们使用来自10x Genomics的PBMC（外周血单核细胞）数据集：
- 细胞数量: ~2,700个
- 基因数量: ~32,000个
- 来源: 健康人外周血

## 实战步骤

### 1. 环境设置和数据加载

```python
# 安装必要库
!pip install scanpy leidenalg umap-learn
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 设置Scanpy
sc.settings.verbosity = 3
sc.logging.print_header()

# 加载数据
adata = sc.read_10x_mtx('pbmc3k/filtered_gene_bc_matrices/hg19/',
                       var_names='gene_symbols',
                       cache=True)

# 查看数据基本信息
print(f"细胞数量: {adata.n_obs}")
print(f"基因数量: {adata.n_vars}")
```

### 2. 质量控制

```python
# 计算质量控制指标
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                          log1p=False, inplace=True)

# 可视化质量控制指标
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# 过滤低质量细胞
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 进一步过滤
adata = adata[adata.obs.pct_counts_mt < 5, :]
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
```

### 3. 数据预处理

```python
# 归一化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 识别高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3,
                           min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# 保存原始数据
adata.raw = adata

# 过滤基因
adata = adata[:, adata.var.highly_variable]

# 线性回归和缩放
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```

### 4. 降维和聚类

```python
# PCA降维
sc.tl.pca(adata, svd_solver='arpack')

# 确定主成分数量
sc.pl.pca_variance_ratio(adata, log=True)

# 构建邻居图
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP降维
sc.tl.umap(adata)

# Leiden聚类
sc.tl.leiden(adata, resolution=0.5)

# 可视化聚类结果
sc.pl.umap(adata, color=['leiden'])
```

### 5. 细胞类型注释

```python
# 标记基因列表
marker_genes = {
    'T cells': ['CD3D', 'CD3E', 'CD8A'],
    'B cells': ['CD79A', 'MS4A1'],
    'NK cells': ['GNLY', 'NKG7'],
    'Monocytes': ['CD14', 'LYZ'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Megakaryocytes': ['PPBP']
}

# 可视化标记基因表达
sc.pl.umap(adata, color=marker_genes['T cells'])
sc.pl.umap(adata, color=marker_genes['B cells'])

# 自动细胞类型注释
def assign_cell_types(adata, marker_genes):
    cell_types = []
    for cell in adata.obs_names:
        scores = {}
        for cell_type, genes in marker_genes.items():
            score = np.mean(adata[cell, genes].X.toarray().flatten())
            scores[cell_type] = score

        assigned_type = max(scores, key=scores.get)
        cell_types.append(assigned_type)

    adata.obs['cell_type'] = cell_types
    return adata

adata = assign_cell_types(adata, marker_genes)

# 可视化细胞类型
sc.pl.umap(adata, color='cell_type', legend_loc='on data')
```

### 6. 差异表达分析

```python
# 找到每个细胞类型的标记基因
sc.tl.rank_genes_groups(adata, 'cell_type', method='t-test')

# 可视化差异表达基因
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# 获取特定细胞类型的标记基因
def get_markers(adata, cell_type, n_top=10):
    markers = pd.DataFrame({
        'names': adata.uns['rank_genes_groups']['names'][cell_type],
        'scores': adata.uns['rank_genes_groups']['scores'][cell_type]
    })
    return markers.head(n_top)

# T细胞标记基因
t_markers = get_markers(adata, 'T cells')
print("T细胞标记基因:")
print(t_markers)
```

### 7. 功能富集分析

```python
# 使用g:Profiler进行功能富集分析
!pip install gprofiler-official
from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)

def enrichment_analysis(gene_list, organism='hsapiens'):
    result = gp.profile(organism=organism, query=gene_list)
    return result

# T细胞标记基因富集分析
t_cell_genes = t_markers['names'].tolist()
enrichment = enrichment_analysis(t_cell_genes)

# 显示top通路
print("T细胞标记基因富集分析:")
print(enrichment[['native', 'name', 'p_value', 'source']].head(10))
```

### 8. 轨迹分析

```python
# PAGA轨迹分析
sc.tl.paga(adata, groups='cell_type')
sc.pl.paga(adata, plot=True)

# diffusion map
sc.tl.diffmap(adata)
sc.pl.diffmap(adata, color='cell_type')

# RNA速度分析（需要spliced/unspliced数据）
# sc.tl.velocity(adata, mode='stochastic')
# sc.pl.velocity_embedding_stream(adata, basis='umap')
```

## 项目成果

完成本项目后，你将获得：

1. **完整分析流程**：从原始数据到生物学解释的全套代码
2. **可视化图表**：多种单细胞数据可视化结果
3. **分析报告**：包含关键发现和生物学解释
4. **可复现代码**：Jupyter Notebook格式的完整代码

## 进阶挑战

### 挑战1: 批次效应校正
使用多个样本数据，练习批次效应校正技术：
- Harmony
- BBKNN
- Scanpy的批次校正方法

### 挑战2: 细胞通讯分析
使用CellPhoneDB或NicheNet分析细胞间的通讯网络。

### 挑战3: 时序分析
如果有时间序列数据，进行拟时序分析。

## 评估标准

项目完成质量将从以下方面评估：

- **代码质量** (30%): 代码结构、注释、可读性
- **分析深度** (40%): 分析的完整性和生物学解释
- **可视化效果** (20%): 图表的清晰度和信息量
- **创新性** (10%): 是否有额外的分析或改进

## 参考资料

### 必读文献
1. Stuart T, Butler A, et al. "Comprehensive Integration of Single-Cell Data." Cell, 2019.
2. Wolf FA, Angerer P, Theis FJ. "SCANPY: large-scale single-cell gene expression data analysis." Genome Biology, 2018.
3. Butler A, Hoffman P, et al. "Integrating single-cell transcriptomic data across different conditions, technologies, and species." Nature Biotechnology, 2018.

### 实用资源
- Scanpy官方教程: https://scanpy-tutorials.readthedocs.io/
- Seurat官方教程: https://satijalab.org/seurat/
- 单细胞数据分析最佳实践: https://www.scrna-tools.org/

---

*准备好挑战这个令人兴奋的项目了吗？开始你的单细胞分析之旅吧！* 🚀