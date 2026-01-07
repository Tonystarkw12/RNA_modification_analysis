# RNA修饰比较分析项目

## 项目概述
比较Pseudouridine (Ψ) 和 m6A 两种RNA修饰在HEK293T细胞中的分布特征

## 目录结构
```
RNA_modification_analysis/
├── data/              # 存放原始数据和GTF文件
├── results/           # 存放分析结果（CSV、中间文件）
├── figures/           # 存放生成的图表
├── scripts/           # Python和R脚本
└── README.md          # 项目说明
```

## 安装环境

```bash
# 创建conda环境
conda env create -f environment.yml
conda activate rna_modif_analysis
```

## 使用说明

1. **数据获取**：运行 `1_data_fetching.ipynb`
2. **完整分析**：运行 `complete_analysis.ipynb`
3. **R可视化**：运行 `visualizations.R`

## 数据来源

- **Pseudouridine (Ψ)**: CeU-seq (GSE102476) 或 RMBase v2.0
- **m6A**: REPIC 或 m6A-Atlas (HEK293T)
- **基因组注释**: GENCODE hg38

## 作者
生成于 2026-01-07
