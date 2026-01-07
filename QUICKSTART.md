# 快速开始指南 (Quick Start Guide)

## 🚀 5分钟上手

本项目已完成，可以直接运行分析！

---

## 📋 前置要求

- **操作系统**：Linux/macOS/Windows (WSL推荐)
- **软件**：Conda (Miniconda或Anaconda)
- **Python版本**：3.10+
- **R版本**：4.3+ (可选，用于高级可视化)

---

## 📦 第一步：安装环境

### 方法1：使用Conda（推荐）

```bash
# 进入项目目录
cd RNA_modification_analysis

# 创建环境并安装依赖
conda env create -f environment.yml

# 激活环境
conda activate rna_modif_analysis
```

### 方法2：手动安装（如果conda失败）

```bash
# Python依赖
pip install pandas numpy matplotlib seaborn biopython logomaker pyfaidx pysam jupyter notebook requests matplotlib-venn

# R依赖（可选）
Rscript -e "install.packages(c('ggplot2', 'dplyr', 'readr', 'gridExtra', 'RColorBrewer', 'scales'))"
Rscript -e "if (!require('BiocManager')) install.packages('BiocManager'); BiocManager::install(c('ChIPseeker', 'VennDiagram'))"
```

---

## 🎯 第二步：运行分析

### 选项A：使用Jupyter Notebook（最简单）

```bash
# 启动Jupyter
jupyter notebook complete_analysis.ipynb

# 然后按照Notebook中的步骤，逐个运行代码单元格
```

**优点**：
- 交互式运行，可以看到每一步的输出
- 图表直接在浏览器中显示
- 适合学习和调试

### 选项B：运行Python脚本（推荐）

```bash
# 1. 数据获取（生成模拟数据）
python scripts/1_data_fetching.py

# 2. Metagene Profile分析
python scripts/3_metagene_profile.py

# 3. 韦恩图分析
python scripts/4_venn_analysis.py

# 如果有真实基因组，还可以运行：
# python scripts/2_motif_analysis.py
```

### 选项C：使用R可视化（高级）

```bash
# 确保已完成Python脚本分析（生成CSV数据）
Rscript scripts/visualizations.R
```

---

## 📊 第三步：查看结果

分析完成后，查看以下目录：

### `results/` - 分析结果
```
results/
├── psi_sites_annotated.csv        # Ψ位点完整信息
├── m6a_sites_annotated.csv        # m6A位点完整信息
├── psi_metagene_profile.csv       # Ψ的metagene profile数据
├── m6a_metagene_profile.csv       # m6A的metagene profile数据
├── colocalized_sites.csv          # 共定位位点信息
├── colocalization_report.txt      # 共定位分析报告
└── ANALYSIS_REPORT.txt            # 完整分析报告
```

### `figures/` - 生成的图表
```
figures/
├── metagene_profile_comparison.png      # Metagene profile对比图
├── feature_distribution_comparison.png  # 特征分布对比图
├── feature_distribution_piecharts.png   # 饼图
├── venn_diagram.png                     # 韦恩图
├── colocalization_features.png          # 共定位特征分析
└── R_*.pdf / R_*.png                    # R生成的高质量图表
```

---

## 🎨 主要图表说明

### 1. Metagene Profile对比图
- **文件**：`metagene_profile_comparison.png`
- **内容**：红线(m6A)和蓝线(Ψ)在基因上的分布曲线
- **关键发现**：m6A在3'UTR有明显峰值

### 2. 特征分布对比图
- **文件**：`feature_distribution_comparison.png`
- **内容**：柱状图展示两种修饰在5'UTR/CDS/3'UTR的分布
- **关键发现**：m6A强烈偏好3'UTR(60%+)，Ψ分布较均匀

### 3. 韦恩图
- **文件**：`venn_diagram.png`
- **内容**：展示Ψ和m6A的共定位情况
- **关键发现**：共定位率约5%

---

## ⚙️ 高级选项

### 使用真实数据

如果想使用真实的RNA修饰数据（而非模拟数据）：

#### 下载数据

```bash
# 进入data目录
cd data

# 下载Pseudouridine数据（示例）
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102476/suppl/GSE102476_CeU-seq_HEK293T.bed.gz
gunzip GSE102476_CeU-seq_HEK293T.bed.gz

# 下载m6A数据（需要从REPIC手动下载）
# 访问: https://repic.idrb.cas.cz/download

# 下载基因组注释和FASTA
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.annotation.gtf.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz
gunzip GRCh38.p14.genome.annotation.gtf.gz
gunzip GRCh38.p14.genome.fa.gz

cd ..
```

#### 修改脚本路径

编辑 `scripts/1_data_fetching.py`，注释掉模拟数据生成部分，改为读取真实BED文件。

---

## 🐛 常见问题

### Q1: Conda环境创建失败
**A**: 尝试使用Mamba（更快的Conda替代品）：
```bash
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

### Q2: ImportError: No module named 'xxx'
**A**: 激活环境后手动安装：
```bash
pip install xxx
```

### Q3: R包安装失败
**A**: 使用BiocManager：
```r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("包名")
```

### Q4: 中文显示乱码
**A**: 修改matplotlib字体设置，或在代码开头添加：
```python
plt.rcParams['font.sans-serif'] = ['SimHei']  # Windows
# 或
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']  # macOS
```

### Q5: 内存不足
**A**: 减少分析的位点数量，在脚本中修改：
```python
# 例如，只分析1000个位点
psi_subset = psi_df.sample(1000, random_state=42)
```

---

## 📈 性能优化建议

- **使用真实数据时**：位点数可能达到数万，建议使用采样分析
- **Motif分析**：提取序列耗时较长，可以先对少量位点测试
- **并行计算**：可以使用Python的multiprocessing加速
- **内存管理**：分批处理大数据集，避免一次性加载

---

## 📧 获取帮助

如遇到问题：
1. 检查 `scripts/` 目录下的代码注释
2. 查看 `PROJECT_ABSTRACT.md` 了解项目细节
3. 查看每个Python脚本顶部的docstring

---

## ✅ 检查清单

运行分析前，确保：
- [ ] 已安装Conda并激活环境
- [ ] 已进入项目目录
- [ ] `data/`, `results/`, `figures/` 目录已创建
- [ ] Python版本 >= 3.10

运行分析后，检查：
- [ ] `results/` 目录中有CSV文件
- [ ] `figures/` 目录中有PNG/PDF文件
- [ ] 统计报告输出正常
- [ ] 图表清晰可读

---

## 🎓 下一步建议

1. **熟悉代码**：阅读 `scripts/` 下的Python脚本，理解分析流程
2. **修改参数**：调整窗口大小、bin数量等参数，观察结果变化
3. **真实数据**：下载真实的CeU-seq和m6A-seq数据，重新运行
4. **扩展功能**：添加新的分析模块（如进化保守性、结构预测等）

---

**祝你分析顺利！如有问题，欢迎提问。** 🎉

---

*更新日期：2026-01-07*
