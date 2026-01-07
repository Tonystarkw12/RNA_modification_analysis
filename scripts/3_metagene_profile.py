#!/usr/bin/env python3
"""
Metagene Profile分析模块 - RNA修饰比较分析
功能：
1. 绘制修饰在mRNA上的分布曲线（5'UTR-CDS-3'UTR）
2. 展示Ψ和m6A的空间分布差异
3. 绘制基因组特征分布图（饼图/柱状图）

关键生物学发现：
- m6A强烈富集在终止密码子附近和3'UTR
- Ψ分布更均匀，CDS和UTR中都有

作者：AI Assistant
日期：2026-01-07
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from scipy import stats
from scipy.signal import savgol_filter
import warnings
warnings.filterwarnings('ignore')

# 配置
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置绘图风格
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 10

# 路径设置
PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures"


# ==============================================================================
# Metagene Profile分析
# ==============================================================================

class MetageneProfiler:
    """
    Metagene分析器：绘制修饰沿基因的分布曲线
    """

    def __init__(self, sites_df: pd.DataFrame, mod_name: str):
        """
        初始化

        参数：
        - sites_df: 位点信息DataFrame，必须包含'relative_pos'列（0-1之间的相对位置）
        - mod_name: 修饰名称（如'Pseudouridine (Ψ)'或'm6A'）
        """
        self.sites_df = sites_df
        self.mod_name = mod_name
        self.bin_counts = None
        self.bins = None

    def compute_metagene_profile(self, n_bins: int = 100,
                                 smooth: bool = True,
                                 window: int = 9) -> pd.DataFrame:
        """
        计算metagene profile（基因被分成100个bin）

        参数：
        - n_bins: bin数量（通常100）
        - smooth: 是否平滑曲线
        - window: 平滑窗口大小（必须为奇数）

        返回：
        - DataFrame包含bin位置和对应的修饰密度
        """
        # 将基因分成100个bin (0% -> 100%)
        bins = np.linspace(0, 1, n_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # 统计每个bin中的位点数
        hist, _ = np.histogram(self.sites_df['relative_pos'], bins=bins)

        # 归一化（总位点数=1）
        density = hist / hist.sum()

        # 应用平滑（Savitzky-Golay滤波器）
        if smooth and len(density) > window:
            try:
                density_smooth = savgol_filter(density, window, polyorder=2)
                # 确保平滑后非负
                density_smooth = np.maximum(density_smooth, 0)
                # 重新归一化
                density_smooth = density_smooth / density_smooth.sum()
            except:
                density_smooth = density
        else:
            density_smooth = density

        # 保存结果
        self.bin_counts = density_smooth
        self.bins = bin_centers

        return pd.DataFrame({
            'bin_center': bin_centers,
            'density': density_smooth
        })

    def get_feature_distribution(self) -> pd.Series:
        """
        获取修饰在不同基因组特征中的分布（CDS, 5'UTR, 3'UTR）

        返回：
        - Series，index为特征名，values为计数
        """
        if 'feature' not in self.sites_df.columns:
            logger.warning("数据中缺少'feature'列，无法计算特征分布")
            return pd.Series()

        feature_counts = self.sites_df['feature'].value_counts()
        return feature_counts


# ==============================================================================
# 可视化函数
# ==============================================================================

def plot_metagene_profile_comparison(psi_profile: pd.DataFrame,
                                     m6a_profile: pd.DataFrame,
                                     output_file: str = None):
    """
    绘制Ψ和m6A的metagene profile对比图

    参数：
    - psi_profile: Ψ的profile DataFrame
    - m6a_profile: m6A的profile DataFrame
    - output_file: 输出文件路径
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    # 绘制曲线
    ax.plot(psi_profile['bin_center'], psi_profile['density'],
           color='#3498db', linewidth=2.5, label='Pseudouridine (Ψ)', marker='o', markersize=3)
    ax.plot(m6a_profile['bin_center'], m6a_profile['density'],
           color='#e74c3c', linewidth=2.5, label='m6A', marker='s', markersize=3)

    # 添加区域标注
    ax.axvspan(0, 0.2, alpha=0.1, color='gray', label="5'UTR")
    ax.axvspan(0.2, 0.8, alpha=0.1, color='green', label='CDS')
    ax.axvspan(0.8, 1.0, alpha=0.1, color='orange', label="3'UTR")

    # 添加终止密码子位置（CDS末端，约80%处）
    ax.axvline(x=0.8, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.text(0.805, ax.get_ylim()[1]*0.95, 'Stop Codon',
           fontsize=10, style='italic', color='black')

    # 添加起始密码子位置（约20%处）
    ax.axvline(x=0.2, color='black', linestyle='--', linewidth=1.5, alpha=0.6)
    ax.text(0.205, ax.get_ylim()[1]*0.95, 'Start Codon',
           fontsize=10, style='italic', color='black')

    # 图形设置
    ax.set_xlabel('Relative position on transcript (5\' → 3\')', fontsize=12, fontweight='bold')
    ax.set_ylabel('Normalized density', fontsize=12, fontweight='bold')
    ax.set_title('Metagene Profile: Pseudouridine (Ψ) vs m6A',
                fontsize=14, fontweight='bold', pad=20)

    ax.legend(loc='upper right', frameon=True, shadow=True)
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"保存Metagene Profile图: {output_file}")

    plt.show()


def plot_feature_distribution_comparison(psi_features: pd.Series,
                                        m6a_features: pd.Series,
                                        output_file: str = None):
    """
    绘制修饰在不同基因组特征中分布的对比柱状图

    参数：
    - psi_features: Ψ的特征分布Series
    - m6a_features: m6A的特征分布Series
    - output_file: 输出文件路径
    """
    # 确保两者包含相同的特征
    all_features = ['5UTR', 'CDS', '3UTR']

    psi_counts = [psi_features.get(f, 0) for f in all_features]
    m6a_counts = [m6a_features.get(f, 0) for f in all_features]

    # 转换为百分比
    psi_percent = np.array(psi_counts) / sum(psi_counts) * 100
    m6a_percent = np.array(m6a_counts) / sum(m6a_counts) * 100

    # 绘图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # 子图1: 柱状图对比
    x = np.arange(len(all_features))
    width = 0.35

    bars1 = ax1.bar(x - width/2, psi_percent, width,
                   label='Pseudouridine (Ψ)', color='#3498db', edgecolor='black')
    bars2 = ax1.bar(x + width/2, m6a_percent, width,
                   label='m6A', color='#e74c3c', edgecolor='black')

    # 添加数值标签
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}%',
                    ha='center', va='bottom', fontsize=9)

    ax1.set_xlabel('Genomic Feature', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution across Genomic Features',
                 fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(["5'UTR", 'CDS', "3'UTR"])
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # 子图2: 堆叠柱状图（绝对计数）
    feature_labels = ["5'UTR", 'CDS', "3'UTR"]
    x_pos = np.arange(len(feature_labels))

    ax2.bar(feature_labels, psi_counts, label='Pseudouridine (Ψ)',
           color='#3498db', alpha=0.8, edgecolor='black')
    ax2.bar(feature_labels, m6a_counts, bottom=psi_counts,
           label='m6A', color='#e74c3c', alpha=0.8, edgecolor='black')

    # 添加数值标签
    for i, (psi, m6a) in enumerate(zip(psi_counts, m6a_counts)):
        ax2.text(i, psi/2, f'{psi:,}', ha='center', va='center',
                color='white', fontweight='bold', fontsize=10)
        ax2.text(i, psi + m6a/2, f'{m6a:,}', ha='center', va='center',
                color='white', fontweight='bold', fontsize=10)

    ax2.set_xlabel('Genomic Feature', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Sites', fontsize=12, fontweight='bold')
    ax2.set_title('Stacked Distribution (Absolute Counts)',
                 fontsize=13, fontweight='bold')
    ax2.legend()

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"保存特征分布对比图: {output_file}")

    plt.show()


def plot_pie_charts(psi_features: pd.Series,
                   m6a_features: pd.Series,
                   output_file: str = None):
    """
    绘制饼图展示两种修饰的基因组特征分布

    参数：
    - psi_features: Ψ的特征分布Series
    - m6a_features: m6A的特征分布Series
    - output_file: 输出文件路径
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # 准备数据
    def prepare_data(features):
        data = [features.get('5UTR', 0),
                features.get('CDS', 0),
                features.get('3UTR', 0)]
        labels = ["5'UTR", 'CDS', "3'UTR"]
        colors = ['#95a5a6', '#2ecc71', '#f39c12']
        return data, labels, colors

    psi_data, psi_labels, psi_colors = prepare_data(psi_features)
    m6a_data, m6a_labels, m6a_colors = prepare_data(m6a_features)

    # 绘制饼图
    wedges1, texts1, autotexts1 = ax1.pie(psi_data, labels=psi_labels, colors=psi_colors,
                                          autopct='%1.1f%%', startangle=90,
                                          explode=(0.05, 0.05, 0.05),
                                          shadow=True, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax1.set_title('Pseudouridine (Ψ)\nDistribution', fontsize=13, fontweight='bold', pad=15)

    wedges2, texts2, autotexts2 = ax2.pie(m6a_data, labels=m6a_labels, colors=m6a_colors,
                                          autopct='%1.1f%%', startangle=90,
                                          explode=(0.05, 0.05, 0.05),
                                          shadow=True, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax2.set_title('m6A\nDistribution', fontsize=13, fontweight='bold', pad=15)

    # 设置百分比文字样式
    for autotext in autotexts1 + autotexts2:
        autotext.set_color('white')
        autotext.set_fontsize(12)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"保存饼图: {output_file}")

    plt.show()


# ==============================================================================
# 统计分析
# ==============================================================================

def compare_distributions(psi_df: pd.DataFrame, m6a_df: pd.DataFrame):
    """
    统计分析两种修饰的分布差异

    参数：
    - psi_df: Ψ位点DataFrame
    - m6a_df: m6A位点DataFrame
    """
    logger.info("\n" + "=" * 80)
    logger.info("统计分析：Ψ vs m6A 分布差异")
    logger.info("=" * 80 + "\n")

    # 1. 比较在不同基因组特征中的分布（卡方检验）
    logger.info("[1] 基因组特征分布差异分析（卡方检验）")

    psi_features = psi_df['feature'].value_counts()
    m6a_features = m6a_df['feature'].value_counts()

    # 构建列联表
    contingency_table = pd.DataFrame({
        'Ψ': [psi_features.get('5UTR', 0), psi_features.get('CDS', 0), psi_features.get('3UTR', 0)],
        'm6A': [m6a_features.get('5UTR', 0), m6a_features.get('CDS', 0), m6a_features.get('3UTR', 0)]
    }, index=['5UTR', 'CDS', '3UTR'])

    logger.info("\n列联表:")
    logger.info(contingency_table)

    # 卡方检验
    chi2, p_value, dof, expected = stats.chi2_contingency(contingency_table)
    logger.info(f"\n卡方检验结果:")
    logger.info(f"  Chi2 = {chi2:.2f}")
    logger.info(f"  P-value = {p_value:.2e}")
    if p_value < 0.001:
        logger.info(f"  结论: Ψ和m6A的分布存在极显著差异 *** (P < 0.001)")
    elif p_value < 0.01:
        logger.info(f"  结论: Ψ和m6A的分布存在显著差异 ** (P < 0.01)")
    elif p_value < 0.05:
        logger.info(f"  结论: Ψ和m6A的分布存在差异 * (P < 0.05)")
    else:
        logger.info(f"  结论: 差异不显著 (P >= 0.05)")

    # 2. 比较在基因上的相对位置分布（Kolmogorov-Smirnov检验）
    logger.info("\n[2] 相对位置分布差异分析（KS检验）")

    psi_positions = psi_df['relative_pos'].values
    m6a_positions = m6a_df['relative_pos'].values

    ks_stat, ks_p = stats.ks_2samp(psi_positions, m6a_positions)
    logger.info(f"\nKS检验结果:")
    logger.info(f"  KS statistic = {ks_stat:.4f}")
    logger.info(f"  P-value = {ks_p:.2e}")

    if ks_p < 0.001:
        logger.info(f"  结论: 两种修饰在基因上的位置分布存在极显著差异 ***")
    elif ks_p < 0.01:
        logger.info(f"  结论: 两种修饰在基因上的位置分布存在显著差异 **")
    else:
        logger.info(f"  结论: 差异不显著")

    # 3. 描述性统计
    logger.info("\n[3] 描述性统计")
    logger.info(f"\nPseudouridine (Ψ):")
    logger.info(f"  总位点数: {len(psi_df)}")
    logger.info(f"  平均相对位置: {psi_positions.mean():.3f} ± {psi_positions.std():.3f}")
    logger.info(f"  中位数相对位置: {np.median(psi_positions):.3f}")

    logger.info(f"\nm6A:")
    logger.info(f"  总位点数: {len(m6a_df)}")
    logger.info(f"  平均相对位置: {m6a_positions.mean():.3f} ± {m6a_positions.std():.3f}")
    logger.info(f"  中位数相对位置: {np.median(m6a_positions):.3f}")

    logger.info("\n" + "=" * 80)


# ==============================================================================
# 主函数
# ==============================================================================

def main():
    """
    主执行函数
    """
    logger.info("\n" + "=" * 80)
    logger.info("RNA修饰分析 - Metagene Profile分析模块")
    logger.info("=" * 80)

    # 加载数据
    psi_file = RESULTS_DIR / "psi_sites_annotated.csv"
    m6a_file = RESULTS_DIR / "m6a_sites_annotated.csv"

    if not psi_file.exists() or not m6a_file.exists():
        logger.error("数据文件不存在！请先运行 1_data_fetching.py")
        return

    psi_df = pd.read_csv(psi_file)
    m6a_df = pd.read_csv(m6a_file)

    logger.info(f"\n数据加载:")
    logger.info(f"  Ψ位点: {len(psi_df)}")
    logger.info(f"  m6A位点: {len(m6a_df)}")

    # 1. 计算Metagene Profile
    logger.info("\n[步骤1] 计算Metagene Profile...")
    psi_profiler = MetageneProfiler(psi_df, 'Pseudouridine (Ψ)')
    m6a_profiler = MetageneProfiler(m6a_df, 'm6A')

    psi_profile = psi_profiler.compute_metagene_profile(n_bins=100, smooth=True, window=9)
    m6a_profile = m6a_profiler.compute_metagene_profile(n_bins=100, smooth=True, window=9)

    # 保存profile数据
    psi_profile.to_csv(RESULTS_DIR / "psi_metagene_profile.csv", index=False)
    m6a_profile.to_csv(RESULTS_DIR / "m6a_metagene_profile.csv", index=False)

    logger.info("  ✓ Profile计算完成")

    # 2. 绘制对比图
    logger.info("\n[步骤2] 绘制Metagene Profile对比图...")
    plot_metagene_profile_comparison(
        psi_profile, m6a_profile,
        output_file=FIGURES_DIR / "metagene_profile_comparison.png"
    )

    # 3. 分析基因组特征分布
    logger.info("\n[步骤3] 分析基因组特征分布...")
    psi_features = psi_profiler.get_feature_distribution()
    m6a_features = m6a_profiler.get_feature_distribution()

    logger.info(f"\nΨ特征分布:")
    logger.info(psi_features)
    logger.info(f"\nm6A特征分布:")
    logger.info(m6a_features)

    # 4. 绘制特征分布对比图
    logger.info("\n[步骤4] 绘制特征分布对比图...")
    plot_feature_distribution_comparison(
        psi_features, m6a_features,
        output_file=FIGURES_DIR / "feature_distribution_comparison.png"
    )

    # 5. 绘制饼图
    logger.info("\n[步骤5] 绘制饼图...")
    plot_pie_charts(
        psi_features, m6a_features,
        output_file=FIGURES_DIR / "feature_distribution_piecharts.png"
    )

    # 6. 统计检验
    logger.info("\n[步骤6] 统计分析...")
    compare_distributions(psi_df, m6a_df)

    logger.info("\n" + "=" * 80)
    logger.info("Metagene Profile分析完成！")
    logger.info(f"图表保存在: {FIGURES_DIR}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
