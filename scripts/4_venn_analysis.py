#!/usr/bin/env python3
"""
韦恩图分析模块 - RNA修饰比较分析
功能：
1. 分析Ψ和m6A的共定位情况（Colocalization）
2. 绘制韦恩图展示交集和独有位点
3. 分析共定位位点的特征分布

生物学意义：
- 共定位可能表示两种修饰在相同RNA位点上存在互作或竞争关系
- 独特位点反映了各修饰的特异性功能

作者：周子航
日期：2026-01-07
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from matplotlib_venn import venn2, venn2_circles
import warnings
warnings.filterwarnings('ignore')

# 配置
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 路径设置
PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures"


# ==============================================================================
# 共定位分析
# ==============================================================================

class ColocalizationAnalyzer:
    """
    共定位分析器：识别和分析两种修饰的重叠位点
    """

    def __init__(self, psi_df: pd.DataFrame, m6a_df: pd.DataFrame,
                 window: int = 50):
        """
        初始化

        参数：
        - psi_df: Ψ位点DataFrame
        - m6a_df: m6A位点DataFrame
        - window: 共定位窗口大小（bp）
                  如果Ψ和m6A位点距离小于window，认为它们共定位
        """
        self.psi_df = psi_df
        self.m6a_df = m6a_df
        self.window = window

        # 创建位点ID（用于韦恩图）
        self.psi_df['site_id'] = ['PSI_' + str(i) for i in range(len(psi_df))]
        self.m6a_df['site_id'] = ['M6A_' + str(i) for i in range(len(m6a_df))]

        # 结果存储
        self.colocalized_sites = None
        self.psi_only = None
        self.m6a_only = None

    def find_colocalized_sites(self) -> pd.DataFrame:
        """
        寻找共定位位点

        返回：
        - colocalized_df: 共定位位点对DataFrame
        """
        logger.info(f"正在寻找共定位位点（窗口±{self.window}bp）...")

        colocalized_pairs = []

        # 遍历每个Ψ位点
        for psi_idx, psi_row in self.psi_df.iterrows():
            psi_chrom = psi_row['chrom']
            psi_center = (psi_row['start'] + psi_row['end']) / 2

            # 在m6A位点中寻找同一染色体上附近的位点
            nearby_m6a = self.m6a_df[
                (self.m6a_df['chrom'] == psi_chrom) &
                (abs((self.m6a_df['start'] + self.m6a_df['end'])/2 - psi_center) <= self.window)
            ]

            for m6a_idx, m6a_row in nearby_m6a.iterrows():
                distance = abs((psi_row['start'] + psi_row['end'])/2 -
                             (m6a_row['start'] + m6a_row['end'])/2)

                pair_info = {
                    'psi_id': psi_row['site_id'],
                    'm6a_id': m6a_row['site_id'],
                    'chrom': psi_chrom,
                    'psi_pos': int(psi_center),
                    'm6a_pos': int((m6a_row['start'] + m6a_row['end'])/2),
                    'distance': int(distance),
                    'psi_feature': psi_row.get('feature', 'Unknown'),
                    'm6a_feature': m6a_row.get('feature', 'Unknown'),
                    'psi_strand': psi_row.get('strand', '+'),
                    'm6a_strand': m6a_row.get('strand', '+')
                }
                colocalized_pairs.append(pair_info)

        if colocalized_pairs:
            self.colocalized_sites = pd.DataFrame(colocalized_pairs)
            logger.info(f"  找到 {len(self.colocalized_sites)} 对共定位位点")
        else:
            self.colocalized_sites = pd.DataFrame()
            logger.info("  未找到共定位位点")

        return self.colocalized_sites

    def categorize_sites(self):
        """
        将位点分类为：仅Ψ、仅m6A、共定位
        """
        if self.colocalized_sites is None:
            self.find_colocalized_sites()

        # 共定位的位点ID
        colocalized_psi_ids = set(self.colocalized_sites['psi_id']) if len(self.colocalized_sites) > 0 else set()
        colocalized_m6a_ids = set(self.colocalized_sites['m6a_id']) if len(self.colocalized_sites) > 0 else set()

        # 仅Ψ
        self.psi_only = self.psi_df[~self.psi_df['site_id'].isin(colocalized_psi_ids)]

        # 仅m6A
        self.m6a_only = self.m6a_df[~self.m6a_df['site_id'].isin(colocalized_m6a_ids)]

        logger.info("\n位点分类结果:")
        logger.info(f"  仅Ψ: {len(self.psi_only)} 个位点")
        logger.info(f"  仅m6A: {len(self.m6a_only)} 个位点")
        logger.info(f"  共定位: {len(colocalized_psi_ids)} 个Ψ位点 + {len(colocalized_m6a_ids)} 个m6A位点")

    def get_venn_counts(self):
        """
        获取韦恩图的计数

        返回：
        - (psi_only_count, m6a_only_count, overlap_count)
        """
        if self.colocalized_sites is None:
            self.find_colocalized_sites()

        overlap_count = len(set(self.colocalized_sites['psi_id'])) if len(self.colocalized_sites) > 0 else 0

        return (len(self.psi_only), len(self.m6a_only), overlap_count)


# ==============================================================================
# 可视化函数
# ==============================================================================

def plot_venn_diagram(psi_only: int, m6a_only: int, overlap: int,
                     output_file: str = None):
    """
    绘制韦恩图

    参数：
    - psi_only: 仅Ψ的位点数
    - m6a_only: 仅m6A的位点数
    - overlap: 共定位位点数
    - output_file: 输出文件路径
    """
    plt.figure(figsize=(10, 8))

    # 创建子图
    subset = (psi_only, m6a_only, overlap)

    # 绘制韦恩图
    v = venn2(subsets=subset, set_labels=('Pseudouridine (Ψ)', 'm6A'),
             set_colors=('#3498db', '#e74c3c'), alpha=0.6)

    # 自定义样式
    if v:
        # 设置标签字体
        for text in v.set_labels:
            text.set_fontsize(14)
            text.set_fontweight('bold')

        # 设置区域数字
        for text in v.subset_labels:
            if text:
                text.set_fontsize(16)
                text.set_fontweight('bold')

    # 添加圆圈边框
    venn2_circles(subsets=subset, linestyle='dashed', linewidth=2)

    # 添加标题和统计信息
    total_psi = psi_only + overlap
    total_m6a = m6a_only + overlap

    plt.title(f'Colocalization Analysis: Ψ vs m6A\n'
             f'Total Ψ: {total_psi:,} | Total m6A: {total_m6a:,} | Overlap: {overlap:,} '
             f'({overlap/max(total_psi, total_m6a)*100:.1f}% of Ψ sites)',
             fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"保存韦恩图: {output_file}")

    plt.show()


def plot_colocalization_features(colocalized_df: pd.DataFrame,
                                 output_file: str = None):
    """
    分析共定位位点的特征分布

    参数：
    - colocalized_df: 共定位位点DataFrame
    - output_file: 输出文件路径
    """
    if len(colocalized_df) == 0:
        logger.warning("没有共定位位点，跳过特征分析")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. 距离分布
    ax1 = axes[0, 0]
    ax1.hist(colocalized_df['distance'], bins=30, color='#9b59b6',
            edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Distance between sites (bp)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax1.set_title('Distribution of Distances\nbetween Colocalized Sites',
                 fontsize=12, fontweight='bold')
    ax1.axvline(colocalized_df['distance'].median(), color='red',
               linestyle='--', linewidth=2, label=f"Median: {colocalized_df['distance'].median():.0f}bp")
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # 2. 特征组合分布
    ax2 = axes[0, 1]
    feature_combinations = colocalized_df['psi_feature'] + ' (Ψ) + ' + colocalized_df['m6a_feature'] + ' (m6A)'
    feature_counts = feature_combinations.value_counts()

    colors = plt.cm.Set3(np.linspace(0, 1, len(feature_counts)))
    feature_counts.plot(kind='barh', ax=ax2, color=colors, edgecolor='black')
    ax2.set_xlabel('Count', fontsize=11, fontweight='bold')
    ax2.set_title('Feature Combinations at Colocalized Sites',
                 fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    # 3. 链方向一致性
    ax3 = axes[1, 0]
    strand_consistency = (colocalized_df['psi_strand'] == colocalized_df['m6a_strand']).value_counts()
    strand_labels = ['Same strand', 'Different strands']
    strand_colors = ['#2ecc71', '#e74c3c']

    ax3.pie(strand_consistency, labels=strand_labels, colors=strand_colors,
           autopct='%1.1f%%', startangle=90, explode=(0.05, 0.05),
           shadow=True, textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax3.set_title('Strand Consistency', fontsize=12, fontweight='bold', pad=15)

    # 4. 染色体分布
    ax4 = axes[1, 1]
    chrom_counts = colocalized_df['chrom'].value_counts().sort_index()
    ax4.bar(range(len(chrom_counts)), chrom_counts.values, color='#3498db',
           edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Chromosome', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Colocalized site count', fontsize=11, fontweight='bold')
    ax4.set_title('Chromosomal Distribution of Colocalized Sites',
                 fontsize=12, fontweight='bold')
    ax4.set_xticks(range(len(chrom_counts)))
    ax4.set_xticklabels([c.replace('chr', '') for c in chrom_counts.index],
                       rotation=45, ha='right')
    ax4.grid(axis='y', alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"保存共定位特征分析图: {output_file}")

    plt.show()


def analyze_colocalization_significance(psi_count: int, m6a_count: int,
                                       overlap_count: int,
                                       genome_size: int = 3e9):
    """
    统计检验：共定位是否显著高于随机期望

    使用超几何分布检验

    参数：
    - psi_count: Ψ位点总数
    - m6a_count: m6A位点总数
    - overlap_count: 实际观察到共定位位点数
    - genome_size: 基因组大小（bp）
    """
    logger.info("\n" + "=" * 80)
    logger.info("共定位显著性分析")
    logger.info("=" * 80 + "\n")

    # 计算随机期望值
    # 假设位点随机分布在基因组上
    psi_density = psi_count / genome_size
    expected_overlap = m6a_count * psi_density

    logger.info(f"基于随机分布的期望值:")
    logger.info(f"  Ψ位点密度: {psi_density:.2e} per bp")
    logger.info(f"  期望共定位位点数: {expected_overlap:.2f}")
    logger.info(f"  实际观察共定位位点数: {overlap_count}")

    # 计算富集倍数
    if expected_overlap > 0:
        enrichment = overlap_count / expected_overlap
        logger.info(f"  富集倍数: {enrichment:.2f}x")

    # 简单的卡方检验
    observed = [overlap_count, psi_count - overlap_count]
    expected = [expected_overlap, psi_count - expected_overlap]

    from scipy import stats
    chi2, p_value = stats.chisquare(f_obs=observed, f_exp=expected)

    logger.info(f"\n统计检验（卡方检验）:")
    logger.info(f"  Chi2 = {chi2:.4f}")
    logger.info(f"  P-value = {p_value:.2e}")

    if p_value < 0.001:
        logger.info(f"  结论: 共定位显著高于随机期望 *** (P < 0.001)")
    elif p_value < 0.01:
        logger.info(f"  结论: 共定位显著高于随机期望 ** (P < 0.01)")
    elif p_value < 0.05:
        logger.info(f"  结论: 共定位高于随机期望 * (P < 0.05)")
    else:
        logger.info(f"  结论: 共定位未达到显著水平")

    logger.info("\n" + "=" * 80)


# ==============================================================================
# 主函数
# ==============================================================================

def main():
    """
    主执行函数
    """
    logger.info("\n" + "=" * 80)
    logger.info("RNA修饰分析 - 韦恩图分析模块")
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
    logger.info(f"  Ψ位点: {len(psi_df):,}")
    logger.info(f"  m6A位点: {len(m6a_df):,}")

    # 1. 寻找共定位位点
    logger.info("\n[步骤1] 寻找共定位位点...")
    analyzer = ColocalizationAnalyzer(psi_df, m6a_df, window=50)
    colocalized_df = analyzer.find_colocalized_sites()

    # 2. 位点分类
    logger.info("\n[步骤2] 位点分类...")
    analyzer.categorize_sites()

    # 3. 绘制韦恩图
    logger.info("\n[步骤3] 绘制韦恩图...")
    psi_only, m6a_only, overlap = analyzer.get_venn_counts()
    plot_venn_diagram(psi_only, m6a_only, overlap,
                     output_file=FIGURES_DIR / "venn_diagram.png")

    # 4. 如果有共定位位点，分析其特征
    if len(colocalized_df) > 0:
        logger.info("\n[步骤4] 分析共定位位点特征...")
        plot_colocalization_features(colocalized_df,
                                    output_file=FIGURES_DIR / "colocalization_features.png")

        # 保存共定位位点信息
        colocalized_file = RESULTS_DIR / "colocalized_sites.csv"
        colocalized_df.to_csv(colocalized_file, index=False)
        logger.info(f"  保存共定位位点: {colocalized_file}")

    # 5. 统计显著性分析
    logger.info("\n[步骤5] 显著性分析...")
    analyze_colocalization_significance(
        len(psi_df), len(m6a_df), overlap,
        genome_size=3e9
    )

    # 6. 生成分析报告
    logger.info("\n[步骤6] 生成分析报告...")
    report = f"""
{'='*80}
共定位分析报告
{'='*80}

位点统计:
  Ψ位点总数: {len(psi_df):,}
  m6A位点总数: {len(m6a_df):,}
  仅Ψ: {psi_only:,} ({psi_only/len(psi_df)*100:.1f}%)
  仅m6A: {m6a_only:,} ({m6a_only/len(m6a_df)*100:.1f}%)
  共定位: {overlap:,} ({overlap/len(psi_df)*100:.1f}% of Ψ, {overlap/len(m6a_df)*100:.1f}% of m6A)

共定位位点特征:
"""

    if len(colocalized_df) > 0:
        report += f"""
  - 平均距离: {colocalized_df['distance'].mean():.1f} ± {colocalized_df['distance'].std():.1f} bp
  - 中位数距离: {colocalized_df['distance'].median():.0f} bp
  - 最大距离: {colocalized_df['distance'].max():.0f} bp
  - 最小距离: {colocalized_df['distance'].min():.0f} bp
  - 链方向一致: {(colocalized_df['psi_strand'] == colocalized_df['m6a_strand']).sum()}/{len(colocalized_df)} ({(colocalized_df['psi_strand'] == colocalized_df['m6a_strand']).sum()/len(colocalized_df)*100:.1f}%)
"""
    else:
        report += "\n  未检测到共定位位点\n"

    report += f"\n{'='*80}\n"

    logger.info(report)

    # 保存报告
    report_file = RESULTS_DIR / "colocalization_report.txt"
    with open(report_file, 'w') as f:
        f.write(report)
    logger.info(f"保存报告: {report_file}")

    logger.info("\n" + "=" * 80)
    logger.info("韦恩图分析完成！")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
