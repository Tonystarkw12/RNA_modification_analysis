#!/usr/bin/env python3
"""
Motif分析模块 - RNA修饰比较分析
功能：
1. 从基因组FASTA中提取Ψ和m6A位点周围的序列
2. 绘制Sequence Logo
3. 验证是否存在已知的motif特征（Ψ: GUUC/UGU; m6A: GGACU/DRACH）

依赖：
- pyfaidx: 用于快速访问FASTA文件
- logomaker: 用于绘制sequence logo
- pandas, numpy, matplotlib

作者：周子航
日期：2026-01-07
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from collections import Counter
from typing import List, Tuple
import warnings
warnings.filterwarnings('ignore')

# 尝试导入专业包
try:
    from pyfaidx import Fasta
    HAS_PYFAIDX = True
except ImportError:
    HAS_PYFAIDX = False
    print("警告: pyfaidx未安装，将使用模拟序列")

try:
    import logomaker
    HAS_LOGOMAKER = True
except ImportError:
    HAS_LOGOMAKER = False
    print("警告: logomaker未安装，将使用matplotlib绘制简单的motif图")

# 配置日志和绘图
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False

# 路径设置
PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
FIGURES_DIR = PROJECT_DIR / "figures"

FIGURES_DIR.mkdir(exist_ok=True)


# ==============================================================================
# 序列提取功能
# ==============================================================================

class SequenceExtractor:
    """
    从基因组中提取指定位置的序列
    """

    def __init__(self, fasta_path: str):
        """
        初始化

        参数：
        - fasta_path: 基因组FASTA文件路径
        """
        self.fasta_path = fasta_path

        if HAS_PYFAIDX and os.path.exists(fasta_path):
            logger.info(f"加载基因组: {fasta_path}")
            self.genome = Fasta(fasta_path)
            self.has_genome = True
        else:
            logger.warning(f"无法加载基因组 {fasta_path}，将使用模拟序列")
            self.has_genome = False
            self.genome = None

    def extract_sequence(self, chrom: str, start: int, end: int,
                        flank: int = 10, strand: str = '+') -> str:
        """
        提取指定位点周围序列

        参数：
        - chrom: 染色体（如'chr1'）
        - start: 位点起始位置（0-based）
        - end: 位点结束位置
        - flank: 侧翼区域长度（bp）
        - strand: 链方向（'+'或'-'）

        返回：
        - 序列字符串（长度为 2*flank + 1）
        """
        if self.has_genome:
            try:
                # 提取序列（包含侧翼区域）
                seq_start = max(0, start - flank)
                seq_end = end + flank

                # pyfaidx使用1-based坐标
                sequence = self.genome[chrom][seq_start:seq_end].seq

                # 如果是负链，需要反向互补
                if strand == '-':
                    sequence = self._reverse_complement(sequence)

                # 取中心区域（确保长度正确）
                center = flank
                seq_extracted = sequence[center-flank:center+flank+1]

                return seq_extracted.upper()

            except Exception as e:
                logger.warning(f"提取序列失败: {chrom}:{start}-{end}, 错误: {e}")
                return self._generate_random_sequence(2*flank + 1)

        else:
            # 使用模拟序列
            return self._generate_random_sequence(2*flank + 1)

    def extract_sequences_batch(self, sites_df: pd.DataFrame,
                               flank: int = 10) -> List[str]:
        """
        批量提取序列

        参数：
        - sites_df: 位点信息DataFrame，必须包含列: chrom, start, end, strand
        - flank: 侧翼区域长度

        返回：
        - 序列列表
        """
        sequences = []

        for idx, row in sites_df.iterrows():
            seq = self.extract_sequence(
                chrom=row['chrom'],
                start=int(row['start']),
                end=int(row['end']),
                flank=flank,
                strand=row.get('strand', '+')
            )
            sequences.append(seq)

            if (idx + 1) % 1000 == 0:
                logger.info(f"已提取 {idx + 1}/{len(sites_df)} 条序列")

        return sequences

    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """
        生成反向互补序列
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                     'N': 'N', 'U': 'A', 'a': 't', 't': 'a',
                     'c': 'g', 'g': 'c', 'n': 'n', 'u': 'a'}

        return ''.join([complement.get(base, base) for base in reversed(sequence)])

    @staticmethod
    def _generate_random_sequence(length: int) -> str:
        """
        生成随机序列（模拟用）
        """
        bases = ['A', 'U', 'C', 'G']
        return ''.join(np.random.choice(bases, length))


# ==============================================================================
# Motif分析和可视化
# ==============================================================================

class MotifAnalyzer:
    """
    Motif分析器：提取序列并绘制Sequence Logo
    """

    def __init__(self, sequences: List[str], motif_name: str = "Motif"):
        """
        初始化

        参数：
        - sequences: 序列列表
        - motif_name: motif名称（用于图标题）
        """
        self.sequences = [seq.upper().replace('T', 'U') for seq in sequences]  # 统一为U
        self.motif_name = motif_name
        self.matrix = None

    def compute_pwm(self) -> pd.DataFrame:
        """
        计算位置权重矩阵（Position Weight Matrix）

        返回：
        - PWM DataFrame (位置 x 碱基)
        """
        if not self.sequences:
            logger.error("序列列表为空！")
            return pd.DataFrame()

        seq_length = len(self.sequences[0])
        bases = ['A', 'U', 'C', 'G']

        # 初始化计数矩阵
        count_matrix = np.zeros((seq_length, len(bases)))

        # 统计每个位置各碱基的出现次数
        for seq in self.sequences:
            if len(seq) != seq_length:
                logger.warning(f"序列长度不一致: {len(seq)} vs {seq_length}")
                continue

            for pos, base in enumerate(seq):
                if base in bases:
                    col_idx = bases.index(base)
                    count_matrix[pos, col_idx] += 1
                else:
                    # 非标准碱基（如N），平均分配
                    count_matrix[pos, :] += 0.25

        # 添加伪计数（避免log(0)）
        count_matrix += 0.5

        # 转换为频率矩阵
        freq_matrix = count_matrix / count_matrix.sum(axis=1, keepdims=True)

        # 转换为DataFrame
        self.matrix = pd.DataFrame(freq_matrix, columns=bases)

        return self.matrix

    def find_enriched_kmers(self, k: int = 4, top_n: int = 10) -> List[Tuple[str, float]]:
        """
        寻找富集的k-mer

        参数：
        - k: k-mer长度
        - top_n: 返回前N个富集的k-mer

        返回：
        - [(kmer, frequency), ...]
        """
        kmer_counts = Counter()

        for seq in self.sequences:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                if 'N' not in kmer:
                    kmer_counts[kmer] += 1

        total = sum(kmer_counts.values())
        kmer_freq = [(kmer, count/total) for kmer, count in kmer_counts.most_common(top_n)]

        return kmer_freq

    def plot_logo(self, output_file: str = None, figsize: Tuple[int, int] = (12, 4)):
        """
        绘制Sequence Logo

        参数：
        - output_file: 输出文件路径（可选）
        - figsize: 图形尺寸
        """
        if self.matrix is None:
            self.compute_pwm()

        fig, ax = plt.subplots(figsize=figsize)

        if HAS_LOGOMAKER:
            # 使用logomaker绘制专业logo
            try:
                # 确保行索引从1开始
                matrix = self.matrix.copy()
                matrix.index = np.arange(1, len(matrix) + 1)

                # 创建logo
                logo = logomaker.Logo(matrix, ax=ax)

                # 样式设置
                logo.style_spines(visible=False)
                logo.style_spines(spines=['left', 'bottom'], visible=True)
                logo.style_xticks(rotation=0, fmt='%d', anchor=0)

                ax.set_ylabel('Bits', fontsize=12, fontweight='bold')
                ax.set_xlabel('Position relative to modification site', fontsize=12, fontweight='bold')
                ax.set_title(f'{self.motif_name} Sequence Logo (n={len(self.sequences)})',
                            fontsize=14, fontweight='bold')

            except Exception as e:
                logger.warning(f"logomaker绘图失败: {e}，使用matplotlib")
                HAS_LOGOMAKER = False

        if not HAS_LOGOMAKER:
            # 使用matplotlib绘制简单版
            self._plot_logo_simple(ax)

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"保存Logo图: {output_file}")

        plt.show()

    def _plot_logo_simple(self, ax):
        """
        使用matplotlib绘制简单的sequence logo（字母高度表示频率）
        """
        matrix = self.matrix.values
        bases = ['A', 'U', 'C', 'G']
        colors = {'A': '#00CC00', 'U': '#FF0000', 'C': '#0000CC', 'G': '#FFB300'}
        seq_length = len(matrix)

        ax.set_xlim(0, seq_length)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax.set_xlabel('Position', fontsize=12, fontweight='bold')
        ax.set_title(f'{self.motif_name} Motif (n={len(self.sequences)})',
                    fontsize=14, fontweight='bold')

        # 在每个位置绘制碱基
        for pos in range(seq_length):
            x = pos
            y_bottom = 0

            # 按频率从高到低排序
            freqs = [(bases[i], matrix[pos, i]) for i in range(4)]
            freqs.sort(key=lambda x: x[1], reverse=True)

            for base, freq in freqs:
                if freq > 0.01:  # 只显示频率>1%的碱基
                    # 绘制矩形背景
                    rect = plt.Rectangle((x - 0.4, y_bottom), 0.8, freq,
                                        facecolor=colors[base], edgecolor='black', linewidth=0.5)
                    ax.add_patch(rect)

                    # 添加字母
                    ax.text(x, y_bottom + freq/2, base,
                           ha='center', va='center', fontsize=10,
                           fontweight='bold', color='white')

                    y_bottom += freq

        # 添加中心线（修饰位点）
        center_pos = seq_length // 2
        ax.axvline(x=center_pos + 0.5, color='black', linestyle='--', alpha=0.5)
        ax.text(center_pos + 0.5, 1.02, 'Modified site', ha='center', fontsize=9)


# ==============================================================================
# 主函数
# ==============================================================================

def analyze_psi_motif(psi_df, fasta_path, flank=10):
    """
    分析Ψ位点的motif
    """
    logger.info("\n" + "=" * 80)
    logger.info("分析Pseudouridine (Ψ) Motif")
    logger.info("=" * 80 + "\n")

    # 提取序列
    extractor = SequenceExtractor(fasta_path)
    sequences = extractor.extract_sequences_batch(psi_df, flank=flank)

    logger.info(f"成功提取 {len(sequences)} 条序列（侧翼±{flank}bp）")

    # 分析motif
    analyzer = MotifAnalyzer(sequences, motif_name="Pseudouridine (Ψ)")

    # 计算PWM
    pwm = analyzer.compute_pwm()
    logger.info("\nPWM矩阵形状:", pwm.shape)

    # 寻找富集的4-mer
    logger.info("\n富集的4-mer motifs:")
    top_kmers = analyzer.find_enriched_kmers(k=4, top_n=10)
    for kmer, freq in top_kmers:
        logger.info(f"  {kmer}: {freq:.4f}")

    # 验证已知motif
    logger.info("\n验证已知Ψ motifs:")
    known_motifs = ['GUUC', 'UGU', 'UUCG', 'GUC']
    for motif in known_motifs:
        count = sum(1 for seq in sequences if motif in seq)
        logger.info(f"  {motif}: {count}/{len(sequences)} ({count/len(sequences)*100:.2f}%)")

    # 绘制logo
    output_file = FIGURES_DIR / "psi_motif_logo.png"
    analyzer.plot_logo(output_file=output_file)

    # 保存PWM矩阵
    pwm_file = RESULTS_DIR / "psi_pwm.csv"
    pwm.to_csv(pwm_file)
    logger.info(f"\n保存PWM矩阵: {pwm_file}")

    return analyzer


def analyze_m6a_motif(m6a_df, fasta_path, flank=10):
    """
    分析m6A位点的motif
    """
    logger.info("\n" + "=" * 80)
    logger.info("分析m6A Motif")
    logger.info("=" * 80 + "\n")

    # 提取序列
    extractor = SequenceExtractor(fasta_path)
    sequences = extractor.extract_sequences_batch(m6a_df, flank=flank)

    logger.info(f"成功提取 {len(sequences)} 条序列（侧翼±{flank}bp）")

    # 分析motif
    analyzer = MotifAnalyzer(sequences, motif_name="m6A")

    # 计算PWM
    pwm = analyzer.compute_pwm()
    logger.info("\nPWM矩阵形状:", pwm.shape)

    # 寻找富集的4-mer
    logger.info("\n富集的4-mer motifs:")
    top_kmers = analyzer.find_enriched_kmers(k=4, top_n=10)
    for kmer, freq in top_kmers:
        logger.info(f"  {kmer}: {freq:.4f}")

    # 验证已知motif（DRACH motif）
    logger.info("\n验证已知m6A motifs (DRACH):")
    known_motifs = ['GGACU', 'GGACA', 'GGACT', 'AGACU', 'GAACU']
    for motif in known_motifs:
        count = sum(1 for seq in sequences if motif in seq)
        logger.info(f"  {motif}: {count}/{len(sequences)} ({count/len(sequences)*100:.2f}%)")

    # 绘制logo
    output_file = FIGURES_DIR / "m6a_motif_logo.png"
    analyzer.plot_logo(output_file=output_file)

    # 保存PWM矩阵
    pwm_file = RESULTS_DIR / "m6a_pwm.csv"
    pwm.to_csv(pwm_file)
    logger.info(f"\n保存PWM矩阵: {pwm_file}")

    return analyzer


def main():
    """
    主执行函数
    """
    logger.info("\n" + "=" * 80)
    logger.info("RNA修饰分析 - Motif分析模块")
    logger.info("=" * 80)

    # 加载数据
    psi_file = RESULTS_DIR / "psi_sites_annotated.csv"
    m6a_file = RESULTS_DIR / "m6a_sites_annotated.csv"

    if not psi_file.exists() or not m6a_file.exists():
        logger.error("数据文件不存在！请先运行 1_data_fetching.py")
        return

    psi_df = pd.read_csv(psi_file)
    m6a_df = pd.read_csv(m6a_file)

    # 为了加快速度，只分析部分位点
    n_sites = 2000
    psi_subset = psi_df.sample(n=min(n_sites, len(psi_df)), random_state=42)
    m6a_subset = m6a_df.sample(n=min(n_sites, len(m6a_df)), random_state=42)

    logger.info(f"Ψ位点数量: {len(psi_subset)}")
    logger.info(f"m6A位点数量: {len(m6a_subset)}")

    # 基因组FASTA文件
    fasta_path = DATA_DIR / "GRCh38.p14.genome.fa"

    # 如果没有真实基因组，使用模拟数据
    if not fasta_path.exists():
        logger.warning("\n真实基因组FASTA文件不存在！")
        logger.info("建议下载:")
        logger.info("  wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz")
        logger.info("  gunzip GRCh38.p14.genome.fa.gz")
        logger.info("\n使用模拟序列继续分析...")

        # 创建一个虚拟的FASTA路径
        fasta_path = DATA_DIR / "mock_hg38.fa"

    # 分析Ψ motif
    analyze_psi_motif(psi_subset, fasta_path, flank=10)

    # 分析m6A motif
    analyze_m6a_motif(m6a_subset, fasta_path, flank=10)

    logger.info("\n" + "=" * 80)
    logger.info("Motif分析完成！")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
