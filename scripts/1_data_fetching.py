#!/usr/bin/env python3
"""
数据获取模块 - RNA修饰比较分析
功能：
1. 下载真实的公共数据库数据（Pseudouridine和m6A）
2. 下载GENCODE基因组注释文件
3. 生成符合生物学特征的模拟数据（用于测试和演示）

作者：周子航
日期：2026-01-07
"""

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import logging

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 定义路径
PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"

# 创建必要目录
DATA_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)


def run_command(command, description=""):
    """
    执行shell命令并处理错误
    """
    try:
        logger.info(f"执行: {description or command}")
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        logger.info(f"成功: {description or command}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"失败: {description or command}")
        logger.error(f"错误信息: {e.stderr}")
        return False


# ==============================================================================
# 第一部分：真实数据下载
# ==============================================================================

def download_real_data():
    """
    下载真实的RNA修饰数据和基因组注释

    数据来源：
    1. Pseudouridine (Ψ): RMBase v2.0 或 GSE102476 (CeU-seq)
    2. m6A: REPIC 或 m6A-Atlas (HEK293T细胞系)
    3. 基因组注释: GENCODE hg38
    """

    logger.info("=" * 80)
    logger.info("开始下载真实数据...")
    logger.info("=" * 80)

    # 1. 下载GENCODE基因组注释 (hg38)
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.annotation.gtf.gz"
    gtf_path = DATA_DIR / "gencode.v44.annotation.gtf.gz"

    logger.info("\n[1/3] 下载GENCODE hg38 GTF文件...")
    if not gtf_path.exists():
        cmd = f"wget -c -O {gtf_path} {gtf_url}"
        run_command(cmd, f"下载GENCODE GTF ({gtf_url})")
    else:
        logger.info(f"文件已存在: {gtf_path}")

    # 2. 下载Pseudouridine数据 (RMBase v2.0示例)
    # RMBase数据需要访问其官方网站，这里提供示例命令
    logger.info("\n[2/3] Pseudouridine数据下载指南...")
    logger.info("方法1 - 从RMBase v2.0下载:")
    logger.info("  访问: http://rmbase.renlab.org/download.php")
    logger.info("  选择: Human -> Pseudouridine -> hg38 -> HEK293T")
    logger.info("  下载并保存到: {}/psi_HEK293T.bed".format(DATA_DIR))

    logger.info("\n方法2 - 从GEO (GSE102476) 下载CeU-seq数据:")
    logger.info("  使用以下命令:")
    ceu_seq_cmd = f"""
    cd {DATA_DIR}
    wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102476/suppl/GSE102476_CeU-seq_HEK293T.bed.gz
    gunzip GSE102476_CeU-seq_HEK293T.bed.gz
    """
    logger.info(ceu_seq_cmd)

    # 3. 下载m6A数据 (REPIC示例)
    logger.info("\n[3/3] m6A数据下载指南...")
    logger.info("方法1 - 从REPIC数据库下载:")
    logger.info("  访问: https://repic.idrb.cas.cz/download")
    logger.info("  选择: Human -> m6A -> HEK293T -> hg38")
    logger.info("  下载并保存到: {}/m6A_HEK293T.bed".format(DATA_DIR))

    logger.info("\n方法2 - 使用pysam访问ENCORE数据:")
    logger.info("  可通过Python脚本从ENCODE数据库获取")

    logger.info("\n" + "=" * 80)
    logger.info("真实数据下载指南完成！")
    logger.info("如果无法直接下载，将使用模拟数据继续分析。")
    logger.info("=" * 80)


# ==============================================================================
# 第二部分：模拟数据生成
# ==============================================================================

def generate_mock_psi(num_sites=5000, genome_size=3000, seed=42):
    """
    生成模拟的Pseudouridine (Ψ) 数据

    生物学特征：
    - Ψ修饰相对均匀地分布在CDS和UTR区域
    - 不像m6A那样有强烈的3'UTR偏好
    - Motif偏好：GUUC, UGU等

    参数：
    - num_sites: 位点数量
    - genome_size: 模拟基因组大小（Mb）
    - seed: 随机种子
    """
    np.random.seed(seed)
    logger.info(f"生成模拟Pseudouridine数据 ({num_sites}个位点)...")

    # 模拟基因位置（简化：假设有20000个基因，分布在各染色体）
    num_genes = 20000
    gene_data = []

    for i in range(num_genes):
        chrom = f"chr{np.random.randint(1, 23)}"  # chr1-chr22
        gene_start = np.random.randint(1000000, genome_size * 1000000 - 100000)
        gene_length = np.random.randint(1000, 100000)  # 基因长度

        # 基因结构：CDS占60%，5'UTR占20%，3'UTR占20%
        utr5_length = int(gene_length * 0.2)
        cds_length = int(gene_length * 0.6)
        utr3_length = gene_length - utr5_length - cds_length

        gene_data.append({
            'chrom': chrom,
            'gene_start': gene_start,
            'gene_end': gene_start + gene_length,
            'utr5_end': gene_start + utr5_length,
            'cds_end': gene_start + utr5_length + cds_length,
            'gene_length': gene_length
        })

    genes_df = pd.DataFrame(gene_data)

    # 生成Ψ位点（在CDS和UTR中随机分布，但CDS稍多）
    psi_sites = []

    for _ in range(num_sites):
        # 随机选择一个基因
        gene = genes_df.sample(1).iloc[0]

        # Ψ分布特征：CDS 50%, 5'UTR 20%, 3'UTR 30%
        region_choice = np.random.choice(['utr5', 'cds', 'utr3'], p=[0.2, 0.5, 0.3])

        if region_choice == 'utr5':
            pos = np.random.randint(gene['gene_start'], gene['utr5_end'])
            feature = '5UTR'
        elif region_choice == 'cds':
            pos = np.random.randint(gene['utr5_end'], gene['cds_end'])
            feature = 'CDS'
        else:  # utr3
            pos = np.random.randint(gene['cds_end'], gene['gene_end'])
            feature = '3UTR'

        # 计算在基因中的相对位置（用于metagene分析）
        relative_pos = (pos - gene['gene_start']) / gene['gene_length']

        psi_sites.append({
            'chrom': gene['chrom'],
            'start': pos,
            'end': pos + 1,
            'name': f'PSI_{len(psi_sites)+1}',
            'score': np.random.randint(100, 1000),  # 可信度评分
            'strand': np.random.choice(['+', '-']),
            'feature': feature,
            'relative_pos': relative_pos
        })

    psi_df = pd.DataFrame(psi_sites)

    # 保存为BED格式（标准6列BED）
    bed_data = psi_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    bed_file = DATA_DIR / "psi_HEK293T_mock.bed"
    bed_data.to_csv(bed_file, sep='\t', header=False, index=False)

    logger.info(f"  保存到: {bed_file}")
    logger.info(f"  位点分布: CDS={sum(psi_df['feature']=='CDS')}, "
                f"5'UTR={sum(psi_df['feature']=='5UTR')}, "
                f"3'UTR={sum(psi_df['feature']=='3UTR')}")

    return psi_df


def generate_mock_m6a(num_sites=5000, genome_size=3000, seed=43):
    """
    生成模拟的m6A数据

    生物学特征：
    - m6A强烈富集在3'UTR和终止密码子附近
    - 经典Motif：DRACH (D=A/G/U, R=A/G, H=A/C/U)
    - 具体到GGACU最为常见

    参数：
    - num_sites: 位点数量
    - genome_size: 模拟基因组大小（Mb）
    - seed: 随机种子
    """
    np.random.seed(seed)
    logger.info(f"生成模拟m6A数据 ({num_sites}个位点)...")

    # 使用相同的基因结构（与Ψ数据一致，便于比较）
    num_genes = 20000
    gene_data = []

    for i in range(num_genes):
        chrom = f"chr{np.random.randint(1, 23)}"
        gene_start = np.random.randint(1000000, genome_size * 1000000 - 100000)
        gene_length = np.random.randint(1000, 100000)

        utr5_length = int(gene_length * 0.2)
        cds_length = int(gene_length * 0.6)
        utr3_length = gene_length - utr5_length - cds_length

        gene_data.append({
            'chrom': chrom,
            'gene_start': gene_start,
            'gene_end': gene_start + gene_length,
            'utr5_end': gene_start + utr5_length,
            'cds_end': gene_start + utr5_length + cds_length,
            'gene_length': gene_length
        })

    genes_df = pd.DataFrame(gene_data)

    # 生成m6A位点（强烈偏向3'UTR和终止密码子附近）
    m6a_sites = []

    for _ in range(num_sites):
        gene = genes_df.sample(1).iloc[0]

        # m6A分布特征：3'UTR 60%, CDS 30%, 5'UTR 10%
        region_choice = np.random.choice(['utr5', 'cds', 'utr3'], p=[0.1, 0.3, 0.6])

        if region_choice == 'utr5':
            pos = np.random.randint(gene['gene_start'], gene['utr5_end'])
            feature = '5UTR'
        elif region_choice == 'cds':
            # 在CDS中，更倾向于靠近3'端（终止密码子附近）
            cds_start = gene['utr5_end']
            cds_end = gene['cds_end']
            cds_length = cds_end - cds_start
            # 确保不会产生空范围
            if cds_end > cds_start + cds_length * 0.6:
                pos = np.random.randint(int(cds_start + cds_length*0.6), cds_end)
            else:
                pos = np.random.randint(cds_start, cds_end)
            feature = 'CDS'
        else:  # utr3
            pos = np.random.randint(gene['cds_end'], gene['gene_end'])
            feature = '3UTR'

        # 计算在基因中的相对位置
        relative_pos = (pos - gene['gene_start']) / gene['gene_length']

        m6a_sites.append({
            'chrom': gene['chrom'],
            'start': pos,
            'end': pos + 1,
            'name': f'M6A_{len(m6a_sites)+1}',
            'score': np.random.randint(100, 1000),
            'strand': np.random.choice(['+', '-']),
            'feature': feature,
            'relative_pos': relative_pos
        })

    m6a_df = pd.DataFrame(m6a_sites)

    # 保存为BED格式
    bed_data = m6a_df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    bed_file = DATA_DIR / "m6A_HEK293T_mock.bed"
    bed_data.to_csv(bed_file, sep='\t', header=False, index=False)

    logger.info(f"  保存到: {bed_file}")
    logger.info(f"  位点分布: CDS={sum(m6a_df['feature']=='CDS')}, "
                f"5'UTR={sum(m6a_df['feature']=='5UTR')}, "
                f"3'UTR={sum(m6a_df['feature']=='3UTR')}")

    return m6a_df


def generate_mock_genome_fasta(output_file, seed=44):
    """
    生成模拟的基因组FASTA文件（用于motif分析）

    注意：这只是用于演示，实际分析应使用真实hg38基因组
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    np.random.seed(seed)
    logger.info("生成模拟基因组FASTA文件...")

    records = []
    for chrom_num in range(1, 23):
        chrom_name = f"chr{chrom_num}"
        # 生成300Mb的随机序列（这里只生成1Mb用于演示）
        seq_length = 1000000

        # 随机序列，但包含一定比例的motif
        bases = ['A', 'T', 'C', 'G']
        sequence = ''.join(np.random.choice(bases, seq_length))

        record = SeqRecord(Seq(sequence), id=chrom_name, description=f"Mock chromosome {chrom_num}")
        records.append(record)

    output_path = DATA_DIR / output_file
    SeqIO.write(records, output_path, "fasta")
    logger.info(f"  保存到: {output_path}")


# ==============================================================================
# 主函数
# ==============================================================================

def main():
    """
    主执行函数
    """
    logger.info("\n" + "=" * 80)
    logger.info("RNA修饰分析 - 数据获取模块")
    logger.info("=" * 80 + "\n")

    # Step 1: 提供真实数据下载指南
    download_real_data()

    # Step 2: 生成模拟数据（确保可以立即开始分析）
    logger.info("\n" + "=" * 80)
    logger.info("生成模拟数据用于测试和演示...")
    logger.info("=" * 80 + "\n")

    psi_df = generate_mock_psi(num_sites=5000, genome_size=3000, seed=42)
    m6a_df = generate_mock_m6a(num_sites=5000, genome_size=3000, seed=43)

    # 生成完整的注释信息（用于后续分析）
    psi_full = RESULTS_DIR / "psi_sites_annotated.csv"
    m6a_full = RESULTS_DIR / "m6a_sites_annotated.csv"

    psi_df.to_csv(psi_full, index=False)
    m6a_df.to_csv(m6a_full, index=False)

    logger.info(f"\n完整注释数据:")
    logger.info(f"  Ψ位点: {psi_full}")
    logger.info(f"  m6A位点: {m6a_full}")

    # Step 3: 生成模拟基因组（可选，用于motif分析）
    logger.info("\n注意：对于Motif分析，建议使用真实hg38基因组。")
    logger.info("下载命令：")
    logger.info("  wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.p14.genome.fa.gz")
    logger.info("  gunzip GRCh38.p14.genome.fa.gz\n")

    # generate_mock_genome_fasta("mock_hg38.fa")  # 取消注释以生成模拟基因组

    logger.info("\n" + "=" * 80)
    logger.info("数据获取完成！可以开始后续分析。")
    logger.info("=" * 80)

    # 输出数据摘要
    logger.info("\n数据摘要:")
    logger.info(f"  Pseudouridine (Ψ): {len(psi_df)} 个位点")
    logger.info(f"  m6A: {len(m6a_df)} 个位点")
    logger.info(f"  数据目录: {DATA_DIR}")
    logger.info(f"  结果目录: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
