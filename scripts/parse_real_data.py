#!/usr/bin/env python3
"""
解析真实RNA修饰数据
- m6A: REPIC数据库 HEK293T细胞数据（725,587个位点）
- Ψ: GSE102476仅包含rRNA数据，将使用高质量模拟数据
"""

import pandas as pd
import numpy as np
import gzip
import os
from pathlib import Path

class RealDataParser:
    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"
        self.data_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)

    def parse_m6a_repic(self, gz_file: str,
                       min_reads: int = 10,
                       fdr_threshold: float = 0.05) -> pd.DataFrame:
        """
        解析REPIC m6A数据

        格式：
        pos,pvalue,fdr,fold_enrichment,region,dataset,sample,tool,cell/tissue,species,assembly,geneid
        chr1:629848-630037[+]	1.00e-268	1.00e-265	29.6	exon	SRP007335	SRX146452_SRX146453	exomePeak	HEK293T	human	hg38	ENSG00000225630.1
        """
        print(f"\n=== 解析m6A REPIC数据 ===")
        print(f"输入文件: {gz_file}")

        # 读取并解析数据
        sites = []
        skipped = 0

        with gzip.open(gz_file, 'rt') as f:
            header = f.readline().strip().split('\t')
            print(f"列名: {header}")

            pos_idx = header.index('pos')
            fdr_idx = header.index('fdr')
            region_idx = header.index('region')
            gene_idx = header.index('geneid')

            for line_num, line in enumerate(f, 2):
                try:
                    cols = line.strip().split('\t')
                    if len(cols) < len(header):
                        skipped += 1
                        continue

                    pos_str = cols[pos_idx]
                    fdr = float(cols[fdr_idx])
                    region = cols[region_idx]
                    gene_id = cols[gene_idx]

                    # FDR过滤
                    if fdr > fdr_threshold:
                        skipped += 1
                        continue

                    # 解析染色体位置: chr1:629848-630037[+]
                    # 提取中间点作为peak位置
                    try:
                        chr_part, rest = pos_str.split(':')
                        range_part, strand = rest.split('[')
                        strand = strand.rstrip(']')

                        start, end = map(int, range_part.split('-'))
                        peak_pos = (start + end) // 2  # 使用中点

                        sites.append({
                            'chrom': chr_part,
                            'start': start,
                            'end': end,
                            'peak_pos': peak_pos,
                            'strand': strand,
                            'fdr': fdr,
                            'region': region,
                            'gene_id': gene_id,
                            'score': 1 - fdr  # 转换为得分
                        })
                    except Exception as e:
                        skipped += 1
                        continue

                except Exception as e:
                    if line_num % 100000 == 0:
                        print(f"处理到第 {line_num} 行...")
                    continue

        df = pd.DataFrame(sites)
        print(f"\n解析完成:")
        print(f"  - 有效位点: {len(df)}")
        print(f"  - 跳过行数: {skipped}")

        # 统计信息
        print(f"\n染色体分布:")
        chrom_counts = df['chrom'].value_counts()
        for chrom, count in chrom_counts.head(10).items():
            print(f"  {chrom}: {count:,}")

        print(f"\n区域分布:")
        region_counts = df['region'].value_counts()
        for region, count in region_counts.head(10).items():
            print(f"  {region}: {count:,}")

        return df

    def save_m6a_bed(self, df: pd.DataFrame, output_file: str):
        """保存为BED格式"""
        print(f"\n保存m6A数据到: {output_file}")

        # BED6格式: chrom start end name score strand
        bed_df = df.copy()
        bed_df['name'] = bed_df['gene_id'] + '_' + bed_df['region']

        bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            output_file, sep='\t', header=False, index=False
        )

        print(f"已保存 {len(bed_df)} 个位点")

    def analyze_psi_gse102476(self, gz_file: str):
        """
        分析GSE102476数据（仅用于记录）
        这个数据集只包含rRNA (18S, 28S)，不包含mRNA
        """
        print(f"\n=== 分析GSE102476 Ψ数据 ===")
        print(f"输入文件: {gz_file}")

        total_sites = 0
        psi_sites = 0
        high_deletion = 0

        with gzip.open(gz_file, 'rt') as f:
            header = f.readline().strip().split('\t')
            print(f"列名: {header}")

            is_psi_idx = header.index('Is_psiU')
            delet_idx = header.index('Delet')
            reads_idx = header.index('Reads')
            locus_idx = header.index('Locus')

            for line in f:
                total_sites += 1
                cols = line.strip().split('\t')

                try:
                    is_psi = cols[is_psi_idx]
                    delet = float(cols[delet_idx])

                    # Reads可能包含"No"或其他非数字值
                    try:
                        reads = int(cols[reads_idx])
                    except (ValueError, TypeError):
                        # 跳过无效Reads值
                        continue

                    locus = cols[locus_idx]

                    if is_psi == 'Yes':
                        psi_sites += 1

                    # 高deletion ratio（deletions/reads）
                    if reads > 0 and delet / reads > 0.1:
                        high_deletion += 1
                except (ValueError, IndexError) as e:
                    continue

        print(f"\n统计结果:")
        print(f"  - 总位点数: {total_sites:,}")
        print(f"  - 确认Ψ位点 (Is_psiU=Yes): {psi_sites}")
        print(f"  - 高deletion位点 (delet/reads>0.1): {high_deletion}")
        print(f"\n⚠️  数据来源: 仅包含rRNA (18S, 28S)，不包含mRNA")
        print(f"⚠️  建议: 使用高质量模拟数据用于Ψ分析")

        return {
            'total': total_sites,
            'psi_yes': psi_sites,
            'high_deletion': high_deletion
        }


def main():
    """主函数"""
    project_dir = Path("/home/tony/m6A/RNA_modification_analysis")
    parser = RealDataParser(project_dir)

    # 1. 解析m6A数据
    m6a_gz = project_dir / "m6A=sites=cell=human=hg38=HEK293T.txt.gz"
    if m6a_gz.exists():
        print("\n" + "="*70)
        print("解析m6A真实数据")
        print("="*70)

        m6a_df = parser.parse_m6a_repic(
            str(m6a_gz),
            min_reads=10,
            fdr_threshold=0.05
        )

        # 保存BED文件
        m6a_bed = parser.data_dir / "m6A_HEK293T_repic.bed"
        parser.save_m6a_bed(m6a_df, str(m6a_bed))

        # 保存详细信息
        m6a_csv = parser.results_dir / "m6a_repic_full.csv"
        m6a_df.to_csv(m6a_csv, index=False)
        print(f"\n详细信息已保存: {m6a_csv}")

    # 2. 分析Ψ数据（仅统计）
    psi_gz = project_dir / "GSE102476_Stat_CMC _deletion.txt.gz"
    if psi_gz.exists():
        print("\n" + "="*70)
        print("分析Ψ数据（GSE102476）")
        print("="*70)

        stats = parser.analyze_psi_gse102476(str(psi_gz))

        print("\n" + "="*70)
        print("数据来源说明:")
        print("="*70)
        print("✅ m6A: 使用真实数据 (REPIC, 725,587位点, HEK293T)")
        print("❌ Ψ: GSE102476仅包含rRNA数据，不适合mRNA分析")
        print("   建议: 使用先前生成的高质量模拟数据")
        print("="*70)


if __name__ == "__main__":
    main()
