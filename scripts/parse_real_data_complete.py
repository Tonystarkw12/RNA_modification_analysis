#!/usr/bin/env python3
"""
完整解析真实RNA修饰数据
- Ψ: RMBase v2.0 (hg19, protein-coding基因)
- m6A: REPIC (hg38, HEK293T)
- 处理hg19/hg38差异
"""

import pandas as pd
import numpy as np
import gzip
import re
from pathlib import Path
from collections import Counter

class RealDataParser:
    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"
        self.data_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)

    def parse_psi_rmbase(self, gz_file: str,
                         min_support: int = 1,
                         max_sites: int = 10000) -> pd.DataFrame:
        """
        解析RMBase Ψ数据

        格式：BED15+
        chromosome, modStart, modEnd, modId, score, strand, modName,
        modType, supportNum, supportList, pubmedIds, geneName, geneType, region, sequence
        """
        print(f"\n=== 解析RMBase Ψ数据 ===")
        print(f"输入文件: {gz_file}")

        sites = []
        total_lines = 0
        skipped_non_coding = 0

        with gzip.open(gz_file, 'rt') as f:
            # 跳过注释行
            for line in f:
                if line.startswith('##'):
                    continue
                elif line.startswith('#chromosome'):
                    # 表头行
                    headers = line.strip().split('\t')
                    print(f"列名: {headers}")
                    continue
                else:
                    break

            # 解析数据行
            for line in f:
                total_lines += 1
                cols = line.strip().split('\t')

                if len(cols) < 14:
                    continue

                try:
                    chrom = cols[0]
                    start = int(cols[1])
                    end = int(cols[2])
                    site_id = cols[3]
                    score = float(cols[4]) if cols[4] != '.' else 0
                    strand = cols[5]
                    gene_name = cols[11]
                    gene_type = cols[12]
                    region = cols[13]

                    # 只保留protein-coding基因
                    if gene_type != 'protein_coding':
                        skipped_non_coding += 1
                        continue

                    # 过滤支持数
                    support_num = int(cols[8]) if cols[8].isdigit() else 1
                    if support_num < min_support:
                        continue

                    sites.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'site_id': site_id,
                        'score': score,
                        'strand': strand,
                        'gene_name': gene_name,
                        'gene_type': gene_type,
                        'region': region,
                        'support_num': support_num
                    })

                except Exception as e:
                    if total_lines % 1000 == 0:
                        print(f"处理到第 {total_lines} 行...")
                    continue

        df = pd.DataFrame(sites)
        print(f"\n解析完成:")
        print(f"  - 总行数: {total_lines:,}")
        print(f"  - protein-coding位点: {len(df):,}")
        print(f"  - 跳过非编码位点: {skipped_non_coding:,}")

        # 随机采样（如果位点数太多）
        if len(df) > max_sites:
            print(f"\n位点数过多，进行随机采样:")
            print(f"  - 当前: {len(df):,}")
            print(f"  - 目标: {max_sites:,}")
            df = df.sample(n=max_sites, random_state=42)
            print(f"  - 采样后: {len(df):,}")

        # 统计信息
        print(f"\n染色体分布:")
        for chrom, count in df['chrom'].value_counts().head(10).items():
            print(f"  {chrom}: {count:,}")

        print(f"\n区域分布:")
        for region, count in df['region'].value_counts().head(10).items():
            print(f"  {region}: {count:,}")

        return df

    def parse_m6a_repic(self, gz_file: str,
                       fdr_threshold: float = 0.01,
                       fold_enrichment_threshold: float = 10.0,
                       max_sites: int = 10000) -> pd.DataFrame:
        """
        解析REPIC m6A数据

        格式：
        pos,pvalue,fdr,fold_enrichment,region,dataset,sample,tool,cell/tissue,species,assembly,geneid
        chr1:629848-630037[+]	1.00e-268	1.00e-265	29.6	exon	...
        """
        print(f"\n=== 解析REPIC m6A数据 ===")
        print(f"输入文件: {gz_file}")

        sites = []
        skipped_fdr = 0
        skipped_enrichment = 0

        with gzip.open(gz_file, 'rt') as f:
            header = f.readline().strip().split('\t')

            pos_idx = header.index('pos')
            pvalue_idx = header.index('pvalue')
            fdr_idx = header.index('fdr')
            fold_idx = header.index('fold_enrichment')
            region_idx = header.index('region')
            gene_idx = header.index('geneid')

            line_num = 1
            for line in f:
                line_num += 1
                try:
                    cols = line.strip().split('\t')
                    if len(cols) < len(header):
                        continue

                    # 解析位置: chr1:629848-630037[+]
                    pos_str = cols[pos_idx]
                    match = re.match(r'(chr\w+):(\d+)-(\d+)\[(\+|\-)\]', pos_str)
                    if not match:
                        continue

                    chrom = match.group(1)
                    start = int(match.group(2))
                    end = int(match.group(3))
                    strand = match.group(4)
                    peak_pos = (start + end) // 2

                    pvalue = float(cols[pvalue_idx])
                    fdr = float(cols[fdr_idx])
                    fold_enrichment = float(cols[fold_idx])
                    region = cols[region_idx]
                    gene_id = cols[gene_idx]

                    # 质量过滤
                    if fdr > fdr_threshold:
                        skipped_fdr += 1
                        continue

                    if fold_enrichment < fold_enrichment_threshold:
                        skipped_enrichment += 1
                        continue

                    sites.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'peak_pos': peak_pos,
                        'strand': strand,
                        'pvalue': pvalue,
                        'fdr': fdr,
                        'fold_enrichment': fold_enrichment,
                        'region': region,
                        'gene_id': gene_id,
                        'score': fold_enrichment
                    })

                except Exception as e:
                    if line_num % 100000 == 0:
                        print(f"处理到第 {line_num} 行...")
                    continue

        df = pd.DataFrame(sites)
        print(f"\n解析完成:")
        print(f"  - 有效位点: {len(df):,}")
        print(f"  - FDR过滤: {skipped_fdr:,}")
        print(f"  - Fold enrichment过滤: {skipped_enrichment:,}")

        # 随机采样
        if len(df) > max_sites:
            print(f"\n位点数过多，进行随机采样:")
            print(f"  - 当前: {len(df):,}")
            print(f"  - 目标: {max_sites:,}")
            df = df.sample(n=max_sites, random_state=42)
            print(f"  - 采样后: {len(df):,}")

        # 统计信息
        print(f"\n染色体分布:")
        for chrom, count in df['chrom'].value_counts().head(10).items():
            print(f"  {chrom}: {count:,}")

        print(f"\n区域分布:")
        for region, count in df['region'].value_counts().head(10).items():
            print(f"  {region}: {count:,}")

        return df

    def normalize_chrom_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """标准化染色体名称"""
        df = df.copy()
        df['chrom'] = df['chrom'].str.replace('chr', '')
        return df

    def save_bed_format(self, df: pd.DataFrame, output_file: str,
                       mod_type: str):
        """
        保存为BED6格式
        """
        print(f"\n保存{mod_type}数据到: {output_file}")

        bed_df = df.copy()

        # 创建name列
        if 'gene_name' in bed_df.columns:
            bed_df['name'] = bed_df['gene_name'] + '_' + bed_df['region']
        elif 'gene_id' in bed_df.columns:
            bed_df['name'] = bed_df['gene_id'] + '_' + bed_df['region']
        else:
            bed_df['name'] = bed_df.get('site_id', [f'{i}' for i in range(len(bed_df))])

        # BED6格式: chrom start end name score strand
        bed_df[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            output_file, sep='\t', header=False, index=False
        )

        print(f"已保存 {len(bed_df)} 个位点")


def main():
    """主函数"""
    project_dir = Path("/home/tony/m6A/RNA_modification_analysis")
    parser = RealDataParser(project_dir)

    print("="*70)
    print("解析真实RNA修饰数据")
    print("="*70)

    # 1. 解析Ψ数据 (RMBase, hg19)
    psi_gz = project_dir / "RMBase_hg19_all_PseudoU_site.txt.gz"
    if psi_gz.exists():
        psi_df = parser.parse_psi_rmbase(
            str(psi_gz),
            min_support=1,
            max_sites=999999  # 获取所有位点
        )

        # 保存
        psi_bed = parser.data_dir / "psi_real_hg19.bed"
        parser.save_bed_format(psi_df, str(psi_bed), "Ψ")

        psi_csv = parser.results_dir / "psi_real_full.csv"
        psi_df.to_csv(psi_csv, index=False)
        print(f"详细信息: {psi_csv}")

    # 2. 解析m6A数据 (REPIC, hg38)
    m6a_gz = project_dir / "m6A=sites=cell=human=hg38=HEK293T.txt.gz"
    if m6a_gz.exists():
        m6a_df = parser.parse_m6a_repic(
            str(m6a_gz),
            fdr_threshold=0.01,
            fold_enrichment_threshold=10.0,
            max_sites=10000
        )

        # 保存
        m6a_bed = parser.data_dir / "m6A_real_hg38.bed"
        parser.save_bed_format(m6a_df, str(m6a_bed), "m6A")

        m6a_csv = parser.results_dir / "m6a_real_full.csv"
        m6a_df.to_csv(m6a_csv, index=False)
        print(f"详细信息: {m6a_csv}")

    print("\n" + "="*70)
    print("数据解析完成！")
    print("="*70)
    print("\n注意:")
    print("  - Ψ数据: hg19基因组版本")
    print("  - m6A数据: hg38基因组版本")
    print("  - 后续分析将使用共同的染色体进行比较")
    print("="*70)


if __name__ == "__main__":
    main()
