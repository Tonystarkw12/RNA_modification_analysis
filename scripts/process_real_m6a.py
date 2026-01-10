#!/usr/bin/env python3
"""
处理和过滤真实的m6A数据
- 质量过滤（FDR, fold enrichment）
- 随机采样到合理数量（5000-10000）
- 准备用于后续分析
"""

import pandas as pd
import numpy as np
from pathlib import Path

class M6ADataProcessor:
    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"

    def load_and_filter(self,
                       bed_file: str,
                       fdr_threshold: float = 0.01,
                       fold_enrichment_threshold: float = 10.0,
                       max_sites: int = 10000,
                       seed: int = 42) -> pd.DataFrame:
        """
        加载、过滤和采样m6A数据

        参数：
        - fdr_threshold: FDR阈值（默认0.01）
        - fold_enrichment_threshold: fold enrichment阈值（默认10）
        - max_sites: 最大位点数（默认10000）
        - seed: 随机种子（用于可重复性）
        """
        print(f"\n=== 处理真实m6A数据 ===")
        print(f"输入文件: {bed_file}")

        # 读取BED文件
        print(f"\n读取数据...")
        df = pd.read_csv(bed_file, sep='\t', header=None,
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

        print(f"原始位点数: {len(df):,}")

        # 读取详细信息（如果有）
        full_csv = self.results_dir / "m6a_repic_full.csv"
        if full_csv.exists():
            full_df = pd.read_csv(full_csv)
            print(f"详细信息文件: {full_csv}")

            # 合并详细信息
            df = df.merge(full_df[['chrom', 'start', 'end', 'fdr', 'region', 'gene_id']],
                         on=['chrom', 'start', 'end'], how='left')

            # 质量过滤
            print(f"\n应用质量过滤:")
            print(f"  原始: {len(df):,} 位点")

            # FDR过滤
            df_filtered = df[df['fdr'] < fdr_threshold].copy()
            print(f"  FDR < {fdr_threshold}: {len(df_filtered):,} 位点")

            # Fold enrichment过滤（如果有的话）
            if 'fold_enrichment' in df.columns:
                df_filtered = df_filtered[df_filtered['fold_enrichment'] > fold_enrichment_threshold].copy()
                print(f"  Fold enrichment > {fold_enrichment_threshold}: {len(df_filtered):,} 位点")

            df = df_filtered

        # 随机采样（如果位点数太多）
        if len(df) > max_sites:
            print(f"\n位点数过多，进行随机采样:")
            print(f"  当前: {len(df):,} 位点")
            print(f"  目标: {max_sites:,} 位点")

            df = df.sample(n=max_sites, random_state=seed)
            print(f"  采样后: {len(df):,} 位点")

        # 统计信息
        print(f"\n最终数据集统计:")
        print(f"  位点数: {len(df):,}")
        print(f"\n染色体分布:")
        for chrom, count in df['chrom'].value_counts().head(10).items():
            print(f"    {chrom}: {count:,}")

        if 'region' in df.columns:
            print(f"\n区域分布:")
            for region, count in df['region'].value_counts().head(10).items():
                print(f"    {region}: {count:,}")

        return df

    def save_processed_data(self, df: pd.DataFrame, output_file: str):
        """保存处理后的数据"""
        print(f"\n保存处理后的数据到: {output_file}")

        # 保存BED6格式
        bed_df = df[['chrom', 'start', 'end', 'name', 'score', 'strand']].copy()
        bed_df.to_csv(output_file, sep='\t', header=False, index=False)

        print(f"已保存 {len(bed_df)} 个位点")

    def prepare_for_analysis(self, df: pd.DataFrame):
        """
        准备用于分析的数据
        - 添加peak位置
        - 创建BED格式
        """
        # 计算peak位置（中点）
        df['peak_pos'] = (df['start'] + df['end']) // 2

        return df


def main():
    """主函数"""
    project_dir = Path("/home/tony/m6A/RNA_modification_analysis")
    processor = M6ADataProcessor(project_dir)

    # 读取并过滤真实m6A数据
    m6a_bed = processor.data_dir / "m6A_HEK293T_repic.bed"

    if m6a_bed.exists():
        # 加载并过滤
        m6a_df = processor.load_and_filter(
            str(m6a_bed),
            fdr_threshold=0.01,
            fold_enrichment_threshold=10.0,
            max_sites=10000,  # 采样到10,000个位点
            seed=42
        )

        # 准备分析数据
        m6a_df = processor.prepare_for_analysis(m6a_df)

        # 保存处理后的数据
        output_bed = processor.data_dir / "m6A_HEK293T_processed.bed"
        processor.save_processed_data(m6a_df, str(output_bed))

        # 保存详细信息
        output_csv = processor.results_dir / "m6a_processed.csv"
        m6a_df.to_csv(output_csv, index=False)
        print(f"\n详细信息已保存: {output_csv}")

        print("\n" + "="*70)
        print("m6A数据处理完成！")
        print("="*70)
        print(f"处理后的数据: {output_bed}")
        print(f"详细信息: {output_csv}")
        print(f"\n下一步: 运行分析pipeline")
        print("="*70)


if __name__ == "__main__":
    main()
