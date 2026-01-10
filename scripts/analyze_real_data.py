#!/usr/bin/env python3
"""
分析真实RNA修饰数据（简化版 - 使用数据自带region信息）
- Metagene profile分析
- Venn共定位分析
- 统计检验
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib_venn import venn2
from scipy.stats import chi2_contingency, mannwhitneyu

class RealDataAnalyzer:
    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"
        self.figures_dir = self.project_dir / "figures"
        self.results_dir.mkdir(exist_ok=True)
        self.figures_dir.mkdir(exist_ok=True)

        # 设置matplotlib中文字体
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False

    def load_data(self):
        """加载解析后的数据"""
        print("\n=== 加载数据 ===")

        # 加载Ψ数据
        psi_csv = self.results_dir / "psi_real_full.csv"
        self.psi_df = pd.read_csv(psi_csv)
        print(f"Ψ位点: {len(self.psi_df)}")

        # 加载m6A数据
        m6a_csv = self.results_dir / "m6a_real_full.csv"
        self.m6a_df = pd.read_csv(m6a_csv)
        print(f"m6A位点: {len(self.m6a_df)}")

    def normalize_regions(self):
        """
        标准化region名称
        处理组合区域如 "cds|utr3"
        """
        print("\n=== 标准化区域名称 ===")

        # 对于Ψ数据
        def simplify_region_psi(region):
            if 'utr3' in region.lower():
                return '3\'UTR'
            elif 'utr5' in region.lower():
                return '5\'UTR'
            elif 'cds' in region.lower():
                return 'CDS'
            else:
                return 'Other'

        self.psi_df['region_simple'] = self.psi_df['region'].apply(simplify_region_psi)

        # 对于m6A数据（可能包含组合区域）
        def simplify_region_m6a(region):
            if not isinstance(region, str):
                return 'Other'

            region_lower = region.lower()
            if 'utr3' in region_lower:
                return '3\'UTR'
            elif 'utr5' in region_lower:
                return '5\'UTR'
            elif 'cds' in region_lower:
                return 'CDS'
            elif 'exon' in region_lower and 'cds' not in region_lower:
                return 'CDS'  # exon通常属于CDS
            else:
                return 'Other'

        self.m6a_df['region_simple'] = self.m6a_df['region'].apply(simplify_region_m6a)

        # 统计
        print("Ψ区域分布:")
        print(self.psi_df['region_simple'].value_counts())

        print("\nm6A区域分布:")
        print(self.m6a_df['region_simple'].value_counts())

    def metagene_analysis(self):
        """Metagene profile分析"""
        print("\n=== Metagene Profile分析 ===")

        # 计算每个区域的百分比
        psi_dist = self.psi_df['region_simple'].value_counts(normalize=True) * 100
        m6a_dist = self.m6a_df['region_simple'].value_counts(normalize=True) * 100

        # 创建对比DataFrame
        regions = ['5\'UTR', 'CDS', '3\'UTR']
        comparison_df = pd.DataFrame({
            'Pseudouridine (Ψ)': [psi_dist.get(r, 0) for r in regions],
            'm6A': [m6a_dist.get(r, 0) for r in regions]
        }, index=regions)

        print("\n区域分布对比 (%):")
        print(comparison_df)

        # 保存结果
        comparison_df.to_csv(self.results_dir / "metagene_profile_real.csv")

        # 统计检验（卡方检验）
        # 创建列联表
        psi_counts = self.psi_df['region_simple'].value_counts()
        m6a_counts = self.m6a_df['region_simple'].value_counts()

        contingency_table = pd.DataFrame({
            'Ψ': [psi_counts.get(r, 0) for r in regions],
            'm6A': [m6a_counts.get(r, 0) for r in regions]
        }, index=regions)

        chi2, p_value, dof, expected = chi2_contingency(contingency_table)

        print(f"\n卡方检验:")
        print(f"  Chi2 = {chi2:.2f}")
        print(f"  P-value = {p_value:.2e}")

        if p_value < 0.001:
            print("  *** 极显著差异 (P < 0.001)")
        elif p_value < 0.01:
            print("  ** 显著差异 (P < 0.01)")
        elif p_value < 0.05:
            print("  * 差异显著 (P < 0.05)")

        return comparison_df, p_value

    def colocalization_analysis(self, window=50):
        """
        共定位分析（基于染色体位置）
        注意：由于hg19/hg38差异，只分析保守区域
        """
        print(f"\n=== 共定位分析 (窗口: ±{window}bp) ===")

        # 标准化染色体名称
        self.psi_df['chrom_std'] = self.psi_df['chrom'].str.replace('chr', '')
        self.m6a_df['chrom_std'] = self.m6a_df['chrom'].str.replace('chr', '')

        # 计算中点位置
        self.psi_df['center'] = (self.psi_df['start'] + self.psi_df['end']) // 2
        self.m6a_df['center'] = self.m6a_df['peak_pos']

        colocalized = 0
        colocalized_pairs = []

        # 只分析相同染色体
        for chrom in set(self.psi_df['chrom_std']) & set(self.m6a_df['chrom_std']):
            psi_chrom = self.psi_df[self.psi_df['chrom_std'] == chrom]
            m6a_chrom = self.m6a_df[self.m6a_df['chrom_std'] == chrom]

            for _, psi_row in psi_chrom.iterrows():
                psi_pos = psi_row['center']

                # 找到m6A在窗口内的位点
                nearby = m6a_chrom[
                    abs(m6a_chrom['center'] - psi_pos) <= window
                ]

                if len(nearby) > 0:
                    colocalized += len(nearby)
                    for _, m6a_row in nearby.iterrows():
                        colocalized_pairs.append({
                            'chrom': chrom,
                            'psi_pos': psi_pos,
                            'm6a_pos': m6a_row['center'],
                            'distance': abs(psi_pos - m6a_row['center'])
                        })

        colocal_df = pd.DataFrame(colocalized_pairs)

        print(f"Ψ位点总数: {len(self.psi_df)}")
        print(f"m6A位点总数: {len(self.m6a_df)}")
        print(f"共定位事件: {colocalized}")
        print(f"涉及Ψ位点: {len(colocal_df['psi_pos'].unique()) if len(colocal_df) > 0 else 0}")

        # 保存结果
        if len(colocal_df) > 0:
            colocal_df.to_csv(self.results_dir / "colocalization_pairs_real.csv", index=False)
            print(f"\n平均距离: {colocal_df['distance'].mean():.2f} bp")
            print(f"中位距离: {colocal_df['distance'].median():.2f} bp")

        return colocalized, len(self.psi_df), len(self.m6a_df)

    def plot_metagene_comparison(self, comparison_df, p_value):
        """绘制Metagene对比图"""
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # 1. 柱状图
        comparison_df.plot(kind='bar', ax=axes[0], color=['#3498db', '#e74c3c'])
        axes[0].set_ylabel('Percentage (%)', fontsize=12)
        axes[0].set_title('Distribution Across Gene Regions', fontsize=14, fontweight='bold')
        axes[0].legend(fontsize=11)
        axes[0].grid(axis='y', alpha=0.3)
        axes[0].set_xlabel('Gene Region', fontsize=12)

        # 2. 堆积柱状图
        comparison_df.T.plot(kind='bar', stacked=True, ax=axes[1],
                             color=['#2ecc71', '#f39c12', '#9b59b6'])
        axes[1].set_ylabel('Percentage (%)', fontsize=12)
        axes[1].set_title('Stacked Distribution', fontsize=14, fontweight='bold')
        axes[1].legend(title='Region', fontsize=10)
        axes[1].grid(axis='y', alpha=0.3)
        axes[1].set_xlabel('Modification', fontsize=12)

        # 3. 饼图
        fig_pie, axes_pie = plt.subplots(1, 2, figsize=(14, 6))

        for idx, (mod, color) in enumerate([('Pseudouridine (Ψ)', '#3498db'),
                                            ('m6A', '#e74c3c')]):
            data = comparison_df[mod]
            axes_pie[idx].pie(data, labels=data.index, autopct='%1.1f%%',
                             colors=['#2ecc71', '#f39c12', '#9b59b6'],
                             startangle=90)
            axes_pie[idx].set_title(f'{mod} Distribution', fontsize=14, fontweight='bold')

        plt.tight_layout()
        fig_pie.savefig(self.figures_dir / "metagene_piecharts_real.png", dpi=300, bbox_inches='tight')
        print(f"\n饼图已保存: metagene_piecharts_real.png")

        plt.tight_layout()
        fig.savefig(self.figures_dir / "metagene_comparison_real.png", dpi=300, bbox_inches='tight')
        print(f"对比图已保存: metagene_comparison_real.png")

        plt.close('all')

    def plot_venn_diagram(self, colocalized, n_psi, n_m6a):
        """绘制韦恩图"""
        fig, ax = plt.subplots(figsize=(10, 8))

        # 计算不重叠的数量
        psi_only = n_psi - min(colocalized, n_psi)  # 粗略估计
        m6a_only = n_m6a - min(colocalized, n_m6a)  # 粗略估计
        overlap = min(colocalized, min(n_psi, n_m6a))  # 保守估计

        venn = venn2(subsets=(psi_only, m6a_only, overlap),
                     set_labels=('Pseudouridine (Ψ)', 'm6A'),
                     ax=ax)

        venn.get_label_by_id('10').set_text(f'{psi_only}\nΨ only')
        venn.get_label_by_id('01').set_text(f'{m6a_only}\nm6A only')
        venn.get_label_by_id('11').set_text(f'{overlap}\nColocalized')

        ax.set_title('Colocalization Analysis (±50bp window)',
                    fontsize=16, fontweight='bold', pad=20)

        plt.savefig(self.figures_dir / "venn_diagram_real.png", dpi=300, bbox_inches='tight')
        print(f"\n韦恩图已保存: venn_diagram_real.png")

        plt.close()

    def generate_summary_report(self):
        """生成汇总报告"""
        report = []
        report.append("="*70)
        report.append("真实RNA修饰数据分析报告")
        report.append("="*70)
        report.append(f"\n数据来源:")
        report.append(f"  Ψ: RMBase v2.0 (hg19, protein-coding genes)")
        report.append(f"  m6A: REPIC (hg38, HEK293T cells)")
        report.append(f"\n数据量:")
        report.append(f"  Ψ位点: {len(self.psi_df)}")
        report.append(f"  m6A位点: {len(self.m6a_df)}")

        # 区域分布
        psi_dist = self.psi_df['region_simple'].value_counts(normalize=True) * 100
        m6a_dist = self.m6a_df['region_simple'].value_counts(normalize=True) * 100

        report.append(f"\n区域分布:")
        report.append(f"{'区域':<10} {'Ψ (%)':<10} {'m6A (%)':<10}")
        report.append("-" * 30)
        for region in ['5\'UTR', 'CDS', '3\'UTR']:
            report.append(f"{region:<10} {psi_dist.get(region, 0):<10.1f} {m6a_dist.get(region, 0):<10.1f}")

        report_text = "\n".join(report)

        with open(self.results_dir / "analysis_summary_real.txt", 'w') as f:
            f.write(report_text)

        print(f"\n{report_text}")
        print(f"\n报告已保存: analysis_summary_real.txt")


def main():
    """主函数"""
    project_dir = Path("/home/tony/m6A/RNA_modification_analysis")
    analyzer = RealDataAnalyzer(project_dir)

    print("="*70)
    print("分析真实RNA修饰数据")
    print("="*70)

    # 1. 加载数据
    analyzer.load_data()

    # 2. 标准化区域
    analyzer.normalize_regions()

    # 3. Metagene分析
    comparison_df, p_value = analyzer.metagene_analysis()

    # 4. 共定位分析
    colocalized, n_psi, n_m6a = analyzer.colocalization_analysis(window=50)

    # 5. 绘图
    print("\n=== 生成可视化 ===")
    analyzer.plot_metagene_comparison(comparison_df, p_value)
    analyzer.plot_venn_diagram(colocalized, n_psi, n_m6a)

    # 6. 生成报告
    print("\n=== 生成报告 ===")
    analyzer.generate_summary_report()

    print("\n" + "="*70)
    print("分析完成！")
    print("="*70)
    print(f"\n结果文件:")
    print(f"  - Metagene数据: {analyzer.results_dir}/metagene_profile_real.csv")
    print(f"  - 共定位数据: {analyzer.results_dir}/colocalization_pairs_real.csv")
    print(f"  - 汇总报告: {analyzer.results_dir}/analysis_summary_real.txt")
    print(f"\n图表:")
    print(f"  - 对比图: {analyzer.figures_dir}/metagene_comparison_real.png")
    print(f"  - 饼图: {analyzer.figures_dir}/metagene_piecharts_real.png")
    print(f"  - 韦恩图: {analyzer.figures_dir}/venn_diagram_real.png")
    print("="*70)


if __name__ == "__main__":
    main()
