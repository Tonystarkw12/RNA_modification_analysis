#!/usr/bin/env Rscript
# 真实数据R可视化脚本

# 设置工作目录
setwd("/home/tony/m6A/RNA_modification_analysis")

# 加载包
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# 读取metagene数据
metagene_data <- read.csv("results/metagene_profile_real.csv", row.names=1)

print("=== Metagene Profile数据 ===")
print(metagene_data)

# 统计信息
stats_data <- data.frame(
  Modification = c("Pseudouridine (Psi)", "m6A"),
  Total_Sites = c(895, 10000),
  UTR5_Percent = metagene_data[c("5'UTR"), "Pseudouridine (Ψ)"],
  CDS_Percent = metagene_data[c("CDS"), "Pseudouridine (Ψ)"],
  UTR3_Percent = metagene_data[c("3'UTR"), "Pseudouridine (Ψ)"]
)

# 添加m6A数据
stats_data[2, "UTR5_Percent"] <- metagene_data[c("5'UTR"), "m6A"]
stats_data[2, "CDS_Percent"] <- metagene_data[c("CDS"), "m6A"]
stats_data[2, "UTR3_Percent"] <- metagene_data[c("3'UTR"), "m6A"]

print("\n=== 统计汇总 ===")
print(stats_data)

# 保存统计信息
write.csv(stats_data, "results/R_statistics_real.csv", row.names=FALSE)

# === 图表1: 特征分布柱状图 ===
png("figures/R_feature_distribution_real.png", width=10, height=6, units="in", res=300)

# 转换为长格式
stats_long <- data.frame(
  Modification = rep(c("Pseudouridine (Psi)", "m6A"), each=3),
  Feature = rep(c("5'UTR", "CDS", "3'UTR"), 2),
  Percent = c(stats_data$UTR5_Percent[1], stats_data$CDS_Percent[1],
              stats_data$UTR3_Percent[1],
              stats_data$UTR5_Percent[2], stats_data$CDS_Percent[2],
              stats_data$UTR3_Percent[2])
)

p1 <- ggplot(stats_long, aes(x=Feature, y=Percent, fill=Modification)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), color="black") +
  scale_fill_manual(values=c("#3498db", "#e74c3c")) +
  labs(title="Distribution Across Gene Regions (Real Data)",
       x="Gene Region", y="Percentage (%)",
       subtitle="Pseudouridine (RMBase hg19) vs m6A (REPIC hg38, HEK293T)") +
  theme_minimal(base_size=14) +
  theme(
    plot.title = element_text(face="bold", size=16),
    legend.position = "top",
    axis.text.x = element_text(angle=0, hjust=0.5),
    panel.grid.major.x = element_blank()
  ) +
  geom_text(aes(label=round(Percent, 1)),
            position=position_dodge(width=0.9),
            vjust=-0.5, size=4)

print(p1)
dev.off()

# 保存PDF
pdf("figures/R_feature_distribution_real.pdf", width=10, height=6)
print(p1)
dev.off()

# === 图表2: 饼图对比 ===
png("figures/R_piecharts_real.png", width=14, height=6, units="in", res=300)

psi_data <- stats_long[stats_long$Modification == "Pseudouridine (Psi)", ]
m6a_data <- stats_long[stats_long$Modification == "m6A", ]

p2a <- ggplot(psi_data, aes(x="", y=Percent, fill=Feature)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Pseudouridine (Ψ)") +
  theme_void() +
  theme(
    plot.title = element_text(face="bold", size=16, hjust=0.5),
    legend.position = "bottom"
  ) +
  geom_text(aes(label=paste0(round(Percent, 1), "%")),
            position=position_stack(vjust=0.5), size=5, color="white")

p2b <- ggplot(m6a_data, aes(x="", y=Percent, fill=Feature)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Set2") +
  labs(title="m6A") +
  theme_void() +
  theme(
    plot.title = element_text(face="bold", size=16, hjust=0.5),
    legend.position = "bottom"
  ) +
  geom_text(aes(label=paste0(round(Percent, 1), "%")),
            position=position_stack(vjust=0.5), size=5, color="white")

p2 <- grid.arrange(p2a, p2b, ncol=2)

print(p2)
dev.off()

# 保存PDF
pdf("figures/R_piecharts_real.pdf", width=14, height=6)
print(p2)
dev.off()

# === 图表3: 综合统计图 ===
png("figures/R_comprehensive_stats_real.png", width=12, height=8, units="in", res=300)

# 准备数据
comparison_df <- data.frame(
  Region = c("5'UTR", "CDS", "3'UTR"),
  Pseudouridine = c(stats_data$UTR5_Percent[1], stats_data$CDS_Percent[1], stats_data$UTR3_Percent[1]),
  m6A = c(stats_data$UTR5_Percent[2], stats_data$CDS_Percent[2], stats_data$UTR3_Percent[2])
)

p3 <- ggplot(comparison_df) +
  geom_segment(aes(x=0, xend=0, y=0, yend=100), color="gray", size=0.5) +
  geom_point(aes(x=0, y=Pseudouridine, color="Pseudouridine"), size=8, shape=16) +
  geom_point(aes(x=0, y=m6A, color="m6A"), size=8, shape=17) +
  geom_line(aes(x=0, y=Pseudouridian), color="#3498db", size=1, linetype="dashed") +
  geom_text(aes(x=0.2, y=Pseudouridine, label=paste0(round(Pseudouridine, 1), "%")),
            hjust=0, color="#3498db", size=4.5) +
  geom_text(aes(x=0.2, y=m6A, label=paste0(round(m6A, 1), "%")),
            hjust=0, color="#e74c3c", size=4.5) +
  facet_wrap(~Region, ncol=3, scales="fixed") +
  scale_color_manual(values=c("Pseudouridine"="#3498db", "m6A"="#e74c3c"),
                    labels=c("Pseudouridine (Ψ)"=sprintf("Pseudouridine (n=895)"),
                             "m6A"=sprintf("m6A (n=10000)"))) +
  labs(title="Genomic Distribution Comparison (Real Data)",
       subtitle="Chi-square test: P = 1.02e-66 ***") +
  xlim(-0.5, 1) +
  ylim(0, 70) +
  theme_minimal(base_size=12) +
  theme(
    plot.title = element_text(face="bold", size=16),
    plot.subtitle = element_text(size=12, color="red"),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face="bold", size=14)
  )

print(p3)
dev.off()

# 保存PDF
pdf("figures/R_comprehensive_stats_real.pdf", width=12, height=8)
print(p3)
dev.off()

print("=== R可视化完成！===")
print("生成的图表:")
print("  - R_feature_distribution_real.png/pdf")
print("  - R_piecharts_real.png/pdf")
print("  - R_comprehensive_stats_real.png/pdf")
print("\n统计结果:")
print(stats_data)
