#!/usr/bin/env Rscript
# 简化版R可视化 - 不依赖复杂包
# 只使用基础R和ggplot2

# 设置用户库路径
lib_path <- "~/R/library"
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

# 设置工作目录
setwd("/home/tony/m6A/RNA_modification_analysis")

cat("=== RNA修饰分析 - R可视化（简化版）===\n\n")

# 加载数据
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# 读取数据
psi_data <- read.csv("results/psi_sites_annotated.csv")
m6a_data <- read.csv("results/m6a_sites_annotated.csv")

cat(sprintf("数据加载:\n"))
cat(sprintf("  Ψ位点: %d\n", nrow(psi_data)))
cat(sprintf("  m6A位点: %d\n", nrow(m6a_data)))

# 准备特征分布数据
psi_features <- psi_data %>%
  group_by(feature) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percent = count / sum(count) * 100,
         modification = "Pseudouridine (Ψ)")

m6a_features <- m6a_data %>%
  group_by(feature) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percent = count / sum(count) * 100,
         modification = "m6A")

feature_data <- bind_rows(psi_features, m6a_features)
feature_data$feature <- factor(feature_data$feature,
                              levels = c("5UTR", "CDS", "3UTR"),
                              labels = c("5'UTR", "CDS", "3'UTR"))

# 图1: 特征分布柱状图
cat("\n[1/3] 绘制特征分布柱状图...\n")
p1 <- ggplot(feature_data, aes(x = feature, y = percent, fill = modification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Pseudouridine (Ψ)" = "#3498db",
                             "m6A" = "#e74c3c")) +
  labs(x = "Genomic Feature",
       y = "Percentage (%)",
       fill = "Modification",
       title = "Distribution across Genomic Features (R version)") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

ggsave("figures/R_feature_distribution.png", p1, width = 8, height = 5, dpi = 300)
ggsave("figures/R_feature_distribution.pdf", p1, width = 8, height = 5, dpi = 300)
cat("  保存: R_feature_distribution.png/pdf\n")

# 图2: 饼图对比
cat("\n[2/3] 绘制饼图对比...\n")

# Ψ饼图
psi_pie_data <- psi_features
psi_pie_data$feature <- factor(psi_pie_data$feature,
                               levels = c("5UTR", "CDS", "3UTR"),
                               labels = c("5'UTR", "CDS", "3'UTR"))

p2a <- ggplot(psi_pie_data, aes(x = "", y = percent, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Pseudouridine (Ψ)", fill = "Feature") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "right")

# m6A饼图
m6a_pie_data <- m6a_features
m6a_pie_data$feature <- factor(m6a_pie_data$feature,
                               levels = c("5UTR", "CDS", "3UTR"),
                               labels = c("5'UTR", "CDS", "3'UTR"))

p2b <- ggplot(m6a_pie_data, aes(x = "", y = percent, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "m6A", fill = "Feature") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "right")

# 组合饼图
library(gridExtra)
p2 <- grid.arrange(p2a, p2b, ncol = 2,
                   top = "Genomic Feature Distribution")

ggsave("figures/R_piecharts.png", p2, width = 10, height = 5, dpi = 300)
ggsave("figures/R_piecharts.pdf", p2, width = 10, height = 5, dpi = 300)
cat("  保存: R_piecharts.png/pdf\n")

# 图3: 综合统计图
cat("\n[3/3] 绘制综合统计图...\n")

# 计算统计摘要
stats_data <- data.frame(
  Modification = c("Pseudouridine (Ψ)", "m6A"),
  Total_Sites = c(nrow(psi_data), nrow(m6a_data)),
  UTR5_Percent = c(sum(psi_data$feature == "5UTR") / nrow(psi_data) * 100,
                   sum(m6a_data$feature == "5UTR") / nrow(m6a_data) * 100),
  CDS_Percent = c(sum(psi_data$feature == "CDS") / nrow(psi_data) * 100,
                 sum(m6a_data$feature == "CDS") / nrow(m6a_data) * 100),
  UTR3_Percent = c(sum(psi_data$feature == "3UTR") / nrow(psi_data) * 100,
                   sum(m6a_data$feature == "3UTR") / nrow(m6a_data) * 100),
  Mean_RelPos = c(mean(psi_data$relative_pos),
                 mean(m6a_data$relative_pos))
)

# 转换为长格式（手动方式，避免依赖tidyr）
stats_long <- data.frame(
  Modification = rep(c("Pseudouridine (Ψ)", "m6A"), each = 3),
  Feature = rep(c("5'UTR", "CDS", "3'UTR"), 2),
  Percent = c(stats_data$UTR5_Percent[1], stats_data$CDS_Percent[1], stats_data$UTR3_Percent[1],
             stats_data$UTR5_Percent[2], stats_data$CDS_Percent[2], stats_data$UTR3_Percent[2])
)

p3 <- ggplot(stats_long, aes(x = Modification, y = Percent, fill = Feature)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("5'UTR" = "#95a5a6",
                             "CDS" = "#2ecc71",
                             "3'UTR" = "#f39c12")) +
  labs(x = "Modification Type",
       y = "Percentage (%)",
       fill = "Genomic Feature",
       title = "Comprehensive Feature Distribution (R version)") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

ggsave("figures/R_comprehensive_stats.png", p3, width = 8, height = 6, dpi = 300)
ggsave("figures/R_comprehensive_stats.pdf", p3, width = 8, height = 6, dpi = 300)
cat("  保存: R_comprehensive_stats.png/pdf\n")

# 保存统计摘要
cat("\n保存统计摘要...\n")
write.csv(stats_data, "results/R_statistics_summary.csv", row.names = FALSE)

cat("\n=== R可视化完成！===\n")
cat("生成的图表:\n")
cat("  - R_feature_distribution.png/pdf\n")
cat("  - R_piecharts.png/pdf\n")
cat("  - R_comprehensive_stats.png/pdf\n")
cat("\n统计结果:\n")
print(stats_data)
