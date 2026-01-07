#!/usr/bin/env Rscript
# ============================================================================
# RNA修饰比较分析 - R可视化脚本
# ============================================================================
#
# 功能：使用R语言的ggplot2、ChIPseeker等包进行高级可视化
# 主要图表：
#   1. Metagene Profile（使用ChIPseeker）
#   2. 基因组特征分布（饼图、柱状图）
#   3. 韦恩图（使用VennDiagram包）
#   4. 综合对比图（多面板）
#
# 作者：周子航
# 日期：2026-01-07
# ============================================================================

# 加载必要的包
required_packages <- c("ggplot2", "ChIPseeker", "VennDiagram",
                      "gridExtra", "RColorBrewer", "scales",
                      "dplyr", "readr", "grid")

# 安装缺失的包
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("ChIPseeker")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
  }
}

# 设置全局主题
theme_set(theme_bw() +
          theme(
            text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(face = "bold"),
            legend.position = "right",
            panel.grid.minor = element_blank()
          ))

# ============================================================================
# 1. 数据加载
# ============================================================================

cat("========================================\n")
cat("数据加载\n")
cat("========================================\n\n")

# 定义路径
PROJECT_DIR <- normalizePath("../")
DATA_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIGURES_DIR <- file.path(PROJECT_DIR, "figures")

# 加载Python脚本生成的数据
psi_data <- read_csv(file.path(RESULTS_DIR, "psi_sites_annotated.csv"),
                     show_col_types = FALSE)
m6a_data <- read_csv(file.path(RESULTS_DIR, "m6a_sites_annotated.csv"),
                     show_col_types = FALSE)

cat(sprintf("Ψ位点数: %d\n", nrow(psi_data)))
cat(sprintf("m6A位点数: %d\n", nrow(m6a_data)))
cat("\n✓ 数据加载完成\n")

# ============================================================================
# 2. Metagene Profile分析（ChIPseeker风格）
# ============================================================================

cat("\n========================================\n")
cat("Metagene Profile分析\n")
cat("========================================\n\n")

# 准备数据：计算相对位置的密度
calculate_metagene_density <- function(data, n_bins = 100) {
  # 分bin
  data$bin <- cut(data$relative_pos,
                  breaks = seq(0, 1, length.out = n_bins + 1),
                  labels = FALSE)

  # 计算密度
  density_table <- table(data$bin) / sum(table(data$bin))

  # 转换为数据框
  density_df <- data.frame(
    bin = as.numeric(names(density_table)),
    density = as.numeric(density_table)
  )

  # 添加平滑（移动平均）
  density_df$density_smooth <- stats::filter(density_df$density,
                                            rep(1/9, 9), sides = 2)
  density_df$density_smooth[is.na(density_df$density_smooth)] <- density_df$density[is.na(density_df$density_smooth)]

  # 归一化
  density_df$density_smooth <- density_df$density_smooth / sum(density_df$density_smooth, na.rm = TRUE)

  # 添加位置百分比
  density_df$position <- density_df$bin / n_bins * 100

  return(density_df)
}

# 计算density
psi_density <- calculate_metagene_density(psi_data, n_bins = 100)
m6a_density <- calculate_metagene_density(m6a_data, n_bins = 100)

# 绘图
p_metagene <- ggplot() +
  # 添加区域背景
  geom_rect(aes(xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.1) +
  geom_rect(aes(xmin = 20, xmax = 80, ymin = -Inf, ymax = Inf),
            fill = "green", alpha = 0.1) +
  geom_rect(aes(xmin = 80, xmax = 100, ymin = -Inf, ymax = Inf),
            fill = "orange", alpha = 0.1) +

  # 绘制曲线
  geom_line(data = psi_density, aes(x = position, y = density_smooth,
                                     color = "Pseudouridine (Ψ)"),
            size = 1.2) +
  geom_line(data = m6a_density, aes(x = position, y = density_smooth,
                                     color = "m6A"),
            size = 1.2, linetype = "solid") +

  # 添加起始和终止密码子线
  geom_vline(xintercept = 20, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = 80, linetype = "dashed", color = "black", alpha = 0.5) +

  # 添加文本标注
  annotate("text", x = 10, y = max(c(psi_density$density_smooth, m6a_density$density_smooth)) * 0.95,
           label = "5'UTR", size = 4, color = "gray40") +
  annotate("text", x = 50, y = max(c(psi_density$density_smooth, m6a_density$density_smooth)) * 0.95,
           label = "CDS", size = 4, color = "darkgreen") +
  annotate("text", x = 90, y = max(c(psi_density$density_smooth, m6a_density$density_smooth)) * 0.95,
           label = "3'UTR", size = 4, color = "darkorange") +
  annotate("text", x = 82, y = max(c(psi_density$density_smooth, m6a_density$density_smooth)) * 0.90,
           label = "Stop", size = 3, angle = 90, color = "black") +
  annotate("text", x = 22, y = max(c(psi_density$density_smooth, m6a_density$density_smooth)) * 0.90,
           label = "Start", size = 3, angle = 90, color = "black") +

  # 颜色和标签
  scale_color_manual(values = c("Pseudouridine (Ψ)" = "#3498db",
                                "m6A" = "#e74c3c")) +
  labs(x = "Relative position on transcript (5' → 3')",
       y = "Normalized density",
       color = "Modification",
       title = "Metagene Profile: Pseudouridine (Ψ) vs m6A") +

  xlim(0, 100) +
  theme_light()

# 保存
ggsave(file.path(FIGURES_DIR, "R_metagene_profile.pdf"),
       p_metagene, width = 10, height = 5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "R_metagene_profile.png"),
       p_metagene, width = 10, height = 5, dpi = 300)

cat("✓ Metagene Profile图已保存\n")

# ============================================================================
# 3. 基因组特征分布图
# ============================================================================

cat("\n========================================\n")
cat("基因组特征分布分析\n")
cat("========================================\n\n")

# 准备数据
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

# 合并数据
feature_data <- bind_rows(psi_features, m6a_features)
feature_data$feature <- factor(feature_data$feature,
                              levels = c("5UTR", "CDS", "3UTR"),
                              labels = c("5'UTR", "CDS", "3'UTR"))

# 图1：分组柱状图（百分比）
p_bar_percent <- ggplot(feature_data, aes(x = feature, y = percent, fill = modification)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", size = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Pseudouridine (Ψ)" = "#3498db",
                               "m6A" = "#e74c3c")) +
  labs(x = "Genomic Feature",
       y = "Percentage (%)",
       fill = "Modification",
       title = "Distribution across Genomic Features") +
  ylim(0, max(feature_data$percent) * 1.2) +
  theme_minimal()

# 图2：堆叠柱状图（绝对值）
feature_data_total <- feature_data %>%
  group_by(modification) %>%
  mutate(total = sum(count)) %>%
  ungroup()

p_bar_count <- ggplot(feature_data, aes(x = feature, y = count, fill = modification)) +
  geom_bar(stat = "identity", position = "stack",
           color = "black", size = 0.3) +
  geom_text(aes(label = comma(count)),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Pseudouridine (Ψ)" = "#3498db",
                               "m6A" = "#e74c3c")) +
  labs(x = "Genomic Feature",
       y = "Number of Sites",
       fill = "Modification",
       title = "Stacked Distribution (Absolute Counts)") +
  scale_y_continuous(labels = comma_format()) +
  theme_minimal()

# 保存
ggsave(file.path(FIGURES_DIR, "R_feature_distribution.pdf"),
       p_bar_percent, width = 8, height = 5, dpi = 300)
ggsave(file.path(FIGURES_DIR, "R_feature_distribution_stacked.pdf"),
       p_bar_count, width = 8, height = 5, dpi = 300)

cat("✓ 特征分布图已保存\n")

# ============================================================================
# 4. 韦恩图（使用VennDiagram包）
# ============================================================================

cat("\n========================================\n")
cat("韦恩图分析\n")
cat("========================================\n\n")

# 计算韦恩图数据
# 简化：假设共定位窗口50bp（与Python脚本一致）
# 这里我们使用Python脚本的结果
if (file.exists(file.path(RESULTS_DIR, "colocalized_sites.csv"))) {
  colocal_data <- read_csv(file.path(RESULTS_DIR, "colocalized_sites.csv"),
                           show_col_types = FALSE)

  psi_unique <- nrow(psi_data) - length(unique(colocal_data$psi_id))
  m6a_unique <- nrow(m6a_data) - length(unique(colocal_data$m6a_id))
  overlap <- length(unique(colocal_data$psi_id))
} else {
  # 如果没有共定位数据，使用估算
  psi_unique <- nrow(psi_data) * 0.95
  m6a_unique <- nrow(m6a_data) * 0.95
  overlap <- min(nrow(psi_data), nrow(m6a_data)) * 0.05
}

# 绘制韦恩图
pdf(file.path(FIGURES_DIR, "R_venn_diagram.pdf"), width = 8, height = 8)
venn.diagram(
  x = list(psi_data$name, m6a_data$name),
  category.names = c("Pseudouridine (Ψ)", "m6A"),
  filename = NULL,  # 不直接保存文件
  output = TRUE,
  scaled = TRUE,
  lwd = 3,
  col = c("#3498db", "#e74c3c"),
  fill = c("#3498db", "#e74c3c"),
  alpha = 0.5,
  cat.cex = 2,
  cat.fontface = "bold",
  margin = 0.1,
  cex = 1.5,
  fontface = "bold",
  main = "Colocalization: Ψ vs m6A",
  sub = sprintf("Total Ψ: %s | Total m6A: %s",
                comma(nrow(psi_data)), comma(nrow(m6a_data)))
)
dev.off()

cat("✓ 韦恩图已保存\n")

# ============================================================================
# 5. 综合多面板图（论文级）
# ============================================================================

cat("\n========================================\n")
cat("生成综合多面板图\n")
cat("========================================\n\n")

# 准备4个子图的数据
# 子图1：Metagene Profile
p1 <- p_metagene +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# 子图2：特征分布柱状图
p2 <- p_bar_percent +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

# 子图3：饼图对比
# 准备饼图数据
create_pie_data <- function(data, mod_name) {
  data %>%
    group_by(feature) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percent = count / sum(count) * 100,
           modification = mod_name,
           label = sprintf("%s\n%.1f%%", feature, percent))
}

psi_pie <- create_pie_data(psi_data, "Ψ")
m6a_pie <- create_pie_data(m6a_data, "m6A")

# 绘制饼图（使用极坐标barplot）
p3 <- ggplot(psi_pie, aes(x = "", y = percent, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Pseudouridine (Ψ)", fill = "Feature") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right")

p4 <- ggplot(m6a_pie, aes(x = "", y = percent, fill = feature)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "m6A", fill = "Feature") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right")

# 组合图形
p_combined <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                           top = "RNA Modification Comparison: Ψ vs m6A")

# 保存
ggsave(file.path(FIGURES_DIR, "R_combined_figure.pdf"),
       p_combined, width = 12, height = 10, dpi = 300)
ggsave(file.path(FIGURES_DIR, "R_combined_figure.png"),
       p_combined, width = 12, height = 10, dpi = 300)

cat("✓ 综合多面板图已保存\n")

# ============================================================================
# 6. 生成统计报告
# ============================================================================

cat("\n========================================\n")
cat("生成统计报告\n")
cat("========================================\n\n")

# 计算统计量
psi_stats <- list(
  total = nrow(psi_data),
  utr5 = sum(psi_data$feature == "5UTR"),
  cds = sum(psi_data$feature == "CDS"),
  utr3 = sum(psi_data$feature == "3UTR"),
  mean_rel_pos = mean(psi_data$relative_pos),
  median_rel_pos = median(psi_data$relative_pos)
)

m6a_stats <- list(
  total = nrow(m6a_data),
  utr5 = sum(m6a_data$feature == "5UTR"),
  cds = sum(m6a_data$feature == "CDS"),
  utr3 = sum(m6a_data$feature == "3UTR"),
  mean_rel_pos = mean(m6a_data$relative_pos),
  median_rel_pos = median(m6a_data$relative_pos)
)

# 打印报告
cat(sprintf("\n========== 统计摘要 ==========\n\n"))

cat(sprintf("Pseudouridine (Ψ):\n"))
cat(sprintf("  总位点数: %s\n", comma(psi_stats$total)))
cat(sprintf("  5'UTR: %s (%.1f%%)\n", comma(psi_stats$utr5),
            psi_stats$utr5/psi_stats$total*100))
cat(sprintf("  CDS: %s (%.1f%%)\n", comma(psi_stats$cds),
            psi_stats$cds/psi_stats$total*100))
cat(sprintf("  3'UTR: %s (%.1f%%)\n", comma(psi_stats$utr3),
            psi_stats$utr3/psi_stats$total*100))
cat(sprintf("  平均相对位置: %.3f\n", psi_stats$mean_rel_pos))
cat(sprintf("  中位数相对位置: %.3f\n", psi_stats$median_rel_pos))

cat(sprintf("\nm6A:\n"))
cat(sprintf("  总位点数: %s\n", comma(m6a_stats$total)))
cat(sprintf("  5'UTR: %s (%.1f%%)\n", comma(m6a_stats$utr5),
            m6a_stats$utr5/m6a_stats$total*100))
cat(sprintf("  CDS: %s (%.1f%%)\n", comma(m6a_stats$cds),
            m6a_stats$cds/m6a_stats$total*100))
cat(sprintf("  3'UTR: %s (%.1f%%)\n", comma(m6a_stats$utr3),
            m6a_stats$utr3/m6a_stats$total*100))
cat(sprintf("  平均相对位置: %.3f\n", m6a_stats$mean_rel_pos))
cat(sprintf("  中位数相对位置: %.3f\n", m6a_stats$median_rel_pos))

cat(sprintf("\n========================================\n"))
cat("分析完成！\n")
cat(sprintf("所有图表已保存到: %s\n", FIGURES_DIR))
cat("========================================\n")

# ============================================================================
# 结束
# ============================================================================
