#!/usr/bin/env Rscript
# 真实数据R可视化脚本（使用基础R，无需额外包）

# 设置工作目录
setwd("/home/tony/m6A/RNA_modification_analysis")

# 读取metagene数据
metagene_data <- read.csv("results/metagene_profile_real.csv", row.names=1, check.names=FALSE)

print("=== Metagene Profile数据 ===")
print(metagene_data)

# 统计数据（按列索引访问，避免特殊字符问题）
psi_5utr <- metagene_data[1, 1]  # 5'UTR, Psi
psi_cds <- metagene_data[2, 1]   # CDS, Psi
psi_3utr <- metagene_data[3, 1]  # 3'UTR, Psi

m6a_5utr <- metagene_data[1, 2]  # 5'UTR, m6A
m6a_cds <- metagene_data[2, 2]   # CDS, m6A
m6a_3utr <- metagene_data[3, 2]  # 3'UTR, m6A

# === 图表1: 柱状图对比 ===
png("figures/R_feature_distribution_real.png", width=10, height=6, units="in", res=300)

par(mar=c(5, 5, 5, 2))

# 准备数据
regions <- c("5'UTR", "CDS", "3'UTR")
psi_values <- c(psi_5utr, psi_cds, psi_3utr)
m6a_values <- c(m6a_5utr, m6a_cds, m6a_3utr)

# 设置颜色
colors_psi <- "#3498db"
colors_m6a <- "#e74c3c"

# 绘制柱状图
x_positions <- barplot(rbind(psi_values, m6a_values),
                      beside=TRUE,
                      col=c(colors_psi, colors_m6a),
                      names.arg=regions,
                      ylim=c(0, 70),
                      las=1,
                      ylab="Percentage (%)",
                      xlab="Gene Region",
                      main="Distribution Across Gene Regions (Real Data)",
                      cex.main=1.5,
                      cex.lab=1.2,
                      cex.names=1.1,
                      border="black")

# 添加数值标签
text(x_positions, c(psi_values, m6a_values),
     labels=round(c(psi_values, m6a_values), 1),
     pos=3, cex=0.9, col="black")

# 添加图例
legend("topright",
       legend=c("Pseudouridine (Ψ)", "m6A"),
       fill=c(colors_psi, colors_m6a),
       border="black",
       cex=1.1)

# 添加副标题
mtext("Pseudouridine (RMBase hg19, n=895) vs m6A (REPIC hg38, HEK293T, n=10000)",
      side=3, line=0.5, cex=0.9)

# 添加统计显著性
mtext("Chi-square test: P = 1.02e-66 ***",
      side=1, line=3.5, cex=0.9, col="red", font=2)

dev.off()

# 保存PDF
pdf("figures/R_feature_distribution_real.pdf", width=10, height=6)
par(mar=c(5, 5, 5, 2))
x_positions <- barplot(rbind(psi_values, m6a_values),
                      beside=TRUE,
                      col=c(colors_psi, colors_m6a),
                      names.arg=regions,
                      ylim=c(0, 70),
                      las=1,
                      ylab="Percentage (%)",
                      xlab="Gene Region",
                      main="Distribution Across Gene Regions (Real Data)",
                      cex.main=1.5,
                      cex.lab=1.2,
                      cex.names=1.1,
                      border="black")
text(x_positions, c(psi_values, m6a_values),
     labels=round(c(psi_values, m6a_values), 1),
     pos=3, cex=0.9, col="black")
legend("topright",
       legend=c("Pseudouridine (Ψ)", "m6A"),
       fill=c(colors_psi, colors_m6a),
       border="black",
       cex=1.1)
mtext("Pseudouridine (RMBase hg19, n=895) vs m6A (REPIC hg38, HEK293T, n=10000)",
      side=3, line=0.5, cex=0.9)
mtext("Chi-square test: P = 1.02e-66 ***",
      side=1, line=3.5, cex=0.9, col="red", font=2)
dev.off()

# === 图表2: 饼图对比 ===
png("figures/R_piecharts_real.png", width=14, height=6, units="in", res=300)

par(mfrow=c(1, 2), mar=c(2, 2, 3, 2))

# Ψ饼图
pie(c(psi_5utr, psi_cds, psi_3utr),
    labels=c(paste0("5'UTR\n", round(psi_5utr, 1), "%"),
             paste0("CDS\n", round(psi_cds, 1), "%"),
             paste0("3'UTR\n", round(psi_3utr, 1), "%")),
    col=c("#2ecc71", "#f39c12", "#9b59b6"),
    main="Pseudouridine (Ψ)\nRMBase hg19",
    cex.main=1.4,
    border="black",
    init.clockwise=TRUE)

# m6A饼图
pie(c(m6a_5utr, m6a_cds, m6a_3utr),
    labels=c(paste0("5'UTR\n", round(m6a_5utr, 1), "%"),
             paste0("CDS\n", round(m6a_cds, 1), "%"),
             paste0("3'UTR\n", round(m6a_3utr, 1), "%")),
    col=c("#2ecc71", "#f39c12", "#9b59b6"),
    main="m6A\nREPIC hg38, HEK293T",
    cex.main=1.4,
    border="black",
    init.clockwise=TRUE)

dev.off()

# 保存PDF
pdf("figures/R_piecharts_real.pdf", width=14, height=6)
par(mfrow=c(1, 2), mar=c(2, 2, 3, 2))
pie(c(psi_5utr, psi_cds, psi_3utr),
    labels=c(paste0("5'UTR\n", round(psi_5utr, 1), "%"),
             paste0("CDS\n", round(psi_cds, 1), "%"),
             paste0("3'UTR\n", round(psi_3utr, 1), "%")),
    col=c("#2ecc71", "#f39c12", "#9b59b6"),
    main="Pseudouridine (Ψ)\nRMBase hg19",
    cex.main=1.4,
    border="black",
    init.clockwise=TRUE)
pie(c(m6a_5utr, m6a_cds, m6a_3utr),
    labels=c(paste0("5'UTR\n", round(m6a_5utr, 1), "%"),
             paste0("CDS\n", round(m6a_cds, 1), "%"),
             paste0("3'UTR\n", round(m6a_3utr, 1), "%")),
    col=c("#2ecc71", "#f39c12", "#9b59b6"),
    main="m6A\nREPIC hg38, HEK293T",
    cex.main=1.4,
    border="black",
    init.clockwise=TRUE)
dev.off()

# === 图表3: 点图对比 ===
png("figures/R_dotplot_real.png", width=10, height=8, units="in", res=300)

par(mfrow=c(3, 1), mar=c(4, 6, 3, 2), oma=c(2, 0, 3, 0))

regions_cn <- c("5'UTR", "CDS", "3'UTR")
psi_vec <- c(psi_5utr, psi_cds, psi_3utr)
m6a_vec <- c(m6a_5utr, m6a_cds, m6a_3utr)

for(i in 1:3) {
  # 绘制连接线
  segments(0.7, psi_vec[i], 0.7, m6a_vec[i], col="gray", lwd=1, lty=2)

  # 绘制点
  points(0.7, psi_vec[i], pch=16, col=colors_psi, cex=3)
  points(0.7, m6a_vec[i], pch=17, col=colors_m6a, cex=3)

  # 添加数值标签
  text(0.85, psi_vec[i], labels=paste0(round(psi_vec[i], 1), "%"),
       pos=4, cex=1.1, col=colors_psi, font=2)
  text(0.85, m6a_vec[i], labels=paste0(round(m6a_vec[i], 1), "%"),
       pos=4, cex=1.1, col=colors_m6a, font=2)

  # 设置坐标轴
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 70))

  axis(2, at=seq(0, 70, 10), las=1)
  mtext(regions_cn[i], side=2, line=4, cex=1.2, font=2)

  if(i == 1) {
    legend("topright",
           legend=c("Pseudouridine (Ψ)", "m6A"),
           pch=c(16, 17),
           col=c(colors_psi, colors_m6a),
           cex=1.1,
           bty="n")
  }
}

mtext("Percentage (%)", side=2, outer=TRUE, line=-0.5, cex=1.2)
mtext("Genomic Distribution Comparison (Real Data)", outer=TRUE, cex=1.5, font=2)
mtext("Chi-square test: P = 1.02e-66 ***", outer=TRUE, line=1, cex=1.1, col="red", font=2)

dev.off()

# 保存PDF
pdf("figures/R_dotplot_real.pdf", width=10, height=8)
par(mfrow=c(3, 1), mar=c(4, 6, 3, 2), oma=c(2, 0, 3, 0))
for(i in 1:3) {
  segments(0.7, psi_vec[i], 0.7, m6a_vec[i], col="gray", lwd=1, lty=2)
  points(0.7, psi_vec[i], pch=16, col=colors_psi, cex=3)
  points(0.7, m6a_vec[i], pch=17, col=colors_m6a, cex=3)
  text(0.85, psi_vec[i], labels=paste0(round(psi_vec[i], 1), "%"),
       pos=4, cex=1.1, col=colors_psi, font=2)
  text(0.85, m6a_vec[i], labels=paste0(round(m6a_vec[i], 1), "%"),
       pos=4, cex=1.1, col=colors_m6a, font=2)
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 70))
  axis(2, at=seq(0, 70, 10), las=1)
  mtext(regions_cn[i], side=2, line=4, cex=1.2, font=2)
  if(i == 1) {
    legend("topright",
           legend=c("Pseudouridine (Ψ)", "m6A"),
           pch=c(16, 17),
           col=c(colors_psi, colors_m6a),
           cex=1.1,
           bty="n")
  }
}
mtext("Percentage (%)", side=2, outer=TRUE, line=-0.5, cex=1.2)
mtext("Genomic Distribution Comparison (Real Data)", outer=TRUE, cex=1.5, font=2)
mtext("Chi-square test: P = 1.02e-66 ***", outer=TRUE, line=1, cex=1.1, col="red", font=2)
dev.off()

# === 保存统计信息 ===
stats_data <- data.frame(
  Modification = c("Pseudouridine (Psi)", "m6A"),
  Total_Sites = c(895, 10000),
  UTR5_Percent = c(psi_5utr, m6a_5utr),
  CDS_Percent = c(psi_cds, m6a_cds),
  UTR3_Percent = c(psi_3utr, m6a_3utr)
)

write.csv(stats_data, "results/R_statistics_real.csv", row.names=FALSE)

print("=== R可视化完成！===")
print("生成的图表:")
print("  - R_feature_distribution_real.png/pdf")
print("  - R_piecharts_real.png/pdf")
print("  - R_dotplot_real.png/pdf")
print("\n统计结果:")
print(stats_data)
