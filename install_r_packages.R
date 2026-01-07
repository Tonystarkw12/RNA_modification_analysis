#!/usr/bin/env Rscript
# 安装R包依赖

# 设置用户库路径
lib_path <- "~/R/library"
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

cat("=== R包安装路径:", lib_path, "===\n\n")

cat("=== 安装CRAN包 ===\n")

# CRAN包
cran_packages <- c("ggplot2", "dplyr", "readr", "gridExtra", "RColorBrewer", "scales")

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("安装 %s...\n", pkg))
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/", quiet = TRUE)
  } else {
    cat(sprintf("✓ %s 已安装\n", pkg))
  }
}

cat("\n=== 安装Bioconductor包 ===\n")

# Bioconductor包
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}

bioc_packages <- c("ChIPseeker", "VennDiagram")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("安装 %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("✓ %s 已安装\n", pkg))
  }
}

cat("\n=== 所有R包安装完成！===\n")
