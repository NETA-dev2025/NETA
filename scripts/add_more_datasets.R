#!/usr/bin/env Rscript
# 快速添加更多真实的神经内分泌肿瘤bulk RNA-seq数据集

library(GEOquery)
library(jsonlite)

# 目标：找到并验证更多真实的bulk RNA-seq数据集
# 这些是经过文献验证的神经内分泌肿瘤相关数据集

cat("\n========================================\n")
cat("添加真实Bulk RNA-seq数据集\n")
cat("========================================\n\n")

# 设置工作目录
setwd("/Users/Apple/Desktop/Modeling/NETA-new/raw_data")

# 高质量的bulk RNA-seq数据集列表（经过文献验证）
candidate_datasets <- c(
  "GSE65286",   # Small cell lung cancer
  "GSE60052",   # Pancreatic NET
  "GSE73338",   # 之前检查过但要重新验证
  "GSE98894"    # 已有的
)

results <- list()

for (gse_id in candidate_datasets) {
  cat(sprintf("\n[检查] %s\n", gse_id))
  cat(paste(rep("-", 50), collapse=""), "\n")

  tryCatch({
    # 获取GEO metadata
    gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

    if (length(gse) > 0) {
      pdata <- pData(gse[[1]])

      cat(sprintf("✓ 样本数: %d\n", nrow(pdata)))
      cat(sprintf("✓ 平台: %s\n", annotation(gse[[1]])))

      # 检查是否有supplementary files
      cat("\n尝试下载supplementary files...\n")

      gse_dir <- file.path(getwd(), gse_id)
      dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)

      supp_files <- getGEOSuppFiles(gse_id, baseDir = getwd())

      files <- list.files(gse_dir, full.names = TRUE)
      cat(sprintf("✓ 下载 %d 个文件\n", length(files)))

      # 查找count文件
      count_keywords <- c("count", "raw", "htseq", "featurecounts", "gene_count")
      count_files <- c()

      for (keyword in count_keywords) {
        matches <- grep(keyword, files, ignore.case = TRUE, value = TRUE)
        count_files <- c(count_files, matches)
      }
      count_files <- unique(count_files)

      if (length(count_files) > 0) {
        cat(sprintf("✓ 找到 %d 个可能的count文件\n", length(count_files)))

        results[[gse_id]] <- list(
          status = "potential",
          samples = nrow(pdata),
          platform = annotation(gse[[1]]),
          count_files = basename(count_files)
        )
      } else {
        cat("⚠ 未找到count文件\n")
        results[[gse_id]] <- list(
          status = "no_counts",
          samples = nrow(pdata)
        )
      }
    }

  }, error = function(e) {
    cat(sprintf("✗ 错误: %s\n", e$message))
    results[[gse_id]] <<- list(status = "error", error = e$message)
  })

  Sys.sleep(2)  # 避免请求过快
}

# 保存结果
write_json(results, "candidate_datasets_check.json",
          pretty = TRUE, auto_unbox = TRUE)

cat("\n\n========================================\n")
cat("检查完成\n")
cat("========================================\n")

# 显示总结
success_count <- sum(sapply(results, function(x) x$status == "potential"))
cat(sprintf("\n✓ 找到 %d 个有count文件的数据集\n", success_count))
cat("\n详细结果已保存到: candidate_datasets_check.json\n\n")
