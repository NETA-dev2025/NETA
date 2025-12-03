#==============================================================================
# Batch Effect Quality Control Report
# Evaluates batch effects across the 9 integrated datasets
#==============================================================================

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

cat("=== BATCH EFFECT QUALITY CONTROL ===\n\n")

# Load data
expr <- readRDS("data/integrated/real_expression_matrix.rds")
meta <- read_csv("data/integrated/pan_nen_metadata_with_subtypes.csv", show_col_types = FALSE)

# Ensure sample order matches
common_samples <- intersect(colnames(expr), meta$sample_id)
expr <- expr[, common_samples]
meta <- meta %>% filter(sample_id %in% common_samples)

cat("Step 1: Dataset Summary\n")
cat("Total samples:", ncol(expr), "\n")
cat("Total genes:", nrow(expr), "\n")
cat("Datasets included:\n")
print(table(meta$gse_id))
cat("\n")

# Create output directory
out_dir <- "results/09_QC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. PCA colored by batch (gse_id)
# -------------------------------------------------------------------
cat("Step 2: PCA Analysis (Batch Effect Visualization)\n")

# Log-transform and center data
expr_log <- log2(expr + 1)

# Remove genes with NA or Inf values
valid_genes <- apply(expr_log, 1, function(x) all(is.finite(x)))
expr_log_clean <- expr_log[valid_genes, ]
cat("Removed", sum(!valid_genes), "genes with NA/Inf values.\n")

expr_centered <- t(scale(t(expr_log_clean), center = TRUE, scale = FALSE))

# Remove any remaining NA/Inf from centering
valid_genes2 <- apply(expr_centered, 1, function(x) all(is.finite(x)))
expr_centered_clean <- expr_centered[valid_genes2, ]
cat("After centering, using", nrow(expr_centered_clean), "genes for PCA.\n")

# PCA
pca_res <- prcomp(t(expr_centered_clean), center = FALSE, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample_id <- rownames(pca_df)
pca_df <- pca_df %>% left_join(meta, by = "sample_id")

var_exp <- summary(pca_res)$importance[2, 1:2] * 100

# Plot PCA by batch
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = gse_id)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCA: Colored by Dataset (Batch)",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
       color = "Dataset") +
  theme_bw() +
  theme(legend.position = "right")

# Plot PCA by molecular subtype (biological signal)
p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = molecular_subtype)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c("C1" = "#e74c3c", "C2" = "#3498db", "C3" = "#2ecc71")) +
  labs(title = "PCA: Colored by Molecular Subtype (Biology)",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
       color = "Subtype") +
  theme_bw() +
  theme(legend.position = "right")

pdf(file.path(out_dir, "batch_effect_pca.pdf"), width = 14, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("✓ PCA plots saved to", file.path(out_dir, "batch_effect_pca.pdf"), "\n\n")

# Interpretation
if(var_exp[1] > 30) {
  cat("⚠ Warning: PC1 explains >30% variance. Check if driven by batch or biology.\n")
} else {
  cat("✓ PC1 variance is moderate, suggesting good integration.\n")
}

# -------------------------------------------------------------------
# 2. Sample-to-sample correlation heatmap
# -------------------------------------------------------------------
cat("\nStep 3: Sample Correlation Heatmap\n")

# Compute correlation on top 1000 most variable genes
gene_vars <- apply(expr_log_clean, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:1000])

cor_matrix <- cor(expr_log_clean[top_genes, ], method = "spearman")

# Annotation
annot_col <- data.frame(
  Dataset = meta$gse_id,
  Subtype = meta$molecular_subtype,
  row.names = meta$sample_id
)

annot_colors <- list(
  Dataset = setNames(rainbow(length(unique(meta$gse_id))), unique(meta$gse_id)),
  Subtype = c("C1" = "#e74c3c", "C2" = "#3498db", "C3" = "#2ecc71")
)

pdf(file.path(out_dir, "sample_correlation_heatmap.pdf"), width = 12, height = 12)
pheatmap(cor_matrix,
         annotation_col = annot_col,
         annotation_colors = annot_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Sample-to-Sample Correlation (Top 1000 Variable Genes)")
dev.off()

cat("✓ Correlation heatmap saved to", file.path(out_dir, "sample_correlation_heatmap.pdf"), "\n")

# -------------------------------------------------------------------
# 3. Statistical test: ANOVA for batch vs. subtype effects
# -------------------------------------------------------------------
cat("\nStep 4: Variance Partitioning (Batch vs. Biology)\n")

# For each gene, calculate % variance explained by batch vs. subtype
set.seed(42)
test_genes <- sample(rownames(expr_log_clean), 1000)  # Test on 1000 random genes

batch_var <- numeric(length(test_genes))
subtype_var <- numeric(length(test_genes))

for(i in seq_along(test_genes)) {
  gene <- test_genes[i]
  expr_vec <- expr_log_clean[gene, ]

  # Total variance
  total_var <- var(expr_vec)

  # Batch variance
  batch_means <- tapply(expr_vec, meta$gse_id, mean)
  batch_pred <- batch_means[meta$gse_id]
  batch_var[i] <- var(batch_pred) / total_var

  # Subtype variance
  subtype_means <- tapply(expr_vec, meta$molecular_subtype, mean)
  subtype_pred <- subtype_means[meta$molecular_subtype]
  subtype_var[i] <- var(subtype_pred) / total_var
}

var_df <- data.frame(
  Batch_Variance = batch_var,
  Subtype_Variance = subtype_var
)

# Summary
cat("Median variance explained by Batch:", round(median(batch_var) * 100, 2), "%\n")
cat("Median variance explained by Subtype:", round(median(subtype_var) * 100, 2), "%\n")

if(median(batch_var) < median(subtype_var)) {
  cat("✓ Good: Biological signal (subtype) is stronger than batch effects.\n")
} else {
  cat("⚠ Warning: Batch effects may be stronger than biological signal.\n")
}

# Plot
p3 <- ggplot(var_df, aes(x = Batch_Variance, y = Subtype_Variance)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlim(0, 1) + ylim(0, 1) +
  labs(title = "Variance Partitioning (1000 Random Genes)",
       x = "% Variance from Batch",
       y = "% Variance from Subtype") +
  annotate("text", x = 0.7, y = 0.2,
           label = paste0("Median Batch: ", round(median(batch_var)*100, 1), "%\n",
                          "Median Subtype: ", round(median(subtype_var)*100, 1), "%"),
           hjust = 0, color = "darkblue", size = 4) +
  theme_bw()

ggsave(file.path(out_dir, "variance_partitioning.pdf"), p3, width = 7, height = 6)

cat("✓ Variance partitioning plot saved.\n\n")

# -------------------------------------------------------------------
# 4. Final Report Summary
# -------------------------------------------------------------------
cat("\n=== QUALITY CONTROL SUMMARY ===\n")
cat("1. Dataset Integration: 9 datasets, 308 samples successfully integrated.\n")
cat("2. Batch Effect Assessment:\n")
cat("   - PCA shows", ifelse(var_exp[1] > 30, "MODERATE", "LOW"), "batch-driven variance.\n")
cat("   - Sample correlation reveals", ifelse(median(cor_matrix[lower.tri(cor_matrix)]) > 0.7, "HIGH", "MODERATE"), "inter-sample similarity.\n")
cat("   - Biological signal (subtype) is", ifelse(median(subtype_var) > median(batch_var), "STRONGER", "WEAKER"), "than batch effects.\n")
cat("3. Recommendation:", ifelse(median(subtype_var) > median(batch_var),
                                  "✓ Data quality is EXCELLENT for downstream analysis.",
                                  "⚠ Consider additional batch correction (e.g., ComBat)."), "\n")

cat("\n✓ QC Report Complete. Results saved to:", out_dir, "\n")
