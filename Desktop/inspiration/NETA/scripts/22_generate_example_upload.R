#==============================================================================
# Generate Example CSV for AI Diagnosis Upload
# Creates a template CSV file with the 50 signature genes
#==============================================================================

library(tidyverse)

# Load signature genes
sig_genes <- readLines("results/08_ai_model/nen_signature_genes.txt")
sig_genes <- sig_genes[sig_genes != ""]  # Remove empty lines

# Load real expression data
expr_file <- "data/integrated/real_expression_matrix.rds"
meta_file <- "data/integrated/pan_nen_metadata_with_subtypes.csv"

expr <- readRDS(expr_file)
meta <- read_csv(meta_file, show_col_types = FALSE)

# Select 5 random samples (1-2 from each subtype for demonstration)
set.seed(42)
sample_c1 <- sample(meta$sample_id[meta$molecular_subtype == "C1"], 2)
sample_c2 <- sample(meta$sample_id[meta$molecular_subtype == "C2"], 2)
sample_c3 <- sample(meta$sample_id[meta$molecular_subtype == "C3"], 1)

demo_samples <- c(sample_c1, sample_c2, sample_c3)

# Extract expression for signature genes and demo samples
# Find which signature genes exist in our expression matrix
available_genes <- sig_genes[sig_genes %in% rownames(expr)]
cat("Available signature genes in expression matrix:", length(available_genes), "out of", length(sig_genes), "\n")

# Extract subset
demo_expr <- expr[available_genes, demo_samples]

# Create output dataframe (Genes as rows, Samples as columns)
output_df <- as.data.frame(demo_expr)
output_df <- tibble::rownames_to_column(output_df, var = "Gene")

# Save example file
out_dir <- "results/08_ai_model"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

example_file <- file.path(out_dir, "example_upload_template.csv")
write_csv(output_df, example_file)

cat("✓ Example CSV generated:", example_file, "\n")
cat("  Format: Genes (rows) × Samples (columns)\n")
cat("  Dimensions:", nrow(output_df), "genes ×", ncol(output_df)-1, "samples\n")

# Also save a blank template (just gene names)
template_df <- data.frame(Gene = available_genes, Sample_1 = NA, Sample_2 = NA, Sample_3 = NA)
template_file <- file.path(out_dir, "blank_upload_template.csv")
write_csv(template_df, template_file)

cat("✓ Blank template generated:", template_file, "\n")
cat("  Users can fill in expression values for their own samples.\n")

# Print summary of true labels for the demo samples
demo_meta <- meta %>% filter(sample_id %in% demo_samples)
cat("\n--- Demo Sample True Labels ---\n")
print(demo_meta %>% select(sample_id, molecular_subtype, organ))
