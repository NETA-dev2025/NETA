#!/bin/bash
#==============================================================================
# NETA Project Release Package Creator
# Creates a clean zip package of the NETA project for distribution
#==============================================================================

echo "=== NETA Project Packaging Tool ==="
echo ""

# Set variables
PROJECT_NAME="NETA_v3.5_Pan-NEN_Atlas"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="$HOME/Desktop"
OUTPUT_FILE="${OUTPUT_DIR}/${PROJECT_NAME}_${TIMESTAMP}.zip"

# Create temporary directory for clean copy
TEMP_DIR=$(mktemp -d)
PACKAGE_DIR="${TEMP_DIR}/NETA"

echo "Step 1: Creating clean copy of project..."
mkdir -p "$PACKAGE_DIR"

# Copy essential files
echo "Step 2: Copying core files..."
cp app.R "$PACKAGE_DIR/"
cp README.md "$PACKAGE_DIR/"

# Copy directories
echo "Step 3: Copying data and results..."
cp -r scripts "$PACKAGE_DIR/"
cp -r data "$PACKAGE_DIR/"
cp -r results "$PACKAGE_DIR/"

# Create a MANIFEST file
echo "Step 4: Creating MANIFEST..."
cat > "$PACKAGE_DIR/MANIFEST.txt" << 'EOF'
NETA v3.5 - Pan-Neuroendocrine Tumor Atlas
==========================================

Package Contents:
- app.R: Shiny web application (main entry point)
- README.md: Project documentation
- scripts/: All R analysis scripts (01-23)
- data/: Raw and integrated expression data
  - data/raw/geo/: GEO datasets (ExpressionSet objects)
  - data/raw/tcga/: TCGA-PCPG validation data
  - data/integrated/: Integrated expression matrix and metadata
- results/: Analysis outputs
  - results/04_subtyping/: Molecular subtyping results
  - results/05_immune/: Immune infiltration scores
  - results/06_drug_response/: Drug sensitivity predictions
  - results/08_ai_model/: XGBoost AI model and signature genes
  - results/09_QC/: Quality control reports

Quick Start:
1. Install R packages: shiny, shinydashboard, tidyverse, plotly, DT,
   pheatmap, survival, survminer, xgboost, SHAPforxgboost
2. Run: Rscript app.R
3. Open browser at displayed URL (http://127.0.0.1:xxxx)

Data Sources:
All data from NCBI Gene Expression Omnibus (GEO):
- GSE118014, GSE151904, GSE208446, GSE213504, GSE229701
- GSE235092, GSE244945, GSE246690, GSE60052

Citation:
[Your Name], et al. (2024). Pan-NEN Atlas: A comprehensive transcriptomic
resource and AI diagnostic tool for neuroendocrine neoplasms.

Generated: $(date)
EOF

# Clean up unnecessary files from the package
echo "Step 5: Cleaning up temporary and log files..."
cd "$PACKAGE_DIR"
find . -name "*.log" -delete
find . -name ".DS_Store" -delete
find . -name "nohup.out" -delete
find . -name ".Rhistory" -delete
find . -name "*.tmp" -delete

# Create zip package
echo "Step 6: Creating zip archive..."
cd "$TEMP_DIR"
zip -r "$OUTPUT_FILE" NETA/ -q

# Clean up temp directory
rm -rf "$TEMP_DIR"

# Report
FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
echo ""
echo "=== Package Created Successfully ==="
echo "Output: $OUTPUT_FILE"
echo "Size: $FILE_SIZE"
echo ""
echo "Contents:"
unzip -l "$OUTPUT_FILE" | tail -20
echo ""
echo "Package is ready for distribution!"
