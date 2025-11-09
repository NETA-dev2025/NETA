# NETA - Pan-Neuroendocrine Tumor Atlas

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Website](https://img.shields.io/badge/Website-neta--atlas.org-green.svg)](https://neta-dev2025.github.io/NETA/)
[![Data](https://img.shields.io/badge/Datasets-75+-orange.svg)](https://neta-dev2025.github.io/NETA/#datasets)

> A comprehensive transcriptomic atlas for neuroendocrine tumors integrating bulk RNA-seq, single-cell RNA-seq, and microarray data.

## 🎯 Overview

NETA (Pan-Neuroendocrine Tumor Atlas) is a curated database and analysis platform designed to facilitate research on neuroendocrine tumors (NETs) across diverse tissue origins. Our resource integrates multi-platform transcriptomic data with clinical annotations, providing researchers with a comprehensive tool for biomarker discovery, molecular subtyping, and therapeutic target identification.

### Key Features

- **📊 Multi-Platform Integration**: Bulk RNA-seq, scRNA-seq, and microarray data
- **🔬 Comprehensive Coverage**: 75+ datasets spanning multiple NET subtypes
- **💾 High-Quality Data**: Standardized processing with rigorous QC
- **🧬 Rich Annotations**: Clinical, genomic, and drug response data
- **🛠️ Analysis Tools**: Built-in differential expression and pathway analysis
- **📖 Open Access**: All data freely available with comprehensive documentation

## 📈 Database Statistics

| Metric | Count |
|--------|-------|
| **Total Datasets** | 75+ |
| **Total Samples** | 2,500+ |
| **Genes Profiled** | 22,593 |
| **Data Platforms** | 3 (RNA-seq, scRNA-seq, Microarray) |
| **Tissue Types** | 10+ |
| **Expression Records** | 4.7M+ |

## 🗂️ Data Types

### Bulk RNA Sequencing

High-throughput RNA sequencing data providing comprehensive gene expression profiles:

- **DESeq2-ready count matrices** for differential expression analysis
- **Normalized expression data** (TPM, FPKM)
- **Quality control metrics** and batch information
- **Clinical annotations** including staging, treatment, and outcomes

### Single-Cell RNA-seq

Single-cell resolution transcriptomic data revealing cellular heterogeneity:

- **Cell-type annotations** and cluster assignments
- **UMAP/t-SNE coordinates** for visualization
- **Marker gene lists** for each cell type
- **Trajectory analysis** results where available

### Microarray Data

Legacy microarray platforms providing historical context:

- **Normalized intensity values**
- **Platform-specific annotations**
- **Quality control reports**
- **Cross-platform comparable data**

## 🚀 Quick Start

### Browse the Web Interface

Visit our online platform to explore datasets interactively:

```
https://neta-dev2025.github.io/NETA/
```

### Download Data

All datasets are available for download in standard formats:

1. Navigate to the [Datasets](https://neta-dev2025.github.io/NETA/#datasets) section
2. Filter by data type, tissue, or tumor characteristics
3. Download count matrices, metadata, and analysis results

### Run Analysis

Use our pre-configured analysis scripts:

```bash
# Clone the repository
git clone https://github.com/NETA-dev2025/NETA.git
cd NETA

# Install dependencies
Rscript install_dependencies.R

# Run example analysis
Rscript examples/deseq2_analysis.R
```

## 📚 Documentation

Comprehensive guides and tutorials are available:

- **[User Guide](docs/user-guide.md)** - Getting started with NETA
- **[Data Dictionary](docs/data-dictionary.md)** - Understanding data formats
- **[API Reference](docs/api.md)** - Programmatic access
- **[Tutorials](docs/tutorials.md)** - Step-by-step analysis examples
- **[FAQ](docs/faq.md)** - Frequently asked questions

## 🛠️ Analysis Workflows

NETA provides ready-to-use analysis pipelines for common tasks:

### Differential Expression Analysis

```r
library(DESeq2)

# Load NETA dataset
counts <- read.csv("data/GSE98894/GSE98894_raw_counts.csv", row.names=1)
metadata <- read.csv("data/GSE98894/sample_metadata.csv", row.names=1)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# Run analysis
dds <- DESeq(dds)
results <- results(dds)
```

### Pathway Enrichment

```r
library(clusterProfiler)

# Extract significant genes
sig_genes <- rownames(subset(results, padj < 0.05))

# GO enrichment
ego <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05
)
```

### Single-Cell Analysis

```r
library(Seurat)

# Load single-cell data
seurat_obj <- readRDS("data/GSE126030/seurat_object.rds")

# Standard workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
```

## 📊 Supported Tumor Types

NETA includes data from diverse neuroendocrine tumor origins:

- **Pancreatic NET** (pNET)
- **Small Intestine NET** (SI-NET)
- **Lung Carcinoid** (TC, AC, SCLC, LCNEC)
- **Gastrointestinal NET**
- **Pheochromocytoma** (PHEO)
- **Paraganglioma** (PGL)
- **Medullary Thyroid Carcinoma** (MTC)
- **Merkel Cell Carcinoma** (MCC)
- **And more...**

## 🤝 Contributing

We welcome contributions from the research community!

### How to Contribute

1. **Submit New Datasets**: Email us with GEO/SRA accession numbers
2. **Report Issues**: Use GitHub Issues for bugs or suggestions
3. **Improve Documentation**: Submit pull requests for documentation
4. **Share Analysis Scripts**: Contribute analysis workflows

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## 📄 Data Usage and Citation

### License

All data and code are released under the **MIT License**, allowing free use for academic and commercial purposes.

### Citation

If you use NETA in your research, please cite:

```bibtex
@article{NETA2025,
  title={NETA: A Pan-Neuroendocrine Tumor Transcriptomic Atlas},
  author={NETA Research Team},
  journal={In Preparation},
  year={2025}
}
```

### Data Sources

NETA aggregates data from multiple public repositories:

- **GEO** (Gene Expression Omnibus)
- **TCGA** (The Cancer Genome Atlas)
- **Published Studies** with available raw data

Please cite original data sources when publishing results.

## 🔗 Links

- **Website**: https://neta-dev2025.github.io/NETA/
- **GitHub**: https://github.com/NETA-dev2025/NETA
- **Documentation**: https://neta-dev2025.github.io/NETA/docs/
- **Contact**: contact@neta-atlas.org

## 🌟 Acknowledgments

We thank the research community for making their data publicly available and the developers of key bioinformatics tools:

- DESeq2, edgeR, limma (differential expression)
- Seurat, Scanpy (single-cell analysis)
- clusterProfiler, GSEA (pathway analysis)
- GEOquery, TCGAbiolinks (data access)

## 📅 Updates

- **2025-11**: Initial release with 75 datasets
- **2025-11**: Added bulk RNA-seq analysis workflows
- **Future**: Single-cell integration, survival analysis

---

**NETA Research Team** | 2025

For questions or feedback, please [open an issue](https://github.com/NETA-dev2025/NETA/issues) or contact us at contact@neta-atlas.org.
