  
# Human Kidney Single-Cell RNA-seq – Injury and Cellular Heterogeneity Analysis

This project explores cellular heterogeneity and injury-related transcriptional programs in human kidney tissue using single-cell RNA sequencing data from **GSE131685**.

The objective was to build a structured and reproducible Seurat v5 pipeline, focusing on biologically meaningful clustering and marker gene identification within a clinically relevant nephrology context.

![UMAP Clustering](results/figures/umap_final_publication.png)
---

## Clinical Context

Kidney disease is characterized by complex interactions between epithelial, immune, and stromal cell populations.

Understanding these interactions at the single-cell level is critical for:

- identifying early injury signatures  
- detecting cell-type–specific vulnerability  
- understanding inflammatory infiltration  
- exploring mechanisms of disease progression  

This analysis focuses on **tubular function**, **immune infiltration**, and **cellular stress responses**, all central to kidney injury.

---

## Dataset

- **Source:** Gene Expression Omnibus (GEO)  
- **Accession:** GSE131685  
- **Data type:** Single-cell RNA sequencing  
- **Framework:** Seurat v5 (R)

---

## Analytical Workflow

The analysis was performed using a reproducible stepwise pipeline:

1. Download GEO supplementary files  
2. Construct Seurat objects from raw count matrices  
3. Merge samples into a unified dataset  
4. Perform quality control filtering  
5. Normalize expression data  
6. Identify variable features  
7. Run PCA  
8. Perform graph-based clustering  
9. Generate UMAP visualizations  
10. Identify cluster-specific marker genes  
11. Perform preliminary biological annotation  

---

## Quality Control Strategy

Cells were filtered using:

- **nFeature_RNA ≥ 200**
- **nFeature_RNA ≤ 6000**
- **percent.mt ≤ 15**

This removed low-quality cells, empty droplets, and high-mitochondrial profiles associated with stressed or dying cells.

---

## Key Results

### 1. Clustering

A total of **13 transcriptionally distinct clusters** were identified (resolution = 0.4), reflecting substantial cellular heterogeneity.

Cluster sizes ranged from large epithelial populations to smaller specialized or immune subsets.

---

### 2. Dimensionality Reduction

Dimensionality reduction was successfully performed using:

- Principal Component Analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)

**Outputs:**
- `elbowplot_pca.png`
- `pca_by_sample.png`
- `umap_by_cluster_res_0_4.png`
- `umap_by_sample_res_0_4.png`

---

### 3. Marker Gene Analysis

Marker detection was performed using `FindAllMarkers()`.

A total of **~5300 marker genes** were identified across clusters.

#### Important technical note

During analysis, a Seurat v5-specific issue was encountered:

- Differential expression initially returned no markers  
- Warning indicated that assay layers were not joined  
- This was resolved using:

```r
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

Biological Interpretation

Cluster-level marker analysis revealed distinct renal and immune populations:

Renal epithelial populations
Proximal tubule:
FABP1, S100A1, GPX3
Distal tubule:
PVALB, DUSP9
Collecting duct:
AQP2, CLDN8
Immune populations
T / NK cells:
TRDC, GZMA, CTSW
B cells:
CD79A, MS4A1, IGHM
Monocytes / macrophages:
CSF1R, FPR1, S100A12
Injury and stress signatures
GSTP1, SOD2, APP

These findings suggest:

active epithelial heterogeneity
immune infiltration within kidney tissue
oxidative stress and injury-related transcriptional programs

kidney-singlecell-injury-heterogeneity-analysis/
├── data/
│   ├── raw/
│   ├── processed/
│   └── temp_10x/
├── results/
│   ├── figures/
│   └── tables/
├── scripts/
│   ├── 01_initialize_seurat.R
│   ├── 02_qc_filtering.R
│   ├── 03_pca_clustering.R
│   ├── 04_clustering_umap.R
│   ├── 05_marker_genes_annotation.R
│   └── 06_cluster_annotation_summary.R
└── README.md


Main Outputs :

Processed objects
gse131685_seurat_raw.rds
gse131685_seurat_qc.rds
gse131685_seurat_pca.rds
gse131685_seurat_umap_clusters.rds

Tables :

QC summaries
PCA embeddings
UMAP embeddings
Cluster sizes
Marker gene tables (~5300 genes)
Top 10 markers per cluster
Cluster annotation summary

Figures :

QC violin plots
QC scatter plots
PCA visualization
UMAP visualization

Limitations :

Cell type annotation is preliminary
No external validation dataset
No trajectory or pseudotime analysis
Limited integration with clinical metadata

Future Directions : 

refined cell-type annotation using curated databases
trajectory analysis (injury → repair transitions)
integration with clinical variables
comparative analysis across disease states

Why This Project :

This project demonstrates:
reproducible single-cell RNA-seq analysis
structured Seurat v5 workflow design
real-world debugging of differential expression issues
clinically oriented interpretation in nephrology
integration of medical expertise with data science

Author

Cristian Arias, MD
Nephrologist | Healthcare Data Scientist | Bioinformatics MSc Candidate