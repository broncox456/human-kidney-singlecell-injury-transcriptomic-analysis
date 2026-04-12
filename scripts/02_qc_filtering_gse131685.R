log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting QC and filtering for GSE131685")

required_packages <- c("Seurat", "ggplot2", "patchwork")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required packages: %s. Please install them and rerun.",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

processed_dir <- file.path("data", "processed")
results_fig_dir <- file.path("results", "figures")
results_tab_dir <- file.path("results", "tables")

for (dir_path in c(processed_dir, results_fig_dir, results_tab_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

input_rds <- file.path(processed_dir, "gse131685_seurat_raw.rds")

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: data/processed/gse131685_seurat_raw.rds")
}

log_step("Loading raw Seurat object")
seurat_obj <- readRDS(input_rds)

log_step("Calculating mitochondrial percentage")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

qc_before <- data.frame(
  metric = c("cells", "genes", "median_nCount_RNA", "median_nFeature_RNA", "median_percent_mt"),
  value = c(
    ncol(seurat_obj),
    nrow(seurat_obj),
    median(seurat_obj$nCount_RNA),
    median(seurat_obj$nFeature_RNA),
    median(seurat_obj$percent.mt)
  )
)

write.csv(
  qc_before,
  file.path(results_tab_dir, "qc_summary_before_filtering.csv"),
  row.names = FALSE
)

log_step("Generating violin plots before filtering")
png(
  file.path(results_fig_dir, "qc_violin_before_filtering.png"),
  width = 1800,
  height = 900,
  res = 150
)
print(
  VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
)
dev.off()

log_step("Generating scatter plots before filtering")
p1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

png(
  file.path(results_fig_dir, "qc_scatter_before_filtering.png"),
  width = 1800,
  height = 900,
  res = 150
)
print(p1 + p2)
dev.off()

log_step("Applying filtering thresholds")
min_features <- 200
max_features <- 6000
max_mt <- 15

log_step(sprintf(
  "Thresholds: nFeature_RNA >= %d, nFeature_RNA <= %d, percent.mt <= %d",
  min_features, max_features, max_mt
))

seurat_qc <- subset(
  seurat_obj,
  subset = nFeature_RNA >= min_features &
    nFeature_RNA <= max_features &
    percent.mt <= max_mt
)

if (ncol(seurat_qc) == 0) {
  stop("All cells were removed during QC filtering. Thresholds are too strict for this dataset.")
}

qc_after <- data.frame(
  metric = c("cells", "genes", "median_nCount_RNA", "median_nFeature_RNA", "median_percent_mt"),
  value = c(
    ncol(seurat_qc),
    nrow(seurat_qc),
    median(seurat_qc$nCount_RNA),
    median(seurat_qc$nFeature_RNA),
    median(seurat_qc$percent.mt)
  )
)

write.csv(
  qc_after,
  file.path(results_tab_dir, "qc_summary_after_filtering.csv"),
  row.names = FALSE
)

filtering_summary <- data.frame(
  cells_before = ncol(seurat_obj),
  cells_after = ncol(seurat_qc),
  cells_removed = ncol(seurat_obj) - ncol(seurat_qc),
  percent_removed = round(
    100 * (ncol(seurat_obj) - ncol(seurat_qc)) / ncol(seurat_obj),
    2
  ),
  genes_before = nrow(seurat_obj),
  genes_after = nrow(seurat_qc),
  min_features = min_features,
  max_features = max_features,
  max_mt = max_mt
)

write.csv(
  filtering_summary,
  file.path(results_tab_dir, "filtering_summary.csv"),
  row.names = FALSE
)

log_step("Generating violin plots after filtering")
png(
  file.path(results_fig_dir, "qc_violin_after_filtering.png"),
  width = 1800,
  height = 900,
  res = 150
)
print(
  VlnPlot(
    seurat_qc,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  )
)
dev.off()

log_step("Saving QC-filtered Seurat object")
saveRDS(
  seurat_qc,
  file.path(processed_dir, "gse131685_seurat_qc.rds")
)

log_step(sprintf("Cells before QC: %d", ncol(seurat_obj)))
log_step(sprintf("Cells after QC: %d", ncol(seurat_qc)))
log_step("Saved filtered Seurat object to data/processed/gse131685_seurat_qc.rds")
log_step("QC and filtering completed successfully")