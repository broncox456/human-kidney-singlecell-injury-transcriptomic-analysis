log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting marker gene analysis for GSE131685")

required_packages <- c("Seurat", "SeuratObject", "dplyr")

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
  library(SeuratObject)
  library(dplyr)
})

processed_dir <- file.path("data", "processed")
results_dir <- file.path("results", "tables")

for (dir_path in c(processed_dir, results_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

input_rds <- file.path(processed_dir, "gse131685_seurat_umap_clusters.rds")

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: data/processed/gse131685_seurat_umap_clusters.rds")
}

log_step("Loading clustered Seurat object")
seurat_obj <- readRDS(input_rds)

log_step("Joining RNA assay layers for Seurat v5 differential expression")
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

joined_rds <- file.path(processed_dir, "gse131685_seurat_umap_clusters_joined.rds")
saveRDS(seurat_obj, joined_rds)
log_step("Saved joined Seurat object for marker analysis")

log_step("Running FindAllMarkers")
markers <- FindAllMarkers(
  object = seurat_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

if (nrow(markers) == 0) {
  stop("FindAllMarkers returned an empty table even after JoinLayers.")
}

if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}

cluster_col <- NULL
if ("cluster" %in% colnames(markers)) {
  cluster_col <- "cluster"
} else if ("cluster_id" %in% colnames(markers)) {
  cluster_col <- "cluster_id"
} else if ("seurat_clusters" %in% colnames(markers)) {
  cluster_col <- "seurat_clusters"
} else {
  stop("No cluster column found in markers table.")
}

if (!"avg_log2FC" %in% colnames(markers)) {
  stop("Column avg_log2FC not found in markers table.")
}

log_step(sprintf("Detected %d marker rows", nrow(markers)))
log_step(sprintf("Using cluster column: %s", cluster_col))

write.csv(
  markers,
  file.path(results_dir, "all_markers.csv"),
  row.names = FALSE
)

log_step("Extracting top 10 markers per cluster")
top_markers <- markers %>%
  group_by(.data[[cluster_col]]) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_markers,
  file.path(results_dir, "top10_markers_per_cluster.csv"),
  row.names = FALSE
)

cluster_marker_counts <- top_markers %>%
  count(.data[[cluster_col]], name = "n_top_markers")

write.csv(
  cluster_marker_counts,
  file.path(results_dir, "top_marker_counts_by_cluster.csv"),
  row.names = FALSE
)

log_step(sprintf("Top marker table rows: %d", nrow(top_markers)))
log_step("Saved all_markers.csv")
log_step("Saved top10_markers_per_cluster.csv")
log_step("Saved top_marker_counts_by_cluster.csv")
log_step("Marker gene analysis completed successfully")