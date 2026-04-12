log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting clustering and UMAP for GSE131685")

required_packages <- c("Seurat", "ggplot2")

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
})

processed_dir <- file.path("data", "processed")
results_fig_dir <- file.path("results", "figures")
results_tab_dir <- file.path("results", "tables")

for (dir_path in c(processed_dir, results_fig_dir, results_tab_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# 🔥 INPUT CORRECTO (del paso 03)
input_rds <- file.path(processed_dir, "gse131685_seurat_pca.rds")

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: data/processed/gse131685_seurat_pca.rds")
}

log_step("Loading PCA-processed Seurat object")
seurat_obj <- readRDS(input_rds)

dims_to_use <- 1:20
log_step(sprintf("Using PCs: %s", paste(dims_to_use, collapse = ",")))

log_step("Running FindNeighbors")
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = dims_to_use,
  verbose = FALSE
)

log_step("Running FindClusters at resolution 0.4")
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.4,
  verbose = FALSE
)

log_step("Running UMAP")
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = dims_to_use,
  verbose = FALSE
)

# 📊 Cluster sizes
cluster_counts <- as.data.frame(table(Idents(seurat_obj)))
colnames(cluster_counts) <- c("cluster", "n_cells")

write.csv(
  cluster_counts,
  file.path(results_tab_dir, "cluster_sizes_resolution_0_4.csv"),
  row.names = FALSE
)

# 📊 UMAP embeddings
umap_embeddings <- Embeddings(seurat_obj, reduction = "umap")
write.csv(
  umap_embeddings,
  file.path(results_tab_dir, "umap_embeddings.csv"),
  row.names = TRUE
)

# 🖼️ UMAP plots
log_step("Generating UMAP by cluster")
png(file.path(results_fig_dir, "umap_by_cluster_res_0_4.png"), width = 1400, height = 1000, res = 150)
print(DimPlot(seurat_obj, reduction = "umap", label = TRUE))
dev.off()

log_step("Generating UMAP by sample")
png(file.path(results_fig_dir, "umap_by_sample_res_0_4.png"), width = 1400, height = 1000, res = 150)
print(DimPlot(seurat_obj, reduction = "umap", group.by = "sample_id"))
dev.off()

# 💾 OUTPUT ÚNICO (CLARO Y CONSISTENTE)
output_rds <- file.path(processed_dir, "gse131685_seurat_umap_clusters.rds")

saveRDS(seurat_obj, output_rds)

log_step(sprintf("Total clusters identified: %d", length(unique(Idents(seurat_obj)))))
log_step(paste("Saved clustered Seurat object to", output_rds))

log_step("Clustering and UMAP completed successfully")