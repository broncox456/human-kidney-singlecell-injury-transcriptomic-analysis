log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting normalization and PCA for GSE131685")

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

# 🔥 INPUT DESDE QC
input_rds <- file.path(processed_dir, "gse131685_seurat_qc.rds")

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found: data/processed/gse131685_seurat_qc.rds")
}

log_step("Loading QC-filtered Seurat object")
seurat_obj <- readRDS(input_rds)

# 🔬 NORMALIZATION
log_step("Running NormalizeData")
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 🔬 VARIABLE FEATURES
log_step("Running FindVariableFeatures")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

variable_features <- VariableFeatures(seurat_obj)

write.csv(
  data.frame(variable_feature = variable_features),
  file.path(results_tab_dir, "variable_features_top2000.csv"),
  row.names = FALSE
)

top10_var <- head(variable_features, 10)

write.csv(
  data.frame(top10_variable_feature = top10_var),
  file.path(results_tab_dir, "top10_variable_features.csv"),
  row.names = FALSE
)

# 🔬 SCALING
log_step("Running ScaleData")
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  verbose = FALSE
)

# 🔬 PCA
log_step("Running PCA")
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj),
  verbose = FALSE
)

# 📊 EXPORT PCA
log_step("Saving PCA embeddings")
pca_embeddings <- Embeddings(seurat_obj, reduction = "pca")

write.csv(
  pca_embeddings,
  file.path(results_tab_dir, "pca_embeddings.csv"),
  row.names = TRUE
)

# 📈 PLOTS
log_step("Generating ElbowPlot")
png(file.path(results_fig_dir, "elbowplot_pca.png"), width = 1400, height = 1000, res = 150)
print(ElbowPlot(seurat_obj, ndims = 50))
dev.off()

log_step("Generating PCA plot by sample")
png(file.path(results_fig_dir, "pca_by_sample.png"), width = 1400, height = 1000, res = 150)
print(DimPlot(seurat_obj, reduction = "pca", group.by = "sample_id"))
dev.off()

# 💾 GUARDADO FINAL (ÚNICO Y LIMPIO)
output_rds <- file.path(processed_dir, "gse131685_seurat_pca.rds")

log_step("Saving Seurat object after PCA")

saveRDS(seurat_obj, output_rds)

# 📊 LOGS FINALES
log_step(sprintf("Cells entering PCA: %d", ncol(seurat_obj)))
log_step(sprintf("Genes entering PCA: %d", nrow(seurat_obj)))
log_step(sprintf("Variable features identified: %d", length(variable_features)))

log_step(paste("Saved normalized + PCA Seurat object to", output_rds))
log_step("Normalization and PCA completed successfully")