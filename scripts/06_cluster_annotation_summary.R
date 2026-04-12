log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting cluster annotation summary for GSE131685")

required_packages <- c("dplyr", "readr", "stringr")

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
  library(dplyr)
  library(readr)
  library(stringr)
})

results_dir <- file.path("results", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

all_markers_file <- file.path(results_dir, "all_markers.csv")
top_markers_file <- file.path(results_dir, "top10_markers_per_cluster.csv")
cluster_sizes_file <- file.path(results_dir, "cluster_sizes_resolution_0_4.csv")

if (!file.exists(all_markers_file)) {
  stop("Missing file: results/tables/all_markers.csv")
}

if (!file.exists(top_markers_file)) {
  stop("Missing file: results/tables/top10_markers_per_cluster.csv")
}

if (!file.exists(cluster_sizes_file)) {
  stop("Missing file: results/tables/cluster_sizes_resolution_0_4.csv")
}

log_step("Loading marker and cluster tables")
all_markers <- readr::read_csv(all_markers_file, show_col_types = FALSE)
top_markers <- readr::read_csv(top_markers_file, show_col_types = FALSE)
cluster_sizes <- readr::read_csv(cluster_sizes_file, show_col_types = FALSE)

if (nrow(all_markers) == 0) {
  stop("all_markers.csv is empty")
}

if (nrow(top_markers) == 0) {
  stop("top10_markers_per_cluster.csv is empty")
}

cluster_col <- NULL
if ("cluster" %in% colnames(top_markers)) {
  cluster_col <- "cluster"
} else if ("cluster_id" %in% colnames(top_markers)) {
  cluster_col <- "cluster_id"
} else if ("seurat_clusters" %in% colnames(top_markers)) {
  cluster_col <- "seurat_clusters"
} else {
  stop("No cluster column found in top marker table.")
}

if (!"gene" %in% colnames(top_markers)) {
  stop("No gene column found in top marker table.")
}

if (!"avg_log2FC" %in% colnames(top_markers)) {
  stop("No avg_log2FC column found in top marker table.")
}

log_step(sprintf("Using cluster column: %s", cluster_col))

# -------------------------------
# 1. Marker counts per cluster
# -------------------------------
log_step("Summarizing marker counts per cluster")

marker_counts <- all_markers %>%
  count(.data[[cluster_col]], name = "n_markers_total") %>%
  arrange(as.numeric(.data[[cluster_col]]))

write_csv(
  marker_counts,
  file.path(results_dir, "marker_counts_by_cluster.csv")
)

# -------------------------------
# 2. Top 5 markers per cluster
# -------------------------------
log_step("Extracting top 5 markers per cluster")

top5_markers <- top_markers %>%
  group_by(.data[[cluster_col]]) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(.data[[cluster_col]]), desc(avg_log2FC))

write_csv(
  top5_markers,
  file.path(results_dir, "top5_markers_per_cluster.csv")
)

# -------------------------------
# 3. Gene signature helper
# -------------------------------
marker_reference <- list(
  podocyte = c("NPHS1", "NPHS2", "WT1", "PODXL", "SYNPO", "PLA2R1"),
  proximal_tubule = c("LRP2", "CUBN", "SLC34A1", "ALDOB", "FABP1", "GPX3"),
  distal_tubule = c("SLC12A3", "PVALB", "TRPM6", "FXYD4"),
  collecting_duct = c("AQP2", "KCNJ1", "SCNN1G", "SCNN1A", "AQP3"),
  endothelial = c("KDR", "EMCN", "ESAM", "PECAM1", "VWF"),
  mesangial = c("PDGFRB", "RGS5", "MCAM", "CSPG4", "TAGLN"),
  fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  immune_tcell = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R"),
  immune_bcell = c("MS4A1", "CD79A", "CD79B", "CD74"),
  macrophage_monocyte = c("LYZ", "CTSS", "TYROBP", "FCER1G", "SAT1"),
  proliferating = c("MKI67", "TOP2A", "STMN1", "UBE2C", "BIRC5"),
  stress_injury = c("LCN2", "HAVCR1", "KRT8", "KRT18", "KRT19", "SOD2", "GSTP1")
)

score_cluster_identity <- function(genes, marker_reference) {
  genes_upper <- unique(toupper(genes))
  
  scores <- sapply(names(marker_reference), function(label) {
    ref_genes <- toupper(marker_reference[[label]])
    sum(genes_upper %in% ref_genes)
  })
  
  best_label <- names(scores)[which.max(scores)]
  best_score <- max(scores)
  
  if (best_score == 0) {
    return(c(predicted_identity = "unassigned", identity_score = 0))
  }
  
  c(predicted_identity = best_label, identity_score = best_score)
}

# -------------------------------
# 4. Preliminary annotation
# -------------------------------
log_step("Generating preliminary cluster annotations")

cluster_annotation <- top5_markers %>%
  group_by(.data[[cluster_col]]) %>%
  summarise(
    top_genes = paste(gene, collapse = ", "),
    .groups = "drop"
  )

annotation_scores <- lapply(cluster_annotation$top_genes, function(x) {
  genes <- str_split(x, pattern = ",\\s*")[[1]]
  score_cluster_identity(genes, marker_reference)
})

annotation_scores_df <- do.call(rbind, annotation_scores) %>%
  as.data.frame(stringsAsFactors = FALSE)

cluster_annotation <- bind_cols(cluster_annotation, annotation_scores_df)

cluster_sizes[[cluster_col]] <- as.character(cluster_sizes[[cluster_col]])
marker_counts[[cluster_col]] <- as.character(marker_counts[[cluster_col]])
cluster_annotation[[cluster_col]] <- as.character(cluster_annotation[[cluster_col]])

cluster_summary <- cluster_sizes %>%
  left_join(marker_counts, by = cluster_col) %>%
  left_join(cluster_annotation, by = cluster_col) %>%
  arrange(as.numeric(.data[[cluster_col]]))

write_csv(
  cluster_summary,
  file.path(results_dir, "cluster_annotation_summary.csv")
)

log_step(sprintf("Clusters summarized: %d", nrow(cluster_summary)))
log_step("Saved marker_counts_by_cluster.csv")
log_step("Saved top5_markers_per_cluster.csv")
log_step("Saved cluster_annotation_summary.csv")
log_step("Cluster annotation summary completed successfully")