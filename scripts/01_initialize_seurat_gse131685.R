log_step <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_step("Starting Seurat project initialization for GSE131685")

required_packages <- c("Seurat", "dplyr", "ggplot2", "GEOquery", "Matrix")

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
  library(dplyr)
  library(ggplot2)
  library(GEOquery)
  library(Matrix)
})

raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")
results_tab_dir <- file.path("results", "tables")
temp_dir <- file.path("data", "temp_10x")

for (dir_path in c(raw_dir, processed_dir, results_tab_dir, temp_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

gse_id <- "GSE131685"
geo_dir <- file.path(raw_dir, gse_id)

log_step(sprintf("Downloading GEO supplementary files for %s", gse_id))
GEOquery::getGEOSuppFiles(GEO = gse_id, baseDir = raw_dir, makeDirectory = TRUE)

archive_files <- list.files(
  geo_dir,
  pattern = "\\.(tar|tar\\.gz|tgz)$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(archive_files) > 0) {
  for (archive_file in archive_files) {
    log_step(sprintf("Extracting archive: %s", basename(archive_file)))
    utils::untar(archive_file, exdir = geo_dir)
  }
} else {
  log_step("No archive files found after GEO download")
}

all_files <- list.files(geo_dir, recursive = TRUE, full.names = TRUE)

detected_files <- data.frame(
  file = all_files,
  stringsAsFactors = FALSE
)

write.csv(
  detected_files,
  file.path(results_tab_dir, "detected_geo_files.csv"),
  row.names = FALSE
)

log_step(sprintf("Total files detected after extraction: %d", length(all_files)))

matrix_files <- all_files[
  grepl("_matrix\\.mtx(\\.gz)?$", basename(all_files), ignore.case = TRUE)
]

if (length(matrix_files) == 0) {
  stop("No sample-specific matrix.mtx(.gz) files found in downloaded GEO directory.")
}

sample_prefixes <- sub(
  "_matrix\\.mtx(\\.gz)?$",
  "",
  basename(matrix_files),
  ignore.case = TRUE
)

sample_prefixes <- unique(sample_prefixes)

log_step(sprintf("Detected %d sample(s)", length(sample_prefixes)))
log_step(sprintf("Sample prefixes: %s", paste(sample_prefixes, collapse = ", ")))

find_single_file <- function(directory, pattern_text) {
  hits <- list.files(
    directory,
    pattern = pattern_text,
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(hits) == 1) {
    return(hits)
  }
  character(0)
}

safe_read_lines <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
  } else {
    con <- file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con)
}

read_mtx_triplet <- function(matrix_file, feature_file, barcode_file) {
  mat <- Matrix::readMM(matrix_file)
  
  features <- utils::read.delim(
    feature_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  barcodes <- utils::read.delim(
    barcode_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  if (ncol(features) >= 2) {
    gene_names <- features[[2]]
  } else {
    gene_names <- features[[1]]
  }
  
  gene_names <- make.unique(as.character(gene_names))
  cell_names <- make.unique(as.character(barcodes[[1]]))
  
  if (nrow(mat) != length(gene_names)) {
    stop(sprintf(
      "Row mismatch: matrix has %d rows but features file has %d rows",
      nrow(mat), length(gene_names)
    ))
  }
  
  if (ncol(mat) != length(cell_names)) {
    stop(sprintf(
      "Column mismatch: matrix has %d columns but barcodes file has %d rows",
      ncol(mat), length(cell_names)
    ))
  }
  
  rownames(mat) <- gene_names
  colnames(mat) <- cell_names
  mat
}

seurat_list <- list()
sample_summary_list <- list()
failed_samples <- character(0)

for (sample_id in sample_prefixes) {
  log_step(sprintf("Processing sample: %s", sample_id))
  
  matrix_file <- find_single_file(
    geo_dir,
    paste0("^", sample_id, "_matrix\\.mtx(\\.gz)?$")
  )
  
  barcode_file <- find_single_file(
    geo_dir,
    paste0("^", sample_id, "_barcodes\\.tsv(\\.gz)?$")
  )
  
  feature_file <- find_single_file(
    geo_dir,
    paste0("^", sample_id, "_features\\.tsv(\\.gz)?$|^", sample_id, "_genes\\.tsv(\\.gz)?$")
  )
  
  if (length(matrix_file) != 1 || length(barcode_file) != 1 || length(feature_file) != 1) {
    log_step(sprintf(
      "Skipping sample %s because required files were not uniquely identified",
      sample_id
    ))
    failed_samples <- c(failed_samples, sample_id)
    next
  }
  
  sample_counts <- tryCatch(
    {
      read_mtx_triplet(matrix_file, feature_file, barcode_file)
    },
    error = function(e) {
      log_step(sprintf("Matrix import failed for %s: %s", sample_id, e$message))
      return(NULL)
    }
  )
  
  if (is.null(sample_counts)) {
    failed_samples <- c(failed_samples, sample_id)
    next
  }
  
  sample_obj <- tryCatch(
    {
      CreateSeuratObject(
        counts = sample_counts,
        project = gse_id,
        min.cells = 3,
        min.features = 200
      )
    },
    error = function(e) {
      log_step(sprintf("CreateSeuratObject failed for %s: %s", sample_id, e$message))
      return(NULL)
    }
  )
  
  if (is.null(sample_obj)) {
    failed_samples <- c(failed_samples, sample_id)
    next
  }
  
  sample_obj$sample_id <- sample_id
  counts_dim <- dim(sample_counts)
  
  log_step(sprintf(
    "Sample %s loaded successfully: %d genes x %d cells",
    sample_id,
    counts_dim[1],
    counts_dim[2]
  ))
  
  seurat_list[[sample_id]] <- sample_obj
  sample_summary_list[[sample_id]] <- data.frame(
    sample_id = sample_id,
    n_genes = counts_dim[1],
    n_cells = counts_dim[2],
    stringsAsFactors = FALSE
  )
}

if (length(seurat_list) == 0) {
  stop("No valid samples could be loaded into Seurat objects.")
}

log_step(sprintf("Successfully loaded %d sample(s)", length(seurat_list)))

if (length(failed_samples) > 0) {
  log_step(sprintf(
    "Failed or skipped samples: %s",
    paste(unique(failed_samples), collapse = ", ")
  ))
}

sample_summary_df <- do.call(rbind, sample_summary_list)

write.csv(
  sample_summary_df,
  file.path(results_tab_dir, "sample_level_qc_summary.csv"),
  row.names = FALSE
)

if (length(seurat_list) == 1) {
  seurat_obj <- seurat_list[[1]]
  log_step("Only one sample loaded. Merge step not required.")
} else {
  seurat_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(seurat_list)],
    add.cell.ids = names(seurat_list),
    project = gse_id
  )
  log_step("Merged individual Seurat objects successfully")
}

metadata <- seurat_obj@meta.data

write.csv(
  metadata,
  file.path(results_tab_dir, "metadata_preview.csv"),
  row.names = TRUE
)

capture.output(
  str(metadata),
  file = file.path(results_tab_dir, "metadata_structure.txt")
)

summary_df <- data.frame(
  dataset = gse_id,
  n_genes = nrow(seurat_obj),
  n_cells = ncol(seurat_obj),
  n_metadata_columns = ncol(metadata),
  metadata_columns = paste(colnames(metadata), collapse = ";"),
  generated_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  file.path(results_tab_dir, "initial_qc_summary.csv"),
  row.names = FALSE
)

saveRDS(
  seurat_obj,
  file.path(processed_dir, "gse131685_seurat_raw.rds")
)

log_step("Saved merged Seurat object to data/processed/gse131685_seurat_raw.rds")
log_step("Saved initial QC summary to results/tables/initial_qc_summary.csv")
log_step("Seurat project initialization completed successfully")