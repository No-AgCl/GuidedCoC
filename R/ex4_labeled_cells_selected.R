library(R.matlab)
library(data.table)

# Path settings
base_dir <- "data/processed data"
input_dir <- "data/processed data"
output_dir <- "data/graph data"

# Read barcode lists
all_barcodes <- fread(file.path(base_dir, "ex4_target_data_barcodes.txt"), header = FALSE)[[1]]
valid_barcodes <- fread(file.path(base_dir, "ex4_labeled_cells_barcodes.txt"), header = FALSE)[[1]]

# Ensure consistent order and extract row indices
row_idx <- which(all_barcodes %in% valid_barcodes)
cat("Extracting", length(row_idx), "rows for labeled cells\n")

# Define processing function
process_and_save <- function(mat_path, row_idx, var_name, out_filename) {
  cat(" Processing:", mat_path, "\n")
  
  # Load .mat file
  mat_list <- readMat(mat_path)
  
  if (!(var_name %in% names(mat_list))) {
    stop(paste("Variable", var_name, "not found in", mat_path))
  }
  
  mat <- mat_list[[var_name]]
  
  # Subset rows
  mat_labeled <- mat[row_idx, , drop = FALSE]
  
  # Construct output path
  out_path <- file.path(output_dir, out_filename)
  
  # Save with original variable name
  args <- list()
  args[[var_name]] <- mat_labeled
  do.call(writeMat, c(list(con = out_path), args))
  
  cat("Saved to:", out_path, "\n\n")
}

# Run processing
process_and_save(file.path(input_dir, "ex4_target_data_RNA.mat"), row_idx, var_name = "q", out_filename = "ex4_labeled_target_data_RNA.mat")
process_and_save(file.path(input_dir, "ex4_target_data_ATAC.mat"), row_idx, var_name = "q0", out_filename = "ex4_labeled_target_data_ATAC.mat")
