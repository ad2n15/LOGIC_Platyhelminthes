# This code was developed to assess the expression of 35 genes in schistosoma mansoni,
#these genes are potential ligand-gated ion channel conservation across Platyhelminthes 

#######################################################################################
# Add cell types as meta data to the seurat object from downloaded from GSE146736
# clusters to cell type map file were generated manually using  PMID: 32973030 and  https://www.collinslab.org/schistocyte/
#######################################################################################

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load the Seurat object
schistodata <- readRDS("GSE146736_adult_scseq_seurat.rds")



# Check if the default assay is set to "RNA" or "integrated"
if (DefaultAssay(schistodata) != "RNA") {
  cat("Setting the default assay to 'RNA'\n")
  DefaultAssay(schistodata) <- "RNA"
}



# Count cells before adding metadata
num_cells_before <- nrow(schistodata@meta.data)
cat("Number of cells before adding metadata:", num_cells_before, "\n")

# Load the cell type mapping file
celltypes_map <- read.table("celltypes_cluster_res5_map_PMID_32973030.txt", header = TRUE)

# Ensure the 'integrated_snn_res.5' column is of the same type
schistodata@meta.data$integrated_snn_res.5 <- as.factor(schistodata@meta.data$integrated_snn_res.5)
celltypes_map$integrated_snn_res.5 <- as.factor(celltypes_map$integrated_snn_res.5)

# Perform left join to add 'cell_type_1' and 'cell_type_2' to the Seurat object metadata
schistodata@meta.data <- schistodata@meta.data %>%
  left_join(celltypes_map, by = c("integrated_snn_res.5" = "integrated_snn_res.5"))

# Count cells after adding metadata
num_cells_after <- nrow(schistodata@meta.data)
cat("Number of cells after adding metadata:", num_cells_after, "\n")

# Validate if the number of cells remained the same
if (num_cells_before == num_cells_after) {
  cat("Validation passed: Number of cells remained the same.\n")
} else {
  cat("Warning: Number of cells changed after adding metadata!\n")
}

# Check the added metadata columns
head(schistodata@meta.data[c("cell_type_1", "cell_type_2")])

# Save the updated Seurat object with added metadata
saveRDS(schistodata, "GSE146736_adult_scseq_seurat_with_celltypes.rds")


#################

##########################################################################################################
# Compare the expression LIGC genes (n=35) between each cell type against other cell types
# repeat this analysis for both original cell types (cell_type_1) and combined cell types (cell_type_2) 
##########################################################################################################


# Check if the default assay is set to "RNA" or "integrated"
if (DefaultAssay(schistodata) != "RNA") {
  cat("Setting the default assay to 'RNA'\n")
  DefaultAssay(schistodata) <- "RNA"
}






# Define a function to read gene lists from a file
read_gene_list <- function(file_path) {
  genes <- readLines(file_path)  # Read the file and return the gene names
  return(genes)
}

# Ensure the data is normalized before running FindAllMarkers
schistodata <- NormalizeData(schistodata, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
schistodata <- FindVariableFeatures(schistodata, selection.method = "vst", nfeatures = 2000)

# Scale the data
schistodata <- ScaleData(schistodata)

# List of gene list files
gene_list_files <- list.files(pattern = "*_list.txt")

# Function to perform FindAllMarkers for a given cell type
run_FindAllMarkers <- function(cell_type_column) {
  Idents(object = schistodata) <- schistodata[[cell_type_column]]  # Set identity class to the specified cell type column
  DefaultAssay(schistodata) <- "RNA"  # Ensure the correct assay is used

  for (gene_file in gene_list_files) {
    # Read the gene list
    genes_of_interest <- read_gene_list(gene_file)

    # Perform differential expression analysis
    markers <- FindAllMarkers(
      object = schistodata,
      features = genes_of_interest,
      only.pos = FALSE,  # Finds both upregulated and downregulated genes
      group.by = cell_type_column
    )

    # Create a name for the output file
    output_file <- paste0(sub(".txt", paste0("_all_", cell_type_column, "_markers.csv"), gene_file))

    # Save results to a CSV file
    write.csv(markers, file = output_file)

    # Print a message indicating successful save
    cat("Saved results for", gene_file, "using", cell_type_column, "to", output_file, "\n")
  }
}

# Run FindAllMarkers for both cell type classifications
run_FindAllMarkers("cell_type_2")
run_FindAllMarkers("cell_type_1")

###################################################

##########################################################################################################
# Compare the expression LIGC genes (n=35) between each cell type against other cell types in each sex group (Female, Male, and IM)
# repeat this analysis for both original cell types (cell_type_1) and combined cell types (cell_type_2)
##########################################################################################################


# Ensure Idents are set before subsetting
Idents(schistodata) <- schistodata$Group

# Function to perform FindAllMarkers for a given cell type within a specific group
run_FindAllMarkers_for_group <- function(cell_type_column, group) {
  cat("Processing group:", group, "\n")

  # Check if the group exists
  available_groups <- unique(schistodata$Group)
  if (!(group %in% available_groups)) {
    cat("Skipping group:", group, "because it is not in available groups.\n")
    return(NULL)
  }

  # Subset based on Idents (which is now Group)
  group_data <- subset(schistodata, idents = group)

  # Check if the subset is empty
  if (ncol(group_data) == 0) {
    cat("Skipping group:", group, "because it contains no cells.\n")
    return(NULL)
  }

  # Transfer metadata for cell_type_1 and cell_type_2 from the original object
  group_data$cell_type_1 <- schistodata$cell_type_1[colnames(group_data)]
  group_data$cell_type_2 <- schistodata$cell_type_2[colnames(group_data)]

  # Normalize and process
  group_data <- NormalizeData(group_data, normalization.method = "LogNormalize", scale.factor = 10000)
  group_data <- FindVariableFeatures(group_data, selection.method = "vst", nfeatures = 2000)
  group_data <- ScaleData(group_data)

  # Check if cell type column exists
  if (!(cell_type_column %in% colnames(group_data@meta.data))) {
    cat("Skipping", cell_type_column, "because it is not found in metadata.\n")
    return(NULL)
  }

  # Set cell type as identity
  Idents(group_data) <- group_data[[cell_type_column]]
  DefaultAssay(group_data) <- "RNA"

  # Loop through gene files
  for (gene_file in gene_list_files) {
    genes_of_interest <- read_gene_list(gene_file)
    table(group_data$cell_type_2)
    # Run FindAllMarkers
    markers <- FindAllMarkers(
      object = group_data,
      features = genes_of_interest,
      only.pos = FALSE,
      group.by = cell_type_column
    )

    if (nrow(markers) == 0) {
      cat("No markers found for", group, "and", cell_type_column, "\n")
      next
    }

    # Save results
    output_file <- paste0(sub(".txt", paste0("_", group, "_", cell_type_column, "_markers.csv"), gene_file))
    write.csv(markers, file = output_file)
    cat("Saved results for", gene_file, "using", group, "and", cell_type_column, "to", output_file, "\n")
  }
}

# Run for each group
for (group in unique(schistodata$Group)) {
  run_FindAllMarkers_for_group("cell_type_2", group)
  run_FindAllMarkers_for_group("cell_type_1", group)
}

#################

##########################################################################################################
# Read the expression invistigation files, adjust p value for multible testing, and generated a final combined file
##########################################################################################################


# Define the LGIC groups, sex, and levels
LGIC_groups <- c("ASIC", "P2X", "cysloop_gaba", "iGluR", "cysloop_acetylcholine")
sex_groups <- c("all", "Female", "Male", "IM")
levels <- c(1, 2)

# Create an empty data frame to store combined results
combined_results <- data.frame()

# Loop over each combination of LGIC group, sex, and level
for (lgic in LGIC_groups) {
  for (sex in sex_groups) {
    for (level in levels) {
      # Construct the file name based on the variables
      file_name <- paste0(lgic, "_Smansoni_list_", sex, "_cell_type_", level, "_markers.csv")
      print(file_name)

      # Check if the file exists and is not empty
      if (file.exists(file_name) && file.info(file_name)$size > 0) {
        # Try reading the file and catch any errors
        tryCatch({
          # Read the CSV file
          df <- read.csv(file_name)

          # Perform Bonferroni correction on the p-value column
          df$p_val_adj_all <- p.adjust(df$p_val, method = "bonferroni")

          # Add columns for LGIC group, sex, and level
          df$LGIC_group <- lgic
          df$sex <- sex
          df$level <- level

          # Combine the data into the main data frame
          combined_results <- bind_rows(combined_results, df)

        }, error = function(e) {
          cat("Error reading file:", file_name, "\n", "Error message:", e$message, "\n")
        })
      } else {
        cat("File not found or is empty:", file_name, "\n")
      }
    }
  }
}


# Remove the '.' and the number after it in the 'X' column
combined_results$X <- sub("\\.\\d+$", "", combined_results$X)


# Write the combined and adjusted results to a new CSV file
write.csv(combined_results, "combined_results_with_bonferroni.csv", row.names = FALSE)

# Optional: Check the first few rows of the combined results
head(combined_results)

