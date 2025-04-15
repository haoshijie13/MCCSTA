# 1. snRNA matrix load and QC

~~~R
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(tidyverse)

# Set the library matrix directory
current_dir <- "./Marmoset_cortex_snRNA_library/"

# Read the sample names from a file
tmp <- read.table("./sampleNames.file.txt", header = FALSE)
filename <- tmp$V1

# Set the quality control parameters
nFeature_RNA_min <- 1000
nFeature_RNA_max <- 6000
mt_ratio_cutoff <- 5
nCount_nFeature_ratio_cutoff <- 1.2

# Define mitochondrial genes
MTGenes <- c("ND6", "COX3", "COX1", "ND5", "ND4", "ND2", "ND4L", "ATP8", "CYTB", "COX2", "ND3", "ATP6", "ND1")

# Initialize an empty list to store Seurat objects
result <- list()
# Initialize an empty data frame to store filtering statistics
Filters <- data.frame()

# Loop through each sample
for (i in seq_along(filename)) {
  # Construct the input directory path
  input_dir <- file.path(current_dir, filename[i])
  cat("Processing sample:", filename[i], "\n")
  
  # Read the expression data
  expr <- Read10X(data.dir = input_dir, gene.column = 1, cell.column = 1)
  cat("Dimensions of expression matrix:", dim(expr), "\n")
  
  # Add sample name to cell barcodes
  colnames(expr) <- paste(colnames(expr), filename[i], sep = "_")
  
  # Create a Seurat object
  seurat <- CreateSeuratObject(expr, project = filename[i])
  
  # Record the number of raw cells
  raw_cells <- ncol(seurat)
  
  # Calculate the percentage of mitochondrial genes
  seurat <- PercentageFeatureSet(seurat, features = MTGenes, col.name = "percent.mt")
  
  # Calculate the ratio of nCount_RNA to nFeature_RNA
  seurat$nCount_nFeature_ratio <- seurat$nCount_RNA / seurat$nFeature_RNA
  
  # Perform quality control filtering
  seurat <- subset(seurat, subset = nFeature_RNA >= nFeature_RNA_min & 
                     percent.mt <= mt_ratio_cutoff & 
                     nCount_nFeature_ratio >= nCount_nFeature_ratio_cutoff)
  
  # Record the number of cells after QC
  qc_cells <- ncol(seurat)
  
  # Store the filtering statistics
  Filters <- rbind(Filters, c(filename[i], nrow(expr), raw_cells, qc_cells))
  
  # Store the Seurat object in the result list
  result[[i]] <- seurat
}

# Merge all Seurat objects
merge_data <- merge(result[[1]], y = result[2:length(result)])
cat("Dimensions of merged data:", dim(merge_data), "\n")

# Set column names for the Filters data frame
colnames(Filters) <- c("SampleName", "NGenes", "Raw.NCells", "QC.NCells")

# Add the total statistics
Filters <- rbind(Filters, c("Total", nrow(merge_data), sum(as.numeric(Filters$Raw.NCells)), ncol(merge_data)))

# Write the filtering statistics to a file
write.table(Filters, file = "Filters.cell.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Create a directory for results
dir.create("./result", showWarnings = FALSE)

# Set the merged data as the main Seurat object
seurat <- merge_data

# Generate violin plots and save them to a PDF file
pdf("./result/seurat.nFeature_RNA.orig.ident.qc.pdf", height = 6, width = 20)
VlnPlot(object = seurat, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + NoLegend()
VlnPlot(object = seurat, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + NoLegend()
VlnPlot(object = seurat, features = c("percent.mt"), group.by = "orig.ident", ncol = 1, pt.size = 0) + NoLegend()
dev.off()
saveRDS(seurat,file="./data_file/Marmoset_QC.rds")
~~~



# 2. Clustering and annotation

~~~R
# Set the default assay to RNA
source("./snRNA_functions.r")
seurat <- readRDS("./data_file/Marmoset_QC.rds")
DefaultAssay(seurat) <- "RNA"

# Perform SCTransform normalization and regress out mitochondrial percentage
seurat <- SCTransform(seurat, assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE)

# Check the default assay
print(DefaultAssay(seurat))

# Perform PCA with 50 principal components
seurat <- RunPCA(seurat, npcs = 50)

# Output PCA plots to a PDF file
output <- "./result/seurat.pca.pdf"
pdf(output, height = 10, width = 10)
DimPlot(seurat, reduction = "pca", group.by = "orig.ident", label = TRUE, label.size = 6, pt.size = 0.5) + coord_fixed(ratio = 1)
DimPlot(seurat, reduction = "pca", group.by = "brain_area", label = TRUE, label.size = 6, pt.size = 0.5) + coord_fixed(ratio = 1)
VizDimLoadings(seurat, dims = 1:4, reduction = "pca")
DimHeatmap(seurat, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Calculate the standard deviations and cumulative variance
stdev <- attributes(seurat[['pca']])$stdev
stdev_square <- stdev^2  # Vectorized operation for squaring
total <- cumsum(stdev_square) / sum(stdev_square)

# Create a data frame for cumulative variance plot
dat <- data.frame(value = total, pc = seq_along(total))

# Find the positions where cumulative variance reaches 85%, 90%, and 95%
Pos1 <- find_pos(0.85)
Pos2 <- find_pos(0.9)
Pos3 <- find_pos(0.95)

# Output elbow plot and cumulative variance plot to a PDF file
pdf("./result/seurat.PCA.elbow&cumulative.pdf", height = 6, width = 6)
ElbowPlot(seurat, ndims = 50)
ggplot(dat, aes(x = pc, y = value)) +
  geom_line() +
  xlab("Number of PC") +
  ylab("Cumulative Rate") +
  geom_vline(xintercept = c(Pos1, Pos2, Pos3), colour = "red", linetype = "dotted") +
  geom_hline(yintercept = c(0.85, 0.9, 0.95), colour = "red", linetype = "dotted") +
  theme_classic()
dev.off()

# Print the default assay and the positions
print(DefaultAssay(seurat))
print(c(Pos1, Pos2, Pos3))

# Define parameters for clustering and visualization
PCs <- c(Pos1, Pos2, Pos3)
resolutions <- c(0.3, 0.5)
K <- c(10, 15, 20)
name <- "Marmo"

# Call the function to perform clustering and visualization
seurat <- cluster_and_visualize(seurat, PCs, resolutions, K, name)
print("All done !")

# Select the final parameter combination
PC <- 43
resolution <- 0.3
k <- 15
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PC, k.param = k)
seurat <- FindClusters(object = seurat, resolution = resolution, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:PC)
print(levels(seurat))
~~~



~~~R
# ########### Determine quality control after classification:
# Remove unnecessary columns from the Seurat object
seurat$X <- NULL
seurat$Y <- NULL
seurat$CellID <- NULL

# Extract UMAP embeddings from the Seurat object
dat <- as.data.frame(seurat@reductions$umap@cell.embeddings)
# If you want to use t-SNE embeddings instead, uncomment the following line
# dat <- as.data.frame(seurat@reductions$tsne@cell.embeddings) 
name <- "all.umap"
names(dat) <- c("X", "Y")
dat$CellID <- rownames(dat)

# Combine metadata and embeddings
lab <- cbind(seurat@meta.data, dat[rownames(seurat@meta.data),])
lab$group <- lab$seurat_clusters

# Create a directory for quality control results
dir.create("QC", showWarnings = FALSE)

# Print information about the groups
cat(paste0("The name of group: ", name, "\n"))
cat(paste0("The levels of group: ", length(unique(lab$group)), "\n"))
print(table(lab$group))

# Calculate the number of rows/columns for facet wrapping
len <- ceiling(sqrt(length(unique(lab$group))))

# Plot quality control plots
plot_quality_control(lab, name, len)

# Split the data by group
tmp <- split(lab, lab$group)
# Print the number of cells in each group
print(sapply(tmp, nrow))

# Identify outliers
result <- identify_outliers(tmp)
re <- result$re
outlier_number <- result$outlier_number
total_outliers <- sum(outlier_number)
re$Sub[which(is.na(re$Sub))] <- "sub"
print(table(re$Sub))

# Plot outlier-related plots
plot_outliers(re, len, total_outliers)

# Remove clusters and cells based on criteria
dat <- remove_clusters_and_cells(re)

# Plot the final cells in each cluster
plot_final_cells(dat, len)

# Add a QC column to the original data
lab$QC <- "no"
lab[rownames(dat), "QC"] <- "yes"
print(table(lab$QC))

seurat$QC <- lab$QC
seurat <- subset(seurat,subset=lab=="yes")
~~~



~~~R
# Class and subclass anno

# Define marker genes for different cell classes
class_markerGenes  <- c(
    "SLC17A6", "SLC17A7", "CUX1", "CUX2",  # Markers for excitatory neurons (EX)
    "GAD1", "GAD2", "RELN", # Markers for inhibitory neurons (IN)
    "SLC1A2", "SLC1A3", "GFAP", "AQP4", # Markers for astrocytes (AST)
    "PDGFRA", # Marker for oligodendrocyte precursor cells (OPC)
    "CLDN11", "MOBP", "MBP", "MOG", # Markers for oligodendrocytes (OLG)
    "RGS5", # Marker for endothelial cells (EC)
    "SLC6A12", "SLC6A13", "COL1A1", "COL1A2", "CEMIP", # Markers for vascular smooth muscle-like cells (VLMC)
    "ITGAM", "P2RY12" # Markers for microglia (MG)
)

# Open a PDF file to save the dot plot
pdf("./plot_file/all_clutser.pdf", width = 218 * 0.2, height = 10)
# Create a dot plot to show the expression of marker genes in each cluster
# group.by="v2_Subcluster" indicates grouping by the v2_Subcluster column
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) is used to rotate the x-axis labels to avoid overlap
# scale_color_gradientn(colours = c('#66CC66', '#FFCC33')) sets the color gradient
# coord_flip() flips the axes
DotPlot(seurat, features = class_markerGenes, group.by = "v2_Subcluster") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_gradientn(colours = c('#66CC66', '#FFCC33')) +
    coord_flip()
# Close the PDF file
dev.off()

# Define new cluster names
new.cluster.ids <- c("EX", "IN", "AST", "OPC", "OLG", "VLMC", "MG")
# Ensure that the new cluster names correspond to the current cluster levels one by one
# The original code names(seurat) <- levels(seurat) was incorrect. It might have intended to assign the cluster levels
# as the names of the new cluster names vector. The following code corrects it.
names(new.cluster.ids) <- levels(seurat)
# Use the RenameIdents function to rename the clusters in the Seurat object
seurat <- RenameIdents(seurat, new.cluster.ids)
# Plot the UMAP plot to show the renamed cluster results
# reduction = "umap" indicates using the UMAP dimensionality reduction results
# label = TRUE indicates labeling the cluster names on the plot
# pt.size = 0.5 sets the size of the points
# NoLegend() removes the legend
DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# where subclass marker 
# We used a similar process to perform dimensionality reduction and clustering on the subclasses, and annotated the unsupervised classes based on the marker genes shown below.
markerGenes  <- c(
    "SLC17A6","SLC17A7", #EX
    "GFAP",#L1
    "HPCAL1",#L2
    "NECAB1","MEF2C",#L2-L4
    "GPR83","CUX1","CUX2","RASGRF1","PCDH8",#L2-L3/#L3
    #"RORB",#L3-L4 "RORB",#L3-L6"RORB",
	"PCP4",#L4
    #"RORB",#L4-L5
    "FOXP2",#L4-L6
    "FEZF2",#L5
    "ETV1","NEFM","NEFH","HTR2C",#L5-L6
    "TLE4","SEMA3E","NR4A2","OPRK1","SYT6",#L6
    "GAD1","GAD2", #IN
    "SOX6","LHX6","PVALB","SNCG","SST",#MGE
    "LAMP5","RELN","VIP","ADARB2",#CGE
    "UNC5B",
    "SLC1A2","SLC1A3","GFAP","AQP4",#AST
    "PDGFRA",#OPC
    "CLDN11","MOBP","MBP","MOG", #OLG
    "RGS5",#EC
    "SLC6A12","SLC6A13","COL1A1","COL1A2","CEMIP",#VLMC
    "ITGAM","P2RY12",#MG
    "ACTA2",#SMC
    "CLDN5",#Endo
    "PTPRC"#Immune
  )

~~~