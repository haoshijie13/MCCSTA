# 1. RandomForest classicifer

## 1.1 Traning and test data process

~~~R
#As described in the Methods section, our process for identifying subclusters is as follows: we use a relatively large resolution to generate an overcluster state. Then, for each pair of clusters, we calculate the Jaccard similarity of the top 50 marker genes between them, and use the expression profiles of these clusters to train a RandomForest classifier. Subsequently, we merge the clusters based on the values of the Jaccard similarity and the predictions of the RandomForest classifier, taking into account the magnitudes of these values.

# RandomForest Classifier

# Suppress package loading messages to keep the console clean
suppressMessages(library(tidyverse))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(getopt))
suppressMessages(library(Seurat))

# Define command - line argument specifications
spec <- matrix(c(
  'dir1', 'D1', 2, "character",
  'dir2', 'D2', 2, "character",
  'cluster_name', 'CN', 2, "character",
  'label', 'L', 2, "character",
  'save_path', 'S', 2, "character"
), byrow = TRUE, ncol = 4)

# Parse command - line arguments
opt <- getopt(spec = spec)
dir1 <- opt$dir1
dir2 <- opt$dir2
cluster_name <- opt$cluster_name
label <- opt$label
save_path <- opt$save_path

# Set a seed for reproducibility
set.seed(1234)

# Load the Seurat object from the specified directory
seurat <- readRDS(dir1)
# Extract the meta - data from the Seurat object
meta <- seurat@meta.data

# Add a 'CellID' column if it doesn't exist
if (!"CellID" %in% colnames(meta)) {
  meta$CellID <- rownames(meta)
}

# Create a 'group' column based on the specified cluster name
meta["group"] <- meta[[cluster_name]]

# Filter out groups with less than 200 cells
group_counts <- table(meta$group)
groups_to_keep <- names(group_counts[group_counts >= 200])
sub_meta <- subset(meta, subset = group %in% groups_to_keep)

# Remove levels with zero count in the 'group' factor
freq_table <- table(meta$group)
levels_with_zero_count <- names(freq_table[freq_table == 0])
print(levels_with_zero_count)
meta$group <- droplevels(meta$group, exclude = levels_with_zero_count)

# Start the down - sampling process
print("start downsample")
sub_meta$group <- as.character(sub_meta$group)

# Split the CellIDs by group
tmp <- sub_meta %>% select(group, CellID)
a <- split(tmp$CellID, tmp$group)

# Perform down - sampling for each group
samples <- data.frame()
for (group in names(a)) {
  cells <- a[[group]]
  if (length(cells) >= 200) {
    CellID <- sample(cells, size = 200, replace = FALSE)
    group_column <- rep(group, 200)
    samples <- rbind(samples, cbind(group = group_column, CellID))
  }
}

# Remove duplicate rows
samples1 <- unique(samples)
rownames(samples) <- samples$CellID

# Mark selected cells in the meta - data
meta$select <- "no"
meta[rownames(samples), "select"] <- "yes"

print("Downsample finished!")

# Update the Seurat object's meta - data
seurat@meta.data <- meta
# Subset the Seurat object to include only selected cells
sub <- subset(seurat, subset = select == "yes")
# Extract the SCT data matrix
counts <- sub@assays$SCT@data

# Load marker genes
markers <- readRDS(dir2)
# Select the top 50 marker genes per cluster based on avg_log2FC
n <- 50
topmarker <- markers %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC)

# Extract the gene names from the top marker genes
gene_list <- topmarker$gene

# Convert the counts matrix to a data frame and subset by the gene list
counts <- as.data.frame(counts)
sct_matrix_subset <- counts[gene_list, ]
# Transpose the subset matrix
sct_matrix_subset <- t(sct_matrix_subset)
# Add the 'group' column from the meta - data
test <- sub@meta.data
sct_matrix_subset <- as.data.frame(sct_matrix_subset)
sct_matrix_subset$group <- test$group

# Split the data into training and test sets
library(caTools)
data <- sct_matrix_subset
data$group <- as.factor(data$group)
split <- sample.split(data$group, SplitRatio = 0.8)
data_train <- data[split, ]
data_test <- data[!split, ]

# Save the training and test sets to an R data file
save(data_train, data_test, file = paste0(save_path, label, "_predict.data.rda"))
~~~

## 1.2 model train

~~~R
# Suppress package loading messages to keep the console clean
suppressMessages(library(randomForest))
suppressMessages(library(caret))
suppressMessages(library(e1071))
suppressMessages(library(dplyr))
suppressMessages(library(getopt))

# Define command - line argument specifications
spec <- matrix(c(
  'dir1', 'D1', 2, "character",
  'label', 'L', 2, "character",
  'save_path', 'S', 2, "character"
), byrow = TRUE, ncol = 4)

# Parse command - line arguments
opt <- getopt(spec = spec)
dir1 <- opt$dir1
label <- opt$label
save_path <- opt$save_path

# Load the data from the specified directory
load(dir1)

# Display the first few values of the 'group' variable in the training and test sets
cat("First few values of 'group' in data_train:\n")
head(data_train$group)
cat("First few values of 'group' in data_test:\n")
head(data_test$group)

# Set the name for saving files
name <- label

# Use the training data and remove rows with NA values
dat <- data_train
dat <- dat[complete.cases(dat), ]
cat(paste0("Dimensions of data after removing NA rows: ", dim(dat), "\n"))

# Remove unused levels from the 'group' factor
dat$group <- droplevels(dat$group)

# Separate the features and the target variable
x <- dat
x$group <- NULL
y <- dat$group

# Calculate the 'mtry' parameter for the random forest
mtry <- floor(sqrt(ncol(x)))
# Set the number of trees in the random forest
ntree <- 1000

# Print the values of 'mtry' and 'ntree'
cat(paste0(" mtry: ", mtry, "\n"))
cat(paste0("ntree: ", ntree, "\n"))

# Replace hyphens in column names with underscores to avoid issues
colnames(dat) <- gsub("-", "_", colnames(dat))

# Train a random forest model
dat.rf <- randomForest(group ~ ., data = dat, mtry = mtry, ntree = ntree, 
                       keep.forest = TRUE, importance = TRUE, prox = TRUE)

# Save the trained random forest model
saveRDS(dat.rf, file = paste0(save_path, name, ".randomForest.rds"))

# Extract variable importance from the random forest model
dat.imp <- importance(dat.rf)

# Save the variable importance plot as a PDF
pdf(paste0(save_path, name, ".var.importance.2.pdf"), 
    height = ceiling(ncol(dat) / 8) + 2, width = 10)
varImpPlot(dat.rf) # Plot the importance of each variable
dev.off()

# Save the margin plot and the error rate plot as a PDF
pdf(paste0(save_path, name, ".point&tree.error.2.pdf"), height = 6, width = 6)
plot(margin(dat.rf)) # Plot the margin of each point
plot(dat.rf, log = "y") # Plot the error rate for different numbers of trees
dev.off()

# The following line is commented out as it's not clear about its usage and may need further adjustment
# MDSplot(dat.rf, dat$group, palette = rep(1, 3), pch = as.numeric(dat$group)) 
~~~



# 2. Jaccard similarity calculate

~~~R
# Calculate Jaccard similarity
# Suppress package loading messages to keep the console clean
suppressMessages(library(tidyverse))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(getopt))

# Define command - line argument specifications
spec <- matrix(c(
  'dir1_markers', 'M', 2, "character",
  'label', 'L', 2, "character",
  'save_path', 'S', 2, "character"
), byrow = TRUE, ncol = 4)

# Parse command - line arguments
opt <- getopt(spec = spec)
dir1_markers <- opt$dir1_markers
label <- opt$label
save_path <- opt$save_path

# Read all marker genes data
markers <- readRDS(dir1_markers)
# Select the top 50 marker genes per cluster based on avg_log2FC
n <- 50
topmarker <- markers %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC)

# Extract relevant columns for further analysis
df <- topmarker[c("cluster", "gene")]
# Create a 'group' column which is a copy of 'cluster'
df$group <- df$cluster

# Define a function to calculate Jaccard similarity between two groups
calculate_jaccard_similarity <- function(group1, group2) {
  # Get the genes for group1
  genes_group1 <- df %>%
    filter(group == group1) %>%
    pull(gene)
  # Get the genes for group2
  genes_group2 <- df %>%
    filter(group == group2) %>%
    pull(gene)
  
  # Calculate the size of the intersection of the two gene sets
  intersection <- length(intersect(genes_group1, genes_group2))
  # Calculate the size of the union of the two gene sets
  union <- length(union(genes_group1, genes_group2))
  
  # Calculate the Jaccard similarity
  jaccard_similarity <- intersection / union
  return(jaccard_similarity)
}

# Get the unique groups in the data
unique_groups <- unique(df$group)

# Initialize a data frame to store Jaccard similarity results
jaccard_data <- data.frame(Group1 = character(), Group2 = character(), Jaccard = double())

# Loop through all pairs of groups to calculate Jaccard similarity
for (i in seq_along(unique_groups)) {
  for (j in seq_along(unique_groups)) {
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]
    
    # Calculate Jaccard similarity for the current pair of groups
    jaccard_similarity <- calculate_jaccard_similarity(group1, group2)
    
    # Store the result in the data frame
    jaccard_data <- rbind(jaccard_data, data.frame(Group1 = group1, Group2 = group2, Jaccard = jaccard_similarity))
  }
}

# Reshape the Jaccard similarity data into a matrix
jaccard_matrix <- reshape2::dcast(jaccard_data, Group1 ~ Group2, value.var = "Jaccard", fun.aggregate = max)

# Save the Jaccard similarity matrix as a CSV file
write.table(jaccard_matrix, file = paste0(save_path, label, "_jaccard_matrix.csv"))

# Set row names of the matrix to the group names
rownames(jaccard_matrix) <- jaccard_matrix$Group1
# Remove the redundant 'Group1' column
jaccard_matrix$Group1 <- NULL

# Calculate the number of unique groups
len <- length(unique_groups)
# Define the path for the PDF heatmap file
pdf_dir <- paste0(save_path, label, "_heatmap.pdf")
# Open a PDF file for saving the heatmap
pdf(pdf_dir, width = 0.2 * len, height = 0.2 * len)
# Plot the heatmap of the Jaccard similarity matrix
print(pheatmap(jaccard_matrix, cluster_col = TRUE, cluster_row = TRUE, color = colorRampPalette(c("white", "blue"))(256)))
# Close the PDF file
dev.off()
~~~



