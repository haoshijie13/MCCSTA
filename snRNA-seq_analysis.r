# Load necessary libraries
# Seurat is used for single-cell RNA sequencing data analysis
library(Seurat)
# pheatmap is used to create heatmaps
library(pheatmap)
# ComplexHeatmap provides more advanced heatmap functions
library(ComplexHeatmap)
# paletteer offers a collection of color palettes
library(paletteer )

# Define a color palette from Archr
Archr_col <- c('#D51F26','#272E6A','#208A42','#89288F','#F47D2B','#FEE500','#8A9FD1',
               '#C06CAB','#E6C2DC','#90D5E4','#89C75F','#F37B7D','#9983BD','#D24B27',
               '#3BBCA8','#6E4B9E','#0C727C','#7E1416','#D8A767')

# Function to order a matrix based on row order
# This function scales the input matrix, finds the maximum value in each column,
# and then orders the resulting data frame according to the given row order.
order_matrix <- function(input_matrix, row_order) {
    # Scale the input matrix
    input_matrix <- scale(input_matrix)
    # Replace NA values with 0
    input_matrix[is.na(input_matrix)] <- 0
    
    # Find the maximum value and its corresponding row label in each column
    return_df <- dplyr::bind_rows(lapply(colnames(input_matrix), function(x) {
        max_value <- max(input_matrix[, x])
        max_label <- rownames(input_matrix)[input_matrix[, x] == max_value][1]
        df <- data.frame(var = x, label = max_label, value = max_value)
        return(df)
    }))
    
    # Set row names of the data frame
    rownames(return_df) <- return_df$var
    
    # Order the data frame according to the given row order
    return_df <- dplyr::bind_rows(lapply(row_order, function(x) {
        tmp_df <- return_df[return_df$label == x, ]
        tmp_df <- tmp_df[order(tmp_df$value, decreasing = TRUE), ]
        return(tmp_df)
    }))
    
    return(return_df)
}

# Function to plot gene module scores
# This function creates a scatter plot of gene module scores using ggplot2.
plot_gene_module_scores <- function(plot_obj, feature_var, group_var, group_value, score_name) {
    ggplot() +
        geom_point(
            data = plot_obj@meta.data[plot_obj[[group_var]] == group_value, ],
            aes(x = new_x, y = new_y, color = !!sym(score_name)),
            size = 1
        ) +
        scale_color_gradientn(
            colours = c('gray50', 'gray', 'gray97', 'red', 'darkred'),
            limits = c(min(plot_obj[[score_name]]), max(plot_obj[[score_name]])),
            name = score_name
        ) +
        theme_void() +
        coord_fixed()
}

# Read the Seurat object from an RDS file
obj <- readRDS('snRNA_EX.rds')

# Update the label_transfer column
obj$label_transfer[obj$label_transfer == 'L2'] <- 'L2/3'

# Create a new group column
obj$group <- paste0(obj$species, obj$label_transfer)

# Calculate the average expression for each group
df <- AverageExpression(obj, group.by = 'group')

# Create a data frame for row annotations
anno_row <- data.frame('anno' = c(rep('human', 8), rep('macaque', 8), rep('marmoset', 8)))
rownames(anno_row) <- colnames(cor(as.matrix(df$RNA)))

# Define annotation colors
annotation_colors1 <- list(
    anno = c(human = Archr_col[9], macaque = Archr_col[2], marmoset = Archr_col[1])
)

# Create a heatmap of the correlation matrix
pheatmap::pheatmap(
    cor(as.matrix(df$RNA)),
    cutree_rows = 3,
    annotation_col = anno_row,
    annotation_colors = annotation_colors1,
    filename = 'label_cor_heat.pdf',
    width = 8,
    height = 8
)

# Read marker gene data from CSV files
mar2mac <- read.csv('new_mar2mac_marker.csv', row.names = 1)
human2mac <- read.csv('new_human2mac_marker.csv', row.names = 1)
human2mar <- read.csv('new_human2mar_marker.csv', row.names = 1)
mac2mar <- read.csv('new_mac2mar_marker.csv', row.names = 1)

# Set thresholds for filtering marker genes
threshold <- 0.25

# Filter marker genes based on thresholds
t1 <- mar2mac[mar2mac$avg_log2FC > threshold & abs(mar2mac$pct.1 - mar2mac$pct.2) > 0.15 & mar2mac$p_val_adj < 0.05, ]
t2 <- human2mac[human2mac$avg_log2FC > threshold & abs(human2mac$pct.1 - human2mac$pct.2) > 0.15 & human2mac$p_val_adj < 0.05, ]
t3 <- human2mar[human2mar$avg_log2FC > threshold & abs(human2mar$pct.1 - human2mar$pct.2) > 0.15 & human2mar$p_val_adj < 0.05, ]
t4 <- mac2mar[mac2mar$avg_log2FC > threshold & abs(mac2mar$pct.1 - mac2mar$pct.2) > 0.15 & mac2mar$p_val_adj < 0.05, ]
t5 <- human2mac[human2mac$avg_log2FC < -threshold & abs(human2mac$pct.1 - human2mac$pct.2) > 0.15 & human2mac$p_val_adj < 0.05, ]
t6 <- human2mar[human2mar$avg_log2FC < -threshold & abs(human2mar$pct.1 - human2mar$pct.2) > 0.15 & human2mar$p_val_adj < 0.05, ]

# Find the intersection of marker genes
human_mar <- intersect(rownames(t1), rownames(t2))
human_mac <- intersect(rownames(t3), rownames(t4))
mar_mac <- intersect(rownames(t5), rownames(t6))
C <- intersect(rownames(t1), rownames(t6))
H <- intersect(rownames(t2), rownames(t3))
M <- intersect(rownames(t4), rownames(t5))

# Scale the expression matrices
tot_matrix <- ScaleData(na.omit(as.matrix(df$RNA[human_mar, ])))
mar_matrix <- ScaleData(na.omit(as.matrix(df$RNA[human_mar, grep('marmoset', colnames(df$RNA))])))
mac_matrix <- ScaleData(na.omit(as.matrix(df$RNA[human_mar, grep('macaque', colnames(df$RNA))])))
human_matrix <- ScaleData(na.omit(as.matrix(df$RNA[human_mar, grep('human', colnames(df$RNA))])))

# Calculate the mean matrix
mean_matrix <- mar_matrix + human_matrix

# Order the genes
gene_order <- order_matrix(t(mean_matrix), colnames(mean_matrix))
gene_order <- gene_order[gene_order$value > 1, ]

# Create a heatmap using ComplexHeatmap
Heatmap(
    t(tot_matrix)[, gene_order$var],
    show_row_names = TRUE,
    show_heatmap_legend = TRUE,
    use_raster = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    col = paletteer_c("viridis::inferno", n = 100)
)

# Read another Seurat object from an RDS file
st_obj <- readRDS('/mnt/SSD16Ta/shenglong/data/inter_6_data/fit_st_obj.rds')

# Subset the Seurat object
plot_obj <- st_obj[, st_obj$slice %in% c('T899', 'T32', 'T454') & st_obj$layer != 'White']
plot_obj@active.assay <- 'RNA'

# Normalize the data
plot_obj <- NormalizeData(plot_obj)

# Add module scores
plot_obj <- AddModuleScore(plot_obj, features = list(gene_order[gene_order$label == 'marmosetL2/3', 'var']), name = 'L2.3')
plot_obj <- AddModuleScore(plot_obj, features = list(gene_order[gene_order$label == 'marmosetL4', 'var']), name = 'L4')

# Define a function to plot scores for different combinations
plot_scores <- function(plot_obj, score_name, group_var, group_values) {
    for (value in group_values) {
        print(plot_gene_module_scores(plot_obj, score_name, group_var, value, score_name))
    }
}

# Plot L2.3 scores
#plot_scores(plot_obj, 'L2.31', 'new_species', c('human'))
plot_scores(plot_obj, 'L2.31', 'new_slice', c('T899','T454', 'T32'))

# Plot L4 scores
#plot_scores(plot_obj, 'L41', 'new_species', c('human'))
plot_scores(plot_obj, 'L41', 'new_slice', c('T899','T454', 'T32'))