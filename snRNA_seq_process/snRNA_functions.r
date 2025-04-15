# Function to perform clustering and visualization with different parameter combinations
cluster_and_visualize <- function(seurat, PCs, resolutions, K, name) {
  count <- 0
  for (PC in PCs) {
    seurat <- RunUMAP(seurat, dims = 1:PC)
    for (resolution in resolutions) {
      for (k in K) {
        count <- count + 1
        print(paste0("Number: ", count, "; PC: ", PC, "; resolution: ", resolution, "; K: ", k))
        filename <- paste(name, count, PC, resolution, k, sep = "_")
        
        # Find neighbors and clusters
        seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PC, k.param = k)
        print("FindNeighbors done !")
        seurat <- FindClusters(object = seurat, resolution = resolution, verbose = FALSE)
        print("FindClusters done !")
        
        # Print cluster levels
        cluster_levels <- levels(seurat)
        print(cluster_levels)
        
        # Store cluster information in meta.data
        seurat@meta.data[[filename]] <- seurat@meta.data$seurat_clusters
        
        # Create a pie chart of cluster frequencies
        cluster_table <- table(Idents(object = seurat))
        cluster_names <- paste(paste("cls", names(cluster_table), sep = ""), cluster_table, sep = ":")
        cluster_df <- data.frame(Var1 = cluster_names, Freq = as.numeric(cluster_table))
        pie_plot <- ggplot(cluster_df, aes(x = "", y = Freq, fill = Var1)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar(theta = "y") +
          labs(x = "", y = "", title = "") +
          theme(axis.ticks = element_blank()) +
          theme(legend.title = element_blank(), legend.position = "top")
        
        # Create the cluster directory if it doesn't exist
        if (!dir.exists("./cluster")) {
          dir.create("./cluster")
        }
        
        # Save the pie chart as a TIFF file
        tiff(paste0("./cluster/", filename, ".pie.tiff"), height = 600, width = 600)
        print(pie_plot)
        dev.off()
        
        # Generate and save UMAP plot
        print("########:UMAP plotting!")
        tiff(paste0("./cluster/", filename, ".umap.tiff"), height = 600, width = 800)
        umap_plot <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 6, pt.size = 2) +
          ggtitle(filename) +
          coord_fixed(ratio = 1)
        print(umap_plot)
        dev.off()
        
        # Generate and save quality plots
        pdf(paste0("./cluster/", filename, ".quality.pdf"), height = 4, width = 12)
        print(VlnPlot(object = seurat, features = c("nFeature_RNA"), group.by = filename, ncol = 1, pt.size = 0) + NoLegend())
        print(VlnPlot(object = seurat, features = c("nCount_RNA"), group.by = filename, ncol = 1, pt.size = 0) + NoLegend())
        print(VlnPlot(object = seurat, features = c("percent.mt"), group.by = filename, ncol = 1, pt.size = 0) + NoLegend())
        dev.off()
      }
    }
  }
  return(seurat)
}

# Find the positions where cumulative variance reaches 85%, 90%, and 95%
find_pos <- function(target) {
  which.min(abs(total - target))
}


# Function to plot UMAP and violin plots for quality control
plot_quality_control <- function(lab, name, len) {
  # Plot UMAP colored by group
  tiff(paste0("./QC/umap.group.tiff"), width = 1000, height = 800)
  p <- ggplot(lab, aes(x = X, y = Y, color = group)) + 
    geom_point(size = 0.5) +
    theme_classic()
  print(p)
  dev.off()
  
  # Plot violin plots for nFeature_RNA and percent.mt
  p1 <- ggplot(lab, aes(x = group, y = nFeature_RNA, fill = group)) + 
    geom_violin() +
    theme_classic()
  p2 <- ggplot(lab, aes(x = group, y = percent.mt, fill = group)) + 
    geom_violin() +
    theme_classic()
  
  # Combine the two violin plots
  library(cowplot)
  p <- plot_grid(p1, p2, ncol = 1, align = "hv", rel_heights = c(0.5, 0.5))
  
  # Save the combined violin plot
  tiff(paste0("./QC/qulity.group.tiff"), width = 3000, height = 1000)
  print(p)
  dev.off()
}

# Function to identify outliers in each cluster
identify_outliers <- function(tmp) {
  outlier_number <- c()
  re <- data.frame()
  
  # Function to remove outliers based on boxplot rules
  remove_outliers <- function(x, na.rm = TRUE) {
    qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    pos1 <- which(x < (qnt[1] - H))
    pos2 <- which(x > (qnt[2] + H))
    c(pos1, pos2)
  }
  
  # Loop through each group to identify outliers
  for (i in seq_along(tmp)) {
    cat(paste0("Processing group: ", names(tmp[i]), "\n"))
    a <- tmp[[i]]
    median_X <- median(a[,"X"])
    median_Y <- median(a[,"Y"])
    
    # Calculate the distance from each cell to the cluster center
    distance_original_point <- apply(a, 1, function(x) 
      (as.numeric(x["X"]) - median_X)^2 + (as.numeric(x["Y"]) - median_Y)^2 
    )
    a$distance_original_point <- distance_original_point
    a$distance_original_point1 <- sqrt(distance_original_point)
    
    # Identify outliers based on the distance
    pos <- remove_outliers(a$distance_original_point1)
    cat(paste0("Number of outliers in group ", names(tmp[i]), ": ", length(pos), "\n"))
    a[pos,"Sub"] <- "outlier"
    re <- rbind(re, a)
    outlier_number <- c(outlier_number, length(pos))
  }
  
  names(outlier_number) <- names(tmp)
  return(list(re = re, outlier_number = outlier_number))
}

# Function to plot outlier-related plots
plot_outliers <- function(re, len, total_outliers) {
  # Plot boxplot of distances with outliers colored red
  p2 <- ggplot(re, aes(x = group, y = distance_original_point1, fill = group)) + 
    geom_boxplot(outlier.colour = "red") +
    theme_classic()
  tiff(paste0("./QC/ggboxplot.outlier.tiff"), width = ceiling(length(unique(re$group))/2) * 300, height = 1000)
  print(p2)
  dev.off()
  
  # Plot UMAP with outliers labeled
  tiff(paste0("./QC/umap.outlier.all.tiff"), width = 400 * len + 200, height = 400 * len)
  p <- ggplot(re, aes(x = X, y = Y, color = Sub)) + 
    geom_point(size = 0.5) +
    theme_classic()
  p1 <- p + facet_wrap(~group, ncol = len, scales = "fixed")
  print(p1)
  dev.off()
  
  # Plot only outlier cells
  dat1 <- subset(re, subset = re$Sub == "outlier")
  dim(dat1)
  tiff(paste0("./QC/umap.outlier.", total_outliers, ".tiff"), width = 400 * len + 200, height = 400 * len)
  p <- ggplot(dat1, aes(x = X, y = Y, color = Sub)) + 
    geom_point(size = 0.5) +
    theme_classic()
  p1 <- p + facet_wrap(~group, ncol = len, scales = "fixed") + coord_fixed(ratio = 1)
  print(p1)
  dev.off()
  
  # Plot non-outlier cells
  dat2 <- subset(re, subset = re$Sub == "sub")
  tiff(paste0("./QC/umap.not.outlier.tiff"), width = 400 * len + 200, height = 400 * len)
  p <- ggplot(dat2, aes(x = X, y = Y, color = group)) + 
    geom_point(size = 0.5) +
    theme_classic()
  p1 <- p + facet_wrap(~group, ncol = len, scales = "fixed") + coord_fixed(ratio = 1)
  print(p1)
  dev.off()
}

# Function to remove clusters and cells based on criteria
remove_clusters_and_cells <- function(re) {
  re$remove_clusters <- "sub"
  
  # Remove clusters with high mitochondrial gene content
  tmp <- split(re$percent.mt, re$group)
  tmp <- sapply(tmp, summary)
  a <- which(tmp["Median",] > 2 & tmp["Mean",] > 2) %>% names %>% as.numeric
  print(a)
  re$remove_clusters[re$group %in% a] <- "remove"
  dat <- subset(re, subset = re$remove_clusters == "sub")
  
  # Remove outliers from each cluster
  dat <- subset(dat, subset = dat$Sub == "sub")
  
  # Remove clusters with less than 100 cells
  tmp <- table(dat$group)
  b <- which(tmp < 100) %>% names %>% as.numeric
  print(b)
  dat$remove_clusters[dat$group %in% b] <- "remove"
  dat <- subset(dat, subset = dat$remove_clusters == "sub")
  
  return(dat)
}

# Function to plot the final cells in each cluster
plot_final_cells <- function(dat, len) {
  level <- unique(dat$group)
  dat$group <- factor(dat$group, levels = level)
  levels(dat$group) <- 1:length(level)
  p <- ggplot(dat, aes(x = X, y = Y, color = group)) + 
    geom_point(size = 0.5) +
    theme_classic()
  
  # Save the first plot of final cells
  tiff(paste0("./QC/umap.final.group.1.tiff"), width = 1000, height = 800)
  p1 <- p + coord_fixed(ratio = 1)
  library(Seurat)
  LabelClusters(plot = p1, id = 'group', colour = "black", size = 5)
  dev.off()
  
  # Save the second plot of final cells
  tiff(paste0("./QC/umap.final.group.2.tiff"), width = 400 * len + 200, height = 400 * len)
  p1 <- p + facet_wrap(~group, ncol = len, scales = "fixed") + coord_fixed(ratio = 1)
  print(p1)
  dev.off()
}