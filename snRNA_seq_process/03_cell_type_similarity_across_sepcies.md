# 1. Calculate Meta cells

~~~python
###MetaCells Pipeline: Merging metacells from different species, converting to RDS, and integrating across species.
# Computing Metacells - One-pass Process
# 1. Setup
import anndata as ad             # For reading/writing AnnData files
import matplotlib.pyplot as plt  # For plotting
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sb             # For plotting
import scipy.sparse as sp        # For sparse matrices
import shutil                    # for filesystem operations
from math import hypot           # For plotting
from typing import *
import argparse
import scanpy as sc

# Use SVG for scalable low-element-count diagrams.
# Set the seaborn style to white for a cleaner look
sb.set_style("white")

# Disable inefficient layout in metacells
mc.ut.allow_inefficient_layout(False)

# Change the current working directory
os.chdir("/home/yangqq/projectold/Marmoset_brain/Script/snRNA_process/MetaCell_process/")

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cluster', type=str, help='Cluster information')
parser.add_argument('--data_path', type=str, help='Path to the data file')
parser.add_argument('--species', type=str, help='Species name')
parser.add_argument('--Class', type=str, help='Class information')

args = parser.parse_args()
cluster = args.cluster
data_path = args.data_path
species = args.species
Class = args.Class
savepath = f'/home/yangqq/projectold/Marmoset_brain/Script/snRNA_process/MetaCell_process/result/{species}/{Class}'

# 2. Reading the data
def create_output_directories():
    """
    Create necessary output directories for the analysis.
    """
    directories = [
        f'{savepath}/output/iterative/iteration-1/figures',
        f'{savepath}/output/iterative/iteration-2',
        f'{savepath}/output/iterative/iteration-3',
        f'{savepath}/output/iterative/iteration-4',
        f'{savepath}/output/iterative/final'
    ]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)

create_output_directories()

# Read the data using scanpy
full = sc.read(data_path)
# Set the data to the top level in metacells
mc.ut.top_level(full)
# Set a name for the data
mc.ut.set_name(full, "hca_bm.full")
print(f"Full: {full.n_obs} cells, {full.n_vars} genes")

# 3. Cleaning the data
# 3.2.1 Excluding cells by UMIs count
PROPERLY_SAMPLED_MIN_CELL_TOTAL = 800
PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000

def analyze_umi_distribution():
    """
    Analyze the distribution of total UMIs per cell and plot it.
    Identify and report cells with too few or too many UMIs.
    """
    total_umis_per_cell = mc.ut.get_o_numpy(full, "__x__", sum=True)
    plot = sb.displot(total_umis_per_cell, log_scale=(10, None))
    plot.set(xlabel="UMIs", ylabel="Density", yticks=[])
    plot.refline(x=PROPERLY_SAMPLED_MIN_CELL_TOTAL, color="darkgreen")
    plot.refline(x=PROPERLY_SAMPLED_MAX_CELL_TOTAL, color="crimson")
    plt.savefig(f'{savepath}/output/iterative/iteration-1/figures/cell_total_umis.svg')

    too_small_cells_count = np.sum(total_umis_per_cell < PROPERLY_SAMPLED_MIN_CELL_TOTAL)
    too_large_cells_count = np.sum(total_umis_per_cell > PROPERLY_SAMPLED_MAX_CELL_TOTAL)

    too_small_cells_percent = 100.0 * too_small_cells_count / full.n_obs
    # Bug fix: Use full.n_obs instead of full.n_vars
    too_large_cells_percent = 100.0 * too_large_cells_count / full.n_obs

    print(
        f"Will exclude {too_small_cells_count} ({too_small_cells_percent:.2f}%) cells "
        f"with less than {PROPERLY_SAMPLED_MIN_CELL_TOTAL} UMIs"
    )
    print(
        f"Will exclude {too_large_cells_count} ({too_large_cells_percent:.2f}%) cells "
        f"with more than {PROPERLY_SAMPLED_MAX_CELL_TOTAL} UMIs"
    )

analyze_umi_distribution()

def exclude_genes_based_on_species():
    """
    Exclude genes based on the species.
    """
    if species == "mouse":
        excluded_patterns = ["mt-.*"]
    elif species == "human":
        excluded_patterns = ["MT-.*"]
    else:
        excluded_names = ["ND6", "COX3", "COX1", "ND5", "ND4", "ND2", "ND4L", "ATP8", "CYTB", "COX2", "ND3", "ATP6", "ND1"]
        mc.pl.exclude_genes(
            full,
            excluded_gene_names=excluded_names,
            random_seed=123456
        )
        return
    mc.pl.exclude_genes(
        full,
        excluded_gene_patterns=excluded_patterns,
        random_seed=123456
    )

exclude_genes_based_on_species()

# Compute the UMIs of excluded genes
mc.tl.compute_excluded_gene_umis(full)
PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION = 0.25
excluded_umis_fraction_regularization = 1e-3  # Avoid 0 values in log scale plot.

def analyze_excluded_umi_fraction():
    """
    Analyze the fraction of excluded UMIs per cell and plot it.
    Identify and report cells with a high fraction of excluded UMIs.
    """
    total_umis_per_cell = mc.ut.get_o_numpy(full, name="__x__", sum=True)
    excluded_umis_per_cell = mc.ut.get_o_numpy(full, "excluded_umis")
    excluded_umis_fraction_per_cell = excluded_umis_per_cell / total_umis_per_cell

    excluded_umis_fraction_per_cell += excluded_umis_fraction_regularization
    plot = sb.displot(excluded_umis_fraction_per_cell, log_scale=(10, None))
    excluded_umis_fraction_per_cell -= excluded_umis_fraction_regularization

    plot.set(xlabel="Fraction of excluded gene UMIs", ylabel="Density", yticks=[])
    plot.refline(x=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION, color="crimson")
    plt.savefig(f'{savepath}/output/iterative/iteration-1/figures/cell_excluded_umis_fraction.svg')

    too_excluded_cells_count = np.sum(
        excluded_umis_fraction_per_cell > PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION
    )
    too_excluded_cells_fraction = too_excluded_cells_count / len(total_umis_per_cell)

    print(
        f"Will exclude {too_excluded_cells_count} ({100 * too_excluded_cells_fraction:.2f}%) cells "
        f"with more than {100 * PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION:.2f}% excluded gene UMIs"
    )

analyze_excluded_umi_fraction()

# Exclude cells based on UMI count and fraction of excluded genes
mc.pl.exclude_cells(
    full,
    properly_sampled_min_cell_total=PROPERLY_SAMPLED_MIN_CELL_TOTAL,
    properly_sampled_max_cell_total=PROPERLY_SAMPLED_MAX_CELL_TOTAL,
    properly_sampled_max_excluded_genes_fraction=PROPERLY_SAMPLED_MAX_EXCLUDED_GENES_FRACTION,
)

# 3.2.3 Extract the clean data
clean = mc.pl.extract_clean_data(full, name="hca_bm.iteration-1.clean")
mc.ut.top_level(clean)
print(f"Clean: {clean.n_obs} cells, {clean.n_vars} genes")

# 3.4 Save the data
full.write_h5ad(f'{savepath}/output/iterative/iteration-1/hca_bm.full.h5ad')
# Release memory by setting full to None
full = None

clean.write_h5ad(f'{savepath}/output/iterative/iteration-1/hca_bm.clean.h5ad')

# 4. Compute the metacells
cells = clean
# Release memory by setting clean to None
clean = None
mc.ut.set_name(cells, "hca_bm.iteration-1.cells")
print(f"Iteration 1: {cells.n_obs} cells, {cells.n_vars} genes")

# Read the lateral gene data
LATERAL_GENE = pd.read_csv("/home/yangqq/projectold/Marmoset_brain/Script/snRNA_process/MetaCell_process/Cellcircle_gene_orth.csv", index_col=0)
LATERAL_GENE.columns = ["human", "marmoset", "macaque", "mouse"]
LATERAL_GENE_NAMES = list(LATERAL_GENE[species])

# Mark lateral genes in the data
mc.pl.mark_lateral_genes(
    cells,
    lateral_gene_names=LATERAL_GENE_NAMES
)

# Get the mask for lateral genes
lateral_gene_mask = mc.ut.get_v_numpy(cells, "lateral_gene")
lateral_gene_names = set(cells.var_names[lateral_gene_mask])
print(sorted([
    name for name in lateral_gene_names
    if not name.startswith("RPL") and not name.startswith("RPS")
]))
print(f"""and {len([
    name for name in lateral_gene_names if name.startswith("RPL") or name.startswith("RPS")
])} RP[LS].* genes""")

# 4.1.3 Parallelization
# Either use the guesstimator:
max_parallel_piles = mc.pl.guess_max_parallel_piles(cells)
# Or, if running out of memory manually override:
# max_parallel_piles = ...
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)

# 4.2 Computation
# 4.2.2 Assigning cells to metacells
with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(cells, random_seed=123456, target_metacell_size=13)

# 4.2.3 Collecting the metacells
metacells = \
    mc.pl.collect_metacells(cells, name="hca_bm.one-pass.preliminary.metacells", random_seed=123456)
print(f"Preliminary: {metacells.n_obs} metacells, {metacells.n_vars} genes")

# Convey cluster information from cells to metacells
mc.tl.convey_obs_to_group(
    adata=cells, gdata=metacells,
    property_name=cluster, to_property_name=cluster,
    method=mc.ut.most_frequent  # This is the default, for categorical data
)

# Convey cluster fraction information from cells to metacells
mc.tl.convey_obs_fractions_to_group(
    adata=cells, gdata=metacells,
    property_name=cluster, to_property_name=cluster
)

# Save the cells and metacells data
cells.write(f'{savepath}/output/iterative/final/cells_save.h5ad')
metacells.write(f'{savepath}/output/iterative/final/metacells_save.h5ad')
~~~

# 2. Cross species integrate

~~~R
# Suppress package loading messages to keep the console clean
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(Matrix))
suppressMessages(library(viridis))
suppressMessages(library(cowplot))
suppressMessages(library(ggsci))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
# suppressMessages(library(LSD))
suppressMessages(library(readxl))
suppressMessages(library(ggrepel))
suppressMessages(library(harmony))
suppressMessages(library(scrattch.hicat))
suppressMessages(library(getopt))

# Set the maximum size of global variables in the future package
options(future.globals.maxSize = 100 * 1024^3)

# Define a function to select marker genes
# norm.dat: Normalized data matrix
# cl: Cluster labels
# n.markers: Number of top marker genes to select for each pair
# de.genes: Pre - calculated differential expression statistics
select_markers = function(norm.dat, cl, n.markers = 20, de.genes = NULL, ...) {
  if (is.null(de.genes)) {
    de.genes = de_stats_all_pairs(norm.dat, cl, ...)
  }
  pairs = names(de.genes)
  select.pairs = pairs
  de.markers = sapply(select.pairs, function(s) {
    tmp = de.genes[[s]]
    c(head(tmp$up.genes, n.markers), head(tmp$down.genes, n.markers))
  }, simplify = F)
  markers = intersect(unlist(de.markers), row.names(norm.dat))
  return(list(markers = markers, de.genes = de.genes[select.pairs]))
}

# Define a function to prepare a Seurat object with homologous gene names
# seurat: A Seurat object
# originalGeneName: Original gene names in the Seurat object
# homoGeneName: Homologous gene names to replace the original ones
# keepMetaCol: Metadata columns to keep
# commonMetaCol: Common metadata columns across species
# speciesPrefix: Prefix to add to cell names and metadata columns
# clusterCol: Column name for cluster information
prepareHomoSeurat = function(seurat, originalGeneName, homoGeneName, keepMetaCol, commonMetaCol, speciesPrefix = "macaque_", clusterCol = "anno3") {
  homoCountMx = seurat[["RNA"]]@counts[originalGeneName, ]
  rownames(homoCountMx) = homoGeneName
  colnames(homoCountMx) = paste0(speciesPrefix, colnames(homoCountMx))
  homoSeurat = CreateSeuratObject(counts = homoCountMx)
  tmpDf = seurat@meta.data[sub(speciesPrefix, "", rownames(homoSeurat@meta.data)), keepMetaCol]
  tmpDf$species = sub("(_|-)", "", speciesPrefix)
  tmpDf$clusterName = paste0(speciesPrefix, tmpDf[[clusterCol]])
  commonMetaCol = c(commonMetaCol, "species", "clusterName")
  colnames(tmpDf) = ifelse(
    colnames(tmpDf) %in% commonMetaCol, colnames(tmpDf),
    paste0(speciesPrefix, colnames(tmpDf))
  )
  rownames(tmpDf) = paste0(speciesPrefix, rownames(tmpDf))
  homoSeurat = AddMetaData(homoSeurat, metadata = tmpDf)
  return(homoSeurat)
}

# Define file paths for Seurat objects of different species
seurat_human = "./MetaCell_process/result_cluster/data/rds/rowcounts/human_120percluster_total22835.rds"
seurat_macaque = "./MetaCell_process/result_cluster/data/rds/rowcounts/macaque_EX.rds"
seurat_marmoset = "./MetaCell_process/result_cluster/data/rds/rowcounts/marmoset_EX.rds"
seurat_mouse = "./MetaCell_process/result_cluster/mouse_zhuangxiaowei/data/rds/rowcounts/all_EX_200percluster.rds"

# Define label and output directory
label = "four_species_EX_metacell_0929v1_"
out_put_dir = "./MetaCell_process/result_cluster/mouse_zhuangxiaowei/data/rds/cross/result/"

# Read the table of one - to - one orthologous genes
oneToOneOrthGeneTb = read.table("./mart_export.humanMacaqeMarmosetMouse.oneToOneOrth.ensembl91.20220428.csv", sep = ",", header = TRUE)

# Calculate the number of intersecting genes between orthologous gene list and each species' Seurat object
# Note: There is a typo in 'mamrosetSeurat', it should be 'marmosetSeurat'
# For now, we assume these Seurat objects are loaded correctly later
print(length(intersect(oneToOneOrthGeneTb$marmosetGene, rownames(marmosetSeurat))))
print(length(intersect(oneToOneOrthGeneTb$macaqueGene, rownames(macaqueSeurat))))
print(length(intersect(oneToOneOrthGeneTb$mouseGene, rownames(mouseSeurat))))
print(length(intersect(oneToOneOrthGeneTb$humanGene, rownames(humanSeurat))))

# Select relevant columns from the orthologous gene table
comGeneTb = oneToOneOrthGeneTb[c("humanGene", "marmosetGene", "macaqueGene", "mouseGene")]
comGeneTb$homoGeneSymbol = comGeneTb$macaqueGene

# Subset the table to include only genes present in all species' Seurat objects
comGeneTb = subset(
  comGeneTb,
  humanGene %in% rownames(humanSeurat[["RNA"]]@counts) &
    marmosetGene %in% rownames(marmosetSeurat[["RNA"]]@counts) &
    macaqueGene %in% rownames(macaqueSeurat[["RNA"]]@counts) &
    mouseGene %in% rownames(mouseSeurat[["RNA"]]@counts)
)
print(comGeneTb)

# Prepare Seurat objects with homologous gene names for each species
marmosetHomoSeurat = prepareHomoSeurat(
  marmosetSeurat, originalGeneName = comGeneTb$marmosetGene, homoGeneName = comGeneTb$homoGeneSymbol,
  keepMetaCol = c("cluster", "total_umis"),
  commonMetaCol = c("total_umis"),
  speciesPrefix = "marmoset_", clusterCol = "cluster"
)

macaqueHomoSeurat = prepareHomoSeurat(
  macaqueSeurat, originalGeneName = comGeneTb$macaqueGene, homoGeneName = comGeneTb$homoGeneSymbol,
  keepMetaCol = c("cluster", "total_umis"),
  commonMetaCol = c("total_umis"),
  speciesPrefix = "macaque_", clusterCol = "cluster"
)

mouseHomoSeurat = prepareHomoSeurat(
  mouseSeurat, originalGeneName = comGeneTb$mouseGene, homoGeneName = comGeneTb$homoGeneSymbol,
  keepMetaCol = c("cluster", "total_umis"),
  commonMetaCol = c("total_umis"),
  speciesPrefix = "mouse_", clusterCol = "cluster"
)

humanHomoSeurat = prepareHomoSeurat(
  humanSeurat, originalGeneName = comGeneTb$humanGene, homoGeneName = comGeneTb$homoGeneSymbol,
  keepMetaCol = c("cluster", "total_umis"),
  commonMetaCol = c("total_umis"),
  speciesPrefix = "human_", clusterCol = "cluster"
)

# Combine the prepared Seurat objects into a list
combineSeuratList = list(
  human = humanHomoSeurat,
  macaque = macaqueHomoSeurat,
  marmoset = marmosetHomoSeurat,
  mouse = mouseHomoSeurat
)

# Process each Seurat object in the list
combineSeuratList = lapply(combineSeuratList, function(seurat) {
  seurat = SCTransform(
    seurat, assay = "RNA",
    ncells = ncol(seurat[["RNA"]]),
    variable.features.n = 3000,
    return.only.var.genes = F,
    method = "glmGamPoi"
  ) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:50) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:50)
  return(seurat)
})

# Modify the cluster names in each Seurat object to remove special characters
combineSeuratList$human$clusterNameMod = gsub("(-| |/)", "_", combineSeuratList$human$clusterName)
combineSeuratList$marmoset$clusterNameMod = gsub("(-| |/)", "_", combineSeuratList$marmoset$clusterName)
combineSeuratList$macaque$clusterNameMod = gsub("(-| |/)", "_", combineSeuratList$macaque$clusterName)
combineSeuratList$mouse$clusterNameMod = gsub("(-| |/)", "_", combineSeuratList$mouse$clusterName)

# Save the combined Seurat object list
saveRDS(combineSeuratList, paste0(out_put_dir, label, "combineSeuratList.rds"))

# Extract variable genes from each Seurat object
marmosetVarGene = combineSeuratList$marmoset@assays$SCT@var.features
macaqueVarGene = combineSeuratList$macaque@assays$SCT@var.features
mouseVarGene = combineSeuratList$mouse@assays$SCT@var.features
humanVarGene = combineSeuratList$human@assays$SCT@var.features

# Find the union of variable genes across all species
hicatMarker = Reduce(union, list(humanVarGene, macaqueVarGene, marmosetVarGene, mouseVarGene))
print(length(hicatMarker))

# Filter hicatMarker to include only genes present in all species' SCT count matrices
hicatMarker = hicatMarker[which(hicatMarker %in% rownames(combineSeuratList$human@assays$SCT@counts))]
print(length(hicatMarker))
hicatMarker = hicatMarker[which(hicatMarker %in% rownames(combineSeuratList$macaque@assays$SCT@counts))]
print(length(hicatMarker))
hicatMarker = hicatMarker[which(hicatMarker %in% rownames(combineSeuratList$marmoset@assays$SCT@counts))]
print(length(hicatMarker))
hicatMarker = hicatMarker[which(hicatMarker %in% rownames(combineSeuratList$mouse@assays$SCT@counts))]
print(length(hicatMarker))

# Save the hicatMarker
saveRDS(hicatMarker, file = paste0(out_put_dir, label, "top3000Marker.rds"))

# Function to perform integration and downstream analysis
perform_integration_analysis = function(combineSeuratList, hicatMarker, out_put_dir, label) {
  print("start integrate!")
  combineSeuratListPrep = PrepSCTIntegration(object.list = combineSeuratList, anchor.features = hicatMarker)
  integrateAnchors = FindIntegrationAnchors(
    object.list = combineSeuratListPrep, normalization.method = "SCT",
    anchor.features = hicatMarker
  )
  integrateSeurat = IntegrateData(anchorset = integrateAnchors, normalization.method = "SCT")
  DefaultAssay(integrateSeurat)
  integrateSeurat = RunPCA(integrateSeurat, features = hicatMarker, npcs = 100, verbose = FALSE)
  
  options(repr.plot.width = 7, repr.plot.height = 4)
  pdf(paste0(out_put_dir, label, "integrateSeurat_ElbowPlot.pdf"), width = 8, height = 8)
  print(ElbowPlot(integrateSeurat, ndims = 100))
  dev.off()
  
  print("finished integrate!")
  
  # Define a function to perform clustering, t - SNE, and UMAP
  perform_clustering = function(integrateSeurat, dims, resolution) {
    integrateSeurat = FindNeighbors(integrateSeurat, reduction = "pca", dims = dims) %>%
      FindClusters(resolution = resolution, n.start = 10, algorithm = 1) %>%
      RunTSNE(dims = dims) %>%
      RunUMAP(dims = dims)
    saveRDS(integrateSeurat, file = paste0(out_put_dir, label, "integrateSeurat.dim", length(dims), ".re", resolution, ".rds"))
    return(integrateSeurat)
  }
  
  # Define parameter combinations for clustering
  param_combinations = list(
    list(dims = 1:100, resolution = 1),
    list(dims = 1:100, resolution = 3),
    list(dims = 1:50, resolution = 1),
    list(dims = 1:50, resolution = 3),
    list(dims = 1:50, resolution = 1, tsne_umap_dims = 1:30),
    list(dims = 1:50, resolution = 3, tsne_umap_dims = 1:30)
  )
  
  for (param in param_combinations) {
    dims = param$dims
    resolution = param$resolution
    tsne_umap_dims = ifelse(is.null(param$tsne_umap_dims), dims, param$tsne_umap_dims)
    integrateSeurat = perform_clustering(integrateSeurat, dims, resolution)
  }
  
  print("finished!")
}

# Call the function to perform integration and analysis
perform_integration_analysis(combineSeuratList, hicatMarker, out_put_dir, label)
~~~

# 3. Calculate cluster similarity across species using Metaneighbor

~~~R
# Suppress package loading messages to keep the console clean
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
# suppressMessages(library(MetaNeighbor))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
# Duplicate library call, can be removed
# suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(Matrix))
suppressMessages(library(viridis))
suppressMessages(library(cowplot))
suppressMessages(library(ggsci))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
# suppressMessages(library(LSD))
suppressMessages(library(readxl))
suppressMessages(library(ggrepel))
suppressMessages(library(harmony))
suppressMessages(library(scrattch.hicat))
suppressMessages(library(getopt))

# Source an external R script for additional functions
source("./Metaneighbor_top5.r")

# Read the integrated Seurat object with resolution 3
EX_integrate_re3 <- readRDS("./MetaCell_process/result_cluster/mouse_zhuangxiaowei/data/rds/cross/result/four_species_EX_metacell_0929v1_integrateSeurat.dim100.re3.rds")

# Bug fix: 'EX_integrate_sub' is not defined, assume it should be 'EX_integrate_re3'
data_EX_sce = as.SingleCellExperiment(EX_integrate_re3)

# Read the list of genes
gene <- readRDS("./MetaCell_process/result_cluster/mouse_zhuangxiaowei/data/rds/cross/result/four_species_EX_metacell_0929v1_top3000Marker.rds")
# Print the length of the gene list
print(length(gene))

# Perform one-vs-best comparison using MetaNeighborUS
# var_genes: Variable genes for the analysis
# dat: SingleCellExperiment object containing the data
# study_id: Identifier for the study of each cell
# cell_type: Cell type annotation for each cell
# fast_version: Use a faster version of the algorithm
# one_vs_best: Perform one-vs-best comparison
# symmetric_output: Whether to output a symmetric matrix
OneVs_EX <- MetaNeighborUS(var_genes = gene, dat = data_EX_sce, study_id = data_EX_sce$orig.ident, cell_type = data_EX_sce$"clusterNameMod", fast_version = TRUE, one_vs_best = TRUE, symmetric_output = FALSE)

# Write the result to a CSV file
write.csv(OneVs_EX, file = "./save_file/Metacluster_EX_all.fourspecies.csv")
~~~

