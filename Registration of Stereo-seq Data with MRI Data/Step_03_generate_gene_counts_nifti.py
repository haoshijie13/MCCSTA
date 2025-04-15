#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Core data processing libraries
import numpy as np  # Numerical computing
import pandas as pd  # Data manipulation
from pandas.api.types import CategoricalDtype  # Categorical data handling
import os  # File system operations

# Visualization
import matplotlib.pyplot as plt
import matplotlib.colors as clr  # Color manipulation

# Single-cell genomics analysis
import scanpy as sc  # Main analysis toolkit
import anndata as ad  # Data structure for genomics data

# Medical imaging processing
import nibabel as nib  # NIfTI file I/O

# Scientific computing
import scipy.io as sio  # MATLAB file support

# Machine learning preprocessing
from sklearn import preprocessing

# Performance optimization
import time
import multiprocessing as mp  # Parallel processing
import numba
from numba import jit  # Just-In-Time compiler for speed

@jit(nopython=True)
def totalRNA(xy, MIDCount, MIDCount_Mat):
    """
    Numba-accelerated function to map gene expression data to spatial coordinates
    
    Parameters:
    xy : array - Spatial coordinates matrix
    MIDCount : array - Molecular Identifier counts
    MIDCount_Mat : array - Target 3D matrix to fill
    
    Returns:
    Updated MIDCount_Mat with spatial expression data
    """
    icount = 0
    for irow in xy:
        # Accumulate UMI counts at corresponding spatial coordinates
        MIDCount_Mat[int(irow[0]), int(irow[1])] += MIDCount[icount]
        icount += 1
    return MIDCount_Mat

# Configuration for top 5000 genes analysis
h5adfolder = '/marmoset_V2_20220627_nifti/coronal_h5ad'  # Input directory
genecsv = '/coronal_allslice_common_genes_replicate.csv'  # Gene list path
savepath = '/coronal_h5ad_gene_nifti_first5000'  # Output directory

# Load and process top 5000 highly expressed genes
genelist = pd.read_csv(genecsv, delimiter=',', header=0, index_col=0)
genelist = genelist.gene.values.tolist()[:5000]  # Select top 5000 genes

# Save gene list for reference
genelistpd = pd.DataFrame({'name': genelist})
genelistpd.to_csv(f'{savepath}/first5000_genelist.csv', sep=',', index=False)

# Main processing pipeline
for i in range(0, len(h5adlist)):
    curh5ad = os.path.join(h5adfolder, h5adlist[i])
    curID = h5adlist[i].split('.')[0]
    print(f'Processing {curID}...')

    # Load spatial genomics data
    adata = sc.read_h5ad(curh5ad)
    
    # Extract spatial coordinates from observation names
    obs_str = list(adata.obs_names)
    coord_x = [int(x.split('_')[0]) for x in obs_str]
    coord_y = [int(x.split('_')[1]) for x in obs_str]
    
    # Store coordinates in AnnData object
    adata.obsm['spatial'] = np.vstack([coord_x, coord_y]).T
    
    # Get spatial dimensions
    x_trans_max = np.max(coord_x)
    y_trans_max = np.max(coord_y)
    
    # Initialize 3D matrix (X x Y x Genes)
    MIDCount_nifti = np.zeros([x_trans_max+1, y_trans_max+1, len(genelist)])

    # Populate matrix for each gene
    for count, curgene in enumerate(genelist):
        # Extract gene expression values
        curgene_exp = np.squeeze(np.array(adata[:, curgene].X.todense()))
        
        # Initialize layer for current gene
        MIDCount_Matrix = np.zeros([x_trans_max+1, y_trans_max+1])
        
        # Map expression to spatial coordinates
        MIDCount_Mat = totalRNA(adata.obsm['spatial'], curgene_exp, MIDCount_Matrix)
        
        # Add to 3D stack
        MIDCount_nifti[:, :, count] = MIDCount_Mat

    # Create and save NIfTI image
    nifti_img = nib.Nifti1Image(MIDCount_nifti, np.eye(4))
    nib.save(nifti_img, os.path.join(savepath, f'{curID}.nii.gz'))

    print(f'Completed {curID}')
