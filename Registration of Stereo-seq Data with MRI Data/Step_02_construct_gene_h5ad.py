#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Core scientific computing libraries
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype

# System operations and visualization
import os
import matplotlib.pyplot as plt

# Single-cell/spatial genomics tools
import scanpy as sc  # Main analysis library
import anndata as ad  # Data structure for genomics data

# Neuroimaging and advanced visualization
import nibabel as nib  # MRI data handling
from mpl_toolkits.axes_grid1 import make_axes_locatable  # Plot layout control

# Color processing and math operations
import matplotlib.colors as colors
import scipy.io as sio
from scipy import sparse  # Sparse matrix support

# Performance optimization
import time
import multiprocessing as mp
import numba
from numba import jit  # Just-In-Time compilation accelerator

#----------------------------------------------------------
# Core Processing Functions
#----------------------------------------------------------

def extractseq(row):
    """Parse genomic coordinate data from tab-delimited strings"""
    tmp = row.split('\t')
    return tmp[0], int(tmp[1]), int(tmp[2]), int(tmp[3])

def xy_coord_trans(xydf, scale):
    """Convert raw spatial coordinates to binned grid positions"""
    x_trans = np.round(np.array(xydf.x) / scale)
    y_trans = np.round(np.array(xydf.y) / scale)
    return np.vstack([x_trans, y_trans]).T, np.max(x_trans), np.max(y_trans)

@jit(nopython=True)
def totalRNA(xy, MIDCount, MIDCount_Mat):
    """Numba-accelerated UMI count matrix population"""
    for icount in range(len(xy)):
        x, y = xy[icount]
        MIDCount_Mat[x, y] += MIDCount[icount]
    return MIDCount_Mat

@jit(nopython=True)
def xycoordlabel(xy, label, seqimg):
    """Label assignment based on spatial barcode positions"""
    for icount in range(len(xy)):
        x, y = xy[icount]
        label[icount] = seqimg[x, y] if seqimg[x, y] != 0 else label[icount]
    return label

#----------------------------------------------------------
# Main Processing Pipeline
#----------------------------------------------------------

# Configuration
seqpath = '/marmoset_V2_20220627/coronal'
seqlist = sorted(os.listdir(seqpath))  # Get sorted list of input files

# Process each sequencing sample
for i_count, curfile in enumerate(seqlist):
    # File handling
    curID = curfile.split('-')[0]
    curseq = os.path.join(seqpath, f"{curID}-all")
    
    # Data loading
    print(f"Processing {curfile}...")
    seqdf = pd.read_csv(curseq, delimiter='\t', comment='#', header=0)
    
    # Spatial coordinate processing
    xybin = 50  # Micron-to-pixel scale factor
    xmin, ymin = seqdf.x.min(), seqdf.y.min()
    x_trans = np.round((seqdf.x - xmin)/xyastype(int)
    y_trans = np.round((seqdf.y - ymin)/xybin).astype(int)
    
    # Create binned data structure
    bindata = pd.DataFrame({
        'xbin': x_trans,
        'ybin': y_trans,
        'umi_count': seqdf.iloc[:, 3],
        'gene_name': seqdf.iloc[:, 0],
        'loci': x_trans.astype(str) + '_' + y_trans.astype(str)
    })
    
    # Aggregate UMI counts
    bininfo = bindata.groupby(['gene_name', 'loci'])['umi_count'].sum().reset_index()
    
    # Create AnnData object (standard genomics data container)
    obs = bininfo.loci.unique()
    var = bininfo.gene_name.unique()
    obs_cat = CategoricalDtype(sorted(obs), ordered=True)
    var_cat = CategoricalDtype(sorted(var), ordered=True)
    
    # Sparse matrix conversion
    obs_index = bininfo["loci"].astype(obs_cat).cat.codes
    var_index = bininfo["gene_name"].astype(var_cat).cat.codes
    coo = sparse.coo_matrix(
        (bininfo["count"], (obs_index, var_index)), 
        shape=(len(obs), len(var))
    )
    
    # Build final AnnData structure
    adata = ad.AnnData(coo.tocsr())
    adata.obs_names = obs_cat.categories
    adata.var_names = var_cat.categories
    
    # Add spatial coordinates
    spatial_coords = adata.obs_names.str.split('_', expand=True).astype(int)
    adata.obsm["spatial"] = spatial_coords.values
    
    # Calculate total RNA content
    adata.obs['totalRNA'] = adata.X.sum(axis=1)
    
    # Visualization and output
    sc.pl.embedding(
        adata, 
        basis="spatial", 
        color='totalRNA',
        s=2,  # Point size
        show=False,
        save=f"{curID}.png"  # Auto-save to path
    )
    
    # Save processed data
    adata.write_h5ad(
        f"/marmoset_V2_20220627_nifti/png/{curID}.h5ad",
        compression='gzip'
    )
