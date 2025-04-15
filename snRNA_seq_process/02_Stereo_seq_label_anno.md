# 1. DestVI_snRNA_model_training

~~~python
# Import necessary libraries
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import ks_2samp, ttest_ind
from statsmodels.stats.multitest import multipletests
import scvi
from scvi.model import CondSCVI, DestVI
import logging
import umap
import torch
from torch.distributions import Gamma
import base64
from io import BytesIO
import argparse
import time

# Initialize logging for scvi
logger = logging.getLogger("scvi")

# Record the start time of the process
start_time = time.time()

# Set up command-line argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--species', type=str, help='Species information')
parser.add_argument('--sc_clustername', type=str, help='Name of the single-cell cluster')
parser.add_argument('--sc_epochs', type=int, help='Number of epochs for single-cell model training')
parser.add_argument('--sc_path', type=str, help='Path to the single-cell data file')
parser.add_argument('--save_path', type=str, help='Path to save the results')

# Parse the command-line arguments
args = parser.parse_args()
species = args.species
sc_clustername = args.sc_clustername
sc_epochs = args.sc_epochs
sc_path = args.sc_path
save_path = args.save_path

# Read the single-cell data
sc_adata = sc.read(sc_path)

# Filter genes based on minimum counts
sc.pp.filter_genes(sc_adata, min_counts=10)

# Create a copy of the raw counts in a separate layer
sc_adata.layers["counts"] = sc_adata.X.copy()

# Normalize the total counts of each cell to a target sum
sc.pp.normalize_total(sc_adata, target_sum=10e4)

# Apply log transformation to the normalized data
sc.pp.log1p(sc_adata)

# Save the raw data
sc_adata.raw = sc_adata

# Load the list of all spatial transcriptomics genes
all_stgene = np.load('/cluster/home/yangqianqian/Marmoset_brain/mouse_steore-seq/destvi_transfer/mouse_allslice.npy')

# Find the intersection of single-cell genes and spatial transcriptomics genes
intersect = np.intersect1d(sc_adata.var_names, all_stgene)

# Save the intersection of genes to a numpy file
np.save(f'{save_path}/mouse_sn_st.npy', np.array(intersect))
print(f"Length of training genes: {len(intersect)}")

# Subset the single-cell data to include only the intersecting genes
sc_adata = sc_adata[:, intersect].copy()

# Set up the AnnData object for the CondSCVI model
scvi.model.CondSCVI.setup_anndata(sc_adata, labels_key=sc_clustername, layer="counts")

# Initialize the CondSCVI model
sc_model = CondSCVI(sc_adata, weight_obs=True, n_latent=4, n_layers=2, n_hidden=128)

# Train the CondSCVI model using GPU if available
sc_model.train(max_epochs=sc_epochs, accelerator="gpu")

# Plot the training history of the model
plt.figure(figsize=(8, 6))
sc_model.history["elbo_train"].plot()
sc_model.history["elbo_train"].iloc[5:].plot()

# Save the training history plot
plt.savefig(f'{save_path}/scmodel_epoch.pdf', dpi=300)

# Save the trained model
sc_model.save(f'{save_path}/sc_model_{species}_snRNA_all', overwrite=True)

# Get the latent representation of the single-cell data from the model
sc_adata.obsm["X_CondSCVI"] = sc_model.get_latent_representation()

# Save the processed single-cell data
sc_adata.write(f'{save_path}/SN_ST_model.h5ad')

# Record the end time of the process
end_time = time.time()

# Calculate the running time in seconds and hours
run_time_seconds = end_time - start_time
run_time_hours = run_time_seconds / 3600

# Print the running time of the process
print(f"Process running time: {run_time_hours:.2f} hours")
~~~

# 2. DestVI_snRNA_label_transfer_to_Stereo_seq_data

~~~python
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import ks_2samp, ttest_ind
from statsmodels.stats.multitest import multipletests
import scvi
from scvi.model import CondSCVI, DestVI
import logging
import umap
import torch
from torch.distributions import Gamma
import base64
from io import BytesIO
import argparse
import time
import anndata as ad

# Initialize logging for scvi
logger = logging.getLogger("scvi")

# Record the start time of the process
start_time = time.time()

# Set up command - line argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--slice', type=str, help='Slice information')
parser.add_argument('--species', type=str, help='Species information')
parser.add_argument('--sc_model', type=str, help='Path to the single - cell model')
parser.add_argument('--st_epochs', type=int, help='Number of epochs for spatial transcriptomics model training')
parser.add_argument('--sc_path', type=str, help='Path to the single - cell data file')
parser.add_argument('--st_path', type=str, help='Path to the spatial transcriptomics data file')
parser.add_argument('--save_path', type=str, help='Path to save the results')

# Parse the command - line arguments
args = parser.parse_args()
species = args.species
slice_info = args.slice
sc_model_path = args.sc_model
st_epochs = args.st_epochs
sc_path = args.sc_path
st_path = args.st_path
save_path = args.save_path

# Load the intersection of genes
intersect = np.load('./destvi_transfer/marmoset_sn_st.npy', allow_pickle=True)

# Read and preprocess single - cell data
sc_adata = sc.read(sc_path)
# Filter genes with less than 10 counts
sc.pp.filter_genes(sc_adata, min_counts=10)
# Save the raw counts in a separate layer
sc_adata.layers["counts"] = sc_adata.X.copy()
# Normalize the total counts of each cell to 100,000
sc.pp.normalize_total(sc_adata, target_sum=10e4)
# Apply log transformation to the normalized data
sc.pp.log1p(sc_adata)
# Save the raw data state
sc_adata.raw = sc_adata
# Add a 'feature' column to the gene metadata
sc_adata.var["feature"] = sc_adata.var.index
# Subset the data to include only genes in the intersection
sc_adata = sc_adata[:, sc_adata.var["feature"].isin(intersect)].copy()

# Load the pre - trained single - cell model
sc_model = CondSCVI.load(sc_model_path, sc_adata)

# Read and preprocess spatial transcriptomics data
st_adata = sc.read(st_path)
# Remove cells from specific areas (background, ROLF, LCNU)
st_adata = st_adata[~st_adata.obs['area_name'].astype(str).str.contains('(background|ROLF|LCNU)', case=False)]

# Find genes in intersect that are not in st_adata
setA = set(intersect)
setB = set(np.intersect1d(intersect, st_adata.var_names))
result = list(setA - setB)
# Create a new DataFrame with zero counts for missing genes
new_genes_df = pd.DataFrame(0, index=st_adata.obs_names, columns=result)
# Concatenate the new DataFrame with the original st_adata
st_adata = ad.concat([st_adata, ad.AnnData(new_genes_df)], axis=1)

# Add a 'features' column to the gene metadata
st_adata.var["features"] = st_adata.var.index
# Subset the data to include only genes in the intersection
st_adata = st_adata[:, st_adata.var["features"].isin(intersect)].copy()
print(f"Spatial transcriptomics data shape: {st_adata.shape}")

# Save the raw counts in a separate layer
st_adata.layers["counts"] = st_adata.X.copy()
print('Before normalization:')
print(st_adata)

# Normalize and transform the spatial transcriptomics data
sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
print(f'After normalization: {sum(st_adata.to_df().iloc[2,])}')

# Save the raw data state
st_adata.raw = st_adata

# Set up the AnnData object for the DestVI model
scvi.model.DestVI.setup_anndata(st_adata, layer="counts")
# Initialize the DestVI model using the pre - trained single - cell model
st_model = DestVI.from_rna_model(st_adata, sc_model, amortization="latent", vamp_prior_p=100)
# Train the DestVI model
st_model.train(max_epochs=st_epochs, plan_kwargs={"lr": 0.001})

# Plot the training history of the model
plt.figure(figsize=(8, 6))
st_model.history["elbo_train"].plot()
st_model.history["elbo_train"].iloc[5:].plot()
# Save the training history plot
plt.savefig(f'{save_path}/{slice_info}_stmodel_epoch.pdf', dpi=300)
# Save the trained spatial transcriptomics model
st_model.save(f'{save_path}/st_model_{species}_{slice_info}_all', overwrite=True)
# Get the cell type proportions from the model
st_adata.obsm["proportions"] = st_model.get_proportions()
# Save the processed spatial transcriptomics data
st_adata.write(f'{save_path}/{species}_{slice_info}_st_adata.h5ad')

# Record the end time of the process
end_time = time.time()
# Calculate the running time in seconds and hours
run_time_seconds = end_time - start_time
run_time_hours = run_time_seconds / 3600
# Print the running time of the process
print(f"Process running time: {run_time_hours:.2f} hours")
~~~

