# ----------------------------
# Module Imports
# ----------------------------
import os
import time
import numpy as np
import nibabel as nib
import pandas as pd
from tqdm import tqdm

# Custom module imports
from gen_slice import (
    gen_slice, rescale, time_elapsed, gen_slice_masked, scale_xy,
    gen_rxy_slice, gen_rxy_slice_ROI, gen_slice_cell, nii_3Dto2D, flip_rotate,
    gen_txt, gen_slice_ROI, make_ROI, gen_slice_ROI_from_parquet, make_ROI_repeat,
    gen_slice_ROI_onecol, gen_txt_onecol, gen_slice_from_parquet
)

# ----------------------------
# Part 1: Transcriptomic Text Data -> NIfTI Conversion
# ----------------------------
# Input paths
path_raw = r'/Transcriptome/marmoset-20230418-v5/txt'
files_txt = sorted([f for f in os.listdir(path_raw) if 'txt' in f])
files_id = [f.split('.')[0] for f in files_txt]

# Output configuration
path_sagittal_slices = r'/Transcriptome/marmoset-20230418-v5/txt_to_nifti'
os.makedirs(path_sagittal_slices, exist_ok=True)

# Batch processing loop
for f in files_txt:
    input_path = os.path.join(path_raw, f)
    output_path = os.path.join(path_sagittal_slices, f'{f}.nii.gz')
    gen_slice(input_path, output_path, 50)

# ----------------------------
# Part 2: Cell Label Data -> NIfTI Conversion
# ----------------------------
# Updated paths for cell labeling data
path_raw = r'/marmoset_V2_20220627/marmoset_cortex_cell_label/SpatialID_result'
files_txt = sorted([f for f in os.listdir(path_raw) if 'T' in f])
files_id = [f.split('.')[0] for f in files_txt]

# Output configuration
path_sagittal_slices = r'/marmoset_V2_20220627/marmoset_cortex_cell_label/tonifti/bin5'
os.makedirs(path_sagittal_slices, exist_ok=True)

# Processing method 1: Cell label slices (5μm)
for f in files_txt:
    input_path = os.path.join(path_raw, f)
    output_path = os.path.join(path_sagittal_slices, f'{f}.nii.gz')
    gen_slice_cell(input_path, output_path, 5)

# Processing method 2: RXY coordinate slices (100μm)
for f in files_txt:
    input_path = os.path.join(path_raw, f)
    output_path = os.path.join(path_sagittal_slices, f'{f}.nii.gz')
    gen_rxy_slice(input_path, output_path, 100)

# ----------------------------
# Part 3: Data Validation & Feature Extraction
# ----------------------------
# Sample data loading (for debugging)
# input_path = <specify actual path>
data_raw_txt = pd.read_csv(input_path, delimiter=',', comment='#', header=0)

# Data structure verification
print(data_raw_txt.head())
print(np.array(data_raw_txt.iloc[:, 2:4]))

# ----------------------------
# Part 4: Mask Relabeling (NIfTI Mask -> Text Annotation)
# ----------------------------
# Path configuration
path_txt_raw = r'/marmoset_V2_20220627/coronal/lgn'
path_nii_mask = r'/marmoset_V2_20220627/coronal_2Dslice_nifti/small/lgn_mask'
path_txt_mask = r'/marmoset_V2_20220627/coronal_2Dslice_nifti/small/lgn_mask_to_txt'

# File path matching
path_txt = sorted([os.path.join(path_txt_raw, f) for f in os.listdir(path_txt_raw)])
path_mask = sorted([os.path.join(path_nii_mask, f) for f in os.listdir(path_nii_mask)])
path_output = sorted([os.path.join(path_txt_mask, f'ROI{f}') for f in os.listdir(path_txt_raw)])

# Batch processing with progress bar
for txt, mask, out in tqdm(zip(path_txt, path_mask, path_output)):
    gen_txt_onecol(txt, mask, out, scale_factor=50)

# ----------------------------
# Part 5: Mask Application Testing
# ----------------------------
# Test configuration
path_raw = r'/marmoset_V2_20220627/coronal_brain_paxROI_txt'
files_txt = sorted(os.listdir(path_raw))
files_id = [f.split('.')[0] for f in files_txt]

# Output configuration
path_sagittal_slices = r'/marmoset_V2_20220627/coronal_brain_paxROI_txt_to_nifti'
os.makedirs(path_sagittal_slices, exist_ok=True)

# Masked slice generation test
for f in files_txt:
    input_path = os.path.join(path_raw, f)
    output_path = os.path.join(path_sagittal_slices, f'{f}.nii.gz')
    gen_slice_masked(input_path, output_path, 50)