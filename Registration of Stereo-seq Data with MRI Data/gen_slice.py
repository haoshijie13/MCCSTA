import os
import time
import numpy as np
import pandas as pd 
import nibabel as nib
from numba import jit
from tqdm import tqdm

# ----------------------------
# Coordinate Scaling Function
# ----------------------------
def scale_xy(xy_raw, scale):
    """Scale spatial coordinates by given factor"""
    return np.round(xy_raw / scale).astype(int)

# ----------------------------
# Numba-accelerated Matrix Generation
# ----------------------------
@jit(nopython=True)
def totalRNA(xy, count, mat):
    """Accumulate RNA counts in spatial matrix"""
    for coor, c in zip(xy, count):
        mat[coor[0], coor[1]] += c
    return mat

# ----------------------------
# Timing Utility Function
# ----------------------------
def time_elapsed(start, str):
    """Track processing time between operations"""
    print(f'{str}: {time.time()-start:.2f}s')
    return time.time()

# ----------------------------
# Core Data Conversion Functions
# ----------------------------
def gen_slice(path_input, path_output, scale_factor=50):
    """Convert tabular data to NIfTI spatial matrix"""
    start = time.time()
    print(f'{path_input}:')
    
    # Load and process text data
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    last_end = time_elapsed(start, 'load_txt')
    
    # Generate spatial matrix
    MIDCount = np.array(data_raw_txt.iloc[:, 3])
    xy_raw = np.array(data_raw_txt.iloc[:, 1:3])
    xy_scaled = scale_xy(xy_raw, scale_factor)
    x_max, y_max = xy_scaled[:, 0].max()+1, xy_scaled[:, 1].max()+1
    
    MIDCount_Matrix = np.zeros((x_max, y_max))
    Mat = totalRNA(xy_scaled, MIDCount, MIDCount_Matrix)
    last_end = time_elapsed(last_end, 'gen_mat')
    
    # Save NIfTI output
    nifti_img = nib.Nifti1Image(Mat, np.eye(4))
    nib.save(nifti_img, path_output)
    time_elapsed(last_end, 'save_nifti')
    print('\n')
    return True

def rescale(path_nii_slice, path_output, edge=1000):
    """Rescale NIfTI image to standardized dimensions"""
    raw = nib.load(path_nii_slice)
    slice = raw.get_fdata()
    new_slice = np.zeros((edge, edge))
    new_slice[:slice.shape[0], :slice.shape[1]] = slice[:1000, :1000]
    nib.save(nib.Nifti1Image(new_slice, np.eye(4)), path_output)

# ----------------------------
# Volume Processing Functions
# ----------------------------
def nii_3Dto2D(path_nii_volume, path_output, slice_dim=2, prefix=''):
    """Convert 3D NIfTI volume to 2D slice sequence"""
    raw = nib.load(path_nii_volume)
    volume = raw.get_fdata()
    new_shape = [volume.shape[slice_dim]] + [s for i,s in enumerate(volume.shape) if i != slice_dim]
    
    for i, slice in tqdm(enumerate(volume.transpose([slice_dim] + [i for i in range(3) if i != slice_dim]))):
        nib.save(nib.Nifti1Image(slice, np.eye(4)), 
                os.path.join(path_output, f'{prefix}_{i:02d}.nii.gz'))

def flip_rotate(path_nii_slice, path_output, manipulation=0):
    """Apply spatial transformations to NIfTI slice"""
    raw = nib.load(path_nii_slice)
    slice = raw.get_fdata()
    
    if manipulation != 0:
        if manipulation < 0:
            slice = np.flip(slice, axis=1)
        slice = np.rot90(slice, abs(manipulation))
    
    nib.save(nib.Nifti1Image(slice, np.eye(4)), path_output)

# ----------------------------
# Mask Application Functions
# ----------------------------
@jit(nopython=True)
def xy_get_mask(xy, mask, xy_new_property):
    """Map mask values to spatial coordinates"""
    for i, coor in enumerate(xy):
        if coor[0] < 1600 and coor[1] < 1600:
            xy_new_property[i] = mask[coor[0], coor[1]]
    return xy_new_property

def gen_txt_onecol(path_raw, path_mask, path_output, new_property='mask', scale_factor=50):
    """Annotate spatial data with mask values"""
    start = time.time()
    print(f'{path_raw}:')
    
    # Load input data
    data = pd.read_csv(path_raw, delimiter='\t', comment='#', header=0)
    mask = nib.load(path_mask).get_fdata()
    last_end = time_elapsed(start, 'load_data')
    
    # Process coordinates
    xy_scaled = scale_xy(np.array(data.iloc[:, 1:3]), scale_factor)
    xy_new_property = np.zeros(xy_scaled.shape[0], dtype=int)
    
    # Apply mask
    xy_new_property = xy_get_mask(xy_scaled, np.squeeze(mask.astype(int)), xy_new_property)
    last_end = time_elapsed(last_end, 'apply_mask')
    
    # Save annotated data
    pd.DataFrame({new_property: xy_new_property}).to_csv(path_output, sep='\t', index=False)
    time_elapsed(last_end, 'save_results')
    print('\n')
    return True

# ----------------------------
# Specialized Processing Functions
# ----------------------------
def gen_slice_masked(path_input, path_output, scale_factor=50):
    """Generate masked spatial matrix from annotated data"""
    start = time.time()
    print(f'{path_input}:')
    
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    MIDCount = np.array(data_raw_txt.iloc[:, 3]) * np.array(data_raw_txt.iloc[:, 4])
    
    xy_scaled = scale_xy(np.array(data_raw_txt.iloc[:, 1:3]), scale_factor)
    x_max, y_max = xy_scaled[:, 0].max()+1, xy_scaled[:, 1].max()+1
    
    Mat = totalRNA(xy_scaled, MIDCount, np.zeros((x_max, y_max)))
    nib.save(nib.Nifti1Image(Mat, np.eye(4)), path_output)
    
    time_elapsed(start, 'processed')
    print('\n')
    return True

# ----------------------------
# ROI Processing Utilities
# ----------------------------
@jit(nopython=True)
def make_ROI(xy, count, mat):
    """Direct ROI assignment without accumulation"""
    for coor, c in zip(xy, count):
        mat[coor[0], coor[1]] = c
    return mat

@jit(nopython=True)
def make_ROI_repeat(xy, count, mat):
    """Accumulative ROI assignment"""
    for coor, c in zip(xy, count):
        mat[coor[0], coor[1]] += c
    return mat

# ----------------------------
# Cell-type Specific Processing
# ----------------------------
def gen_slice_cell(path_input, path_output, scale_factor=50, data_col=3, data_x=6, data_y=7):
    """Process cell-type specific spatial data"""
    start = time.time()
    print(f'{path_input}:')
    
    data_raw_txt = pd.read_csv(path_input, delimiter='\t', comment='#', header=0)
    xy_scaled = scale_xy(np.array(data_raw_txt.iloc[:, data_x:data_y+1]), scale_factor)
    
    x_max, y_max = xy_scaled[:, 0].max()+1, xy_scaled[:, 1].max()+1
    Mat = make_ROI(xy_scaled, np.array(data_raw_txt.iloc[:, data_col]), np.zeros((x_max, y_max)))
    
    nib.save(nib.Nifti1Image(Mat, np.eye(4)), path_output)
    time_elapsed(start, 'processed')
    print('\n')
    return True

# ----------------------------
# Parquet Data Handling
# ----------------------------
def gen_slice_from_parquet(path_input, path_output, x_column='x', y_column='y', 
                          umi_count_column='umi_count', scale_factor=50):
    """Process spatial data from Parquet format"""
    data_raw_txt = pd.read_parquet(path_input)
    xy_scaled = scale_xy(np.array(data_raw_txt[[x_column, y_column]]), scale_factor)
    
    x_max, y_max = xy_scaled[:, 0].max()+1, xy_scaled[:, 1].max()+1
    Mat = totalRNA(xy_scaled, np.array(data_raw_txt[umi_count_column]), np.zeros((x_max, y_max)))
    
    nib.save(nib.Nifti1Image(Mat, np.eye(4)), path_output)
    return True
