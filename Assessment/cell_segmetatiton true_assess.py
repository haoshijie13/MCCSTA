import os
import pandas as pd
import numpy as np
import cv2
from skimage import io
import matplotlib.pyplot as plt
from scipy.spatial import distance
from skimage.measure import label, regionprops
from sklearn.neighbors import KDTree

# Function to convert an image to a DataFrame
def imag2df2(image):
    """
    Convert an image to a DataFrame containing x, y coordinates and label information.

    Args:
        image (np.ndarray): Input image.

    Returns:
        pd.DataFrame: DataFrame with columns 'x', 'y', and 'label'.
    """
    test = label(image[:, :, 0])
    y, x = np.where(test + 1)
    df = pd.DataFrame({"x": x, "y": y, "label": test[y, x]})
    return df

# Function to calculate the minimum distance
def calculate_min_distance(row, label_centers_mean):
    """
    Calculate the minimum distance between a row and the centers of labels.

    Args:
        row (pd.Series): A row from the DataFrame containing 'x' and 'y' coordinates.
        label_centers_mean (pd.DataFrame): DataFrame containing the centers of labels.

    Returns:
        float: Minimum distance.
    """
    distances = np.sqrt((label_centers_mean['x'] - row['x'])**2 + (label_centers_mean['y'] - row['y'])**2)
    return distances.min()

# Function to query the minimum distance using KDTree
def query_distance(ai_center, mask_center):
    """
    Query the minimum distance between AI centers and mask centers using KDTree.

    Args:
        ai_center (pd.DataFrame): DataFrame containing the centers of AI labels.
        mask_center (pd.DataFrame): DataFrame containing the centers of mask labels.

    Returns:
        pd.DataFrame: AI centers DataFrame with an additional 'min_distance' column.
    """
    X_ai = np.array(ai_center[['x', 'y']])
    X_mask = np.array(mask_center[['x', 'y']])
    kdt = KDTree(X_mask, leaf_size=30, metric='euclidean')
    result = kdt.query(X_ai, k=1, return_distance=True)
    ai_center['min_distance'] = result[0]
    return ai_center

# Function to determine if AI labels are close to mask labels
def getTrue2(ai, mask, threshold = 20):
    """
    Determine if AI labels are close to mask labels based on a distance threshold.

    Args:
        ai (pd.DataFrame): DataFrame containing AI label information.
        mask (pd.DataFrame): DataFrame containing mask label information.
        threshold (float, optional): Distance threshold. Defaults to 20.

    Returns:
        pd.DataFrame: AI DataFrame with additional 'is_close' and 'pre' columns.
    """
    label_centers_mean = mask[mask['label'] != 0].groupby('label').agg({'x': 'mean', 'y': 'mean'})
    df_ai_center = ai[ai['label'] != 0].groupby('label').agg({'x': 'mean', 'y': 'mean'}).reset_index()
    label_centers_mean['x'] = label_centers_mean['x'].astype(float)
    label_centers_mean['y'] = label_centers_mean['y'].astype(float)
    df_ai_center['x'] = df_ai_center['x'].astype(float)
    df_ai_center['y'] = df_ai_center['y'].astype(float)
    df_ai_center = query_distance(df_ai_center, label_centers_mean)
    df_ai_center['is_close'] = df_ai_center['min_distance'] < threshold
    ai['is_close'] = ai['label'].map(df_ai_center.set_index('label')['is_close'])
    ai['is_close'] = ai['is_close'].fillna(False)
    ai.loc[(ai['is_close'] == False) & (mask['label'] != 0), 'is_close'] = True
    ai['is_close'] = ai['is_close'].astype(int)
    ai['pre'] = ai['label'].apply(lambda x: 0 if x == 0 else 1)
    return ai

# Directory containing the mask and ai files
data_dir = "/data/fig"
# Lists to store precision and recall values
precision_list = []
recall_list = []
# Lists to store file names
mask_file_list = []
ai_file_list = []

# Loop through all files in the directory
for filename in os.listdir(data_dir):
    if filename.endswith(".mask.tif"):
        # Get the corresponding ai file name
        ai_filename = filename.replace(".mask.tif", ".ai.tif")
        ai_file_path = os.path.join(data_dir, ai_filename)
        mask_file_path = os.path.join(data_dir, filename)
        # Check if the corresponding ai file exists
        if os.path.exists(ai_file_path):
            # Read the mask and ai images
            test = cv2.imread(mask_file_path)
            test1 = cv2.imread(ai_file_path)
            # Convert the images to DataFrames
            mask = imag2df2(test)
            ai = imag2df2(test1)
            # Determine if AI labels are close to mask labels
            ai = getTrue2(ai, mask)
            # Calculate the true positive count
            TP = ((ai['is_close'] == 1) & (ai['pre'] == 1)).sum()
            # Calculate the false negative count
            FN = ((ai['is_close'] == 1) & (ai['pre'] == 0)).sum()
            # Calculate the true negative count
            TN = ((ai['is_close'] == 0) & (ai['pre'] == 0)).sum()
            # Calculate the false positive count
            FP = ((ai['is_close'] == 0) & (ai['pre'] == 1)).sum()
            # Calculate the recall
            R = TP / (TP + FN)
            # Calculate the precision
            P = TP / (TP + FP)
            # Append the precision and recall values to the lists
            precision_list.append(P)
            recall_list.append(R)
            # Append the file names to the lists
            mask_file_list.append(filename)
            ai_file_list.append(ai_filename)

# Create a DataFrame to store the results
results_df = pd.DataFrame({
    'mask_file': mask_file_list,
    'ai_file': ai_file_list,
    'precision': precision_list,
    'recall': recall_list
})

# Save the results to a CSV file
results_df.to_csv('results.csv', index=False)  