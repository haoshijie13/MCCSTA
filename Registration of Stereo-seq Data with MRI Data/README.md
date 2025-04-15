# ​**​Marmoset Spatial Transcriptomics-MRI Integration​**​  
*Combining single-cell resolution gene expression with neuroimaging in 3D space*

---

## 🔬 ​**​Project Overview​**​  
This pipeline integrates ​**​spatial transcriptomics​**​ (50μm resolution) and ​**​high-field MRI​**​ (25μm isotropic) data from marmoset brain through:  
- *Multi-modal image registration*  
- *3D gene expression mapping*  
- *Cellular density modeling*  
- *Cross-modal feature correlation*  

---

## 📂 ​**​Pipeline Architecture​**​

### 🧬 ​**​Core Modules​**​

#### ​**​1. Spatial Genomics Processing​**​  
*Document 1-2*:  
- Processes raw sequencing data → AnnData objects  
- Generates binned spatial matrices (NIfTI format)  
- Visualizes total RNA distributions  

#### ​**​2. MRI-Gene Registration​**​  
*Document 3-7*:  
- ​**​ANTs-based 2D/3D registration​**​  
- Slice-wise affine + diffeomorphic mapping  
- Intensity normalization & QC checks  
- Transforms applied to 5000+ genes  

#### ​**​3. Cellular Modeling​**​  
*Document 8-9*:  
- 3D Gaussian kernel density estimation  
- Cortical layer-specific distributions  
- Missing data imputation  

#### ​**​4. Multivariate Analysis​**​  
*Document 11*:  
- Partial Least Squares (PLS) regression  
- Gene-network correlation mapping  
- Contribution thresholding (r > 0.95)  

---

## 🛠️ ​**​Key Technical Components​**​

### ​**​Data Structures​**​  
| Type | Format | Description |  
|------|--------|-------------|  
| Gene Data | AnnData | Scanpy-compatible matrices |  
| MRI Data | NIfTI | 80μm → 25μm resampled |  
| Masks | NRRD | Cortical ROI annotations |  

### ​**​Registration Workflow​**​  
```mermaid
graph TD
    A[Raw Sequencing] --> B(Spatial Binning)
    B --> C{NIfTI Conversion}
    C --> D[2D Slice Registration]
    D --> E[3D Volume Reconstruction]
    E --> F[Template Alignment]
