# Single-Cell RNA-Seq Label Transfer Pipeline

This repository contains two Jupyter notebooks for annotating single-cell RNA-seq data using reference atlases.

---

## **Notebook 1: HLCA Atlas Annotation**  

### **Purpose**  
Prepares and annotates single-cell RNA-seq data using the Human Lung Cell Atlas (HLCA) as a reference.  

### **Inputs**  
- Query dataset (single-cell RNA-seq counts matrix in `.h5ad` or similar format).  
- Preprocessed HLCA reference atlas (provided in the notebook).  

### **Outputs**  
- Annotated query dataset with cell type labels transferred from the HLCA.  
- Visualizations (UMAP/t-SNE plots with transferred labels).  

### **Workflow**  
1. Load and preprocess the query dataset.  
2. Align query data to the HLCA reference using PCA and harmony integration.  
3. Transfer cell type labels using `scANVI` or `Symphony`.  
4. Generate visualization plots.  

---

## **Notebook 2: scRNA-seq Label Transfer**  

### **Purpose**  
Transfers cell type labels from a reference dataset to a query dataset using Seurat or SCANVI. (not specific to HLCA) 

### **Inputs**  
- Query dataset (counts matrix in `.h5ad`, `.loom`, or `.rds` format).  
- Reference dataset (pre-annotated single-cell data, e.g., HLCA or custom).  

### **Outputs**  
- Query dataset with predicted cell type labels.  
- Confusion matrices or label transfer confidence scores.  
- UMAP plots showing reference and query integration.  

### **Workflow**  
1. Normalize and log-transform query data.  
2. Integrate query and reference using PCA/CCA (Seurat) or variational inference (SCANVI).  
3. Transfer labels and assess confidence.  
4. Visualize results.  


---


For detailed execution, refer to the notebooks' inline documentation. 

### Contact

Hao Cheng, Hao.Cheng@osumc.edu
Cankun, Cankun.Wang@osumc.edu
Qi, Qi.Guo@osumc.edu
