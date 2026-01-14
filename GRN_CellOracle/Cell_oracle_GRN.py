import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co

####load your scRNA-seq (has row count)
adata=sc.read_h5ad('Pancrease_acinar.h5ad')

##normalize your data
adata.raw = adata.copy()  # Preserve raw counts before normalization

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
##extract the HVG
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000
)

adata = adata[:, adata.var.highly_variable]


# Renormalize after filtering
sc.pp.normalize_per_cell(adata)

# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()


# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')


sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)

##optional to run diffmap
sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

sc.pl.diffmap(adata,save='_diffmap_HVG_.pdf',components=['2,3'],color='Sample_2')





adata.write_h5ad('adata_acinar_3000.h5ad')

# Instantiate Oracle object

##use build in base GRN 
base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN()####you could load human and different organism depend on your model


oracle = co.Oracle()

adata.X = adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="Sample_2",###column in your metadata, could be celltype, cluster (to generate cluster specific GRN)
                                   embedding_name="X_umap")####if you want to use diffmap you could replace X_umap with X_diffmap


# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)


# Perform PCA
oracle.perform_PCA()

n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=20)



#cluster-specific GRN for all clusters.
links = oracle.get_links(cluster_name_for_GRN_unit="Sample_2", alpha=10,
                         verbose_level=10)



links.to_hdf5(file_path="acinar_Sample_2_GRN.celloracle.links")##save the file


##GRN processing
#The raw network data is stored in the links_dict attribute,
#  while the filtered network data is stored in the filtered_links attribute.
#Remove weak network edge.  keep the top 5000 edges ranked by edge strength.
#links.filter_links(p=0.001, weight="coef_abs", threshold_number=5000)



##load the GRN
# load files 

links = co.load_hdf5(file_path="acinar_Sample_2_GRN.celloracle.links")
links.filter_links(p=0.05, weight="coef_abs")

#
links.links_dict.keys()#####this is dictionary and you can save it to csv using links[key].to_csv('your_GRN.csv')
