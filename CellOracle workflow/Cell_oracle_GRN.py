
##cicero was peromed to get cis-regulatory networks

#/N/project/OHagan_single_cell/MultiOME/Chrm_404_OHagan_IUB_Multiome2_May2023/Velocity

conda activate celloracle_env
import pandas as pd
import numpy as np
import celloracle as co
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns
from celloracle import motif_analysis as ma
from tqdm.notebook import tqdm

#workflow
# identify the peaks around transcription starting sites (TSS). then merge the Cicero data with the TSS peak information and filter any peaks with weak connections to the TSS peaks.
#As such, the filtered peak data will only include TSS peaks and peaks with strong TSS connections. These will be our active promoter/enhancer elements for our base GRN

# Load scATAC-seq ISX9 peak list.
peaks = pd.read_csv("all_ISX9_peaks.csv", index_col=0)
peaks = peaks.x.values
peaks

# Load Cicero coaccessibility scores.
cicero_connections = pd.read_csv("ISX9_cicero_connections.csv", index_col=0)
cicero_connections.head()

##Annotate transcription start sites (TSSs)
ma.SUPPORTED_REF_GENOME

##my peaks is chr1-100036819-100039115 and it shound be chr1_1000368191_100039115
##do thet converstion
corrected_peaks=[]
for peak in peaks:
    corrected_peaks.append('_'.join(peak.split('-')))

corrected_peaks=np.array(corrected_peaks,dtype=object)
  

tss_annotated = ma.get_tss_info(peak_str_list=corrected_peaks, ref_genome="hg38")##bedtool module has to be loaded

# Check results
tss_annotated.tail()

##Integrate TSS info and cicero connections
chagne cicero_connection as well
corrected_peaks1=[]
corrected_peaks2=[]
for peak in cicero_connections.Peak1:
    corrected_peaks1.append('_'.join(peak.split('-')))

for peak in cicero_connections.Peak2:
    corrected_peaks2.append('_'.join(peak.split('-')))
    
##add the list to dataframe
cicero_connections['Cor_Peak1']=corrected_peaks1
cicero_connections['Cor_Peak2']=corrected_peaks2

#delet the other colmn
del cicero_connections['Peak1']
del cicero_connections['Peak2']

# Or rename the existing DataFrame (rather than creating a copy) 
cicero_connections.rename(columns={'Cor_Peak1': 'Peak1', 'Cor_Peak2': 'Peak2'}, inplace=True)

##reoder it
cicero_connections_ = cicero_connections.reindex(columns=['Peak1','Peak2','coaccess'])
 
##integrate TSS peaks annotated with coaccess cicero   
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections_)
                                               
print(integrated.shape)
integrated.head()
#peak_id” is either the TSS peak or the peaks that have a connection to a TSS peak.
#“gene_short_name” is the gene name that associated with the TSS site.
#“coaccess” is the coaccessibility score between the peak and a TSS peak. If the score is 1, it means that the peak is a TSS itself.  

##save unfilterd peaks
integrated.to_csv('integrated_TSS_Cicero_unfiltered.csv')       

#Remove peaks with weak coaccessibility scores.
peak = integrated[integrated.coaccess >= 0.8]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)      

print(peak.shape)
peak.head()
##save the filtered peak
peak.to_csv("processed_peak_file.csv") 


#####to scan for TF binding motifs. The base GRN will be generated by combining the ATAC-seq peaks and motif information.
# Load annotated peak data.
peaks = pd.read_csv("processed_peak_file.csv", index_col=0)
peaks.head()

#make sure the reference genome is installed
ref_genome = "hg38"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome,
                                             genomes_dir=None)
print(ref_genome, "installation: ", genome_installation)

##install ref. genome
if not genome_installation:
    import genomepy
    genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=None)
else:
    print(ref_genome, "is installed.")

##check peak data formate
peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=None)

##I am using here CellOracle’s default motifs during the motif analysis.
##start the motif analysis using TFinfo

#Converts a peak data into a DNA sequences.
#Scans the DNA sequences searching for TF binding motifs.
#Post-processes the motif scan results.

#If your reference genome file are installed in non-default location, please speficy the location using genomes_dir

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome,
                genomes_dir=None)

# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)

# Save tfinfo object
tfi.to_hdf5("scanned_tf.celloracle.tfinfo")

# Check motif scan results
tfi.scanned_df.head()

# filter the motifs with low scores.

# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

#Get final base GRN
df = tfi.to_dataframe()
df.head()

# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet("base_GRN_dataframe.parquet")

#Parquet can skip reading irrelevant data, resulting in faster query execution times.
# In contrast, CSV files need to read entire rows even if only a subset of columns is needed

##construct GRN using RNA count data from ISX9
##read ISX9 data
ISX9=sc.read_h5ad('ISX9.h5ad')

#To begin, we will instantiate a new Oracle object and input our gene expression data (anndata) and TF info (base GRN).
# Instantiate Oracle object
oracle = co.Oracle()

# Check data in anndata
print("Metadata columns :", list(ISX9.obs.columns))
print("Dimensional reduction: ", list(ISX9.obsm.keys()))

# Instantiate Oracle object
# Get the highly variable genes first.
sc.pp.normalize_total(ISX9, target_sum=1e4)#normalize to the libirary size
#sc.pp.log1p(ISX9)
sc.pp.highly_variable_genes(ISX9, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(ISX9, flavor='cell_ranger', n_top_genes=3000)##for cell orcle

#confirming
print('\n','Number of highly variable genes: {:d}'.format(np.sum(ISX9.var['highly_variable'])))

##select only HVG
ISX9_ = ISX9[:, ISX9.var.highly_variable]

ISX9.raw = ISX9

oracle.import_anndata_as_raw_count(adata=ISX9_,
                                   cluster_column_name='Cluster_ID',
                                   embedding_name='X_umap')

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=df) ##df = tfi.to_dataframe() from our base GRN scRNA-seq

3. KNN imputation

#CellOracle uses the same strategy as velocyto for visualizing cell transitions. This process requires KNN imputation in advance.
# Perform PCA

oracle.perform_PCA()


n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]

print(n_comps)
n_comps = min(n_comps, 50)
#Estimate the optimal number of nearest neighbors for KNN imputation.
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=10)

# Save oracle object.
oracle.to_hdf5("ISX9.celloracle.oracle")

#The next step constructs a cluster-specific GRN for all clusters.

# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="Cluster_ID", alpha=10,
                         verbose_level=10)

# Save Links object.
links.to_hdf5(file_path="UNFILTERED_ISX9_links.celloracle.links")

#Network preprocessing
#Remove uncertain network edges based on the p-value.

#Remove weak network edge. In this tutorial, we keep the top 3000 edges ranked by edge strength

links.filter_links(p=0.001, weight="coef_abs", threshold_number=3000)

# Calculate network scores.
links.get_network_score()

# Save Links object. use this file in in silico TF perturbation analysis.
links.merged_score.head()

links.to_hdf5(file_path="filtered_ISX9_links.celloracle.links")
##export LINK dictionary
ISX9.obs['Cluster_ID']
for cluster in  ISX9.obs['Cluster_ID'].unique():
    links.links_dict[cluster].to_csv(f'raw_GRN_for_{cluster}.csv')

##Network analysis
##read the data
links = co.load_hdf5(file_path="filtered_ISX9_links.celloracle.links")
oracle=co.load_hdf5('ISX9.celloracle.oracle')

# Visualize top 30-th genes with high scores.
links.plot_scores_as_rank(cluster="EEC", n_gene=30, save="ranked_score")
##analyze the difference between 2 cluster using:
#eignvector centerality...meansure the power of the node in the network
#degree of centerality...measure tje number of the edges the nodes has

#compare EEC and SPC
links.plot_score_comparison_2D(value="eigenvector_centrality",
                               cluster1="SPC", cluster2="EEC",
                               percentile=98,
                               save="score_comparison")

links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="SPC", cluster2="EEC",
                               percentile=98,
                               save="score_comparison")

##compate ELC and ENLC
links.plot_score_comparison_2D(value="eigenvector_centrality",
                               cluster1="ELC", cluster2="ENLC",
                               percentile=98,
                               save="score_comparison")

links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="ELC", cluster2="ENLC",
                               percentile=98,
                               save="score_comparison")

#visualize networks scores dynamics for NEUROG3, INSM1, and E2F1, KLF6, PBX3
links.plot_score_per_cluster(goi="PBX3", save="network_score_per_gene")

######################################################Make predictive models for simulation***********************************************************
# Here, we will need to fit the ridge regression models again. This process will take
# less time than the GRN inference in the previous notebook, because we are using the filtered GRN models.

# Make predictive models for simulation
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

# Save oracle object.
oracle.to_hdf5("ISX9_link.celloracle.oracle")
#In silico TF perturbation analysis
# Check gene expression
goi = "PBX3"
sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis",save='ISX9.PBX3.pdf')

#To simulate E2F1 KO, we will set its  expression to 0., OVEREXPRESSION OF E2F1 make it to 1.5
oracle.simulate_shift(perturb_condition={goi: 1.5},
                      n_propagation=3)

# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

#If the vectors are not visible, you can try a smaller scale parameter to magnify the vector length.
##sacale
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale = 20
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} OE")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

plt.savefig('PBX3OE_ISX9_.pdf')


#########KO FOSL1

goi = "FOSL1"
sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis",save='ISX9.Fosl3.pdf')

#To simulate Neurog3 KO, we will set its  expression to 0.
oracle.simulate_shift(perturb_condition={goi: 0.0},
                      n_propagation=3)

# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

#If the vectors are not visible, you can try a smaller scale parameter to magnify the vector length.
##sacale
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])

scale = 10
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")

plt.savefig('Fosl1KO_ISX9_scale10.pdf')
