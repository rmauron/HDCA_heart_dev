"""
title: "11_HDCA_heart_scVelo"
author: "Raphaël Mauron"
"""

"""
In this script, RNA velocity is performed on different combination of High and Deep Level clusters.

"""

# Import libraries
import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import anndata as ad
import os
from scipy.sparse import csr_matrix
print(ad.__version__)
import h5py
import warnings
import pandas as pd
from matplotlib import pyplot as plt


# Setups
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
warnings.simplefilter("ignore", category=UserWarning)
wd = os.getcwd()
save_path = wd + "/hdca_DevHeart/output/HDCA_heart_sc_analysis_docker/trajectory/"


# Import spliced and unspliced objects created from Seurat 
spliced = wd + "/hdca_DevHeart/data/AnnData/data_spliced.h5ad"
spliced_data = ad.read_h5ad(spliced)

unspliced = wd + "/hdca_DevHeart/data/AnnData/data_unspliced.h5ad"
unspliced_data = ad.read_h5ad(unspliced)

# Create AnnData object & add spliced and unspliced objects as layers of the AnnData
adata = spliced_data.copy()
adata.layers['spliced'] = spliced_data.X
adata.layers['unspliced'] = unspliced_data.X

# Return expression matrices
adata.X

######################################################
# High Level
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_HL = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta = pd.read_table(wd + "/hdca_DevHeart/data/AnnData/metadata.tsv", sep='\t', index_col=0)
meta
adata_HL = adata_HL[meta.index]
adata_HL.obs = pd.concat([adata_HL.obs, meta], axis=1)
adata_HL
adata_HL.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_HL = pd.read_csv(wd + "/hdca_DevHeart/data/AnnData/UMAP.csv", delimiter=',',index_col=0)
UMAP_HL = UMAP_HL.to_numpy()
adata_HL.obsm["X_umap"]=UMAP_HL

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_HL, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_HL)
scv.pp.filter_genes_dispersion(adata_HL, n_top_genes=2000)
scv.pp.log1p(adata_HL)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
## adata_HL_temp = adata_HL[(adata_HL.obs.high_level_clusters.isin(["4_35", "5", "9", "11", "17"]))]
adata_HL_temp = adata_HL[(adata_HL.obs.age_groups.isin(["12-14"]))]
adata_HL_temp = adata_HL_temp[(adata_HL_temp.obs.high_level_clusters.isin(["4_35", "5", "9", "11", "17"]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_HL_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_HL_temp)

# Plot the velocity on the UMAP imported
# scv.pl.velocity_embedding_stream(adata_HL_temp, basis='umap',  color=["high_level_clusters"], save=save_path+"HL_4_35_5_9_11_17_umap.svg")
scv.pl.velocity_embedding_stream(adata_HL_temp, basis='pca', color=["high_level_clusters"], save=save_path+"HL_age_12-14_4_35_5_9_11_17_pca.svg")




######################################################
# Fibroblast
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_FB = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_FB = pd.read_table(wd + "/hdca_DevHeart/data/AnnData/metadata_FB.tsv", sep='\t', index_col=0)
meta_FB
adata_FB = adata_FB[meta_FB.index]
adata_FB.obs = pd.concat([adata_FB.obs, meta_FB], axis=1)
adata_FB
adata_FB.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_FB = pd.read_csv(wd + "/hdca_DevHeart/data/AnnData/UMAP_FB.csv", delimiter=',',index_col=0)
UMAP_FB = UMAP_FB.to_numpy()
adata_FB.obsm["X_umap"]=UMAP_FB

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_FB, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_FB)
scv.pp.filter_genes_dispersion(adata_FB, n_top_genes=2000)
scv.pp.log1p(adata_FB)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_FB_temp = adata_FB[(adata_FB.obs.age_groups.isin(["12-14"]))]
adata_FB_temp = adata_FB_temp[(adata_FB_temp.obs.clusters_subset.isin([3, 7, 8, 18]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_FB_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_FB_temp)

# Plot the velocity on the UMAP imported
scv.pl.velocity_embedding_stream(adata_FB_temp, basis='pca', color=["clusters_subset"], save=save_path+"FB_age_12-14_3_7_8_18_pca.svg")



#sc.tl.pca(adata_FB, random_state=0)
#sc.pp.neighbors(adata_FB, random_state=0)

#sc.pl.embedding(adata_FB, basis="umap", color=["cell.types", "clusters_subset"])


######################################################
# Endothelial
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_EN = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_EN = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_EN.tsv", sep='\t', index_col=0)
meta_EN
adata_EN = adata_EN[meta_EN.index]
adata_EN.obs = pd.concat([adata_EN.obs, meta_EN], axis=1)
adata_EN
adata_EN.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_EN = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_EN.csv", delimiter=',',index_col=0)
UMAP_EN = UMAP_EN.to_numpy()
adata_EN.obsm["X_umap"]=UMAP_EN

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_EN, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_EN)
scv.pp.filter_genes_dispersion(adata_EN, n_top_genes=2000)
scv.pp.log1p(adata_EN)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_EN_temp = adata_EN[(adata_EN.obs.age_groups.isin(["12-14"]))]
adata_EN_temp = adata_EN_temp[(adata_EN_temp.obs.clusters_subset.isin([15,13,10,2]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_EN_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_EN_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_EN_temp, basis='pca', color=["clusters_subset"], save=save_path+"EN_age_12-14_2_10_13_15_pca.svg")



######################################################
# Cardiomyocyte
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_CM = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_CM = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_CM.tsv", sep='\t', index_col=0)
meta_CM
adata_CM = adata_CM[meta_CM.index]
adata_CM.obs = pd.concat([adata_CM.obs, meta_CM], axis=1)
adata_CM
adata_CM.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_CM = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_CM.csv", delimiter=',',index_col=0)
UMAP_CM = UMAP_CM.to_numpy()
adata_CM.obsm["X_umap"]=UMAP_CM

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_CM, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_CM)
scv.pp.filter_genes_dispersion(adata_CM, n_top_genes=2000)
scv.pp.log1p(adata_CM)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_CM_temp = adata_CM[(adata_CM.obs.clusters_subset.isin([18,8]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_CM_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_CM_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_CM_temp, basis='umap',  color=["cell.types", "clusters_subset_split"], save=save_path+"CM_22_24_pca.svg")



######################################################
# Innervation
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_IN = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_IN = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_IN.tsv", sep='\t', index_col=0)
meta_IN
adata_IN = adata_IN[meta_IN.index]
adata_IN.obs = pd.concat([adata_IN.obs, meta_IN], axis=1)
adata_IN
adata_IN.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_IN = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_IN.csv", delimiter=',',index_col=0)
UMAP_IN = UMAP_IN.to_numpy()
adata_IN.obsm["X_umap"]=UMAP_IN

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_IN, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_IN)
scv.pp.filter_genes_dispersion(adata_IN, n_top_genes=2000)
scv.pp.log1p(adata_IN)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_IN_temp = adata_IN[(adata_IN.obs.age_groups.isin(["12-14"]))]
adata_IN_temp = adata_IN[(adata_IN.obs.clusters_subset.isin([1, 3, 4, 7, 8, 11, 12, 13, 17, 19]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_IN_temp, random_state=0, use_rep='X')

# Calculate velocity graphs
scv.tl.velocity_graph(adata_IN_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_IN_temp, basis='pca', color=["cell.types", "clusters_subset"], save=save_path+"IN_1_3_4_7_8_11_12_13_17_19_age_12-14_pca.svg")

scv.pl.velocity_embedding_stream(adata_IN_temp, basis='pca', components= '1,2', size = 80, alpha=0.7, color=["clusters_subset", "latent_time"])

scv.tl.recover_dynamics(adata_IN_temp)
scv.tl.latent_time(adata_IN_temp)
scv.pl.scatter(adata_IN_temp, basis='umap', color='latent_time', color_map='gnuplot', size=80)

######################################################
# Innervation SELECTION
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_IN_sub = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_IN_sub = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_IN_selection.tsv", sep='\t', index_col=0)
meta_IN_sub
adata_IN_sub = adata_IN_sub[meta_IN_sub.index]
adata_IN_sub.obs = pd.concat([adata_IN_sub.obs, meta_IN_sub], axis=1)
adata_IN_sub
adata_IN_sub.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_IN_sub = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_IN_selection.csv", delimiter=',',index_col=0)
UMAP_IN_sub = UMAP_IN_sub.to_numpy()
adata_IN_sub.obsm["X_umap"]=UMAP_IN_sub

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_IN_sub, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_IN_sub)
scv.pp.filter_genes_dispersion(adata_IN_sub, n_top_genes=4000) #2000 first, here to try
scv.pp.log1p(adata_IN_sub)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
#adata_IN_sub_temp = adata_IN_sub[(adata_IN_sub.obs.age_groups.isin(["12-14"]))]
# adata_IN_sub_temp = adata_IN_sub[(adata_IN_sub.obs.clusters_subset.isin([3, 13, 17]))]
adata_IN_sub_temp = adata_IN_sub[(adata_IN_sub.obs.clusters_subset.isin([4, 7, 8, 11, 19]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_IN_sub_temp, random_state=0, use_rep='X')

# Calculate velocity graphs
scv.tl.velocity_graph(adata_IN_sub_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_IN_sub_temp, basis='pca' , components= '1,2', size = 300, alpha=0.7, color=["clusters_subset", "latent_time"])
scv.pl.velocity_embedding_stream(adata_IN_sub_temp, basis='umap', size = 300, alpha=0.7, color=["clusters_subset", "latent_time"])


scv.pl.velocity_embedding_stream(adata_IN_sub_temp, basis='umap', size = 500, color=["cell.types", "clusters_subset"], save=save_path+"IN_sub_1_3_13_17_umap_4000feat.svg")


## Latent time
scv.tl.recover_dynamics(adata_IN_sub_temp)
scv.tl.latent_time(adata_IN_sub_temp)
scv.pl.scatter(adata_IN_sub_temp, color='latent_time', color_map='gnuplot', size=80)


######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################

# Based on cell-rank

######################################################
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_IN_sub = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_IN_sub = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_IN_selection.tsv", sep='\t', index_col=0)
meta_IN_sub
adata_IN_sub = adata_IN_sub[meta_IN_sub.index]
adata_IN_sub.obs = pd.concat([adata_IN_sub.obs, meta_IN_sub], axis=1)
adata_IN_sub
adata_IN_sub.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_IN_sub = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_IN_selection.csv", delimiter=',',index_col=0)
UMAP_IN_sub = UMAP_IN_sub.to_numpy()
adata_IN_sub.obsm["X_umap"]=UMAP_IN_sub

#adata_IN_sub = adata_IN_sub[(adata_IN_sub.obs.clusters_subset.isin([4, 7, 8, 11, 12, 19]))]


scv.pp.filter_and_normalize(adata_IN_sub, min_shared_counts=20, n_top_genes=2000, subset_highly_variable=False)

sc.tl.pca(adata_IN_sub)
sc.pp.neighbors(adata_IN_sub, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata_IN_sub, n_pcs=None, n_neighbors=None)

scv.tl.recover_dynamics(adata_IN_sub, n_jobs=8)
scv.tl.velocity(adata_IN_sub, mode="dynamical")


vk = cr.kernels.VelocityKernel(adata_IN_sub)
vk.compute_transition_matrix()

ck = cr.kernels.ConnectivityKernel(adata_IN_sub)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
print(combined_kernel)

combined_kernel.plot_projection()
combined_kernel.plot_random_walks(start_ixs={"clusters_subset": 17}, max_iter=200, seed=0)

# scv.tl.velocity_graph(adata_IN_sub)
# scv.pl.velocity_embedding_stream(adata_IN_sub , basis='umap', size = 500, color=["cell.types", "clusters_subset"])


# scv.tl.latent_time(adata_IN_sub)
# scv.pl.scatter(adata_IN_sub, basis='pca', color='latent_time', color_map='gnuplot', size=80)

######################################################
# Deep Level & High Level
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_comb = adata.copy()

# Append HIGH LEVEL metadata
meta = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata.tsv", sep='\t', index_col=0)
meta
#adata_comb = adata_comb[meta_CM.index] #to keep only cells present in metadata.index
adata_comb.obs = pd.concat([adata_comb.obs, meta], axis=1)
adata_comb
adata_comb.obs


# Append DEEP LEVEL metadata
## CM
meta_CM = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_CM.tsv", sep='\t', index_col=0)
meta_CM["clusters_subset_split"] #what we want to append to adata_comb.obs
meta_CM["clusters_subset_split"] = meta_CM["clusters_subset_split"].apply(lambda x: str(x))
# rename column
meta_CM=meta_CM.rename(columns={"clusters_subset_split": "clusters_subset_split_CM"})
#adata_comb = adata_comb[meta_CM.index] #to keep only cells present in metadata.index
adata_comb.obs = pd.concat([adata_comb.obs, meta_CM["clusters_subset_split_CM"]], axis=1)
adata_comb
adata_comb.obs


## FB
meta_FB = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_FB.tsv", sep='\t', index_col=0)
meta_FB["clusters_subset"] #what we want to append to adata_comb.obs
meta_FB["clusters_subset"] = meta_FB["clusters_subset"].apply(lambda x: str(x))
# rename column
meta_FB=meta_FB.rename(columns={"clusters_subset": "clusters_subset_FB"})
#adata_comb = adata_comb[meta_CM.index] #to keep only cells present in metadata.index
adata_comb.obs = pd.concat([adata_comb.obs, meta_FB["clusters_subset_FB"]], axis=1)
adata_comb
adata_comb.obs


## EN
meta_EN = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_EN.tsv", sep='\t', index_col=0)
meta_EN["clusters_subset"] #what we want to append to adata_comb.obs
meta_EN["clusters_subset"] = meta_EN["clusters_subset"].apply(lambda x: str(x))
# rename column
meta_EN=meta_EN.rename(columns={"clusters_subset": "clusters_subset_EN"})
#adata_comb = adata_comb[meta_CM.index] #to keep only cells present in metadata.index
adata_comb.obs = pd.concat([adata_comb.obs, meta_EN["clusters_subset_EN"]], axis=1)
adata_comb
adata_comb.obs

# Import and append the UMAP from High Level Seurat to the High Level AnnData 
UMAP_HL = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP.csv", delimiter=',',index_col=0)
UMAP_HL = UMAP_HL.to_numpy()
adata_comb.obsm["X_umap"]=UMAP_HL


# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_comb, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_comb)
scv.pp.filter_genes_dispersion(adata_comb, n_top_genes=2000)
scv.pp.log1p(adata_comb)


# Define the clusters you want to keep
to_keep_CM = []
to_keep_EN = []
to_keep_FB = ["3", "5", "7", "8"]
to_keep_HL = ["20"]

# Create the final condition based on your filtering criteria
final_condition = (
    adata_comb.obs['clusters_subset_split_CM'].isin(to_keep_CM) |
    adata_comb.obs['clusters_subset_EN'].isin(to_keep_EN) |
    adata_comb.obs['clusters_subset_FB'].isin(to_keep_FB) |
    adata_comb.obs['high_level_clusters_removed'].isin(to_keep_HL)
)

#adata_subset = adata_comb[adata_comb.obs['clusters_subset_EN'].isin(['2', '3'])]


# Use the final condition to filter and split the AnnData object
adata_comb_temp = adata_comb[final_condition, :]


# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_comb_temp, random_state=0, use_rep='X')

# Calculate velocity graphs
scv.tl.velocity_graph(adata_comb_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_comb_temp, basis='pca',  color=["cell.types", "high_level_clusters_removed", "age_group"], save=save_path+"comb_umap.svg")
