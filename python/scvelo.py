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
spliced = "/Users/raphaelmauron/hdca_DevHeart/data/AnnData/data_spliced.h5ad"
spliced_data = ad.read_h5ad(spliced)

unspliced = "/Users/raphaelmauron/hdca_DevHeart/data/AnnData/data_unspliced.h5ad"
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
meta = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata.tsv", sep='\t', index_col=0)
meta
adata_HL = adata_HL[meta.index]
adata_HL.obs = pd.concat([adata_HL.obs, meta], axis=1)
adata_HL
adata_HL.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_HL = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP.csv", delimiter=',',index_col=0)
UMAP_HL = UMAP_HL.to_numpy()
adata_HL.obsm["X_umap"]=UMAP_HL

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_HL, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_HL)
scv.pp.filter_genes_dispersion(adata_HL, n_top_genes=2000)
scv.pp.log1p(adata_HL)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_HL_temp = adata_HL[(adata_HL.obs.high_level_clusters.isin(["4_35", "5", "9", "11", "17"]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_HL_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_HL_temp)

# Plot the velocity on the UMAP imported
scv.pl.velocity_embedding_stream(adata_HL_temp, basis='umap',  color=["high_level_clusters"], save=save_path+"HL_4_35_5_9_11_17_umap.svg")




######################################################
# Fibroblast
######################################################

# Create copy of adata to use for FB Deep Level analysis
adata_FB = adata.copy()

# Append FB Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_FB = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_FB.tsv", sep='\t', index_col=0)
meta_FB
adata_FB = adata_FB[meta_FB.index]
adata_FB.obs = pd.concat([adata_FB.obs, meta_FB], axis=1)
adata_FB
adata_FB.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_FB = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_FB.csv", delimiter=',',index_col=0)
UMAP_FB = UMAP_FB.to_numpy()
adata_FB.obsm["X_umap"]=UMAP_FB

# Preprocessing steps on the Deep Level AnnData
scv.pp.filter_genes(adata_FB, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_FB)
scv.pp.filter_genes_dispersion(adata_FB, n_top_genes=2000)
scv.pp.log1p(adata_FB)

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_FB_temp = adata_FB[(adata_FB.obs.clusters_subset.isin([13,11]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_FB_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_FB_temp)

# Plot the velocity on the UMAP imported
scv.pl.velocity_embedding_stream(adata_FB_temp, basis='umap',  color=["cell.types", "clusters_subset"])



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
adata_EN_temp = adata_EN[(adata_EN.obs.clusters_subset.isin([15,13,10,2]))]

# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_EN_temp, random_state=0)

# Calculate velocity graphs
scv.tl.velocity_graph(adata_EN_temp)

# Plot the velocity on the UMAP imported
#plt.savefig("myImagePDF.pdf", format="pdf", bbox_inches="tight")
#plt.show()
scv.pl.velocity_embedding_stream(adata_EN_temp, basis='umap',  color=["cell.types", "clusters_subset"])



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
