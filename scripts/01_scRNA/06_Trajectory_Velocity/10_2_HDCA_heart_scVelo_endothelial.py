"""
title: "10_HDCA_heart_scVelo_IN"
author: "Raphaël Mauron"
"""

"""
In this script, RNA velocity is performed on Innervation clusters
"""

# Import libraries
import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import anndata as ad
import os, sys
from scipy.sparse import csr_matrix
print(ad.__version__)
import h5py
import warnings
import pandas as pd
from matplotlib import pyplot as plt


# Setups
sc.set_figure_params()
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
warnings.simplefilter("ignore", category=UserWarning)
wd = os.getcwd()
save_path = wd + "/hdca_DevHeart/output/HDCA_heart_sc_analysis_docker/trajectory/scVelo/"


#############################################
# Build AnnData object from exported Seurat objects, metadata & embeddings
#############################################

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

# Create copy of adata to use for EN Deep Level analysis
adata_EN_sub = adata.copy()

# Append EN Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
meta_EN_sub = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_EN.tsv", sep='\t', index_col=0)
meta_EN_sub
adata_EN_sub = adata_EN_sub[meta_EN_sub.index]
adata_EN_sub.obs = pd.concat([adata_EN_sub.obs, meta_EN_sub], axis=1)
adata_EN_sub
adata_EN_sub.obs

# Import and append the UMAP from Deep Level Seurat to the Deep Level AnnData 
UMAP_EN = pd.read_csv("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/UMAP_EN.csv", delimiter=',',index_col=0)
UMAP_EN = UMAP_EN.to_numpy()
adata_EN_sub.obsm["X_umap"]=UMAP_EN

# Create copy of adata to use for INNERVATION Deep Level analysis
adata_EN = adata_EN_sub.copy()

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_EN = adata_EN[(adata_EN.obs.clusters_subset.isin([1, 12, 3, 6, 18, 17, 7, 19, 15, 13, 10, 2, 5, 21, 16, 20, 11, 14]))]


#############################################
# Preprocessing
#############################################

adata_EN_sub = adata_EN.copy()

scv.pp.filter_genes(adata_EN_sub, min_cells=3)
scv.pp.normalize_per_cell(adata_EN_sub)
scv.pp.filter_genes_dispersion(adata_EN_sub, n_top_genes=5000)
scv.pp.log1p(adata_EN_sub)
adata_EN_sub.raw = adata_EN_sub
sc.pp.scale(adata_EN_sub)
sc.pp.pca(adata_EN_sub)
sc.pp.neighbors(adata_EN_sub, n_neighbors=25)
#sc.tl.umap(adata_EN_sub)
sc.pp.pca(adata_EN_sub)

sc.pp.neighbors(adata_EN_sub, random_state=0)
#adata_EN_sub.obsm['X_umap'][:, 1] = -adata_EN_sub.obsm['X_umap'][:, 1] # flip the y-coordinates to visualized the umap correctly


#############################################
# Estimate RNA velocity
#############################################

scv.tl.velocity(adata_EN_sub)
scv.tl.velocity_graph(adata_EN_sub)

scv.pl.velocity_embedding_stream(adata_EN_sub, basis='pca',  color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_pca.svg")
scv.pl.velocity_embedding_stream(adata_EN_sub, basis='umap', color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_umap.svg")


#############################################
# Pseudotime
#############################################

scv.tl.velocity_pseudotime(adata_EN_sub)
scv.pl.scatter(adata_EN_sub, color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=False, save=save_path+"EN_velocity_pseudotime_umap.svg")
scv.pl.scatter(adata_EN_sub, basis='pca', color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=False, save=save_path+"EN_velocity_pseudotime_pca.svg")

scv.pl.velocity_embedding_stream(adata_EN_sub, basis='pca',  color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_velocity_pseudotime_pca.svg")
scv.pl.velocity_embedding_stream(adata_EN_sub, basis='umap', color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_velocity_pseudotime_umap.svg")





#############################################
# cl 6 and 18
#############################################

# Create copy of adata to use for INNERVATION Deep Level analysis
adata_EN = adata_EN_sub.copy()

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_EN = adata_EN[(adata_EN.obs.clusters_subset.isin([6, 18]))]


#############################################
# Preprocessing
#############################################

adata_EN_sub = adata_EN.copy()

scv.pp.filter_genes(adata_EN_sub, min_cells=3)
scv.pp.normalize_per_cell(adata_EN_sub)
scv.pp.filter_genes_dispersion(adata_EN_sub, n_top_genes=5000)
scv.pp.log1p(adata_EN_sub)
adata_EN_sub.raw = adata_EN_sub
sc.pp.scale(adata_EN_sub)
sc.pp.pca(adata_EN_sub)
sc.pp.neighbors(adata_EN_sub, n_neighbors=25)
#sc.tl.umap(adata_EN_sub)
sc.pp.pca(adata_EN_sub)

sc.pp.neighbors(adata_EN_sub, random_state=0)
#adata_EN_sub.obsm['X_umap'][:, 1] = -adata_EN_sub.obsm['X_umap'][:, 1] # flip the y-coordinates to visualized the umap correctly


#############################################
# Estimate RNA velocity
#############################################

scv.tl.velocity(adata_EN_sub)
scv.tl.velocity_graph(adata_EN_sub)

scv.pl.velocity_embedding_stream(adata_EN_sub, basis='pca',  color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_6_18_velocity_graph_pca.svg")
scv.pl.velocity_embedding_stream(adata_EN_sub, basis='umap', color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_6_18_velocity_graph_umap.svg")


#############################################
# Pseudotime
#############################################

scv.tl.velocity_pseudotime(adata_EN_sub)
scv.pl.scatter(adata_EN_sub, color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=False, save=save_path+"EN_6_18_velocity_pseudotime_umap.svg")
scv.pl.scatter(adata_EN_sub, basis='pca', color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=False, save=save_path+"EN_6_18_velocity_pseudotime_pca.svg")

scv.pl.velocity_embedding_stream(adata_EN_sub, basis='pca',  color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_6_18_velocity_graph_velocity_pseudotime_pca.svg")
scv.pl.velocity_embedding_stream(adata_EN_sub, basis='umap', color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_6_18_velocity_graph_velocity_pseudotime_umap.svg")
