"""
title: "06_4_HDCA_heart_scVelo_innervation"
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

# Create copy of adata to use for INNERVATION Deep Level analysis
adata_IN_sub = adata.copy()

# Append INNERVATION Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
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

# Create copy of adata to use for INNERVATION Deep Level analysis
adata_IN = adata_IN_sub.copy()


#############################################
# Preprocessing
#############################################

adata_IN_sub = adata_IN.copy()

scv.pp.filter_genes(adata_IN_sub, min_cells=3)
scv.pp.normalize_per_cell(adata_IN_sub)
scv.pp.filter_genes_dispersion(adata_IN_sub, n_top_genes=5000)
scv.pp.log1p(adata_IN_sub)
adata_IN_sub.raw = adata_IN_sub
sc.pp.scale(adata_IN_sub)
sc.pp.pca(adata_IN_sub)
sc.pp.neighbors(adata_IN_sub, n_neighbors=25)
sc.tl.umap(adata_IN_sub)
sc.pp.pca(adata_IN_sub)

sc.pp.neighbors(adata_IN_sub, random_state=0)
#adata_IN_sub.obsm['X_umap'][:, 1] = -adata_IN_sub.obsm['X_umap'][:, 1] # flip the y-coordinates to visualized the umap correctly


#############################################
# Estimate RNA velocity (Fig. 4B)
#############################################

scv.tl.velocity(adata_IN_sub)
scv.tl.velocity_graph(adata_IN_sub)

scv.pl.velocity_embedding_stream(adata_IN_sub, basis='pca',  color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"IN_velocity_graph_pca.svg")


#############################################
# Pseudotime (Fig. 4B)
#############################################

scv.tl.velocity_pseudotime(adata_IN_sub)
scv.pl.scatter(adata_IN_sub, basis='pca', color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"IN_velocity_pseudotime_pca.pdf")

scv.pl.velocity_embedding_stream(adata_IN_sub, basis='pca',  color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=True, save=save_path+"IN_velocity_graph_velocity_pseudotime_pca.pdf")
