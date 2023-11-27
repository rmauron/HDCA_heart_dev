"""
title: "13_HDCA_heart_scFates"
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

# Create copy of adata to use for Endothelial Deep Level analysis
adata_EN_sub = adata.copy()

# Append Endothelial Deep Level metadata and subset Deep Level AnnData object based on indexes present in metadata
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

# Selection of clusters from Deep Level Seurat to calculated Velocity on
#adata_IN_sub_temp = adata_IN_sub[(adata_IN_sub.obs.age_groups.isin(["12-14"]))]
adata_EN_sub = adata_EN_sub[(adata_EN_sub.obs.clusters_subset.isin([1, 12, 3, 6, 18, 17, 7, 19, 15, 13, 10, 2, 5, 21, 16, 20, 11, 14]))]

# Create copy of adata to use for INNERVATION Deep Level analysis
adata_EN = adata_EN_sub.copy()


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
adata_EN_sub.obsm['X_umap'][:, 1] = -adata_EN_sub.obsm['X_umap'][:, 1] # flip the y-coordinates to visualized the umap correctly


#############################################
# Estimate RNA velocity
#############################################

scv.tl.velocity(adata_EN_sub)
scv.tl.velocity_graph(adata_EN_sub)

scv.pl.velocity_embedding_stream(adata_EN_sub, basis='pca',  color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_pca.svg")
scv.pl.velocity_embedding_stream(adata_EN_sub, basis='umap', color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"EN_velocity_graph_umap.svg")


#############################################
# Interprete the velocities
#############################################

scv.pl.velocity(adata_EN_sub, ["SOX10", "FOXD3", "PHOX2B", "PHOX2A", "PRPH", "ASCL1", "CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"], basis='umap', ncols=2, add_margin = 0, dpi=100, show=True, save=save_path+"EN_velocity_markers_umap.pdf")
scv.pl.velocity(adata_EN_sub, ["SOX10", "FOXD3", "PHOX2B", "PHOX2A", "PRPH", "ASCL1", "CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"], basis='pca', ncols=2, add_margin = 0, dpi=100, show=True, save=save_path+"EN_velocity_markers_pca.pdf")


scv.pl.umap(adata_EN_sub, color=['age_groups'])
scv.pl.pca(adata_EN_sub, color=['age_groups'])

#############################################
# Identify important genes
#############################################

scv.tl.rank_velocity_genes(adata_EN_sub, groupby='clusters_subset', min_corr=.3)

df = scv.get_df(adata_EN_sub.uns['rank_velocity_genes']['names'])
df.head()


#############################################
# Latent time
#############################################

scv.tl.recover_dynamics(adata_EN_sub, n_jobs =8)
scv.tl.latent_time(adata_EN_sub)
scv.pl.scatter(adata_EN_sub, basis='umap', color='latent_time', color_map='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"EN_latent_time_umap.pdf")
scv.pl.scatter(adata_EN_sub, basis='pca', color='latent_time', color_map='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"EN_latent_time_pca.pdf")

#############################################
# Pseudotime
#############################################

scv.tl.velocity_pseudotime(adata_EN_sub)
scv.pl.scatter(adata_EN_sub, color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"EN_velocity_pseudotime_umap.pdf")
scv.pl.scatter(adata_IN_sub, basis='pca', color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"velocity_pseudotime_pca.pdf")


#############################################
# PAGA velocity graph
#############################################

adata_IN_sub.uns['neighbors']['distances'] = adata_IN_sub.obsp['distances']
adata_IN_sub.uns['neighbors']['connectivities'] = adata_IN_sub.obsp['connectivities']

scv.tl.paga(adata_IN_sub, groups='clusters_subset')
df = scv.get_df(adata_IN_sub, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata_IN_sub, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, show=True, save=save_path+"paga_umap.pdf")
scv.pl.paga(adata_IN_sub, basis='pca', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, show=True, save=save_path+"paga_pca.pdf")





