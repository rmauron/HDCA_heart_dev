"""
title: "10_5_HDCA_heart_scVelo_IN"
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
adata_IN_sub.obsm['X_umap'][:, 1] = -adata_IN_sub.obsm['X_umap'][:, 1] # flip the y-coordinates to visualized the umap correctly


#############################################
# Estimate RNA velocity
#############################################

scv.tl.velocity(adata_IN_sub)
scv.tl.velocity_graph(adata_IN_sub)

scv.pl.velocity_embedding_stream(adata_IN_sub, basis='pca',  color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"IN_velocity_graph_pca.svg")
scv.pl.velocity_embedding_stream(adata_IN_sub, basis='umap', color=["clusters_subset"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"IN_velocity_graph_umap.svg")

#############################################
# Interprete the velocities
#############################################

# scv.pl.velocity(adata_IN_sub, ["SOX10", "FOXD3", "PHOX2B", "PHOX2A", "PRPH", "ASCL1", "CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"], basis='umap', ncols=2, add_margin = 0, dpi=100, show=True, save=save_path+"velocity_markers_umap.pdf")
# scv.pl.velocity(adata_IN_sub, ["SOX10", "FOXD3", "PHOX2B", "PHOX2A", "PRPH", "ASCL1", "CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"], basis='pca', ncols=2, add_margin = 0, dpi=100, show=True, save=save_path+"velocity_markers_pca.pdf")


# scv.pl.umap(adata_IN_sub, color=['age_groups'])
# scv.pl.pca(adata_IN_sub, color=['age_groups'])

#############################################
# Identify important genes
#############################################

# scv.tl.rank_velocity_genes(adata_IN_sub, groupby='clusters_subset', min_corr=.3)

# df = scv.get_df(adata_IN_sub.uns['rank_velocity_genes']['names'])
# df.head()


#############################################
# Latent time
#############################################

# scv.tl.recover_dynamics(adata_IN_sub)
# scv.tl.latent_time(adata_IN_sub)
# scv.pl.scatter(adata_IN_sub, basis='umap', color='latent_time', color_map='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"latent_time_umap.pdf")
# scv.pl.scatter(adata_IN_sub, basis='pca', color='latent_time', color_map='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"latent_time_pca.pdf")

#############################################
# Pseudotime
#############################################

scv.tl.velocity_pseudotime(adata_IN_sub)
scv.pl.scatter(adata_IN_sub, color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=False, save=save_path+"IN_velocity_pseudotime_umap.pdf")
scv.pl.scatter(adata_IN_sub, basis='pca', color='velocity_pseudotime', cmap='gnuplot', size=80, add_margin = 0.1, dpi=300, show=True, save=save_path+"IN_velocity_pseudotime_pca.pdf")

scv.pl.velocity_embedding_stream(adata_IN_sub, basis='pca',  color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=True, save=save_path+"IN_velocity_graph_velocity_pseudotime_pca.pdf")
scv.pl.velocity_embedding_stream(adata_IN_sub, basis='umap', color=["velocity_pseudotime"], size=300, alpha=0.5, add_margin = 0.05, dpi=300, show=False, save=save_path+"IN_velocity_graph_velocity_pseudotime_umap.pdf")


#############################################
# PAGA velocity graph
#############################################

# adata_IN_sub.uns['neighbors']['distances'] = adata_IN_sub.obsp['distances']
# adata_IN_sub.uns['neighbors']['connectivities'] = adata_IN_sub.obsp['connectivities']

# scv.tl.paga(adata_IN_sub, groups='clusters_subset')
# df = scv.get_df(adata_IN_sub, 'paga/transitions_confidence', precision=2).T
# df.style.background_gradient(cmap='Blues').format('{:.2g}')

# scv.pl.paga(adata_IN_sub, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, show=True, save=save_path+"paga_umap.pdf")
# scv.pl.paga(adata_IN_sub, basis='pca', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, dpi=300, show=True, save=save_path+"paga_pca.pdf")
