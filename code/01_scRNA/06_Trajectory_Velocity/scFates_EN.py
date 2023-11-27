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
import scFates as scf
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
save_path = wd + "/hdca_DevHeart/output/HDCA_heart_sc_analysis_docker/trajectory/scFates/"


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

# Selection of clusters from Deep Level Seurat to calculated Velocity on
#adata_IN_sub_temp = adata_IN_sub[(adata_IN_sub.obs.age_groups.isin(["12-14"]))]
adata_EN_temp = adata_EN[(adata_EN.obs.clusters_subset.isin([1, 12, 3, 6, 18, 17, 7, 19, 15, 13, 10, 2, 5, 21, 16, 20, 11, 14]))]



#############################################
#############################################
#############################################
# PCA
#############################################
#############################################
#############################################

adata_EN_sub = adata_EN_temp.copy()

#############################################
# Preprocessing
#############################################

sc.pp.filter_genes(adata_EN_sub, min_cells=3)
sc.pp.normalize_total(adata_EN_sub)
sc.pp.log1p(adata_EN_sub, base=10)
sc.pp.highly_variable_genes(adata_EN_sub, n_top_genes=5000, flavor='cell_ranger')
adata_EN_sub.raw = adata_EN_sub
adata_EN_sub=adata_EN_sub[:,adata_EN_sub.var.highly_variable]
sc.pp.scale(adata_EN_sub)
sc.pp.pca(adata_EN_sub)
sc.pp.neighbors(adata_EN_sub, n_neighbors=25)
#sc.tl.umap(adata_EN_sub)
sc.pp.pca(adata_EN_sub)


#############################################
# Learn curve using ElPiGraph algorithm
#############################################

scf.tl.curve(adata_EN_sub, Nodes=30, use_rep="X_pca", ndims_rep=2,)
adata_EN_sub.obsm['X_R']
scf.pl.graph(adata_EN_sub, basis="pca")

sc.pl.pca(sc.AnnData(adata_EN_sub.obsm["X_R"],obsm=adata_EN_sub.obsm),color="1",cmap="Reds")


#############################################
# Selecting a root and computing pseudotime
#############################################

scf.tl.root(adata_EN_sub,"NPR3")
scf.tl.pseudotime(adata_EN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.pca(adata_EN_sub,color="t")

scf.pl.trajectory(adata_EN_sub, basis="pca", arrows=True, arrow_offset=7)

sc.pl.pca(adata_EN_sub, color=[ "clusters_subset"])






#############################################
#############################################
#############################################
# UMAP
#############################################
#############################################
#############################################

adata_EN_sub = adata_EN_temp.copy()

#############################################
# Preprocessing
#############################################

sc.pp.filter_genes(adata_EN_sub, min_cells=3)
sc.pp.normalize_total(adata_EN_sub)
sc.pp.log1p(adata_EN_sub, base=10)
sc.pp.highly_variable_genes(adata_EN_sub, n_top_genes=5000, flavor='cell_ranger')
adata_EN_sub.raw = adata_EN_sub
adata_EN_sub=adata_EN_sub[:,adata_EN_sub.var.highly_variable]
sc.pp.scale(adata_EN_sub)
sc.pp.pca(adata_EN_sub)
sc.pp.neighbors(adata_EN_sub, n_neighbors=25)
#sc.tl.umap(adata_EN_sub)
sc.pp.pca(adata_EN_sub)


#############################################
# Learn curve using ElPiGraph algorithm
#############################################

scf.tl.curve(adata_EN_sub, Nodes=30, use_rep="X_umap", ndims_rep=2,)
adata_EN_sub.obsm['X_R']
scf.pl.graph(adata_EN_sub, basis="umap")

sc.pl.umap(sc.AnnData(adata_EN_sub.obsm["X_R"],obsm=adata_EN_sub.obsm),color="1",cmap="Reds")


#############################################
# Selecting a root and computing pseudotime
#############################################

scf.tl.root(adata_EN_sub,"NPR3")
scf.tl.pseudotime(adata_EN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.umap(adata_EN_sub,color="t") #this guy is nice

scf.pl.trajectory(adata_EN_sub, basis="umap", arrows=True, arrow_offset=7)

sc.pl.umap(adata_EN_sub, color=['clusters_subset'])


