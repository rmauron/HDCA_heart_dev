"""
title: "10_1_HDCA_heart_scFates"
author: "Raphaël Mauron"
"""

"""
In this script, RNA velocity is performed on innvervation cells.
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



#############################################
#############################################
#############################################
# PCA
#############################################
#############################################
#############################################

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_EN_temp = adata_EN[(adata_EN.obs.clusters_subset.isin([1, 12, 3, 6, 18, 17, 7, 19, 15, 13, 10, 2, 5, 21, 16, 20, 11, 14]))]


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

scf.tl.root(adata_EN_sub,"ENPP1")
scf.tl.pseudotime(adata_EN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.pca(adata_EN_sub,color="t")



#############################################
#############################################
#############################################
# UMAP
#############################################
#############################################
#############################################

# Selection of clusters from Deep Level Seurat to calculated Velocity on
adata_EN_temp = adata_EN[(adata_EN.obs.clusters_subset.isin([1, 12, 3, 6, 18, 17, 7, 19, 15, 13, 10, 2, 5, 21, 16, 20, 11, 14]))]


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

scf.tl.root(adata_EN_sub,"ENPP1")
scf.tl.pseudotime(adata_EN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.umap(adata_EN_sub,color="t")
