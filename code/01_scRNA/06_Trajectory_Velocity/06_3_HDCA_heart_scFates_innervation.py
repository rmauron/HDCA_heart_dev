"""
title: "06_3_HDCA_heart_scFates_innervation"
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
adata_IN_sub = adata_IN.copy()


#############################################
# Preprocessing
#############################################

sc.pp.filter_genes(adata_IN_sub, min_cells=3)
sc.pp.normalize_total(adata_IN_sub)
sc.pp.log1p(adata_IN_sub, base=10)
sc.pp.highly_variable_genes(adata_IN_sub, n_top_genes=5000, flavor='cell_ranger')
adata_IN_sub.raw = adata_IN_sub
adata_IN_sub=adata_IN_sub[:,adata_IN_sub.var.highly_variable]
sc.pp.scale(adata_IN_sub)
sc.pp.pca(adata_IN_sub)
sc.pp.neighbors(adata_IN_sub, n_neighbors=25)
sc.tl.umap(adata_IN_sub)
sc.pp.pca(adata_IN_sub)


#############################################
# Learn curve using ElPiGraph algorithm
#############################################

scf.tl.curve(adata_IN_sub, Nodes=30, use_rep="X_pca", ndims_rep=2,)
adata_IN_sub.obsm['X_R']


#############################################
# Selecting a root and computing pseudotime
#############################################

scf.tl.root(adata_IN_sub,"PENK")
scf.tl.pseudotime(adata_IN_sub,n_jobs=20,n_map=100,seed=42)

# Plot pseudotime
sc.pl.pca(adata_IN_sub,color="t")

# Plot clusters
sc.pl.pca(adata_IN_sub, color="clusters_subset")

# Plot features on embedding (Fig. 4C, 4G, Supp Fig. 4C)
sc.pl.pca(adata_IN_sub,color=["SOX10", "FOXD3", "MBP", "MPZ", "ASCL1", "PRPH"],cmap="RdBu_r") #Fig 4C
sc.pl.pca(adata_IN_sub,color=["EPAS1", "COX4I2", "HIGD1C", "CHGB"],cmap="RdBu_r") #Fig 4G
sc.pl.pca(adata_IN_sub,color=["PMP22", "PHOX2B", "PHOX2A", "HTR3A", "CHGA", "PENK", "SST", "CHAT", "DBH", "TH", "PNMT", "VIP"],cmap="RdBu_r") #Supp Fig. 4C
