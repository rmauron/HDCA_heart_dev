# import libraries
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

# Setups
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
warnings.simplefilter("ignore", category=UserWarning)


# Import my data
os.getcwd()
spliced = "/Users/raphaelmauron/hdca_DevHeart/data/AnnData/data_spliced.h5ad"
spliced_data = ad.read_h5ad(spliced)

unspliced = "/Users/raphaelmauron/hdca_DevHeart/data/AnnData/data_unspliced.h5ad"
unspliced_data = ad.read_h5ad(unspliced)

# Create AnnData object
adata = spliced_data.copy()

adata.layers['spliced'] = spliced_data.X
adata.layers['unspliced'] = unspliced_data.X

# Return expression matrices
adata.X

# Create copy of adata to use for Deep Level analysis
adata_FB = adata.copy()


# Append High Level metadata
meta = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata.tsv", sep='\t', index_col=0)
meta
adata.obs = pd.concat([adata.obs, meta], axis=1)
adata
adata.obs

# Append FB Deep Level metadata
meta_FB = pd.read_table("/Users/raphaelmauron/hdca_DevHeart/data/AnnData/metadata_FB.tsv", sep='\t', index_col=0)
meta_FB
adata_FB.obs = pd.concat([adata_FB.obs, meta_FB], axis=1)
adata_FB
adata_FB.obs
adata_FB = adata_FB[adata_FB.obs['cell.types'].notna()]
adata_FB
adata_FB.obs

# Check proportions
#scv.pl.proportions(adata)

# Subset adata to get less cells


# Preprocessing of the data
scv.pp.filter_and_normalize(adata_FB, min_shared_counts=2500, n_top_genes=2000, subset_highly_variable=False) #AnnData object with n_obs × n_vars = 76991 × 8276
adata_FB

sc.tl.pca(adata_FB)
sc.pp.neighbors(adata_FB, n_pcs=30, n_neighbors=30, random_state=0)
scv.pp.moments(adata_FB, n_pcs=None, n_neighbors=None)


# Run scVelo
scv.tl.recover_dynamics(adata_FB, n_jobs=8) # this guy fails

scv.tl.velocity(adata_FB, mode="dynamical")

vk = cr.kernels.VelocityKernel(adata_FB)
vk.compute_transition_matrix()


ck = cr.kernels.ConnectivityKernel(adata_FB)
ck.compute_transition_matrix()

combined_kernel = 0.8 * vk + 0.2 * ck

print(combined_kernel)


vk.plot_projection(adata.obsm['X_pca']) # The goal

vk.plot_random_walks(start_ixs={"clusters": "Ngn3 low EP"}, max_iter=200, seed=0)