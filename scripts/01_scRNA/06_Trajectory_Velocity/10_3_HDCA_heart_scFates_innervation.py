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
# Parameters optimization (plots clusters and pseudotime on multiple PCA and UMAP)
#############################################

for n_top_genes in [2000, 3000, 4000, 5000, 6000]:
    for n_neighbors in [10, 15, 20, 25, 30, 35, 40, 45, 50]:
        print(n_top_genes, n_neighbors)

        # Your existing code goes here to generate the plots
        adata_IN_sub = adata_IN.copy()
        
        sc.pp.filter_genes(adata_IN_sub, min_cells=3)
        sc.pp.normalize_total(adata_IN_sub)
        sc.pp.log1p(adata_IN_sub, base=10)
        sc.pp.highly_variable_genes(adata_IN_sub, n_top_genes=n_top_genes, flavor='cell_ranger')
        adata_IN_sub.raw = adata_IN_sub
        adata_IN_sub=adata_IN_sub[:,adata_IN_sub.var.highly_variable]
        sc.pp.scale(adata_IN_sub)
        sc.pp.pca(adata_IN_sub)
        sc.pp.neighbors(adata_IN_sub, n_neighbors=n_neighbors)
        sc.tl.umap(adata_IN_sub)
        sc.pp.pca(adata_IN_sub)
        scf.tl.curve(adata_IN_sub, Nodes=30, use_rep="X_pca", ndims_rep=2,)
        adata_IN_sub.obsm['X_R']
        scf.tl.root(adata_IN_sub,"PENK")
        scf.tl.pseudotime(adata_IN_sub,n_jobs=20,n_map=100,seed=42)

        fig, axs = plt.subplots(1, 6, figsize=(30, 5))  # Create a 1x6 grid of subplots

        # Generate and place plots in the grid
        sc.pl.pca(adata_IN_sub, color="t", ax=axs[0], show=False)
        sc.pl.umap(adata_IN_sub, color="t", ax=axs[1], show=False)
        scf.pl.trajectory(adata_IN_sub, basis="pca", arrows=True, arrow_offset=7, ax=axs[2], show=False)
        scf.pl.trajectory(adata_IN_sub, basis="umap", arrows=False, arrow_offset=0, ax=axs[3], show=False)
        sc.pl.umap(adata_IN_sub, color=['clusters_subset'], ax=axs[4], show=False)
        sc.pl.pca(adata_IN_sub, color=["clusters_subset"], ax=axs[5], show=False)

        # Set the combined figure title
        title = f"n_top_genes_{n_top_genes}_n_neighbors_{n_neighbors}"
        fig.suptitle(title)

        # Save the figure with the parameter values in the file name
        plt.tight_layout()
        file_name = f"{title}.png"
        plt.savefig(save_path + file_name)
        plt.close()  # Close the figure


#############################################
#############################################
#############################################
# n_top_genes = 5000 & n_neighbors = 25
#############################################
#############################################
#############################################

#############################################
#############################################
#############################################
# PCA
#############################################
#############################################
#############################################

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
scf.pl.graph(adata_IN_sub, basis="pca")

sc.pl.pca(sc.AnnData(adata_IN_sub.obsm["X_R"],obsm=adata_IN_sub.obsm),color="1",cmap="Reds")


#############################################
# Selecting a root and computing pseudotime
#############################################

scf.tl.root(adata_IN_sub,"PENK")
scf.tl.pseudotime(adata_IN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.pca(adata_IN_sub,color="t")

scf.pl.trajectory(adata_IN_sub, basis="pca", arrows=True, arrow_offset=7)

sc.pl.pca(adata_IN_sub, color=[ "clusters_subset"])

sc.pl.pca(adata_IN_sub,color=["SOX10", "FOXD3", "MBP", "MPZ", "PRPH", "PHOX2B", "PHOX2A", "ASCL1"],cmap="RdBu_r")
sc.pl.pca(adata_IN_sub,color=["CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"],cmap="RdBu_r")
sc.pl.pca(adata_IN_sub,color=["VIP", "SST", "EPAS1", "COX4I2", "CHAT", "DBH", "PMEL", "DCT"],cmap="RdBu_r")
sc.pl.pca(adata_IN_sub,color=["VIP", "SST", "EPAS1", "COX4I2", "NDUFS2", "DBH", "PMEL", "DCT"],cmap="RdBu_r")


#############################################
# Assign and plot milestones
#############################################

sc.pl.pca(adata_IN_sub,color="milestones")

scf.tl.rename_milestones(adata_IN_sub,new={"2":"Neurons", "0": "Radial Glia"})
adata_IN_sub.uns["graph"]["milestones"]
sc.pl.pca(adata_IN_sub,color="milestones")
scf.pl.milestones(adata_IN_sub,basis="pca",annotate=True)


#############################################
# Linearity deviation assessment
#############################################

scf.tl.linearity_deviation(adata_IN_sub,
                           start_milestone="27",
                           end_milestone="Neurons",
                           n_jobs=20,plot=True,basis="pca")
scf.pl.linearity_deviation(adata_IN_sub,
                           start_milestone="27",
                           end_milestone="Neurons")
sc.pl.pca(adata_IN_sub,color=["STMN2","ELAVL4","CD24"],cmap="RdBu_r")



scf.tl.linearity_deviation(adata_IN_sub,
                           start_milestone="27",
                           end_milestone="Radial Glia",
                           n_jobs=20,plot=True,basis="pca")
scf.pl.linearity_deviation(adata_IN_sub,
                           start_milestone="27",
                           end_milestone="Radial Glia")
sc.pl.pca(adata_IN_sub,color=["HSPB2","POSTN","PLAT"],cmap="RdBu_r")



#############################################
# Significantly changing feature along pseudotime test
#############################################

scf.tl.test_association(adata_IN_sub, n_jobs=20)

scf.tl.test_association(adata_IN_sub, reapply_filters=True, A_cut=.5)
scf.pl.test_association(adata_IN_sub)


#############################################
# Fitting & clustering significant features
#############################################

scf.tl.fit(adata_IN_sub,n_jobs=20)

marker_trend_list = ["ASCL1", "STMN2", "ELAVL4", "COL1A2", "ALDH1A1", "XKR4"]

for i in marker_trend_list:
    scf.pl.single_trend(adata_IN_sub, feature=i, basis="pca", color_exp="k")

# ?
scf.pl.single_trend(adata_IN_sub, feature=marker_trend_list, basis="pca", color_exp="k")


scf.pl.single_trend(adata_IN_sub,"ASCL1",basis="pca",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"STMN2",basis="pca",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"ELAVL4",basis="pca",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"COL1A2",basis="pca",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"ALDH1A1",basis="pca",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"XKR4",basis="pca",color_exp="k")



scf.tl.cluster(adata_IN_sub,n_neighbors=100,metric="correlation")
adata_IN_sub.var.clusters.unique()

for c in adata_IN_sub.var["clusters"].unique():
    scf.pl.trends(adata_IN_sub,features=adata_IN_sub.var_names[adata_IN_sub.var.clusters==c],basis="pca")




#############################################
#############################################
#############################################
# UMAP
#############################################
#############################################
#############################################

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

scf.tl.curve(adata_IN_sub, Nodes=30, use_rep="X_umap", ndims_rep=2,)
adata_IN_sub.obsm['X_R']
scf.pl.graph(adata_IN_sub, basis="umap")

sc.pl.umap(sc.AnnData(adata_IN_sub.obsm["X_R"],obsm=adata_IN_sub.obsm),color="1",cmap="Reds")


#############################################
# Selecting a root and computing pseudotime
#############################################

scf.tl.root(adata_IN_sub,"PENK")
scf.tl.pseudotime(adata_IN_sub,n_jobs=20,n_map=100,seed=42)

sc.pl.umap(adata_IN_sub,color="t")

scf.pl.trajectory(adata_IN_sub, basis="umap", arrows=True, arrow_offset=7) 

sc.pl.umap(adata_IN_sub, color=['clusters_subset'])

sc.pl.umap(adata_IN_sub, color=['age_groups'])
sc.pl.pca(adata_IN_sub, color=['age_groups'])

sc.pl.umap(adata_IN_sub,color=["SOX10", "FOXD3", "PHOX2B", "PHOX2A", "PRPH", "ASCL1", "CHGB", "CHGA", "PENK", "MBP", "PMP22", "HTR3A", "TH", "PNMT"],cmap="RdBu_r")


#############################################
# Assign and plot milestones
#############################################

sc.pl.umap(adata_IN_sub,color="milestones")

scf.tl.rename_milestones(adata_IN_sub,new={"2":"Radial Glia", "0": "Neurons"})
adata_IN_sub.uns["graph"]["milestones"]
sc.pl.umap(adata_IN_sub,color="milestones")
scf.pl.milestones(adata_IN_sub,basis="umap",annotate=True)


#############################################
# Linearity deviation assessment
#############################################

scf.tl.linearity_deviation(adata_IN_sub,
                           start_milestone="25",
                           end_milestone="Neurons",
                           n_jobs=20,plot=True,basis="umap")
scf.pl.linearity_deviation(adata_IN_sub,
                           start_milestone="25",
                           end_milestone="Neurons")
sc.pl.umap(adata_IN_sub,color=["STMN2","CENPV", "MLLT11"],cmap="RdBu_r")



scf.tl.linearity_deviation(adata_IN_sub,
                           start_milestone="25",
                           end_milestone="Radial Glia",
                           n_jobs=20,plot=True,basis="umap")
scf.pl.linearity_deviation(adata_IN_sub,
                           start_milestone="25",
                           end_milestone="Radial Glia")
sc.pl.umap(adata_IN_sub,color=["CRABP1","GASK1B","DHRS3"],cmap="RdBu_r")



#############################################
# Significantly changing feature along pseudotime test
#############################################

scf.tl.test_association(adata_IN_sub, n_jobs=20)

scf.tl.test_association(adata_IN_sub, reapply_filters=True, A_cut=.5)
scf.pl.test_association(adata_IN_sub)


#############################################
# Fitting & clustering significant features
#############################################

scf.tl.fit(adata_IN_sub,n_jobs=20)

marker_trend_list = ["ASCL1", "STMN2", "ELAVL4", "COL1A2", "ALDH1A1", "XKR4"]

for i in marker_trend_list:
    scf.pl.single_trend(adata_IN_sub, i, basis="umap", color_exp="k")



scf.pl.single_trend(adata_IN_sub,"ASCL1",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"STMN2",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"ELAVL4",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"COL1A2",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"ALDH1A1",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"XKR4",basis="umap",color_exp="k")
scf.pl.single_trend(adata_IN_sub,"CHAT",basis="umap",color_exp="k")



scf.tl.cluster(adata_IN_sub,n_neighbors=100,metric="correlation")
adata_IN_sub.var.clusters.unique()

for c in adata_IN_sub.var["clusters"].unique():
    scf.pl.trends(adata_IN_sub,features=adata_IN_sub.var_names[adata_IN_sub.var.clusters==c],basis="umap")




# Re-calculate neighborhood between the points from the cluster selection
sc.pp.neighbors(adata_IN_sub, random_state=0)

scv.tl.velocity_graph(adata_IN_sub)

# Plot the velocity on the UMAP imported
# scv.pl.velocity_embedding_stream(adata_HL_temp, basis='umap',  color=["high_level_clusters"], save=save_path+"HL_4_35_5_9_11_17_umap.svg")
scv.pl.velocity_embedding_stream(adata_IN_sub, basis='pca', color=["high_level_clusters"], save=save_path+"HL_age_12-14_4_35_5_9_11_17_pca.svg")