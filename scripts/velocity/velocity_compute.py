#!/usr/bin/env python
#  coding: utf-8
import os
import argparse
import scvelo as scv
import pandas as pd
import numpy as np
import pickle, gzip
import scipy.sparse as sp


hlp = "Compute RNA velocity on a given dimensionality reduction."

parser = argparse.ArgumentParser(description=hlp)
parser.add_argument("-e", "--experiment", dest="experiment", help="Experiment ID", required=True)
parser.add_argument("-L", "--loom", dest="loom", help="Gene expression in .loom format", required=True)
parser.add_argument("-r", "--reduction", dest="reduction", help="Dimensionality reduction to use (pca, cca, scscope)",
                    required=True)
args = parser.parse_args()

# Arguments
reduction = args.reduction
in_loom = args.loom
exp_id = args.experiment

# Paths
top_dir = os.path.join("output", exp_id)
out_dir = os.path.join(top_dir, "velocity")
out_pkl = os.path.join(out_dir, "data.pkl.gz")
out_csr = os.path.join(out_dir, "transitions.npz")

# Read reduction
in_pca = os.path.join(top_dir, "pca", "tsne_coordinates.csv")       # PCA on t-SNE
in_tsne = os.path.join(top_dir, reduction, "tsne_coordinates.csv")  # PCA on chosen dim. reduction
assert os.path.exists(in_pca)
assert os.path.exists(in_tsne)

# Make output
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print("Makedir %s" % out_dir)

# Load PCA
pca_df = pd.read_csv(in_pca)
barcode2coord = dict(zip(pca_df.iloc[:, 0], pca_df.iloc[:, 1:].values))

# Load TSNE
tsne_df = pd.read_csv(in_tsne)
barcode2tsne_coord = dict(zip(tsne_df.iloc[:,0], tsne_df.iloc[:,1:].values))

# Load Loom
print("Load loom")
xenoData = scv.read(in_loom)

cellsNotFound = [x for x in xenoData.obs_names if x not in barcode2coord.keys()]
print(cellsNotFound)
xenoDataAll = xenoData.copy()

# Remove the missing cell
xenoData = xenoData[[i for i,x in enumerate(xenoData.obs_names ) if x  in barcode2coord.keys()], :]

# Show proportions
scv.utils.show_proportions(xenoData)


# ## Velocity Computation and Projection on TSNE embedding

# Basic preprocess
print("Basic preprocess")
scv.pp.filter_genes(xenoData, min_shared_counts=10)
scv.pp.normalize_per_cell(xenoData)
scv.pp.filter_genes_dispersion(xenoData, n_top_genes=3000)
scv.pp.log1p(xenoData)
scv.pp.filter_and_normalize(xenoData, min_shared_counts=10, n_top_genes=3000)
xenoBackup = xenoData.copy()

# Adding PCA and tSNE
print("Adding own tSNE")
xenoData.X_pca = np.array([barcode2coord[x] for x in xenoData.obs_names])
xenoData.obsm['X_pca'] = xenoData.X_pca
xenoData.X_tsne = np.array([barcode2tsne_coord[x] for x in xenoData.obs_names])
xenoData.obsm['X_tsne'] = xenoData.X_tsne

# Velocity
print("Velocity")
scv.pp.moments(xenoData, n_pcs=30, n_neighbors=30)
scv.tl.velocity(xenoData)
scv.tl.velocity_graph(xenoData)
scv.tl.velocity_embedding(xenoData, basis="tsne")

# Add transitions to output
print("Compute transitions")
xenoTransitionProb = scv.tl.transition_matrix(xenoData)
print("Type:", type(xenoTransitionProb))
print("Size:", xenoTransitionProb.shape)
print("Density:", xenoTransitionProb.astype(bool).mean())
sp.save_npz(out_csr, xenoTransitionProb)
print("Written %s" % out_csr)

# Save to pickle
pickle.dump(xenoData, gzip.open(out_pkl, "w"))
print("Written %s" % out_pkl)
