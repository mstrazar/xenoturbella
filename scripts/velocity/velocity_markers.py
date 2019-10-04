#!/usr/bin/env python
#  coding: utf-8
import os
import argparse
import scvelo as scv
import pandas as pd
import pickle, gzip
import matplotlib.pyplot as plt

hlp = "Draw plots with velocities and marker genes."

parser = argparse.ArgumentParser(description=hlp)
parser.add_argument("-e", "--experiment", dest="experiment", help="Experiment ID", required=True)
args = parser.parse_args()

# Arguments
exp_id = args.experiment

# Paths
top_dir = os.path.join("output", exp_id)
in_dir = os.path.join(top_dir, "velocity")
in_meta =os.path.join(top_dir, "metadata.csv")
in_pkl = os.path.join(in_dir, "data.pkl.gz")
out_dir = os.path.join(in_dir, "markers")

# Make output
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print("Makedir %s" % out_dir)

# Load input
data = pickle.load(gzip.open(in_pkl))
print("Loaded file %s" % in_pkl)
print(data)

# Load markers and assign
meta = pd.read_csv(in_meta, index_col=0)

# Add metadata by concatenation
obc = pd.concat([data.obs, meta.loc[data.obs.index.values]], axis=1)
data.obs = obc

# Iterate through columns
for col in meta.columns:
    fname = os.path.join(out_dir, "velocity_%s.tiff" % col)
    scv.pl.velocity_embedding_stream(data, basis="tsne", dpi=300, color=[col], show=False,
                                     figsize=(3, 2.9), alpha=0.9, color_map="cool")
    plt.savefig(fname)
    print("Written %s" % fname)
    plt.close()

    fname = os.path.join(out_dir, "velocity_%s_no_legend.tiff" % col)
    scv.pl.velocity_embedding_stream(data, basis="tsne", dpi=300, color=[col], show=False, title="",
                                     figsize=(3, 2.9), alpha=0.9, color_map="cool", legend_loc="none")
    plt.savefig(fname)
    print("Written %s" % fname)
    plt.close()