#!/usr/bin/env python3
import argparse
from os import path

import anndata
import numpy as np
import scprep
import scscope
import scipy.sparse as sp
from os.path import basename, join
from sklearn.preprocessing import OneHotEncoder

# Set up CLI
help_msg = "Impute missing values using scScope."

parser = argparse.ArgumentParser(description=help_msg)
parser.add_argument(
    "-i", "--input", dest="input", help="Input data matrix (.h5ad)", required=True
)
parser.add_argument(
    "-o", "--output-dir", dest="output_dir", help="Output directory", required=True
)
group = parser.add_argument_group(
    "scScope", "Parameters that will be passed on to the scScope imputation method"
)
group.add_argument("--n-latent", default=30, type=int)
group.add_argument("--depth", default=2, type=int)
group.add_argument("--use-mask", default=True, type=bool)
group.add_argument("--batch-size", default=256, type=int)
group.add_argument("--batch-var", default=None, type=str)
group.add_argument("--max-epoch", default=100, type=int)
group.add_argument("--depth-encoder", default=[64], type=int, nargs="+")
group.add_argument("--depth-decoder", default=[128], type=int, nargs="+")
group.add_argument("--rate", default=0.0001, type=float)
group.add_argument("--beta", default=0.05, type=float)
args = parser.parse_args()

# Make sure that we've actually gotten an H5AD file
fname, ext = path.splitext(basename(args.input))
ext = ext[1:]
assert ext == "h5ad", "The input data must be an H5AD formatted file."

# Varying and hardcoded parameters
scscope_args = dict(
          latent_code_dim=args.n_latent,
          use_mask=args.use_mask,
          batch_size=args.batch_size,
          max_epoch=args.max_epoch,
          T=args.depth,
          encoder_layers=args.depth_encoder,
          decoder_layers=args.depth_decoder,
          learning_rate=args.rate,
          beta1=args.beta,
          epoch_per_check=5,
          num_gpus=1,
)

# Read input
f_in = args.input
adata = anndata.read_h5ad(f_in)

# Refer to these two tutorial notebooks on the official MAGIC GH repo
X = adata.X

# Store scaling parameters
rescale, library_size = scprep.normalize._get_scaled_libsize(X, return_library_size=True)

# Scale
X = scprep.normalize.library_size_normalize(X)
X = scprep.transform.log(X, 1)

# Batch indicator, if any. Passed as one-hot encoding.
B = []
bv = args.batch_var
if bv is not None:
    assert bv in adata.obs
    y = adata.obs[bv].as_matrix().ravel().astype(str)
    B = OneHotEncoder(handle_unknown='ignore').fit_transform(y.reshape((len(y), 1))).toarray()
scscope_args["exp_batch_idx_input"] = B

# Convert data to dense
X = X.toarray() if sp.isspmatrix(X) else X
d = (X > 0).mean()
print("Original density: %f" % d)

# Summarize arguments
print("Running scScope with: ")
for ky, val in scscope_args.items():
    print("\t%s: %s" % (ky, str(val)[:25]))

# Run model
model = scscope.train(X, **scscope_args)

# Predict
Z, I, _ = scscope.predict(X, model, batch_effect=B)
print("New density: %f" % (I > 0).mean())

# Do not store batch indicator
del scscope_args["exp_batch_idx_input"]

# Un-normalize the imputed matrix
X_imputed = adata.copy()
X_imputed.X = (np.exp(I) - 1) * library_size[:, None] / rescale
X_imputed.uns = adata.uns
X_imputed.uns["args"] = scscope_args
X_imputed.uns["latent"] = Z

# Write output
fname = join(args.output_dir, f"scScope__{fname}.{ext}")
X_imputed.write_h5ad(fname)
print("Written %s" % fname)
