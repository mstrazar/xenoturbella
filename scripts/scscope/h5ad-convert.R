# Read pipeline arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  message(sprintf("Wrong input"))
  message(sprintf("Usage: h5ad-convert.R dset [nmax]"))
  q(1, save = "no")
}
dset = args[1]
nmax = as.numeric(args[2])

# Set home dir on server
home = "/ufrc/moroz/mstrazar"
if(file.exists(home)){
  Sys.setenv(HOME=home)  
}

# Start
require(Seurat)
hlp = "Convert data set to H5AD."

# Load auxilliary data
setwd("~/Dev/xeno/scripts/")
source("utils.R")
source("h5ad.R")

# Paths
design = xeno.load.design(dset=dset)

# Load and normalize data
data = xeno.load.data(design, nmax=nmax, preprocess=FALSE)

# Write H5Ad matrix
inxs = which(rowSums(data@raw.data) > 0)
X = as.matrix(t(data@raw.data[inxs,]))   # If conversion fails.
# X = t(data@raw.data[inxs,])
cell_data = data@meta.data
gene_data = data.frame(geneID=colnames(X), geneID2=colnames(X))
fname = file.path(in.h5ad.dir, sprintf("%s.h5ad", dset))
write.h5ad(fname, X, cell_data = cell_data, gene_data = NULL)
message(sprintf("Written %d x %d matrix to %s", nrow(X), ncol(X), fname))
fname = file.path(in.h5ad.dir, sprintf("%s.genes.csv", dset))
write.csv(gene_data, fname, row.names = FALSE)
message(sprintf("Written %d rows to %s", nrow(gene_data),  fname))
