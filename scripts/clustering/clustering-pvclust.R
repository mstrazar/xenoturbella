#!/usr/bin/env Rscript
require(optparse)

# Parameters
nboot = 1000    # Number of bootstrap iterations
alpha = 0.9     # Percentage of bootstrap experiments in which cluster must be found to 
                # be deemed significant

hlp = "Compute a clustering tree with pvclust on a given dimensionality reduction."

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-r", "--reduction"), type="character", help="Data representation/embedding (hvg, pca, cca, scscope)", 
              default="pca"),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", 
              default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
config = opt$config
reduction.use = opt$reduction
stopifnot(reduction.use %in% c("hvg", "pca", "cca", "scscope"))

# Load configturation
source(config)

# Create output
top.dir = file.path(out.base.dir, exp_id)
out.dir = file.path(top.dir, "clustering", reduction.use)
data.file = file.path(top.dir, "data.Robj")
xeno.makedirs(out.dir)

# Read existing Seurat data frame
load(data.file, envir = .GlobalEnv)
data = .GlobalEnv$data

# Cluster assignment matrix
  ident = as.character(data@ident)
  clusters = levels(data@ident)
  C = matrix(0, nrow=length(ident), ncol=length(clusters))
  row.names(C) = data@cell.names
  colnames(C) = clusters
  C[cbind(data@cell.names, ident)] = 1


# Cell embeddings to cluster
if(reduction.use %in% c("scscope", "pca", "cca")){
  width = 30
  height = 10
  Z = data@dr[[reduction.use]]@cell.embeddings
  X = (t(C) %*% Z) / colSums(C)
  tree <- pvclust(t(X), method.dist="euclidean", method.hclust="complete", nboot=nboot)
} else if(reduction.use == "hvg"){
  # Cluster variable genes
  width = 30
  height = 10
  W = t(as.matrix(data@data[data@var.genes,]))
  Y = (t(C) %*% W) / colSums(C)
  tree <- pvclust(t(Y), method.dist="euclidean", method.hclust="complete", nboot=nboot)
}


# Save clustering plot
fname = file.path(out.dir, sprintf("clustering_%d_pvclust.pdf", length(levels(data@ident))))
pdf(fname, width = width, height = height)
plot(tree)
pvrect(tree, alpha=alpha)
dev.off()
message(sprintf("Written %s", fname))

# Save clustering tree to data
fname = file.path(out.dir, "tree.Robj")
save(tree, file=fname, compress = TRUE)
message(sprintf("Written %s", fname))
