#!/usr/bin/env Rscript
require(optparse)

# Parameters
set.seed(42)    # Fix random seed
perplexity = 30 # t-SNE perplexity
nboot = 1000    # Number of bootstrap iterations
alpha = 0.9     # Percentage of bootstrap experiments in which cluster must be found to 
                # be deemed significant

hlp = "Compute hierarchical clustering between neurons with pvclust on a given dimensionality reduction.
       Add a corresponding t-SNE plot colored by library of origin, neurotransmitter and inferred neuron cluster."

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

# Assign points to clusters based on significance given a clustering tree
cluster.assign <- function(tree, alpha){
  
  # Structure is addressed by clustering step;
  # There are 3 possibilities at each step:
  # 1. Two negatives form a new cluster at step i
  # 2. A negative is added to a positive: cluster already exists
  # 3. Two positives merge two existing branches, add them to step i and delete from previous locations
  struct = list()  # Significant hits
  total = list()   # All edges; something might be skipped on the wa
  edges = tree$edges
  merge = tree$hclust$merge
  n = nrow(edges)
  for(i in 1:(n-1)){
    a = min(merge[i,])
    b = max(merge[i,])
    if(a > 0) {A = total[[a]]} else {A = abs(a)} 
    if(b > 0) {B = total[[b]]} else {B = abs(b)} 
    total[[i]] = c(A, B)
    if (edges$au[i] > alpha){
      struct[[i]] = c(A, B)
      if(i == 1) next;
      for(j in 1:(i-1)){
        if(length(intersect(struct[[i]], struct[[j]]))){
          struct[[j]] = 0
        }
      }
    }
  }

  # Pruned trees
  j = 1
  prune = list()
  for(i in 1:length(struct)){
    s = struct[[i]]
    s = s[s>0]
    if(length(s) > 0){
      prune[[j]] = s
      j = j + 1
    } 
  }
  
  # Provide cluster mapping
  result = rep(NA, n)
  for(j in 1:length(prune)){
    result[prune[[j]]] = j
  }
  result
}


# Create output
top.dir = file.path(out.base.dir, exp_id)
out.dir = file.path(top.dir, "clustering", sprintf("%s-neurons", reduction.use))
data.file = file.path(top.dir, "data.Robj")
xeno.makedirs(out.dir)

# Read existing Seurat data frame
load(data.file, envir = .GlobalEnv)
data = .GlobalEnv$data

# Load neurotransmitters for colors
markers = read.csv(in.neuro.markers, stringsAsFactors = FALSE)
data = xeno.set.neurons.specific(data, markers)  

# Set colors
markers = markers[!duplicated(markers$Cell.type),]
row.names(markers) = markers$Cell.type
col.neuron = markers[as.character(data@ident), "Color"]

# Compute clustering tree
if(reduction.use == "hvg"){
  X = t(as.matrix(data@data[data@var.genes,]))
} else {
  X = data@dr[[reduction.use]]@cell.embeddings  
}
tree <- pvclust(t(X), method.dist="euclidean", method.hclust="complete", nboot=nboot)
tree$hclust$labels = as.character(data@ident)
tree$hclust$cells = as.character(data@cell.names)

# Load design for colors
design = xeno.load.design()
col.sample = design[data@meta.data$orig.ident[tree$hclust$order], "color"]

# Assign clusters to points
r = cluster.assign(tree, alpha)
r[tree$hclust$order]
tree$hclust$labels = sprintf("%2d: %s", r, tree$hclust$labels)

# Assign colors to clusters
n <- length(unique(r))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col.cluster = col_vector[r]

# Plot tree
width = 50
height = 20
fname = file.path(out.dir, sprintf("clustering_%d_pvclust.pdf", length(levels(data@ident))))
pdf(fname, width = width, height = height)
plot(tree)
pvrect(tree, alpha=alpha)
points(c(1:233), rep(-1.75, 233), col=col.sample, pch=20, cex=1.2)
dev.off()

# Save clustering tree to data
fname = file.path(out.dir, "tree.Robj")
save(tree, file=fname, compress = TRUE)
message(sprintf("Written %s", fname))

# Compute t-SNE on given dimensionality reduction
Z = Rtsne::Rtsne(X, perplexity=perplexity)$Y

# t-SNE by sample
fname = file.path(out.dir, "combined-tsne-sample.pdf")
pdf(fname, width = 6, height = 6)
plot(Z[,1], Z[, 2], pch=20, col=col.sample, 
     xlab="", cex=1.3,
     ylab="", axes = FALSE,
     main = sprintf("Neurons by library (n=%d)", nrow(Z)))
legend(9, 14, fill = design$color, legend = design$name, cex=0.5, box.lwd = 0)
dev.off()
message(sprintf("Written %s", fname))

# t-SNE by neurotransmitter
fname = file.path(out.dir, "combined-tsne-neurotransmitter.pdf")
pdf(fname, width = 6, height = 6)
plot(Z[,1], Z[, 2], pch=20, col=col.neuron, 
     xlab="", cex=1.3,
     ylab="", axes = FALSE,
     main = sprintf("Neurotransmitters"))
legend(9, 14, fill = markers$Color, legend = markers$Cell.type, cex=0.5, box.lwd = 0)
dev.off()
message(sprintf("Written %s", fname))

# t-SNE by cluster
fname = file.path(out.dir, "combined-tsne-cluster.pdf")
pdf(fname, width = 6, height = 6)
plot(Z[,1], Z[, 2], pch=20, col=col.cluster, 
     xlab="", cex=1.3,
     ylab="", axes = FALSE,
     main = sprintf("Clusters"))
for(ri in unique(r)){
  zm = colMeans(Z[r == ri,], na.rm = TRUE)
  text(zm[1], zm[2], ri)
}
dev.off()
message(sprintf("Written %s", fname))


