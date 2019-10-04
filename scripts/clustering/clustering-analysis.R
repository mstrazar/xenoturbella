#!/usr/bin/env Rscript
require(optparse)

hlp = "Select the number of clusters using two experiments"

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name (main)."),
  make_option(c("-E", "--secondary"), type="character", help="Experiment name (secondary)."),
  make_option(c("-r", "--reduction"), type="character", help="Dimensionality reduction method (pca, cca, scscope)", default="pca"),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)",
              default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id_1 = opt$experiment
exp_id_2 = opt$secondary
redc = opt$reduction
config = opt$config
stopifnot(redc %in% c("pca", "cca", "scscope"))

# Load configturation
source(config)

# Paths
out.dir = file.path("~/Dev/data/xeno/output/")
top.dir1 = file.path(out.base.dir, exp_id_1)
top.dir2 = file.path(out.base.dir, exp_id_2)
in.assig1 = file.path(top.dir1, "clustering",  "cassig.csv")
in.assig2 = file.path(top.dir2, "clustering",  "cassig.csv")

# Create output
res.dir = file.path(top.dir1, "clustering", "analysis")
xeno.makedirs(res.dir)
U = read.csv(in.assig1, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
V = read.csv(in.assig2, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Count clusters and sort by size
count.clusters <- function(x){
  z = rep(0, length(unique(x)))
  for(xi in (x + 1)){
    z[xi] = z[xi] + 1
  }
  z = z / sum(z)
  z
}

# Compare clustering distributinos
compare.dist <- function(u, v){
  sum(abs(u - v))
}

# Number of clusters
Nu = apply(U, 2, function(x) length(unique(x)))
Nv = apply(V, 2, function(x) length(unique(x)))

# Add data frame
data = data.frame()

# Compare matches
for(i in 1:length(Nu)){
  for(j in 1:length(Nv)){
    if(Nu[i] != Nv[j]) next;
    message(sprintf("Match %d/%d (n=%d)", i, j, Nu[i])) 
    u = count.clusters(U[,i])
    v = count.clusters(V[,j])
    duv = compare.dist(u, v)
    df = data.frame(ui=i, vi=j, N=Nu[i], d=duv)
    data = rbind(data, df)
    
    # Plot
    fname = file.path(res.dir, sprintf("d-%d-%d-%d.pdf", Nu[i], i, j))
    pdf(fname, width = 5, height = 3.5)
    plot(u, xlab = "Cluster", ylab = "Fraction of cells", type="l", col="blue",
         main=sprintf("N=%d    d(p, q) = %.3f", Nu[i], duv))
    lines(v, col="orange")
    grid()
    dev.off()
    message(sprintf("Written %s", fname))
    }
}


# Choose best configuration
result = data[data$d == min(data$d),]


# Distance summary
fname = file.path(res.dir, "_summary.pdf")
pdf(fname, width = 5, height = 3.5)
plot(data$N, data$d, pch=20, ylab="d(p, q)", xlab="Clusters")
lines(c(0, data$N[j]), c(0, data$d[j]), col="black")
grid()
dev.off()
message(sprintf("Written %s", fname))

# Dataset 1
clustering_id = result$ui
load(file.path(top.dir1, "clustering", redc, clustering_id, "tree.Robj"), envir = .GlobalEnv)
tree = .GlobalEnv$tree
load(file.path(top.dir1, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data
data@ident = as.factor(U[,clustering_id])
data@meta.data$ident = data@ident
data@cluster.tree = tree

# Clustering assignment 1 (if applicable)
data@meta.data$cluster.assign = NA
in_assig = file.path(top.dir1, "clustering", redc, clustering_id, "assignment.csv")
if(file.exists(in_assig)){
  af = read.csv(in_assig, header=TRUE, row.names = 1, stringsAsFactors = FALSE)
  data@meta.data$cluster.assign = af[as.character(data@ident), "Cell.type"]
}

# Write dataset 1
fname = file.path(top.dir1, "data.Robj")
save(data, file=fname, compress = TRUE)
message(sprintf("Written %d clusters in %s", length(unique(data@ident)), fname))

# Dataset 2
clustering_id = result$vi
load(file.path(top.dir2, "clustering", redc, clustering_id, "tree.Robj"), envir = .GlobalEnv)
tree = .GlobalEnv$tree
load(file.path(top.dir2, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data
data@ident = as.factor(V[,clustering_id])
data@meta.data$ident = data@ident
data@cluster.tree = tree
fname = file.path(top.dir2, "data.Robj")
save(data, file=fname, compress = TRUE)
message(sprintf("Written %d clusters in %s", length(unique(data@ident)), fname))

# Clustering assignment 1 (if applicable)
data@meta.data$cluster.assign = NA
in_assig = file.path(top.dir2, "clustering", redc, clustering_id, "assignment.csv")
if(file.exists(in_assig)){
  af = read.csv(in_assig, header=TRUE, row.names = 1, stringsAsFactors = FALSE)
  data@meta.data$cluster.assign = af[as.character(data@ident), "Cell.type"]
}