#!/usr/bin/env Rscript
require(optparse)

# Parameters
save = TRUE
random.start = 30
res.low = 0.1
res.high = 10.1

hlp = sprintf("Cluster data by Louvain clustering, trying different resolutions in [%.2f, %.2f]", res.low, res.high)

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-r", "--reduction"), type="character", help="Dimensionality reduction method (pca, cca, scscope)", default="pca"),
  make_option(c("-n", "--num-clusterings"), type="integer", help="Number of clustering settings to try (odd integer).", default="pca"),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", 
              default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
n.clusterings = opt$`num-clusterings`
config = opt$config
reduction.use = opt$reduction
stopifnot(reduction.use %in% c("pca", "cca", "scscope"))

# Load configturation
source(config)

# Number of test parameters must be odd to compute median correctly
resolution = linspace(res.low, res.high, n.clusterings)

# Create output
top.dir = file.path(out.base.dir, exp_id)
out.dir = file.path(top.dir, "clustering")
xeno.makedirs(out.dir)

# Read existing Seurat data frame (migh)
load(file.path(top.dir, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data

# Process markers
marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE)
data = xeno.frame.markers(data, marker.frame)

# Create a clustering output
cell.ids = row.names(data@meta.data)
cscores = data.frame(resolution=resolution, nclusters=0, silhouette=0, modularity=0)
cassig = matrix(0, nrow = length(cell.ids),
                   ncol = length(resolution))
row.names(cassig) = cell.ids
colnames(cassig) = row.names(cscores)

# Clustering
for(i in 1:nrow(cscores)){
  # Subdirectory
  res = cscores$resolution[i]
  sub.dir = file.path(out.dir, reduction.use, sprintf("%02d", i))
  if(!file.exists(sub.dir)) dir.create(sub.dir, recursive = TRUE)
  
  # Clustering and sink
  di = 1:ncol(data@dr[[reduction.use]]@cell.embeddings)
  data = FindClusters(object = data, 
                      reduction.type = reduction.use, 
                      dims.use = di, 
                      resolution = res, 
                      n.start = random.start,
                      save.SNN = TRUE)
  
  # Compute silhouette; cluster conversion is 1-based
  cassig[,i] = as.numeric(data@ident) - 1
  cscores$nclusters[i] = length(unique(data@ident))
  
  # Build cluster tree; nasty trick to force Seurat to use correct reduction.
  fname = file.path(sub.dir, "tree.pdf")
  pdf(fname, width=12, height = 10)
  temp = data@dr$pca@cell.embeddings
  data@dr$pca@cell.embeddings = data@dr[[reduction.use]]@cell.embeddings
  data = BuildClusterTree(data, pcs.use=di, do.plot = TRUE)
  data@dr$pca@cell.embeddings = temp
  dev.off()
  message(sprintf("Written %s", fname))
  fname = file.path(sub.dir, "tree.Robj")
  tree = data@cluster.tree
  save(tree, file=fname, compress = TRUE)
  message(sprintf("Written %s", fname))
  
  # Marker heatmaps
  H = xeno.marker.heatmap(data, unique(marker.frame$Cell.type))
  fname = file.path(sub.dir, sprintf("cluster_heatmap.pdf"))
  ph = pheatmap(H, cluster_rows=TRUE, show_rownames=TRUE, 
                cluster_cols=TRUE, scale="column",
                filename = fname, width = 14, height = 8)
  message(sprintf("Written %s", fname))
  
  # Plot clusters
  for (ident in c("orig.ident", "ident")){
    for(leg in c(T, F)){
      fname = file.path(sub.dir, sprintf("tsne_plot_%s_leg-%s.pdf", ident, leg))
      pdf(fname, width = 10, height = 8)
      Seurat::DimPlot(data, group.by = ident,
                       do.label = leg, no.legend = T,
                      reduction.use = sprintf("%s.tsne", reduction.use))
      message(sprintf("Written %s", fname))
      dev.off()  
    }
  }
}

# Plot scores
for(sc in c("modularity")){
  fname = file.path(out.dir, sprintf("scores_%s.pdf", sc))
  pdf(fname, width=6, height = 3)
  plot(cscores$nclusters, cscores[,sc], type = "o", 
       ylab=sc, xlab="Clusters")
  grid()
  dev.off()
  message(sprintf("Written %s", fname))
}

# Save clustering results
fname = file.path(out.dir, "cscores.csv")
write.csv(cscores, fname)
message(sprintf("Written %s", fname))

fname = file.path(out.dir, "cassig.csv")
write.csv(cassig, fname)
message(sprintf("Written %s", fname))

# Set clusters and plot distribution; does not set cluster tree
data = xeno.set.clusters(data, cscores, cassig, out.dir = out.dir)

# Save Seurat object
message("Clustering done!")
if(save){
  fname = file.path(top.dir, "data.Robj")
  save(data, file=fname, compress = TRUE)
  message(sprintf("Written %s", fname))
}
