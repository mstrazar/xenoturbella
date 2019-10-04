#!/usr/bin/env Rscript

hlp = "Pre-processing and dimensionality reduction script."

require(optparse)


# Fixed parameters
dims.use = 1:100      # Number of principal components to use
top.var.genes = 500   # Top variable genes to use
min.genes = 0         # Minimum genes per cell
num.threads = 4       # Number of CPU threads
max_iter = 3000       # Max t-SNE iterations.
perplexity = 50       # t-SNE perplecity
save = TRUE           # Save experiment object

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-m", "--max"), type="integer", help="Maximum number of cells to load per library.", default=Inf),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", 
              default=file.path("scripts", "utils.R")),
  make_option(c("-d", "--data"), type="character", help="Data subset (all/mix)", default="all"),
  make_option(c("-r", "--reduction"), type="character", help="Dimensionality reduction method (pca, cca, scscope)", 
              default="pca")
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
nmax = opt$max
config = opt$config
dset = opt$data
redc = opt$reduction
stopifnot(redc %in% c("pca", "cca", "scscope"))

# Load configturation
source(config)

# Parameters
reductions = c("pca", redc)

# Create output
out.dir = file.path(out.base.dir, exp_id)
xeno.makedirs(out.dir)

# Paths
design = xeno.load.design(dset=dset)

# Load and normalize data
data = xeno.load.data(design, min.genes=min.genes, nmax=nmax,
                      low.thresholds = 1, high.thresholds = Inf,
                      n.genes=top.var.genes)

# Process markers
marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE)
data = xeno.frame.markers(data, marker.frame)

# Store meta data to disk
out.meta.csv = file.path(out.dir, "metadata.csv")
write.csv(data@meta.data, out.meta.csv, quote=FALSE)
message(sprintf("Written %d rows to %s", nrow(data@meta.data), out.meta.csv))

# Subset and use CCA
if("cca" %in% reductions){
  data.list = list()
  for(name in design$name){
    inxs = which(data@meta.data$orig.ident == name)
    data.list[[length(data.list)+1]] = SubsetData(data, cells.use = inxs)
  }
  set.seed(42)
  if(length(data.list) == 2){
    data <- RunCCA(data.list[[1]], data.list[[2]], genes.use = data@var.genes, num.cc = max(dims.use))
  } else {
    data <- RunMultiCCA(object.list = data.list, genes.use = data@var.genes, num.ccs = max(dims.use))  
  }
  data <- AlignSubspace(data, reduction.type = "cca", grouping.var = "orig.ident", dims.align = dims.use)
}

if("pca" %in% reductions){
  # Run PCA  with seeding
  message("Running PCA ...")
  set.seed(42)
  data = RunPCA(data, pcs.print = FALSE, pcs.compute = max(dims.use))
  
  fname = file.path(out.dir, "pca_plot.pdf")
  pdf(fname, width = 10, height = 8)
  Seurat::PCAPlot(data, group.by = "orig.ident", pt.size=0.35)
  message(sprintf("Written %s", fname))
  dev.off()
  
  fname = file.path(out.dir, "pca_elbow.pdf")
  pca.sd = data@dr$pca@sdev
  pdf(fname)
  plot(pca.sd, xlab="PC", ylab="Standard deviation", pch=20)
  grid()
  message(sprintf("Written %s", fname))
  dev.off()
}
if ("scscope" %in% reductions){
  data = xeno.scscope.reduction(data, dset)
}


# Run TSNE for each reduction
for(red in reductions){
  
  # Make dirs
  sub.dir = file.path(out.dir, red)
  xeno.makedirs(sub.dir)
  
  message(sprintf("Running TSNE (%s)... ", red))
  set.seed(42)
  di = 1:ncol(data@dr[[red]]@cell.embeddings)
  data = RunTSNE(data, dims.use = di,  reduction.use = red,
                 max_iter=max_iter, perplexity=perplexity, verbose=TRUE,
                 num_threads=num.threads, reduction.name = sprintf("%s.tsne", red))
  fname = file.path(sub.dir, "tsne_plot.pdf")
  pdf(fname, width = 9, height = 7.5)
  Seurat::DimPlot(data, 
                   group.by = "orig.ident", 
                   pt.size=0.35, 
                   reduction.use = sprintf("%s.tsne", red))
  message(sprintf("Written %s", fname))
  dev.off()
  
  # Store tSNE coordinates to disk
  fname = file.path(sub.dir, "tsne_coordinates.csv")
  write.csv(data@dr[[sprintf("%s.tsne", red)]]@cell.embeddings, fname, quote=FALSE, row.names=TRUE)
  message(sprintf("Written %s", fname))
}

# Save Seurat object
message("Preprocessing done!")
if(save){
  fname = file.path(out.dir, "data.Robj")
  save(data, file=fname, compress = TRUE)
  message(sprintf("Written %s", fname))
}
