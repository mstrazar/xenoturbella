### Basic configuration for running Xenoturbella analysis ###

# Set base paths
out.base.dir = "output"
in.base.dir = "data"
in.script.dir = "scripts"

# Outputs
in.dir = file.path(in.base.dir, "matrix")
in.h5ad.dir = file.path(in.base.dir, "h5ad")

# Marker sheets
in.models = file.path(in.base.dir, "markers", "Xenoturbella_gene_models_all.csv")
in.markers = file.path(in.base.dir, "markers", "Xenoturbella_protein_gene_markers.csv")
in.neuro.markers = file.path(in.base.dir, "markers", "Xenoturbella_neurotransmitters.csv")

# Imports
source(file.path(in.script.dir, "stats.R"))
require(MAST)
require(Matrix)
require(RColorBrewer)
require(Seurat)
require(cluster)
require(pheatmap)
require(pracma)
require(pvclust)
require(ramify)
require(reticulate)


# Define experimental design
xeno.load.design <- function(dset="all"){
  design = data.frame(path=c(
    "Xenoturbella_dorsal_matrix",
    "Xenoturbella_mix1_matrix",
    "Xenoturbella_mix2_matrix",
    "Xenoturbella_skin_matrix",
    "Xenoturbella_ventral_1_5hr_matrix",
    "Xenoturbella_ventral_2hr_matrix"),
      animal = c("animalC", 
                 "animalAB", 
                 "animalAB", 
                 "animalC", 
                 "animalC", 
                 "animalC"),
      color = c("#f8766d", "#b79f00", "#00ba38", "#00bfc4", "#619cff", "#f564e3"),
    stringsAsFactors = FALSE)
  design$name = unlist(lapply(strsplit(design$path, "_"), 
                              FUN = function(x){paste(x[2:3], collapse = "_")}))
  
  if(dset != "all"){
    design = design[grep(dset, design$path),]
  }
  row.names(design) = design$name
  return(design)
}


# Load and normalize Xenoturbella data
xeno.load.data <- function(design, min.genes=50, nmax=NA,
                           low.thresholds = 1, high.thresholds = Inf, n.genes = 1000,
                           preprocess=TRUE){
  # Read files
  data = NA
  for(i in 1:nrow(design)){
    path = design[i,]$path
    name = design[i,]$name
    animal = design[i,]$animal
    message(sprintf("Loading %s", path))
    data.10x = Read10X(file.path(in.dir, path))
    if(!is.na(nmax)) data.10x = data.10x[,1:min(nmax, ncol(data.10x))]
    colnames(data.10x) = sprintf("%s-%d", colnames(data.10x), i) 
    seu = CreateSeuratObject(data.10x, project = name,
                             min.genes = min.genes)
    seu@meta.data$animal = animal
    message(sprintf("Read %d cells from %s", ncol(seu@data), name))
    if(i == 1){
      data = seu
    } else {
      data = MergeSeurat(data, seu, do.normalize = FALSE)
    }
  }
  
  # Filter cells
  if(preprocess){
    data <- FilterCells(object = data, 
                          subset.names = "nGene", 
                          low.thresholds = low.thresholds, 
                          high.thresholds = high.thresholds)
    # Normalize and HVG
    data = NormalizeData(data,
                         normalization.method = "LogNormalize", 
                         scale.factor = 1e4)
    data = FindVariableGenes(data, mean.function = ExpMean,
                             dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125,
                             x.high.cutoff = 10, 
                             y.cutoff = 0.5,
                             do.plot = FALSE)
    
    
    # Select variable genes
    d = data@hvg.info$gene.dispersion
    names(d) = row.names(data@hvg.info)
    data@var.genes = names(rev(sort(d))[1:n.genes])
    message(sprintf("Reduced to %d variable genes", n.genes))
    
    data = ScaleData(data,
                     genes.use = data@var.genes, 
                     vars.to.regress = "nUMI")
  }
  data
}


# Load gene models
xeno.load.models <- function(in.models){
  df = read.csv(in.models, header=TRUE, stringsAsFactors = FALSE)
  df$description = gsub(".t[0-1]*", "", df$description) 
  row.names(df) <- df$description
  return(df)
}


# Aggregate marker data to a seruat object
xeno.agg.markers <- function(obj, marker.list, name){
  aval.markers = row.names(obj@data)
  markers = intersect(aval.markers, marker.list)
  if(length(markers) == 0){
    message(sprintf("No markers for group %s", name))
    return(obj)
  } else if (length(markers) == 1) {
    agg = as.vector(obj@data[markers,])
  } else {
    agg = colSums(as.matrix(obj@data[markers,]))
  }  
  obj@meta.data[,name] = agg
  message(sprintf("Aggregated %d markers for %s", length(markers), name))
  return(obj)
}

# Add columns for each marker group in a data frame
xeno.frame.markers <- function(obj, marker.frame){
  for (cell in unique(marker.frame$Cell.type)){
    marker.list = subset(marker.frame, Cell.type == cell)$Name
    obj = xeno.agg.markers(obj, marker.list, cell)
  }
  return(obj)
}

# Compute a heatmap for clusters and markers
xeno.marker.heatmap <- function(obj, marker.list, FUN=mean, max_clusters=10){
  ident = obj@ident
  markers = marker.list
  ac = aggregate(ident, by=list(cluster=ident), length)
  ac = ac[rev(order(ac$x)),]
  if(!is.null(max_clusters)){
    ac = ac[1:max_clusters,]
  }
  clusters = ac$cluster
  cells = which(ident %in% clusters)
  agg = aggregate(obj@meta.data[cells, markers], 
                  by=list(cluster=ident[cells]), FUN)
  row.names(agg) = agg$cluster
  M = as.matrix(agg[,markers])
  M = M[, colSums(M != 0) > 0]
  return(M)
}


# Aggregate reads for each marker group
xeno.marker.reduction <- function(obj, marker.frame, noise=0.001, reduction.name="gene.groups"){
  marker.frame = marker.frame[!marker.frame$Validation,]
  markers = marker.frame$Name
  groups = unique(marker.frame$Cell.type)
  stopifnot(!("noise" %in% groups))
  M = matrix(0, nrow = length(groups), ncol=ncol(obj@raw.data))
  row.names(M) = groups
  colnames(M) = colnames(obj@raw.data)
  keep = colnames(M)
  for(g in groups){ 
    tmp = subset(marker.frame, Cell.type == g)$Name
    ms = intersect(row.names(obj@raw.data), tmp)
    if(length(tmp) != length(ms)){
      message(sprintf("Markers for %s missing from data (%d): %s ...", 
                      g, length(tmp) - length(ms), setdiff(tmp, ms)[1]))
    }
    if(length(ms) == 0){
      message(sprintf("No markers for %s, skipping ...", g))
      keep = setdiff(keep, g)
    } else if(length(ms) == 1){
      M[g,] = obj@raw.data[ms,]
    } else {
      M[g,] = colSums(as.matrix(obj@raw.data[ms,]))  
    }
  }
  # Perform standard Seurat normalization
  seu = CreateSeuratObject(M)
  seu = NormalizeData(seu,
                     normalization.method = "LogNormalize", 
                     scale.factor = 1e4)
  seu = ScaleData(seu, vars.to.regress = "nUMI")

  # Add a noise column
  emb = as.data.frame(t(seu@scale.data))
  emb$noise = rnorm(nrow(emb), sd=noise)
  
  # Filter to match only existing cells
  emb = emb[obj@cell.names,]
  
  # Fill in the slots of original objec
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "cell.embeddings", 
                  new.data = as.matrix(emb))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "key", 
                  new.data = "gene.groups")
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "gene.loadings", 
                        new.data = as.matrix(rowMeans(M > 0)))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "sdev", 
                        new.data = mean(M > 0))
  return(obj)
}





# Filter reads for a subset of marker genes
xeno.special.reduction <- function(obj, special.frame, col.name="Short", terms.use=c("None"), dims.use=1:15,
                                   reduction.name="noannot", noise=0.001){
  reduced.frame = special.frame[special.frame[,col.name] %in% terms.use,]
  markers = intersect(reduced.frame$Name, row.names(obj@raw.data)) 
  message(sprintf("Found %d markers", length(markers)))
  M = obj@data[markers,]
  
  # Subset data to exclude invalid cells
  total = length(obj@cell.names)
  keep = which(colSums(M) > 0 & !duplicated(as.matrix(t(M)))) 
  M = M[,keep]
  obj = SubsetData(obj, cells.use = keep)
  message(sprintf("Keeping %d/%d cells", length(keep), total))
  
  # Perform standard Seurat normalization; do not throw away data
  seu = CreateSeuratObject(M, min.genes = -1)
  seu = ScaleData(seu, vars.to.regress = "nUMI")
  seu = RunPCA(seu, pcs.compute = max(dims.use), pc.genes = markers, do.print = FALSE)
  emb = as.data.frame(seu@dr$pca@cell.embeddings)
  # emb$noise = rnorm(nrow(emb), sd=noise)
  
  # Fill in the slots of original objec
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "cell.embeddings", 
                        new.data = as.matrix(emb))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "key", 
                        new.data = "special")
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "gene.loadings", 
                        new.data = as.matrix(rowMeans(M > 0)))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "sdev", 
                        new.data = mean(M > 0))
  return(obj)
}

# Join multiple reductions
xeno.hybrid.reduction <- function(obj, reductions, reduction.name){
    emb = data@dr[[reductions[1]]]@cell.embeddings
    for (ri in 2:length(reductions)){
      emb = cbind(emb, data@dr[[reductions[ri]]]@cell.embeddings)
    }
    obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "cell.embeddings", 
                          new.data = as.matrix(emb))
    obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "key", 
                          new.data = reduction.name)
    return(obj)
}

# Load DCA reduction from disk
xeno.dca.reduction <- function(obj, dset, reduction.name = "dca", noise=0.001, type="latent",
                               dims.use=1:15){
  cells = obj@cell.names
  if(type == "latent"){
    dca = read.table(file.path(in.dca.dir, sprintf("%s-imputed", dset), sprintf("%s.tsv", type)),
                     row.names=1, header = FALSE)
    emb = dca[cells,]  
    emb$noise = rnorm(nrow(emb), sd=noise)
  } else if(type == "mean"){
    dca = read.table(file.path(in.dca.dir, sprintf("%s-imputed", dset), sprintf("%s.tsv", type)),
                     row.names=1, header = TRUE)
    expr = t(as.matrix(dca[,gsub("-", ".", cells)]))
    pc = prcomp(expr, rank. = max(dims.use))
    emb = as.data.frame(pcoc$x)
    row.names(emb) = cells
    emb$noise = rnorm(nrow(emb), sd=noise)
  } 
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "cell.embeddings", 
                        new.data = as.matrix(emb))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "key", 
                        new.data = reduction.name)
  return(obj)
}


# Load scScope reduction from disk; use all dimensions
xeno.scscope.reduction <- function(obj, dset, reduction.name = "scscope", noise=0.001, type="latent"){
  source("h5ad.R")
  cells = obj@cell.names
  fname = file.path(in.scscope.dir, sprintf("scScope__%s.h5ad", dset))
  message(sprintf("Reading %s ...", fname))
  df = read.h5ad(fname)
  M  = as.matrix(df$metadata$latent)
  row.names(M) = row.names(df$cell_data)
  emb = data.frame(M[cells,])
  emb$noise = rnorm(nrow(emb), sd=noise)
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "cell.embeddings", 
                        new.data = as.matrix(emb))
  obj = SetDimReduction(object = obj, reduction.type = reduction.name, slot = "key", 
                        new.data = reduction.name)
  return(obj)
}

# Load scScope reduction from disk; use all dimensions
xeno.cca.reduction <- function(obj, reduction.name = "cca.animal", dims.use = 1:10,
                               out.dir = NULL){

  # Split data on animal
  data.list = list()
  data.list[[1]] = SubsetData(obj, cells.use = which(obj@meta.data$animal == "animalAB"))
  data.list[[2]] = SubsetData(obj, cells.use = which(obj@meta.data$animal == "animalC"))
  message(sprintf("Size of dataset AB: %d", length(data.list[[1]]@cell.names)))
  message(sprintf("Size of dataset  C: %d", length(data.list[[2]]@cell.names)))
  
  # Compute CCA and order components
  genes.use = obj@var.genes
  Xt = as.matrix(t(data.list[[1]]@data[genes.use,]))
  Yt = as.matrix(t(data.list[[2]]@data[genes.use,]))
  keep = which((colSums(Xt) > 0) & (colSums(Yt) > 0))
  X = Xt[,keep]
  Y = Yt[,keep]
  message(sprintf("Keeping %d variable genes present in both data sets", length(keep)))
  res <- cca.low.rank(X, Y, rank=max(dims.use))
  A = res$As
  B = res$Bs
  if(!is.null(out.dir)){
    sub.dir = file.path(out.dir, reduction.name)
    xeno.makedirs(sub.dir)
    fname = file.path(sub.dir, "cca_scores.pdf")
    pdf(fname, width = 6, height = 4)
    plot.cca.low.rank.scores(res)
    dev.off()
    message(sprintf("Written %s", fname))
    fname = file.path(sub.dir, "cca_plot.pdf")
    pdf(fname, width = 7, height = 7)
    minus = min(min(A[,1:2]), min(B[,1:2]))
    plus = max(max(A[,1:2]), max(B[,1:2]))
    plot(NA, ylim=c(minus, plus), xlim=c(minus, plus), main="CCA", xlab="CCA 1", ylab="CCA 2")
    points(A[, 1], A[,2], pch=20, col="blue")
    points(B[, 1], B[,2], pch=20, col="orange")
    dev.off()
  }
  dims.use = 1:ncol(A)
  row.names(A) = data.list[[1]]@cell.names
  row.names(B) = data.list[[2]]@cell.names
  embedding = rbind(res$A, res$B)  
  row.names(embedding) = c(row.names(A), row.names(B))
  embedding = embedding[obj@cell.names,]
  obj = SetDimReduction(object = obj, reduction.type = "cca", slot = "cell.embeddings", 
                        new.data = as.matrix(embedding))
  obj = SetDimReduction(object = obj, reduction.type = "cca", slot = "key", 
                        new.data = "cca")
  
  # Align subspace
  # obj <- Seurat::AlignSubspace(obj, reduction.type = "cca", grouping.var = "animal", dims.align = dims.use)
  # obj@dr[[reduction.name]] = obj@dr[["cca.aligned"]]
  obj@dr[[reduction.name]] = obj@dr[["cca"]]
  return(obj)
}


# Choose clustering ID and assign clusters based on median value
xeno.set.clusters<- function(obj, cscores, cassig, out.dir=NULL){
  md = median(cscores$nclusters)
  clustering_id = which(cscores$nclusters == md)[1]
  obj@ident = as.factor(cassig[,clustering_id])
  obj@meta.data$ident = obj@ident
  if(!is.null(out.dir)){
    # Plot clustering distribution
    fname = file.path(out.dir, "cluster_select.pdf")
    pdf(fname, width=4, height = 3)
    hist(cscores$nclusters, 
         xlab="Clusters", 
         main=sprintf("Median: %d (c. id=%d)", md, clustering_id))
    dev.off()
    message(sprintf("Written %s", fname))
  }
  return(obj)
}

# Set neuronal clusters
# Use only cells that express marker uniquely (exactly one)
xeno.set.neurons <- function(obj, marker.frame){
  markers = marker.frame
  row.names(markers) = marker.frame$Name
  cells = which(colSums(as.matrix(obj@data[markers$Name,])) > 0)
  X = as.matrix(obj@data[markers$Name, cells])
  types = rep("other", length(cells))
  for(i in 1:length(cells)){
    types[i] = names(which(X[,i] > 0))
  }
  colors = markers[types,]$Color
  names = markers[types,]$Cell.type
  ident = rep("other", length(obj@cell.names))
  ident[cells] = names
  names(ident) =  obj@cell.names
  obj@ident = as.factor(ident)
  obj
}

# Set neuronal clusters - specific neuron families 
# Use only cells that express marker uniquely (exactly one)
xeno.set.neurons.specific <- function(obj, marker.frame){
  markers = marker.frame
  row.names(markers) = marker.frame$Name
  cells = which(colSums(as.matrix(obj@data[markers$Name,])) > 0)
  X = as.matrix(obj@data[markers$Name, cells])
  types = rep("other", length(cells))
  for(i in 1:length(cells)){
    types[i] = names(which(X[,i] > 0))
  }
  colors = markers[types,]$Color
  names = markers[types,]$Cell.type
  o = SubsetDataRepaired(obj, cells.use = names(cells))
  o@ident = as.factor(names)
  names(o@ident) = names(cells)
  o
}

# Set sample ID as a cluster
xeno.set.parts <- function(obj, parts=NULL){
  ident = obj@meta.data$orig.ident
  if(!is.null(parts)){
    inxs = which(ident %in% parts)
    ident = ident[inxs]
    cells = data@cell.names[inxs]
    obj = SubsetDataRepaired(obj, cells.use = cells)
  }
  obj@ident = as.factor(ident)
  obj
}

# Split data into superclusters
xeno.set.superclusters <- function(obj){
  ident = obj@meta.data$cluster.assign
  obj@ident = as.factor(ident)
  obj
}

# Identify top node off tree tree; This is a node that nothing links to
tree.top <- function(edges){
  S = edges[,1]   # start (sources)
  E = edges[,2]   # end (sinks)
  t = setdiff(unique(S), unique(E))
  stopifnot(length(t) == 1)
  t
}

# Get leaves; these are nodes that link to nothing
tree.leaves <- function(edges){
  S = edges[,1]   # start (sources)
  E = edges[,2]   # end (sinks)
  leaves = setdiff(unique(E), unique(S))
  leaves
}

# Traverse tree using a queue; return leaves only
tree.subset <- function(edges, start){
  S = edges[,1]      # starts (sources)
  E = edges[,2]      # ends (sinks)
  nodes = c()        # Nodes that fall under start
  queue = c(start)   # Newly discovered nodes go to queue
  while(length(queue) > 0){
    t = queue[1]
    nodes = c(nodes, t)
    es = setdiff(E[S == t], nodes)
    queue = unique(c(tail(queue, n=length(queue)-1), es))
  }
  leaves = tree.leaves(edges[E %in% nodes, ])
  leaves
}

# Set clusters as branches in the tree ; if branch point is NA, set topmost split
xeno.set.tree.clusters <- function(obj, start=NA){
  edges = obj@cluster.tree[[1]]$edge
  if(is.na(start)){
    start = edges[1, 2]
  }
  leaves = tree.subset(edges, start)
  clusters = as.character(leaves - 1)
  ident = rep(sprintf("%d-", start), length(obj@cell.names))
  inxs = obj@meta.data$ident %in% clusters
  ident[inxs] = sprintf("%d+", start)
  message(sprintf("Set %d cluster %s cells ", sum(inxs), sprintf("%d+", start)))
  obj@ident = as.factor(ident)
  obj
} 

# Make dirs if they don't exist
xeno.makedirs <- function(dirs){
  for (d in dirs){
    if(!file.exists(d)){
      dir.create(d, recursive = TRUE)
      message(sprintf("Created %s", d))
    }
  }  
}


# Feature plot without too much clutter ; distinguish expression strength if you have a lot of cells
xeno.feature.plot <- function(obj, features = "BHLH", reduction.use = "pca", ncells.alpha=100, colors="orange"){
  # Initial assertions
  if(length(colors) == 1){
    colors = rep(colors, length(features))
  }
  stopifnot(length(features) == length(colors))
  
  # Initial setting
  E = obj@dr[[reduction.use]]@cell.embeddings
  plot(E[,1], E[, 2], pch=20, col=alpha("gray", 0.3), 
       xlab="", 
       ylab="", axes = FALSE, cex=0.3)

  # Process multiple features
  count = 0
  for(i in 1:length(features)){
    feature = features[i]
    col = colors[i]
    if(feature %in% colnames(obj@meta.data)){
      Y = obj@meta.data[,feature]
      cells = which(Y > 0)
      count = count + length(cells)
    } else {
      cells = c()
    }
    if(length(cells) > ncells.alpha){
      rbPal <- colorRampPalette(c('blue', col))
      scores = Y[cells] / max(Y)
      pal <- rbPal(10)[as.numeric(cut(Y, breaks = 10))]  
      cells = cells[order(Y[cells])]
      points(E[cells,1], E[cells, 2], col=alpha(pal[cells], 0.7), pch=20, cex=1.3)
    } else {
      if(is.na(col)){
        col = alpha("orange", 0.7)
      }
      points(E[cells,1], E[cells, 2], col=col, pch=20, cex=1.3)
    }
  }
  if(length(features) == 1){
    tit = sprintf("%s (%d)", features, count)
  } else {
    tit = sprintf("Combined (%d)", count)
  }
  title(main=tit)
}
                  
# Repair broken Seurat function SubsetData
SubsetDataRepaired <- function (object, cells.use = NULL, subset.name = NULL, ident.use = NULL,
          ident.remove = NULL, accept.low = -Inf, accept.high = Inf,
          accept.value = NULL, do.center = FALSE, do.scale = FALSE,
          max.cells.per.ident = Inf, random.seed = 1, do.clean = FALSE,
          subset.raw, ...)
{
  data.use <- NULL
  object@cell.names <- cells.use
  object@data <- object@data[, cells.use]
  gc(verbose = FALSE)
  if (!is.null(x = object@scale.data)) {
    if (length(x = colnames(x = object@scale.data) > 0)) {
      object@scale.data[, cells.use]
      object@scale.data <- object@scale.data[, cells.use]
    }
    gc(verbose = FALSE)
  }
  if (do.scale) {
    object <- ScaleData(object = object, do.scale = do.scale,
                        do.center = do.center)
    gc(verbose = FALSE)
  }
  object@ident <- droplevels(x = object@ident[cells.use])
  if (length(x = object@dr) > 0) {
    for (i in 1:length(object@dr)) {
      if (length(object@dr[[i]]@cell.embeddings) > 0) {
        object@dr[[i]]@cell.embeddings <- object@dr[[i]]@cell.embeddings[cells.use,
                                                                         , drop = FALSE]
      }
    }
    gc(verbose = FALSE)
  }
  if (!.hasSlot(object = object, name = "assay")) {
    object@assay <- list()
  }
  if (length(object@assay) > 0) {
    for (i in 1:length(object@assay)) {
      if ((!is.null(x = object@assay[[i]]@raw.data)) &&
          (ncol(x = object@assay[[i]]@raw.data) > 1)) {
        object@assay[[i]]@raw.data <- object@assay[[i]]@raw.data[,
                                                                 cells.use]
      }
      if ((!is.null(x = object@assay[[i]]@data)) && (ncol(x = object@assay[[i]]@data) >
                                                     1)) {
        object@assay[[i]]@data <- object@assay[[i]]@data[,
                                                         cells.use]
      }
      if ((!is.null(x = object@assay[[i]]@scale.data)) &&
          (ncol(x = object@assay[[i]]@scale.data) > 1)) {
        object@assay[[i]]@scale.data <- object@assay[[i]]@scale.data[,
                                                                     cells.use]
      }
      gc(verbose = FALSE)
    }
  }
  object@meta.data <- data.frame(object@meta.data[cells.use,
                                                  ])
  gc(verbose = FALSE)
  if (do.clean) {
    calcs.to.keep <- c("CreateSeuratObject", "NormalizeData",
                       "ScaleData")
    object@calc.params <- object@calc.params[calcs.to.keep]
    object@var.genes <- vector()
    object@hvg.info <- data.frame()
    object@cluster.tree <- list()
    object@snn <- as(matrix(), "dgCMatrix")
    object@scale.data <- matrix()
    object@misc <- NULL
    object@kmeans <- NULL
    object@dr <- list()
    object@meta.data
    if (missing(x = subset.raw)) {
      subset.raw <- TRUE
    }
    object@meta.data[, sapply(colnames(object@meta.data),
                              function(x) {
                                grepl("res", x)
                              })] <- NULL
    gc(verbose = FALSE)
  }
  if (!missing(x = subset.raw)) {
    if (subset.raw) {
      object@raw.data <- object@raw.data[, cells.use]
      gc(verbose = FALSE)
    }
  }
  return(object)
}
