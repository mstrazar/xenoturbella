#!/usr/bin/env Rscript
require(optparse)

# Parameters
lfc.thresh = 0.6     # Minimal log-fold change
min.cells = 20       # Minimal number of cells expressing a gene
width = 15/107       # With per gene


hlp = "Display differential expression results in a heatmap."


# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-r", "--reduction"), type="character", help="Dimensionality reduction method (pca, cca, scscope)", default="pca"),
  make_option(c("-t", "--type"), type="character", help="DE test type: clusters, neurons, neurons-specific, parts, tree, supercluster", 
              default="clusters"),
  make_option(c("-M", "--Max"), type="integer", help="Genes to display per cluster.", default=1),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
diffexp_id = opt$type
reduction.use = opt$reduction
no_genes = opt$Max
config = opt$config

# Load configturation
source(config)

# Select top n result per cluster ; assume table is ordered
filter.top <- function(tab, n=1){
  tab$rank = 0
  for (clu in unique(tab$cluster)){
    inxs = which(tab$cluster == clu)
    tab[inxs,]$rank = 1:length(inxs)
  }  
  tab[tab$rank<=n,]
}

# Reorder table by cluster
reorder.tab <- function(tab, ord){
  ident = as.character(tab$cluster)
  o = unlist(lapply(ident, function(x) which(x == ord)))
  tab[order(o),]
}

# Paths
top.dir = file.path(out.base.dir, exp_id)
dif.dir = file.path(top.dir, sprintf("diffexp-%s", diffexp_id))
res.dir = file.path(dif.dir, "heatmaps")
xeno.makedirs(res.dir)

# Load data
load(file.path(top.dir, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data

# Set clusters
if(diffexp_id == "supercluster"){
  data = xeno.set.superclusters(data)
} 

# Load tree
tree.file = file.path(top.dir, "clustering", reduction.use, "tree.Robj")
if(file.exists(tree.file)){
  message(sprintf("Reordering clusters according to %s", tree.file))
  load(tree.file, envir = .GlobalEnv)
  tree = .GlobalEnv$tree
  cl.order = levels(data@ident)[tree$hclust$order]  
} else {
  cl.order = levels(data@ident)
}


# Load differential expression
dexp = read.csv(file.path(dif.dir, "diffexp.csv"), header = TRUE, stringsAsFactors = FALSE)
keep = which(dexp$lfc >= lfc.thresh & dexp$nexp_cluster >= min.cells)
de = dexp[keep,]

# Gene names
genes = de$geneID
de$label = unlist(lapply(strsplit(strTrim(de$annotation), " "), function(x) x[[1]]))
de$label[de$annotation == "No annotation"] = "No annot."
if (is.numeric(de$cluster)){
  de$label = sprintf("%03d: %s %s", de$cluster, de$geneID, de$label)  
}

# One gene per cluster (top)
de.top = filter.top(de, no_genes)
de.top = reorder.tab(de.top, cl.order)

# One gene per cluster (remove overlapping) - reorder be intersect
de.noo = de[!duplicated(de$geneID),]
de.noo = filter.top(de.noo, no_genes)
de.noo = reorder.tab(de.noo, cl.order)

################
### Heatmaps ###
################

for (mode in c("top", "nooverlap")){
  if (mode == "top") dx = de.top
  else if (mode == "nooverlap") dx = de.noo

  # Sample by Cluster matrix
  C = cluster.one.hot(data, factors = FALSE)
  C = C[,as.character(dx$cluster)]
  csize = colSums(C)
  Z = t(as.matrix(data@data[dx$geneID,]))  
  
  for (quant in c("percentage", "expression")){
    if(quant == "percentage"){
      scale = "none"
      H = t(Z > 0) %*% C
      H = sweep(H, 2, csize, "/")     # Divide by cluster size
    } 
    else if(quant == "expression"){
      clip  = 0.3
      scale = "none"
      H = t(Z) %*% C
      H = sweep(H, 2, csize, "/")     # Divide by column (cluster size)
      H = H / rowSums(H)              # Scale by row (gene)
      H = pmin(H, clip)
    }
    H = H[,colnames(H)[!duplicated(colnames(H))]]
    h = nrow(H) * width
    w = 3 + ncol(H) * width
    for (ext in c("tiff", "pdf")){
      fname = file.path(res.dir, sprintf("cluster_groups_%s_%s_%d.%s", mode, quant, no_genes, ext))
      ph = pheatmap(H, cluster_rows=FALSE, show_rownames=TRUE,
                    cluster_cols=FALSE, scale=scale, labels_row = dx$label,
                    filename = fname, width = w, height = h,
                    na_col = "darkblue")
      message(sprintf("Written %s (%.2f x %.2f)", fname, w, h))
    }
  }
  
  # Write table
  fname = file.path(res.dir, sprintf("diffexp_%s.csv", mode))
  write.csv(dx, fname, quote=TRUE, row.names=FALSE)
  message(sprintf("Written %s", fname))
}
