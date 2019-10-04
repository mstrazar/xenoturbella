#!/usr/bin/env Rscript
require(optparse)

# Parameters
save = TRUE
lfc.thresh = 0.6
fdr.thresh = 0.05

hlp = "Differential expression analysis between (super-)clusters, neurons or animal parts."

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-M", "--Max"), type="integer", help="Maximum number of genes to test", default=Inf),
  make_option(c("-t", "--type"), type="character", 
                                 help="Test type: clusters, neurons, neurons-specific, parts, tree, supercluster", 
                                 default="clusters"),
  make_option(c("-s", "--tree-start"), type="integer", help="Clustering tree split", default=NA),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", 
              default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
max.genes = opt$Max
type = opt$type
tree.start = opt$`tree-start`
config = opt$config
stopifnot(type %in% c("clusters", "neurons", "neurons-specific", "parts", "tree", "supercluster"))

# Load configturation
source(config)

# Create output
top.dir = file.path(out.base.dir, exp_id)
out.dir = file.path(top.dir, sprintf("diffexp-%s", type))
if(!is.na(tree.start)){
  out.dir = file.path(top.dir, sprintf("diffexp-%s-%d", type, tree.start))
}
temp.dir = file.path(out.dir, "temp")
xeno.makedirs(c(out.dir, temp.dir))

# Read existing Seurat data frame (migh)
load(file.path(top.dir, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data

# Load gene models
models.frame = xeno.load.models(in.models)

# Which genes to test
if(is.finite(max.genes)){
  genes.test = row.names(data@hvg.info[1:max.genes,])  
} else {
  inxs = which(data@hvg.info$gene.dispersion > 0)
  genes.test = row.names(data@hvg.info[inxs,])  
}

# Set predicted neuronal cell types 
if (type == "clusters"){
  test.type = "wilcox"
  marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE, sep=";")
  data = xeno.frame.markers(data, marker.frame)  
} else if(type == "neurons"){
  # need marker frame to exclude genesmarkers from interpretation
  test.type = "prop"
  marker.frame = read.csv(in.neuro.markers, header=TRUE, stringsAsFactors = FALSE)
  data = xeno.set.neurons(data, marker.frame)
} else if (type == "neurons-specific"){
  # need marker frame to exclude genesmarkers from interpretation
  test.type = "prop"
  marker.frame = read.csv(in.neuro.markers, header=TRUE, stringsAsFactors = FALSE)
  data = xeno.set.neurons.specific(data, marker.frame)
} else if (type == "parts"){
  test.type = "wilcox"
  parts = c("dorsal_V3", "skin_V3", "ventral_1", "ventral_2hr")
  marker.frame = read.csv(in.neuro.markers, header=TRUE, stringsAsFactors = FALSE)
  data = xeno.set.parts(data, parts)
} else if (type == "tree"){
  test.type = "wilcox"
  marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE)
  data = xeno.frame.markers(data, marker.frame)
  data = xeno.set.tree.clusters(data, start = tree.start)
} else if (type == "supercluster"){
  test.type = "wilcox"
  marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE)
  data = xeno.frame.markers(data, marker.frame)
  data = xeno.set.superclusters(data)
}

# Run test
results = st.diff(data, genes.use = genes.test, temp.dir = temp.dir, 
                  lfc.thresh = lfc.thresh, fdr.thresh = fdr.thresh, 
                  marker.frame = marker.frame, models.frame = models.frame,
                  test.type = test.type)

# Write output
fname = file.path(out.dir, "diffexp.csv")
write.csv(results, fname)
message(sprintf("Written %d rows to %s", nrow(results), fname))