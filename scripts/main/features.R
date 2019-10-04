#!/usr/bin/env Rscript 

hlp = "Augment dimensionality reductions by marker genes for all dimensionality reductions."

require(optparse)

# Parameters
ext = "tiff"  # Output file extension

# Parse command line arguments
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
config = opt$config

# Load auxilliary data
source(config)

# Create output
top.dir = file.path(out.base.dir, exp_id)

# Read existing Seurat data frame (migh)
load(file.path(top.dir, "data.Robj"), envir = .GlobalEnv)
data = .GlobalEnv$data

# Process markers
marker.frame = read.csv(in.markers, header=TRUE, stringsAsFactors = FALSE)
data = xeno.frame.markers(data, marker.frame)

# Currect dimensionality reductions
reductions = names(data@dr)
reductions = reductions[grep("tsne", reductions, invert = TRUE)]

# Plot markers
data@meta.data$ident.char = as.character(data@ident)
for(red in reductions){

  # Ensure t-SNE exists
  red.tsne = sprintf("%s.tsne", red)
  if(!red.tsne %in% names(data@dr)){
    next
  }
  
  # Make dirs
  sub.dir = file.path(top.dir, red)
  marker.dir = file.path(sub.dir, "markers")
  xeno.makedirs(c(sub.dir, marker.dir))
  
  # Draw features
  for(cell in unique(marker.frame$Cell.type)){
    fname = file.path(marker.dir, sprintf("%s.%s", gsub(" ", "-", cell), ext))
    tiff(fname, width = 800, height = 800)
    xeno.feature.plot(data, feature = cell, reduction.use = red.tsne)
    dev.off()
    message(sprintf("Written %s", fname))
  }

  # Plot clusters
  for (column in c("orig.ident", "ident.char")){
    fname = file.path(sub.dir, sprintf("tsne_plot_%s.%s", column, ext))
    tiff(fname, width = 2000, height = 1800)
    Seurat::DimPlot(data, 
                    group.by = column, 
                    pt.size=1.0, label.size = 12,
                    do.label = T, no.legend = T, no.axes=T,
                    reduction.use = red.tsne)
    message(sprintf("Written %s", fname))
    dev.off()
  }
}
message("Feature plots done!")
