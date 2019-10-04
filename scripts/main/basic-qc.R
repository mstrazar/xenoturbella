#!/usr/bin/env Rscript

hlp = "Basic quality control script. Count the number of genes/cell and gene expression matrix density."

require(optparse)

# Parameters
min.genes = 0   # Minimum genes detected per cell


# Parse hyperparameters
option_list = list(
  make_option(c("-e", "--experiment"), type="character", help="Experiment name."),
  make_option(c("-m", "--max"), type="integer", help="Maximum number of cells to load per library.", default=Inf),
  make_option(c("-c", "--config"), type="character", help="Configuration script path (scripts/utils.R)", 
              default=file.path("scripts", "utils.R"))
);
opt_parser = OptionParser(option_list=option_list, description=hlp);
opt = parse_args(opt_parser);
exp_id = opt$experiment
nmax = opt$max
config = opt$config

# Load configuration
source(config)

# Create output
exp.dir = file.path(out.base.dir, exp_id, "basic_qc")
xeno.makedirs(exp.dir)

# Load design
design = xeno.load.design()

# Iterate through samples
for(i in 1:nrow(design)){
  
  # Custom output
  out.dir = file.path(exp.dir, design[i,]$name)
  xeno.makedirs(out.dir)
  
  # Load and normalize data
  data = xeno.load.data(design[i,], min.genes=min.genes, nmax=nmax)
  bin = data@raw.data > 0
  
  # Histogram (genes)
  fname = file.path(out.dir, "ngene.pdf")
  pdf(fname, width=5, height = 3)
  m = mean(data@meta.data$nGene)
  s = sd(data@meta.data$nGene)
  hist(data@meta.data$nGene, xlab = "Num. genes", main = sprintf("%s (%.2f ± %.2f, density %.2f %%)", 
                                                                 gsub(data@project.name, "matrix", ""),
                                                                 m, s, 100*mean(bin)))
  dev.off()
  message(sprintf("Written %s", fname))
  
  # Histogram (cells)
  fname = file.path(out.dir, "ncell.pdf")
  pdf(fname, width=5, height = 3)
  p = mean(Matrix::rowSums(bin) > 0)
  hist(log2(Matrix::rowSums(bin)), xlab = "log2 (cells)", 
       main = sprintf("%s (%.2f %% genes n.z.)", 
                      data@project.name, 100*p), xaxt="n")
  axis(1, at=0:10, labels=2**(0:10))
  dev.off()
  message(sprintf("Written %s", fname))
}