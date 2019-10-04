# *Xenoturbella* single-cell RNA-seq analysis

The *Xenoturbella* single-cell RNA analysis pipeline is set up as a
series of scripts, forming an experiment.


### Installation

The project makes use of R and Python 3 (for scScope, RNA velocity). To
install required packages, run the scripts with active interpreters or
virtual environments.

R:

    Rscript requirements.R

Python:

    pip install -r requirements.txt



### Setup

An experiment is identified
by a name and is stored in the output folder under a single directory.
Optionally, the maximum number of cells per library can be specified
- experiments on reduced data can be used for debugging.
The data subset used in this study can comprise all libraries
(`DATA=all`) or the two whole animals (`DATA=mix`).


Example values for debug:

    EXP=test_experiment_all
    MAX_CELLS=1000
    DATA=all


The code is meant to be run in the project root directory, but can
be easily modified to do otherwise. Global paths to data
and code are set in `./scripts/utils.R` .

    # Set base paths
    out.base.dir = "output"
    in.base.dir = "data"
    in.script.dir = "scripts"

    # Outputs
    in.dir = file.path(in.base.dir, "matrix")
    in.h5ad.dir = file.path(in.base.dir, "h5ad")

    # Marker sheets
    in.models = file.path(in.base.dir, "markers", "Xenotubella_gene_models_all.csv")
    in.markers = file.path(in.base.dir, "markers", "Xenoturbella_protein_gene_markers.csv")
    in.neuro.markers = file.path(in.base.dir, "markers", "Xenoturbella_neurotransmitters.csv")


Some analyses require to uncompress the following large files inside the `data/` directory.

    - `matrix.tar.gz` (Main analysis)
    - `h5ad.tar.gz ` (scScope imputation and batch effect correction)
    - `loom.tar.gz` (RNA velocity)



### Loading own files

To run the analysis pipeline on custom data sets, modify the
`xeno.load.design` function, which constructs should return
a data frame with columns

- `path`: Name of directory in 10x format in `matrix/`.
- `name`: Human-readable name of the library.
- `color`: Color of library to be used in plots.
- `animal`: Optional "batch" parameter.

Default definition:

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




## Pre-processing, QC and basic visualization

Scripts related to basic pre-processing and data display.

### scScope data imputation and batch effect removal

The data sparsity and batch effect provide the major challenges in this
study.  The <a href="https://github.com/AltschulerWu-Lab/scScope">scScope</a>
autoencoder model is used to correct the gene expression matrix for
these confounding factors. The autoencoder depth and the number of latent
 factors was chosen according to the methods section (see paper).
 The resulting files with imputed gene expression matrix
 are readily available in this repository (`./data/h5ad/scScope__all.h5a`).

    IN="./data/h5ad/all.h5ad"
    OUT_DIR="./data/h5ad/all.h5ad"
    DEPTH=128
    LATENT=50
    ./scScope.py --input $IN --output-dir $OUT_DIR --batch-var orig.ident
    --n-latent $LATENT --depth-encoder $DEPTH --depth-decoder $DEPTH


### Basic quality control

This scripts counts the number of genes per cell, computes distribution
of detected genes per cell and gene expression matrix density. An example
run

    ./scripts/main/basic-qc.R -e $EXP -m $MAX_CELLS

Stores the output in `./output/test_experiment/basic_qc/`.

The parameter specification for each script is obtained with `-h` switch.

    ./scripts/main/basic-qc.R -h

prints:

    Usage: ./scripts/main/basic-qc.R [options]
    Basic quality control script. Count the number of genes/cell and gene expression matrix density.

    Options:
        -e EXPERIMENT, --experiment=EXPERIMENT
            Experiment name.

        -m MAX, --max=MAX
            Maximum number of cells to load per library.

        -c CONFIG, --config=CONFIG
            Configuration script path (scripts/utils.R)

        -h, --help
            Show this help message and exit

### Dimensionality reduction and visualization

This basic script constructs the central Seurat data object and
performs one or more dimensionality reductions: `pca`, `cca`
or `scscope`. Each dimensionality reduction method is the used
by t-SNE to produce a 2D visualization.

    ./scripts/main/preprocess.R -e $EXP -m $MAX_CELLS -r scscope

Produces

    output/test_experiment/pca/tsne_plot.pdf
    output/test_experiment/scscope/tsne_plot.pdf

The R data object to be used by downstream scripts is stored in
`output/test_experiment/data.Robj`.



### Gene expression overlay

Augment 2D visualization with gene expression plots. Gene markers
are read from `data/markers/Xenoturbella_protein_gene_markers.csv`

     ./scripts/main/features.R -e $EXP

Produces marker gene expression plots for all dimensionality reductions
performed previously.


    output/test_experiment/pca/markers/Cell_Adhesion.tiff
    output/test_experiment/pca/markers/Innexins.tiff
    output/test_experiment/pca/markers/Ion_Channels.tiff
    output/test_experiment/pca/markers/Neuropeptides.tiff
    output/test_experiment/pca/markers/Neurotransmitter_Biosynthesis.tiff

Example visualization of Muscle cells:

![Muscle](./output/test_experiment/pca/markers/Muscle.1.tiff)


## Clustering and cluster analysis

In the following, we devise methods for finding clusters of cells in *Xenoturbella*.


### Identifying clusters of cells by Louvain clustering

The potential number of different cell clusters - cell types - are identified by Louvain clustering.  This method automatically find the number of clusters using modularity optimization. It however depends on a resolution hyperparameter, governing the minimal distance between for two cells to be deemed similar. Here, we try *N* uniformly separated values of the resolution between empirically selected 0.1 and 10.1. 

    REDC=pca
    N=10
    ./scripts/clustering/clustering-louvain.R -e $EXP -n $N -r $REDC


For each configuration, the script computes the clustering outcome: clustering assignment and a t-SNE visualization augmented with cluster identity.

    output/test_experiment/clustering/pca/1/tsne_plot_orig.ident_leg-TRUE.pdf
    output/test_experiment/clustering/pca/1/tsne_plot_ident_leg-TRUE.pdf
    ...
    output/test_experiment/clustering/pca/20/tsne_plot_orig.ident_leg-TRUE.pdf
    output/test_experiment/clustering/pca/20/tsne_plot_ident_leg-TRUE.pdf

Clustering assignments for different configurations are stored in:

    output/test_experiment/clustering/pca/
    output/test_experiment/clustering/pca/


### Choosing the number of clusters with a second experiment

We use a biologically-inspired strategy to decide on the number of clusters. Namely, we perform another experiment, consisting of only cells from the two whole animals (libraries labelled *mix*). Intuitively, the two whole animal samples should contain similar cell types as the four individual parts, and the cell types should be present in roughly consistent proportions. This strategy disregards and experimental/technical biases.

Suppose we run another preprocessing and clustering experiment, this time selecting only the cells from the two whole animals.

    EXP2=test_experiment_mix
    DATA2=mix
    MAX_CELLS=2000
    REDC=pca
    ./scripts/main/preprocess.R -e $EXP2 -d $DATA2 -m $MAX_CELLS -r $REDC
    ./scripts/clustering/clustering-louvain.R -e $EXP2 -n $N -r $REDC


This second experiment (`test_experiment_mix`) is then used to select the number of clusters in our main experiment (`test_experiment_all`). Briefly, the two experiments are combined to find a clustering configuration that is consistent in both experiments - having a similar number of clusters and similar cluster size distribution (see paper). The chosen clustering assignment is stored in the R object corresponding to the first experiment (`-e $EXP`).

    ./scripts/clustering/clustering-analysis.R -e $EXP -E $EXP2 -r $REDC


In this artificial case, we decide to select X clusters, supported by both experiments.


### Tree of clusters

The clusters can be represented using hierarchical clustering with significance (obtained by package <a href="">pvclust</a>) to distinguish significant versus non-significant clusters of cells.

The following dimensionality reduction are available.

- **hvg**: Highly variable genes.
- **scscope**, **pca**, **cca**: Output of scScope, PCA or CCA, respectively.

The script is run as follows.

    ./scripts/clustering/clustering-pvclust.R -e $EXP -r $REDC

The clustering tree is stored as follows.

    output/test_experiment_all/clustering/pca/clustering_4_pvclust.pdf
    output/test_experiment_all/clustering/pca/tree.Robj

### Tree and clusters of neurons

A similar visualization is used to show relations between potential neurons in the dataset, stratified by the expressing neurotransmitter molecule. The following script also computes a matching t-SNE visualization that emphasizes the clusters and the library of origin.
    
Dimensionality reductions are the same as in the tree of clusters (above).

    ./scripts/clustering/clustering-pvclust-neurons.R -e $EXP -r $REDC

The output files display three t-SNE plots and one hierarchical clustering tree.

    output/test_experiment_all/clustering/pca-neurons/clustering_5_pvclust.pdf
    output/test_experiment_all/clustering/pca-neurons/combined-tsne-cluster.pdf
    output/test_experiment_all/clustering/pca-neurons/combined-tsne-neurotransmitter.pdf
    output/test_experiment_all/clustering/pca-neurons/combined-tsne-sample.pdf

## Differential gene expression

Differential gene expression between clusters, super-clusters (see paper) and neurons.


### Differential expression test

The test type (`--type $TYPE`) determines that data subset and stratification to be tested:

- **clusters**: Differentially expressed genes between cell clusters.
- **neurons**: Test between 6 neurotramsmitters and all other cells.
- **neurons-specific**: Test specifically for the cells expressing the neurotransmitters.
- **tree**: Test between different branches of clustering tree. The parameter `--tree-start` must be set to determine at which cluster to split the tree.
- **superclusters**: Test between empirically determined cell types: extracellular, primed and undifferentiated.
- **parts**: Test between six libraries.

The test is run as follows:

    TYPE="clusters"
    MAX_GENES=100
    ./scripts/de/diffexp.R -e $EXP -M $MAX_GENES -t $TYPE

The result is a differential gene expression table `output/test_experiment_all/diffexp-clusters/diffexp.csv`, with one gene per row.

- **geneID**: Gene identifier.  
- **n_cluster**: Number of cells in the group (cluster).   
- **n_other**: Number of other cells. 
- **nexp_cluster**: Number cells in the cluster expressing the gene. 
- **nexp_other**: Number of other cells expressing the gene.   
- **pexp_cluster**: Fraction cells in the cluster expressing the gene. 
- **pexp_other**: Fraction of other cells expressing the gene.   
- **lfc**: Log-fold change between cluster and other cells. 
- **pvalue**: Significance of the test result.  
- **fdr**: False discovery rate, Benjamini-Hochberg corrected p-value. 
- **test.type**: Selected test type.   
- **cluster**: Cluster identifier. 
- **marker**: Is the gene a previously known marker (File `in.markers`).  
- **annotation**: Gene annotation.

### Differential expression heatmaps

The results of a given differential expression can be shown in a heatmap. 
Dimesionality reduction is used to load a clustering tree, if available. 

    TYPE="clusters"
    N_GENES=1
    ./scripts/de/heatmap.R -e $EXP -M $N_GENES -t $TYPE -r $REDC


The output files are gene expression heatmaps in `output/test_experiment_all/diffexp-clusters/heatmaps/`.

![heatmap](./output/test_experiment_all/diffexp-clusters/heatmaps/cluster_groups_top_expression_1.tiff)



## RNA velocity

The RNA velocity is computed on a given dimensionality reduction. The velocity estimates
are obtained with the `scvelo` package. A `.loom` file containing all six libraries is 
prepared with this repository.

    LOOM=data/loom/xeno-combined-bc.loom

The script computes RNA velocity estimates and stores them in a pickled object, along with 
the transition matrix. The velocity vectors are embedded on a chosen dimensionality reduction (t-SNE) plot.

    ./scripts/velocity/velocity_compute.py -e $EXP -L $LOOM -r $REDC

The resulting objects are:

    output/test_experiment_all/velocity/transitions.npz
    output/test_experiment_all/velocity/data.pkl.gz

To overlay the velocity plot with the marker genes (stored in `metadata.csv`), run

    ./scripts/velocity/velocity_markers.py -e $EXP

The result are RNA velocity plots embedded into a chosed t-SNE and overlayed with marker genes expression.

    output/test_experiment_all/velocity/markers/velocity_nGene.tiff
    output/test_experiment_all/velocity/markers/velocity_nGene_no_legend.tiff
    output/test_experiment_all/velocity/markers/velocity_nUMI.tiff
    output/test_experiment_all/velocity/markers/velocity_nUMI_no_legend.tiff
    output/test_experiment_all/velocity/markers/velocity_orig.ident.tiff
    output/test_experiment_all/velocity/markers/velocity_orig.ident_no_legend.tiff


Example output (Library of origin, `orig.ident`):
![velocity-example](./output/test_experiment_all/velocity/markers/velocity_orig.ident.tiff)


### R session info

    > sessionInfo()
    R version 3.5.1 (2018-07-02)
    Platform: x86_64-apple-darwin15.6.0 (64-bit)
    Running under: macOS Sierra 10.12.6

    Matrix products: default
    BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
     [1] RColorBrewer_1.1-2          Rtsne_0.15                  ramify_0.3.3                pracma_2.2.2
     [5] pheatmap_1.0.12             MAST_1.8.2                  SingleCellExperiment_1.4.1  SummarizedExperiment_1.12.0
     [9] DelayedArray_0.8.0          BiocParallel_1.16.0         matrixStats_0.54.0          Biobase_2.42.0
    [13] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1         IRanges_2.16.0              S4Vectors_0.20.1
    [17] BiocGenerics_0.28.0         cluster_2.0.7-1             Seurat_2.3.4                Matrix_1.2-15
    [21] cowplot_0.9.3               ggplot2_3.1.0               pvclust_2.0-0

    loaded via a namespace (and not attached):
      [1] colorspace_1.3-2       class_7.3-14           modeltools_0.2-22      ggridges_0.5.1         mclust_5.4.2
      [6] htmlTable_1.12         XVector_0.22.0         base64enc_0.1-3        rstudioapi_0.8         proxy_0.4-22
     [11] npsurv_0.4-0           flexmix_2.3-14         bit64_0.9-7            mvtnorm_1.0-8          codetools_0.2-15
     [16] splines_3.5.1          R.methodsS3_1.7.1      lsei_1.2-0             robustbase_0.93-3      knitr_1.20
     [21] Formula_1.2-3          jsonlite_1.5           ica_1.0-2              kernlab_0.9-27         png_0.1-7
     [26] R.oo_1.22.0            compiler_3.5.1         httr_1.3.1             backports_1.1.2        assertthat_0.2.0
     [31] lazyeval_0.2.1         lars_1.2               acepack_1.4.1          htmltools_0.3.6        tools_3.5.1
     [36] igraph_1.2.2           GenomeInfoDbData_1.2.0 gtable_0.2.0           glue_1.3.0             RANN_2.6
     [41] reshape2_1.4.3         dplyr_0.8.0.1          Rcpp_1.0.1             trimcluster_0.1-2.1    gdata_2.18.0
     [46] ape_5.2                nlme_3.1-137           iterators_1.0.10       fpc_2.1-11.1           gbRd_0.4-11
     [51] lmtest_0.9-36          stringr_1.3.1          irlba_2.3.3            gtools_3.8.1           DEoptimR_1.0-8
     [56] zlibbioc_1.28.0        MASS_7.3-51.1          zoo_1.8-4              scales_1.0.0           doSNOW_1.0.16
     [61] yaml_2.2.0             reticulate_1.10        pbapply_1.3-4          gridExtra_2.3          rpart_4.1-13
     [66] segmented_0.5-3.0      latticeExtra_0.6-28    stringi_1.2.4          foreach_1.4.4          checkmate_1.8.5
     [71] caTools_1.17.1.1       bibtex_0.4.2           Rdpack_0.10-1          SDMTools_1.1-221       rlang_0.3.4.9003
     [76] pkgconfig_2.0.2        dtw_1.20-1             prabclus_2.2-6         bitops_1.0-6           lattice_0.20-38
     [81] ROCR_1.0-7             purrr_0.2.5            htmlwidgets_1.3        bit_1.1-14             tidyselect_0.2.5
     [86] plyr_1.8.4             magrittr_1.5           R6_2.3.0               snow_0.4-3             gplots_3.0.1
     [91] Hmisc_4.1-1            pillar_1.3.1           foreign_0.8-71         withr_2.1.2            fitdistrplus_1.0-11
     [96] mixtools_1.1.0         abind_1.4-5            RCurl_1.95-4.11        survival_2.43-1        nnet_7.3-12
    [101] tibble_2.1.1           tsne_0.1-3             crayon_1.3.4           hdf5r_1.0.1            KernSmooth_2.23-15
    [106] grid_3.5.1             data.table_1.11.8      metap_1.0              digest_0.6.18          diptest_0.75-7
    [111] tidyr_0.8.2            R.utils_2.7.0          munsell_0.5.0