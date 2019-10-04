require(reticulate)
require(Matrix)
use_python("~/Dev/py3/bin/python")

read.h5ad = function (fname) {
  #' Read an H5AD file
  #' 
  #' @param fname str: The path to the H5AD file. Note that sparse matrices MUST
  #' be stored in the CSC format.
  ad = import("anndata", convert = FALSE)
  sp = import("scipy.sparse", convert = FALSE)
  
  # Extract the file name from the file path and use that as the project name
  project_name = tools::file_path_sans_ext(basename(fname))
  
  data = ad$read_h5ad(fname)
  
  if (py_to_r(sp$issparse(data$X))) {
    # Read in the matrix where rows are cells and columns are genes
    x = sparseMatrix(
      i = as.numeric(x = data$X$indices),
      p = as.numeric(x = data$X$indptr),
      x = as.numeric(x = data$X$data),
      dims = c(py_to_r(data$X$shape[0]), py_to_r(data$X$shape[1])),
      index1 = F
    )
  } else {
    x = py_to_r(data$X)
  }
  
  # We've got to assign cell and gene names onto the x matrix
  colnames(x) = py_to_r(data$var_names$values)
  rownames(x) = py_to_r(data$obs_names$values)
  
  # Read in gene/cell metadata
  obs_data = py_to_r(data$obs)
  var_data = py_to_r(data$var)
  
  # Read in unstructed metadata
  uns_data = py_to_r(data$uns)
  
  return(list(name = project_name, x = x, cell_data = obs_data,
              gene_data = var_data, metadata = uns_data))
}

write.h5ad = function(fname, x, cell_data = NULL, gene_data = NULL, metadata = NULL) {
  #' Write an expression matrix into an H5AD file
  #'
  #' @param x matrix: An expression matrix with cells in rows and genes in columns
  #' @param cell_data data.frame: Contains a row for each cell in the expression matrix
  #' @param gene_data data.frame: Contains a row for each gene in the expression matrix
  #' @param metadata list: Additional metadata such as method parameters
  ad = import("anndata", convert = FALSE)
  
  # Constuct a python instance of AnnData
  adata = ad$AnnData(x)
  
  # Add cell data if exists
  if (!is.null(cell_data)) {
    if (dim(cell_data)[1] != dim(x)[1]) {
      stop("Please ensure that `cell_data` contains the same number of rows as there are cells.")
    }
    # If the data frame is empty, we can't add it
    if (dim(cell_data)[2] > 0) {
      adata$obs = cell_data
    }
  }
  
  # Add gene data if exists
  if (!is.null(gene_data)) {
    if (dim(gene_data)[1] != dim(x)[2]) {
      stop("Please ensure that `gene_data` contains the same number of rows as there are genes.")
    }
    # If the data frame is empty, we can't add it
    if (dim(gene_data)[2] > 0) {
      adata$obs = gene_data
    }
  }
  
  # Add metadata if exists
  if (!is.null(metadata)) {
    adata$uns$update(metadata)
  }
  
  adata$write_h5ad(fname)
}
