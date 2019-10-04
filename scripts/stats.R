# Various statistical utilities
require(cluster)
require(MAST)

# Approximate silhouette score stratified by cluster
appr.silhouette <- function(obj, n.samples=NULL, 
                            reduction.use="pca", dims.use=1:15, FUN=median){
  n = ncol(obj@scale.data)
  emb = GetDimReduction(obj, reduction.use, "cell.embeddings")
  if(!is.null(n.samples)){
    inxs = sample(1:n, n.samples, replace = FALSE)
    D = dist(emb[inxs, dims.use])
    x = as.numeric(data@ident[inxs])
  } else {
    D = dist(emb[, dims.use])
    x = as.numeric(data@ident)
  }
  S = silhouette(x, dist = D)
  return(FUN(S[,"sil_width"]))
}

# Create a one-hot encoding for clusters
cluster.one.hot <- function(obj, factors=TRUE){
  clusters = levels(obj@ident)
  cells = obj@cell.names
  M = matrix(0, ncol = length(clusters), nrow=length(cells))
  row.names(M) = cells
  colnames(M) = clusters
  Mf = as.data.frame(M)
  for(c in clusters){
    inxs = which(obj@ident == c)
    M[inxs, c] = 1
    Mf[,c] = as.factor(M[,c])  
  }
  if(factors) return(Mf)
  else return(M)
}
  
# Initialize DE data frame - standard for all tests
st.diff.init.df <- function(genes.use, marker.frame, models.frame, cid=""){
  df = data.frame(geneID=genes.use, n_cluster=0, n_other=0, nexp_cluster=NA, nexp_other=NA, 
                  pexp_cluster=NA, pexp_other=NA, lfc=0,
                  pvalue=1, fdr=0, test.type=NA, cluster=cid, marker=NA, groups=NA, annotation=NA)
  df$marker = genes.use %in% marker.frame$Name
  df$annotation = models.frame[genes.use,]$best_annotation
  df[df$marker,]$groups = unlist(lapply(df$geneID[df$marker], 
                                 function(i) 
                                   paste(marker.frame[marker.frame$Name == i, ]$Cell.type, 
                                         collapse = ";")))
  row.names(df) = genes.use
  return(df)
}


# Wilcox rank-sum test
st.diff <- function(obj, genes.use=NULL, 
                           lfc.thresh=0.6, fdr.thresh=0.05, temp.dir=NULL, marker.frame = NULL,
                           models.frame=NULL, test.type="wilcox"){
  
  # Supported test
  stopifnot(test.type %in% c("wilcox", "prop"))
  
  # Convert to one-hot
  Mf = cluster.one.hot(obj, factors = TRUE)
  
  # Run differential expression
  genes = row.names(obj@data)
  if(!is.null(genes.use)) genes = genes.use
  X = obj@data
  
  # Run test for all genes
  results = data.frame()
  count = 1
  for(cont in colnames(Mf)){
    # Determine candidate genes for cluster
    cells.cluster = which(Mf[,cont] == "1")
    cells.other = which(Mf[,cont] == "0")
    rc = Matrix::rowMeans(X[,cells.cluster])
    ro = Matrix::rowMeans(X[,cells.other])
    genes.cluster = intersect(names(rc[rc > ro]), genes)
    df = st.diff.init.df(genes.cluster, marker.frame, models.frame, cont)
    message(sprintf("Testing %d candidate genes for cluster %s, %d cells, (%d / %d) %s ...", 
                    length(genes.cluster), cont, sum(Mf[,cont] == "1"), 
                    count, length(colnames(Mf)), test.type))
    for (g in genes.cluster){
      x0 = X[g, Mf[,cont] == "0"]  
      x1 = X[g, Mf[,cont] == "1"] 
      m0 = mean(x0)
      m1 = mean(x1)
      n_other = length(x0)
      n_cluster = length(x1)
      nexp_other = sum(x0 > 0)
      nexp_cluster = sum(x1 > 0)
      pexp_other = mean(x0 > 0)
      pexp_cluster = mean(x1 > 0)
      if(test.type == "wilcox"){
        pseudo = 1e-5
        wt = wilcox.test(x1, x0, alternative = "greater")
        pvalue   = wt$p.value  
        lfc = log((m1 + pseudo) / (m0 + pseudo)) / log(2)
      } else if (test.type == "prop"){
        pseudo = 1e-2
        Y = t(matrix(c(nexp_cluster, n_cluster - nexp_cluster,
                     nexp_other, n_other - nexp_other), nrow=2))
        pt = prop.test(Y, alternative = "greater")
        pvalue = pt$p.value  
        lfc = log((pexp_cluster + pseudo) / (pexp_other + pseudo)) / log(2)
      }
      
      # Stats
      df[g, "test.type"] = test.type
      df[g, "lfc"] = lfc
      df[g, "pvalue"] = pvalue  
      df[g, "pexp_cluster"] = pexp_cluster
      df[g, "pexp_other"] = pexp_other
      df[g, "n_cluster"] = n_cluster 
      df[g, "n_other"] = n_other
      df[g, "nexp_cluster"] = nexp_cluster 
      df[g, "nexp_other"] = nexp_other
    }
    
    # Filter
    row.names(df) = NULL
    df$fdr = p.adjust(df$pvalue, method = "BH")
    df = df[(df$lfc >= lfc.thresh) & (df$fdr <= fdr.thresh),]
    df = df[rev(order(df$lfc)),]
    
    # Write temp to disk
    if(!is.null(temp.dir)){
      fname = file.path(temp.dir, sprintf("test_%s.csv", cont))
      write.csv(df, fname)  
    }
    message(sprintf("Cluster %s, found genes: %d, written %s", cont, nrow(df), fname))
    results = rbind(results, df)
    count = count + 1
  }
  return(results)
}

# Differential expression - MAST
st.diff.mast <- function(obj, lfc.thresh=0.6, fdr.thresh=0.05, temp.dir=NULL, marker.frame = NULL){
  
  # Convert to one-hot
  Mf = cluster.one.hot(obj)
  
  # Run differential expression
  sca <- FromMatrix(as.matrix(obj@data[obj@var.genes,]), 
                    Mf, data.frame(primerid=obj@var.genes))
  results = data.frame()
  for(cont in colnames(Mf)){
    df = mast.summary(sca, 
                      cont, 
                      lfc.thresh = lfc.thresh, 
                      fdr.thresh = fdr.thresh)
    cid = as.numeric(gsub("c", "", cont))
    df$cluster = cid

    # Mark which genes where used as markers and append groups
    if(!is.null(marker.frame)){
      df$marker = df$primerid %in% marker.frame$Name
      df$groups = ""
      ids = df$primerid[df$marker]
      groups = c()
      for(i in ids){
        groups = c(groups, paste(marker.frame[marker.frame$Name == i, ]$Cell.type, collapse = ";"))
      }
      df$groups[df$marker] = groups
    }
    # Write temp to disk
    if(!is.null(temp.dir)){
      fname = file.path(temp.dir, sprintf("test_%03d.csv", cid))
      write.csv(df, fname)  
    }
    message(sprintf("Cluster %d, found genes: %d, written %s", cid, nrow(df), fname))
    
    results = rbind(results, df)
  }
  return(results)
}


# Summarize differential expressions in a table
mast.summary <- function(sca, cont, lfc.thresh=0.6, fdr.thresh=0.05){
  zlmCond <- zlm(as.formula(sprintf("~%s", cont)), sca)
  cont = colnames(zlmCond@coefD)[2]
  summaryCond <- summary(zlmCond, doLRT=cont) 
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast == cont & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast == cont & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
                    by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle[fdr < fdr.thresh & abs(coef) > lfc.thresh], as.data.frame(mcols(sca)), by='primerid')
  fcHurdleSig = fcHurdleSig[order(fcHurdleSig$fdr),]
  return(fcHurdleSig)
}


# Number of shared entries between matrix X and Y (genes are in rows)
dist.shared.genes <- function(X, Y){
  Bx = as.matrix(X) > 0
  By = as.matrix(Y) > 0
  G = 1 * t(Bx) %*% By
  t = as.vector(G[upper.tri(G, diag=FALSE)])
  t
}

# Euclidean distances matrix X and Y (genes are in rows)
dist.euclidean <- function(X, Y){
  D = proxy::dist(as.matrix(X), as.matrix(Y), by_rows=FALSE)
  d = as.vector(D[upper.tri(D, diag = FALSE)])
  d
}



