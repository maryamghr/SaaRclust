#' Hard clustering using k-means
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export

hardClust <- function(counts.l=NULL, num.clusters=NULL, nstart=10, iter.max=10, by_chrom=F) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Kmeans clustering for ",num.clusters," clusters") 
  
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
   
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.na(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
   
  ratios.m <- do.call(cbind, ratios.l)
  
  if (by_chrom){
    ratios.m <- abs(ratios.m)
  }

  #hard clustering using kmeans
  km <- suppressWarnings( kmeans(ratios.m, centers = num.clusters, nstart = nstart, iter.max = iter.max) )
  ord <- km$cluster
  names(ord) <- rownames(counts.l[[1]])
  #ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)
  return(ord)
}

#' Getting wfrac matrix
#'
#' @param counts.dt A \code{data.table} of w/c read counts per lib/long read with the column names rname, w, c, lib.name.
#' @inheritParams SaaRclust
#' @return A \code{matrix} of w.frac values (w-c)/(w+c) with rownames=long read names and column names=single-cell library names
#' @author Maryam Ghareghani
#' @export
#' 
get.wfrac.matrix <- function(counts.dt){
  w.frac <- counts.dt
  
  w.frac[, W.frac:=(w-c)/(w+c)]
  #w.frac[, `:=`(w=NULL, c=NULL)]
  
  mat <- dcast(w.frac, rname~lib, value.var='W.frac')
  
  row.names <- mat[, rname]
  mat[, rname:=NULL]
  mat <- as.matrix(mat)
  rownames(mat) <- row.names
  
  mat
}

#' Hard clustering using hclust
#'
#' @param wfrac.matrix A \code{matrix} of of w.frac values (w-c)/(w+c) with rownames=long read names and column names=single-cell library names
#' @param num.clusters Number of output clusters
#' @param by_chrom a logical value (default FALSE) indicating whether to cluster only by chromosome
#' @return A numeric vector including the clusters, names with object (long read/unitig) names
#' @author Maryam Ghareghani
#' @export

hard.hclust <- function(wfrac.matrix=NULL, num.clusters=NULL, by_chrom=F) {
  
  message("Hard clustering")
  ptm <- startTimedMessage("    Hierarchical clustering for ",num.clusters," clusters") 
  
  mat <- wfrac.matrix
  if (by_chrom){
    mat <- abs(wfrac.matrix)
  }
  
  d=dist(mat)
  
  hc = hclust(d)
  hc.clust <- cutree(hc, k=num.clusters)
  ord = hc.clust
  
  stopTimedMessage(ptm)
  return(ord)
}
  
#' Estimate theta values based on hard clustering
#'
#' This function takes results of hard clustering and estimates majority cell types for each Strand-seq library
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @param hard.clust A \code{integer} of cluster assignments for each PacBio read. 
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export

# 
estimateTheta <- function(counts.l=NULL, hard.clust=NULL, alpha=0.1) {
  ptm <- startTimedMessage("Estimate theta values") 
  theta.estim <- list()
  for (j in 1:length(counts.l)) {
    minus.c <- split(counts.l[[j]][,1], hard.clust)
    plus.c <- split(counts.l[[j]][,2], hard.clust)
    
    # adjustment for data.table count format
    if (length(minus.c)>0 & "data.table" %in% class(minus.c[[1]])){
      minus.c <- lapply(minus.c, function(x) x[,w])
      plus.c <- lapply(plus.c, function(x) x[,c])
    }
    #minus.counts <- sapply(minus.c, sum)
    #plus.counts <- sapply(plus.c, sum)
    #probs <- countProb2(minusCounts = minus.counts, plusCounts = plus.counts)
  
    clust.prob <- mapply(function(X,Y) { countProb(X,Y) }, X=minus.c, Y=plus.c)
    clust.prob.norm <- lapply(clust.prob, function(x) colSums(log(x)))
    estimates <- sapply(clust.prob.norm, which.max)
    # FIXME: theta_wc is always the most likely theta!!!
    
    #Assign cell type probs based on the majority type in each cluster
    probs <- list()
    for (i in 1:length(clust.prob)) {
      estim <- estimates[i]
      if (estim == 1) {
        theta <- c(1-alpha, alpha/2, alpha/2)
      } else if (estim == 2) {
        theta <- c(alpha/2, 1-alpha, alpha/2)
      } else {
        theta <- c(alpha/2, alpha/2, 1-alpha)
      }
      probs[[i]] <- theta
    }
    probs <- do.call(rbind, probs)
    theta.estim[[j]] <- probs
    names(theta.estim)[j] <- names(counts.l)[j]
  }
  stopTimedMessage(ptm)
  
  return(theta.estim)
}

                             
#' Hierarchical clustering for merging the kmeans clusters.
#'
#' This function takes as input the kmeans hard clustering output and the initialized thetas and merges the kmeans clusters based on thetas
#'
#' @param theta.l A \code{list} of estimated theta values for each cluster and cell.
#' @param hard.clust The kmeans hard clustering.
#' @param k Desired number of clusters after merging.
#' @param row.clusters the cluster ids for matrix rows to sort the rows by cluster.
#' @inheritParams SaaRclust
#' @return A new hard clustering with the correct number of clusters
#' @author Maryam Ghareghani
#' @export

mergeClusters <- function(hard.clust, theta.l, k=46, row.clusters=NULL)
{ # order the matrix not tested yet
  ptm <- startTimedMessage("Merging clusters")
  theta.all <- do.call(cbind, theta.l)
  
  d <- dist(theta.all)
  hc <- hclust(d)
  hc.clust <- cutree(hc, k=k)
  
  hc.clust.ord <- NULL
  
  for (i in 1:length(hc.clust)){
    hc.clust.ord[row.clusters[i]] <- row.clusters[hc[i]]
  }
  
  stopTimedMessage(ptm)
  return(sapply(hard.clust, function(i) hc.clust.ord[i]))
}
