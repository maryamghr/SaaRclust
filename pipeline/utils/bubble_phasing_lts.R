<<<<<<< HEAD
log <- file(snakemake@log[[1]], open='wt')
sink(file=log, type='message')
sink(file=log, type='output')


=======
>>>>>>> 0cf722152bcf6a7cff658e785779d8cf437b54ad
## Load required libraries
library(dplyr)
library(tidyr)
library(reshape2)
<<<<<<< HEAD
library(data.table)

## FUNCTIONS ##
###############

output_phased_strand_states <- function(bubble.cov.files, clust.pairs, select.libs, output.file){
	phased.strand.states = data.table()

	for (i in 1:nrow(clust.pairs)) {
		cluster.pair <- clust.pairs[i,]
		cluster1 <- as.character(cluster.pair$clust.forward)
		cluster2 <- as.character(cluster.pair$clust.backward)
		
		clust1.file.idx <- grep(paste0("_cluster", cluster1, "_"), bubble.cov.files)
		clust2.file.idx <- grep(paste0("_cluster", cluster2, "_"), bubble.cov.files)

		cluster1.file <- bubble.cov.files[clust1.file.idx]
		cluster2.file <- bubble.cov.files[clust2.file.idx]

		phased.strand.states = rbind(phased.strand.states, phase_strand_states(cluster1.file, cluster2.file, cluster1, cluster2, select.libs))
	}

	fwrite(phased.strand.states, file=output.file, row.names=F, sep='\t')
}

phase_strand_states <- function(cluster1.file, cluster2.file, cluster1, cluster2, select.libs) {
=======

## Get WC regions ##
#WC.regions <- read.table("/home/porubsky/WORK/Phase_bubbles/LTS_results/wc_lib_clust.data", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
WC.regions <- read.table("/home/porubsky/WORK/Phase_bubbles/LTS_results/adjusted_bubble/wc_cells_clusters.data", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
## Get selected library names
select.libs <- unique(WC.regions$lib)
#select.libs <- gsub(select.libs, pattern = "\\.", replacement = "-")

## Get cluster partners ##
clust.pairs <- read.table("/home/porubsky/WORK/Phase_bubbles/LTS_results/clust_partners.txt", header = TRUE)
## Keep only unique pairs
forw <- as.numeric(gsub("[^\\d]+", "", clust.pairs$clust.forward, perl=TRUE))
backw <- as.numeric(gsub("[^\\d]+", "", clust.pairs$clust.backward, perl=TRUE))
pairs <- cbind(pmin(forw, backw), pmax(forw, backw))
clust.pairs <- clust.pairs[!duplicated(pairs),]

## Go over all cluster pairs and perform phasing ## (test: i <- 6; cluster pair V7 & V10)
for (i in 1:nrow(clust.pairs)) {
  cluster.pair <- clust.pairs[i,]
  #cluster1.file <- paste0("/home/porubsky/WORK/Phase_bubbles/LTS_results/bubble_SSlib_cov_cluster", cluster.pair$clust.forward,"_snv_bubbles_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00514.data")
  #cluster2.file <- paste0("/home/porubsky/WORK/Phase_bubbles/LTS_results/bubble_SSlib_cov_cluster", cluster.pair$clust.backward,"_snv_bubbles_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00514.data")
  cluster1.file <- paste0("/home/porubsky/WORK/Phase_bubbles/LTS_results/adjusted_bubble/adjusted_bubble_SSlib_cov_cluster", cluster.pair$clust.forward,"_snv_bubbles_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00514.data")
  cluster2.file <- paste0("/home/porubsky/WORK/Phase_bubbles/LTS_results/adjusted_bubble/adjusted_bubble_SSlib_cov_cluster", cluster.pair$clust.backward,"_snv_bubbles_k63_a3_l23_kminimap15_w1_f0.1_z500_HG00514.data")
>>>>>>> 0cf722152bcf6a7cff658e785779d8cf437b54ad
  ## Load bubbles into matrices
  matrices <- loadClusterBubbles(cluster1.file = cluster1.file,
                                 cluster2.file = cluster2.file)
  
  ## Sort matrices
  srt.matrices <- sortClusterBubbles(matrices = matrices, select.libs = select.libs)
<<<<<<< HEAD

  ## adding the haplo-phased strand states to the 

  clust1.haplo.strand.states <- data.table(lib=srt.matrices$cluster1.libs, cluster=cluster1)
  clust1.haplo.strand.states[,`:=`(lib=lapply(lib, function(x) strsplit(x, "__")[[1]][1]), haplotype=sapply(lib, function(x) as.integer(strsplit(x, "__C")[[1]][2])-1))]
  haplo.strand.states <- clust1.haplo.strand.states

  clust2.haplo.strand.states = clust1.haplo.strand.states
  clust2.haplo.strand.states[, `:=`(cluster=cluster2, haplotype=1-haplotype)]
  haplo.strand.states <- rbind(haplo.strand.states, clust2.haplo.strand.states)
  
=======
>>>>>>> 0cf722152bcf6a7cff658e785779d8cf437b54ad
 
  ## Get consensus
  h1.cons <- exportConsensus(data.matrix = srt.matrices$cluster1.m)
  h2.cons <- exportConsensus(data.matrix = srt.matrices$cluster2.m)
<<<<<<< HEAD

  return(haplo.strand.states)
}

=======
  
}


## FUNCTIONS ##
###############

>>>>>>> 0cf722152bcf6a7cff658e785779d8cf437b54ad
loadClusterBubbles <- function(cluster1.file = NULL, cluster2.file = NULL) {
    
  cluster1 <- read.table(cluster1.file, stringsAsFactors = FALSE, header=TRUE)
  cluster1.libs <- colnames(cluster1)[-1]
  #cluster1.libs <- gsub(cluster1.libs, pattern = "\\.|_", replacement = "-")
  cluster1.bubble.id <- cluster1[,1]
  cluster1 <- cluster1[,-1]
  # Switch missing values to 0 and reference alleles to 1 and alternative alleles to 2
  cluster1.m <- t(apply(cluster1, 2, function(x) recode(x, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)))
  attr(cluster1.m, 'dimnames') <- NULL
    
  cluster2 <- read.table(cluster2.file, stringsAsFactors = FALSE, header=TRUE)
  cluster2.libs <- colnames(cluster2)[-1]
  #cluster2.libs <- gsub(cluster2.libs, pattern = "\\.|_", replacement = "-")
  cluster2.bubble.id <- cluster2[,1]
  cluster2 <- cluster2[,-1]
  # Switch missing values to 0 and reference alleles to 1 and alternative alleles to 2
  cluster2.m <- t(apply(cluster2, 2, function(x) recode(x, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)))
  attr(cluster2.m, 'dimnames') <- NULL
    
  final.object <- list(cluster1.m = cluster1.m, cluster1.libs = cluster1.libs, cluster1.bubble.id = cluster1.bubble.id,
                       cluster2.m = cluster2.m, cluster2.libs = cluster2.libs, cluster2.bubble.id = cluster2.bubble.id)
  return(final.object)
}  
  
sortClusterBubbles <- function(matrices = NULL, select.libs = NULL) {
  ## Load data
  cluster1.m <- matrices[['cluster1.m']]
  cluster1.libs <- matrices[['cluster1.libs']]
  cluster1.bubble.id <- matrices[['cluster1.bubble.id']]
  cluster2.m <- matrices[['cluster2.m']]
  cluster2.libs <- matrices[['cluster2.libs']]
  cluster2.bubble.id <- matrices[['cluster2.bubble.id']]
  
  ## Filter and sort shared bubbles IDs
  shared.bubble.id <- intersect(cluster1.bubble.id, cluster2.bubble.id)
  shared.bubble.id <- sort(shared.bubble.id)
  cluster1.m <- cluster1.m[,match(shared.bubble.id, cluster1.bubble.id)]
  cluster1.bubble.id <- shared.bubble.id
  cluster2.m <- cluster2.m[,match(shared.bubble.id, cluster2.bubble.id)]
  cluster2.bubble.id <- shared.bubble.id
  
  ## Filter shared cells
  shared.libs <- intersect(cluster1.libs, cluster2.libs)
  ## NOTE library names are not consistent: sometimes they include dot and sometimes don't
  if (!is.null(select.libs)) {
    shared.libs <- intersect(shared.libs, select.libs)
  }
  ## Select rows present in shared cells
  mask.rows.cluster1 <- match(shared.libs, cluster1.libs)
  mask.rows.cluster2 <- match(shared.libs, cluster2.libs)
  cluster1.m <- cluster1.m[mask.rows.cluster1,]
  cluster2.m <- cluster2.m[mask.rows.cluster2,]
  cluster1.libs <- cluster1.libs[mask.rows.cluster1]
  cluster2.libs <- cluster2.libs[mask.rows.cluster2]
  cluster1.libs <- paste(cluster1.libs, "C1", sep="__")
  cluster2.libs <- paste(cluster2.libs, "C2", sep="__")
  
  ## Filter homozygous SNV positions
  mask.cl1 <- apply(cluster1.m, 2, function(x) any(x == 3))
  mask.cl2 <- apply(cluster2.m, 2, function(x) any(x == 3))
  mask <- c(which(mask.cl1 == TRUE), which(mask.cl2 == TRUE))
  mask <- unique(mask)
  cluster1.m <- cluster1.m[,-mask]
  cluster2.m <- cluster2.m[,-mask]
  
  ## sort parallel matrices
  for ( i in 1:length(shared.libs)) {
    filename <- shared.libs[i]
    message("Processing ", filename, " ...")
    
    cluster1.pos <- which(cluster1.m[i,] > 0)
    cluster2.pos <- which(cluster2.m[i,] > 0)
    cov.pos <- union(cluster1.pos, cluster2.pos)
    
    ## initialize score
    cluster1.m.score <- calcMatrixScore(matrix = cluster1.m, covered.pos = cov.pos)
    cluster2.m.score <- calcMatrixScore(matrix = cluster2.m, covered.pos = cov.pos)
    
    ## swap rows in matrices
    cluster1.filename <- cluster1.libs[i]
    cluster2.filename <- cluster2.libs[i]
    cluster1.libs[i] <- cluster2.filename
    cluster2.libs[i] <- cluster1.filename
    
    cluster1.row <- cluster1.m[i,]
    cluster2.row <- cluster2.m[i,]
    cluster1.m[i,] <- cluster2.row
    cluster2.m[i,] <- cluster1.row
    
    ## calculate new score
    curr.cluster1.m.score <- calcMatrixScore(matrix = cluster1.m, covered.pos = cov.pos)
    curr.cluster2.m.score <-  calcMatrixScore(matrix = cluster2.m, covered.pos = cov.pos)
    
    ## compare previous matrix score with score after swapping rows
    if ( (cluster1.m.score + cluster2.m.score) < (curr.cluster1.m.score + curr.cluster2.m.score) ) {
      cluster1.libs[i] <- cluster1.filename
      cluster2.libs[i] <- cluster2.filename
      
      cluster1.m[i,] <- cluster1.row
      cluster2.m[i,] <- cluster2.row
    }
  }
  
  ## Export sorted data
  sorted.matrices <- list()
  sorted.matrices[['cluster1.m']] <- cluster1.m
  sorted.matrices[['cluster1.libs']] <- cluster1.libs
  sorted.matrices[['cluster2.m']] <- cluster2.m
  sorted.matrices[['cluster2.libs']] <- cluster2.libs
  sorted.matrices[['cluster.bubble.id']] <- shared.bubble.id
  return(sorted.matrices)
}

calcMatrixScore <- function(matrix, covered.pos) {
  sub.m <- matrix[,covered.pos]
  sub.m <- sub.m[,colSums(sub.m) > 0]
  sub.m.long <- melt(sub.m)
  base.freq <- sub.m.long %>% group_by(Var2, value) %>% summarise(count=n())
  base.freq <- base.freq[base.freq$value > 0,]
  
  scores <- base.freq %>% group_by(Var2) %>% summarise(count = sum(count) - max(count))

  return(sum(scores$count))
}

exportConsensus <- function(data.matrix) {
  #base.freq <- apply(data.bases, 2, function(x) table(x[x>0])) ## TODO avoid warnings here
  indices <- which(data.matrix != 0,arr.ind = TRUE) # get indices of non-zero values
  values <- data.matrix[indices] # get non-zero values based on indices
  col.vals <- split(values, (indices[,2])) # split values based on column indices
  
  base.freq <- lapply(col.vals, table) # get base frequencies for each column
  
  positions <- as.numeric(names(base.freq)) # get position of each snv/column
  
  ## Filter SNV positions (columns) which do not pass set criteria (coverge, ambiguity)
  max.cov <- sapply(base.freq, function(x) max(x)) # get coverage of highest covered base
  score <- sapply(base.freq, function(x) sum(x)-max(x))
  
  if (length(positions) != 0) {
    
    ## Calculate entropy values for each column
    entropy <- c(rep(0, length(positions)))
    entropy <- sapply(col.vals, calcEnt)
    
    ## Calculate Phred score values for each column
    scores <- c(rep(0, length(positions)))
    
    cov = sapply(col.vals, length)
    
    hap <- sapply(base.freq, function(x) as.numeric(names(which.max(x))))
    
    assem.haps <- data.frame(pos=positions, hap=hap, cov=cov, score=scores, ent=entropy)
    rownames(assem.haps) <- NULL
    return(assem.haps)
  } else {
    return(0)
  }  
}

calcEnt <- function(bases) {
  
  # log2 function (modif: return zero for zero values)
  getlog2 <- function(x) {
    if (x == 0 ) { 
      return(0)
    } else {  
      return(log2(x))
    }  
  }
  
  counts <- c(0,0,0,0)
  
  # count frequency of each base in a column 
  for (i in 1:4) {
    counts[i] <- length(bases[bases == i])
  }
  probs <- counts/sum(counts) # get probabilities
  
  ent <- -( probs[1]*getlog2(probs[1]) + probs[2]*getlog2(probs[2]) + probs[3]*getlog2(probs[3]) + probs[4]*getlog2(probs[4]) )
  return(ent)    
}
<<<<<<< HEAD

=======
>>>>>>> 0cf722152bcf6a7cff658e785779d8cf437b54ad
