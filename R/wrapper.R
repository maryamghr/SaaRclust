#' Wrapper function to run saarclust pipeline.
#'
#' @param inputfolder A folder name where minimap files are stored.
#' @param store.bestAlign Store best alignements in RData object.
#' @param HC.only Perform only hard clustering and skip the rest of the pipeline.
#' @param numAlignments Required number of best PBvsSS alignmnets to selest for hard clustering.
#' @param verbose Set to \code{TRUE} to print function messages.
#' @param HC.input Filaname where hard clustering results are stored
#' @param cellNum specifies the number of single cells to be used in clustering
#' @inheritParams SaaRclust
#' @inheritParams EMclust
#' @export
#' @author David Porubsky, Maryam Ghareghani


runSaaRclust <- function(inputfolder=NULL, outputfolder="SaaRclust_results", input_type="bam", 
                         num.clusters=54, EM.iter=100, alpha=0.01, minLib=10, upperQ=0.95, 
                         logL.th=1, theta.constrain=FALSE, hardclustMethod="hclust",
                         hard.theta.estim.method="median", hardclustMinLib = 35, hardclustLowerQ=0.7, 
                         add.garbage.cluster=FALSE, hardclustUpperQ=0.9, min.cluster.frequency=0.002, 
                         store.counts=TRUE, store.bestAlign=TRUE, numAlignments=30000, minSScov=30, 
                         HC.only=TRUE, verbose=TRUE, cellNum=NULL, log.scale=FALSE, 
                         hardclust.file="hard_clusters.RData", softclust.file="soft_clusters.RData",
                         ss.clust.file="ss_clusters.data",
                         ref.aln.bam=NULL, store.chrom.flag=TRUE, chrom.flag=NULL, numCPU=20) {
  
  #=========================#
  ### Create directiories ###
  #=========================#
  
  #Create a master output directory
  outputfolder.destination <- file.path(outputfolder)
  if (!file.exists( outputfolder.destination)) {
    dir.create( outputfolder.destination)
  }
  
  #Directory to store raw read counts and best alignments
  if (store.counts) {
    rawdata.store <- file.path(outputfolder.destination, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }
  
  if (store.bestAlign) {
    rawdata.store <- file.path(outputfolder.destination, 'RawData')
    if (!file.exists(rawdata.store)) {
      dir.create(rawdata.store)
    }
  }

  #Directory to store processed/clustered data
  Clusters.store <- file.path(outputfolder.destination, 'Clusters')
  if (!file.exists(Clusters.store)) {
    dir.create(Clusters.store)
  }
  
  #Directory to store plots
  plots.store <- file.path(outputfolder.destination, 'Plots')
  if (!file.exists(plots.store)) {
    dir.create(plots.store)
  }

  #Directory to store 'difficult' PacBio reads for later processing [TODO]
  trashbin.store <- file.path(outputfolder.destination, 'TrashBin')
  if (!file.exists(trashbin.store)) {
    dir.create(trashbin.store)
  }
  
  # getting Strand-seq read counts
  #reuse existing data if they were already created and saved in a given location
  destination1 <- file.path(rawdata.store, "alignments.RData")
  destination2 <- file.path(rawdata.store, "read_counts.RData")
  destination3 <- file.path(rawdata.store, "read_selected_counts.RData")
  
  if (file.exists(destination1)) {
    alignments <- get(load(destination1))
  }
  
  if (file.exists(destination2)) {
    counts.l.all <- get(load(destination2))
  }
  
  if (file.exists(destination3)) {
    counts.l <- get(load(destination3))
  } else  {
    if (input_type=="bam"){
      # counting w/c in bam files
      bam.files <- list.files(path=inputfolder, pattern='.mdup.bam$', full.names=TRUE)
      if (!is.null(cellNum)) {
        bam.files <- bam.files[1:cellNum]
      }
      
      cat('counting w/c reads...\n')
      counts <- count.wc.bam(bam.files, numCPU=numCPU)
      counts.l.all <- counts[[1]]
      alignments <- counts[[2]]
      # getting representative counts table (alignments with highest ss coverage)
      cat('getting representative alignments...\n')
      counts <- get_representative_counts(counts.l=counts.l.all, 
                                          num.alignments = numAlignments, 
                                          min.ss.cov=minSScov)
      counts.l.all <- counts[[1]]
      counts.l <- counts[[2]]
      
    } else { # input file is minimap
      ### Get representative alignments to estimate theta and pi values ###
      destination <- file.path(rawdata.store, "representativeAligns.RData")
      
      if (!file.exists(destination)) {
        # setting the min number of SS libs as a cutoff to use PB reads in hard clustering
        
        if (!is.null(cellNum)) {
          hardclustMinLib <- round(cellNum/4)
        }
        best.alignments <- getRepresentativeAlignments(inputfolder=inputfolder, 
                                                       numAlignments=numAlignments, 
                                                       quantileSSreads=c(hardclustLowerQ, hardclustUpperQ), 
                                                       minSSlibs=c(hardclustMinLib,Inf))
        if (store.bestAlign) {
          save(file = destination, best.alignments)
        }
      } else {
        best.alignments <- get(load(destination))
      }
      
      #use PB read names as factor in order to export counts for every PB read (also for zero counts)
      best.alignments$PBreadNames <- factor(best.alignments$PBreadNames, 
                                            levels=unique(best.alignments$PBreadNames))
      
      #split data by Strand-seq library
      tab.l <- split(best.alignments, best.alignments$SSlibNames)
      
      ### Count directional reads ###
      counts.l <- countDirectionalReads(tab.l)
    }
    
    if (store.counts) {
      if (exists("counts.l.all")) {
        save(file = destination1, alignments)
        save(file = destination2, counts.l.all)
      }
      save(file = destination3, counts.l)
    }
  }
  
  #subsetting single cell libraries
  if (!is.null(cellNum)) {
    if (exists("counts.l.all")) {
      counts.l.all <- counts.l.all[1:cellNum]
    }
    
    if(class(counts.l)[1]=='list') {
      counts.l <- counts.l[1:cellNum]
    } else {
      selected.libs <- head(counts.l[, unique(lib)], cellNum)
      counts.l <- counts.l[lib %in% selected.libs]
    }
  }
  
  if (!is.null(ref.aln.bam)){
    # getting chrom/flag information for long reads/unitigs
    destination <- file.path(outputfolder.destination, 'chrom_flags.data')
    
    if (file.exists(destination)) {
      chrom.flag <- fread(destination)
    } else {
      chrom.flag <- getChromFlag(ref.aln.bam)
      if (store.chrom.flag) {
        fwrite(chrom.flag, file=destination, sep='\t', quote=FALSE)
      }
    }
  }
  
  #Load Hard clustering results if they were already created
  destination <- file.path(Clusters.store, hardclust.file)
  if (file.exists(destination)) {
    message("Loading Hard clustering results")
    cat('destination =', destination)
    hard.clust <- get(load(destination))
    
  } else {
    message("Hard clustering results not available!!!")
    message("Running Hard clustering")

    ### Perform hard clustering ###
    hardClust.ord <- hardClust(counts.l=counts.l, method=hardclustMethod, 
                               num.clusters=num.clusters, 
                               min.cluster.frequency=min.cluster.frequency, 
                               chrom.flag=chrom.flag)
    
    theta.estim <- estimateTheta(counts.l, hard.clust=hardClust.ord, alpha=alpha, 
                                 method=hard.theta.estim.method)
    
    #Merge splitted clusters after hard clustering
    k <- 46
    if (add.garbage.cluster) {k <- 47}
    hardClust.ord.merged <- mergeClusters(hard.clust=hardClust.ord, theta.l=theta.estim, k=k)

    #Re-estimate theta parameter after cluster merging. 
    #Use the current theta as initial theta parameter for soft clustering.
    theta.param <- estimateTheta(counts.l, hard.clust=hardClust.ord.merged, alpha=alpha, 
                                 method=hard.theta.estim.method)
    
    # assert that all theta params per cell/cluster sum up to 1
    assertthat::assert_that(all(sapply(theta.param, function(theta) rowSums(theta)==1L))) %>% invisible()
    
    #Estimate pi parameter based on number of long reads in each cluster
    readsPerClusts <- table(hardClust.ord.merged)
    pi.param <- readsPerClusts/sum(readsPerClusts)
    
    #save hard clustering results into a file
    hard.clust <- list(ord=hardClust.ord.merged, unmerged.ord=hardClust.ord, 
                       theta.param=theta.param, pi.param=pi.param)
    if (!is.null(ref.aln.bam)){
      # adding ground true info to the clustering objects
      hard.clust <- addChromFlag(clust.obj=hard.clust, chrom.flag=chrom.flag, clust.type = 'hard')
    }
    
    destination <- file.path(Clusters.store, hardclust.file)
    save(file = destination, hard.clust)
    
  }
  
  clust.to.chrom <- NULL
  
  if (!is.null(chrom.flag)){
    print('hard clustering on merged clusters')
    clust.to.chrom <- numFoundClusters(ord=hard.clust$ord, chrom.flag)
  }
  clust.pairs <- findClusterPartners_simple(hard.clust$theta.param, clust.to.chrom)
  
  if(!HC.only) {
    destination <- file.path(Clusters.store, softclust.file)
    if (file.exists(destination)) {
      message("Loading Soft clustering results")
      cat('destination =', destination)
      soft.clust <- get(load(destination))
    } else {
      #Initialize theta parameter
      theta.param <- hard.clust$theta.param
      #Initialize pi parameter
      pi.param <- hard.clust$pi.param
      
      if (input_type=="bam"){
        read.names <- rownames(counts.l.all[[1]])
        counts.l.all <- lapply(counts.l.all, as.data.frame)
        soft.clust <- SaaRclust(counts.l=counts.l.all, outputfolder=outputfolder.destination, 
                                num.clusters=length(pi.param), EM.iter=EM.iter, alpha=alpha, 
                                minLib=minLib, upperQ=upperQ, theta.param=theta.param, pi.param=pi.param, 
                                logL.th=logL.th, theta.constrain=theta.constrain, log.scale=log.scale, 
                                HC.input=hard.clust, read.names=read.names, clust.pairs=clust.pairs,
                                numCPU=numCPU)
      } else {
        #List files to process
        file.list <- list.files(path = inputfolder, pattern = "chunk.+maf.gz$", full.names = TRUE)
      
        ### Main loop to process all files using EM algorithm ###
        for (file in file.list) {
          if (verbose) {
            soft.clust <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, 
                                    num.clusters=length(pi.param), EM.iter=EM.iter, alpha=alpha, 
                                    minLib=minLib, upperQ=upperQ, theta.param=theta.param, 
                                    pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain, 
                                    log.scale=log.scale, HC.input=hard.clust)
          } else {
            suppressMessages( soft.clust <- SaaRclust(minimap.file=file, outputfolder=outputfolder.destination, 
                                                      num.clusters=num.clusters, EM.iter=EM.iter, alpha=alpha, 
                                                      minLib=minLib, upperQ=upperQ, theta.param=theta.param, 
                                                      pi.param=pi.param, logL.th=logL.th, theta.constrain=theta.constrain, 
                                                      log.scale=log.scale, HC.input=hard.clust) )
          }
        }
      }
      if (!is.null(ref.aln.bam)){
        # adding ground true info to the clustering objects
        soft.clust <- addChromFlag(clust.obj=soft.clust, chrom.flag=chrom.flag, clust.type = 'soft')
      }
      
      destination <- file.path(Clusters.store, softclust.file)
      dt.destination <- file.path(Clusters.store, "MLclust.data")
      save(file = destination, soft.clust)
      orig.colnames <- colnames(soft.clust$ML.clust)
      colnames(soft.clust$ML.clust)[1] <- paste0('#', colnames(soft.clust$ML.clust)[1])
      fwrite(soft.clust$ML.clust, dt.destination, sep='\t', quote=F, row.names=F)
      colnames(soft.clust$ML.clust) <- orig.colnames
    }
    
    destination <- file.path(Clusters.store, 'clust_partners.txt')
    if (!file.exists(destination)) {
      fwrite(clust.pairs, file=destination, quote = F, row.names = F)
    }
    
    wc.cells.clusters <- get.wc.cells.clusters(soft.clust)
    destination <- file.path(Clusters.store, 'wc_cells_clusters.data')
    if (!file.exists(destination)) {
      fwrite(wc.cells.clusters, file=destination, quote = F, row.names = F)
    }
    
    if (!is.null(chrom.flag)) {
      # computing the accuracy of ML clustering resulting from EM soft clustering algorithm
      clust.to.chrom <- numFoundClusters(soft.clust$ML.clust, chrom.flag)
      destination <- file.path(Clusters.store, 'clust_to_chrom.data')
      fwrite(clust.to.chrom, file=destination, quote = F, row.names = F)
    }
    
    # clustering ss reads
    destination <- file.path(Clusters.store, ss.clust.file)
    if (file.exists(destination)) {
      message("Loading SS clustering results")
      cat('destination =', destination)
      ss.clust <- fread(destination)
    } else {
      ss.clust <- cluster.ss.reads(alignments, soft.clust$ML.clust, clust.pairs)
      orig.colnames <- colnames(ss.clust)
      colnames(ss.clust)[1] <- paste0('#', colnames(ss.clust)[1])
      ss.clust.sp <- split(ss.clust, by="chrom_clust")
      fwrite(ss.clust, file=destination, sep='\t', quote = F, row.names = F)
      lapply(names(ss.clust.sp), function(clust) fwrite(ss.clust.sp[[clust]], 
                                                        file=paste0(Clusters.store,"/ss_clusters_",clust,".data"),
                                                        sep='\t', quote = F, row.names = F)) %>% invisible()
    }
    
    output.list <- list(hard.clust=hard.clust, soft.clust=soft.clust, 
                        clust.pairs=clust.pairs, ss.clust=ss.clust)
    
    if (!is.null(clust.to.chrom)){
      output.list$clust.to.chrom <- clust.to.chrom
    }
    
    return(output.list)
  
  } else {
    return(hard.clust)
  }

}
