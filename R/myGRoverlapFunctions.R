#' Get a dataframe from findOverlaps from a GRangesList with 2 items where colnames match the names of the items
#'
#' @param my2GRsToOverlap a GRangesList with 2 items to overlap
#' @return a dataframe with 2 columns with colnames identicals to names of the input with the indices of the overlaps
#' @importFrom GenomicRanges findOverlaps
#' @importFrom BiocGenerics as.data.frame
#' @export
#' @examples
#'gr1 <- GenomicRanges::GRanges(
#'seqnames = "chr1",
#'ranges = IRanges::IRanges(start = c(1, 11, 199),
#'                          end = c(10, 12, 200)),
#'score = c(20, 30, 100))
#'gr2 <- GenomicRanges::GRanges(
#'seqnames = "chr1",
#'ranges = IRanges::IRanges(start = c(5, 101),
#'                          end = c(15, 102)),
#'score = c(1, 2))
#'overlap2GR(list(first = gr1, second = gr2))
overlap2GR <- function(my2GRsToOverlap){
  if (length(my2GRsToOverlap) != 2){
    stop("This is not a list with 2 items.")
  }
  mAll <- BiocGenerics::as.data.frame(
    as.matrix(
      GenomicRanges::findOverlaps(my2GRsToOverlap[[1]],
                                  my2GRsToOverlap[[2]])))
  colnames(mAll) <- names(my2GRsToOverlap)
  return(mAll)
}

### THIS IS VERY TIME CONSUMING IT WOULD GO FASTER WITH GRAPHS...

#' Get a dataframe from all possible pairwise overlaps
#' from a GRangesList where multiple overlaps are on the same line when possible.
#' Or a list of dataframes of any combination using any items in the list
#'
#' @param myGRsToOverlap a GRangesList with items to overlap
#' @param allMergedCalculated a List with previously calculated overlaps. Default is NULL.
#' @param verbose logical Default is F.
#' @param returnAllMergedCalculated logical to decide if only the final dataframe with all the merges are given or if the list with all the overlaps should be returned. Default is F.
#' @return Either a dataframe with one column per item in the list with colnames identicals to names of the input with the indices of the overlaps
#' Or a list with this this dataframe + all other which have been calculated in the process.
#' @details This function is really slow because it will make all pair, triple... comparison possible.
#' Because when you merge 3 dataframes with 2 columns with in total 3 column names, you have 3 ways to do the merge:
#' merge 1 with the merge of 2 and 3 or merge 2 with the merge of 1 and 3 or merge 3 with the merge of 1 and 2.
#' The algorithm will do all and report all.
#' This function is recursive.
#' @importFrom combinat combn permn
#' @importFrom usefulLDfunctions subsetByNamesOrIndices
#' @export
#' @examples
#'gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
#'ranges = IRanges::IRanges(start = c(1, 11, 199),
#'                          end = c(10, 12, 200)),
#'score = c(20, 30, 100))
#'gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(5, 101),
#'                                                        end = c(15, 102)),
#'                              score = c(1, 2))
#'gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(1, 100),
#'                                                        end = c(4, 120)),
#'                              score = c(10, 20))
#'grL <- list(first = gr1, second = gr2, third = gr3)
#'findAllOverlap(grL)
#'findAllOverlap(grL, returnAllMergedCalculated = TRUE)
findAllOverlap <- function(myGRsToOverlap,
                           allMergedCalculated = NULL,
                           verbose = F,
                           returnAllMergedCalculated = F){
  # Uses overlap2GR from the same package
  if (verbose){
    cat("Evaluating the overlap of ")
    cat(names(myGRsToOverlap))
    cat("\n")
  }
  n <- length(myGRsToOverlap)
  if (is.null(allMergedCalculated)){
    #This is to avoid to calculate twice the same merge
    allMergedCalculated <- list()
  }
  if (n == 2){
    if (verbose){
      cat("Do a merge\n")
    }
    mAll <- usefulLDfunctionsGR::overlap2GR(myGRsToOverlap)
    if (returnAllMergedCalculated){
      nameOfN <- paste(sort(names(myGRsToOverlap)), collapse = "____")
      allMergedCalculated[[nameOfN]] <- mAll[, names(myGRsToOverlap)]
      return(allMergedCalculated)
    } else {
      return(mAll[, names(myGRsToOverlap)])
    }
  }
  # Will contain all the merge of all except one the (n-1)s merges
  mNm1 <- list()
  allNm1 <- combinat::combn(names(myGRsToOverlap), n - 1)
  for (i in 1:ncol(allNm1)){
    nameOfNm1 <- paste(sort(allNm1[, i]), collapse = "____")
    if (!nameOfNm1 %in% names(allMergedCalculated)){
      allMergedCalculated <- findAllOverlap(
        subsetByNamesOrIndices(myGRsToOverlap,
                               sort(allNm1[, i])),
        allMergedCalculated,
        verbose = verbose,
        returnAllMergedCalculated = T)
      if (verbose){
        cat("back to ")
        cat(names(myGRsToOverlap))
        cat("\n")
      }
    }
    mNm1[[nameOfNm1]] <- allMergedCalculated[[nameOfNm1]]
  }
  mAll <- NULL
  # The permutations will be use to decide the order of the merge.
  allPerm <- matrix(unlist(combinat::permn(names(mNm1))),
                    byrow = T, ncol = n)
  allPerm <- allPerm[allPerm[, 1] < allPerm[, 2], ]
  for (i in 1:nrow(allPerm)){
    if (verbose){
      cat("perm", i, ",")
    }
    v <- allPerm[i, ]
    temp.m <- mNm1[[v[1]]]
    for (j in 2:length(v)){
      temp.m <- merge(temp.m, mNm1[[v[j]]], all = T)
    }
    mAll <- unique(rbind(mAll, temp.m))
  }
  if (verbose){
    cat("\n")
  }
  if (returnAllMergedCalculated){
    nameOfN <- paste(sort(names(myGRsToOverlap)), collapse = "____")
    allMergedCalculated[[nameOfN]] <- mAll[, names(myGRsToOverlap)]
    return(allMergedCalculated)
  } else {
    return(mAll[, names(myGRsToOverlap)])
  }
}

#' Find the overlaps of GRangesList items to have only once each indice of each GRanges
#' maximizing the sum of scores or normed scores.
#'
#' @param myGRsToOverlap a GRangesList with items to overlap
#' @param mAll a dataframe obtained with findAllOverlap with all possible overlaps between items of GRangeList (default is NULL).
#' @param colWithScore the name of the column in the meta columns of the GRanges of the GRList
#' to use to choose the best overlap (default is `"score"``).
#' @param useNormScore a boolean to specify if the score of each GRanges should be normalized
#' so the sum of all scores is 1 before summing between the overlaps to maximize.
#' @param verbose logical (default is `FALSE`).
#' @return a dataframe with column names which correspond to the names of the GRanges items of the GRList.
#' Each indice of each GRanges will be present. When indices from different columns are on the same row,
#' this means that the specified GRanges overlap (at least with another one).
#' @importFrom GenomicRanges mcols
#' @export
#' @examples
#'gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
#'ranges = IRanges::IRanges(start = c(1, 11, 199),
#'                          end = c(10, 12, 200)),
#'score = c(20, 30, 100))
#'gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(5, 101),
#'                                                        end = c(15, 102)),
#'                              score = c(60, 90))
#'gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(1, 100),
#'                                                        end = c(4, 120)),
#'                              score = c(1, 2))
#'grL <- list(first = gr1, second = gr2, third = gr3)
#'filterByScore(grL, useNormScore = TRUE)
#'filterByScore(grL, useNormScore = FALSE)
filterByScore <- function(myGRsToOverlap,
                          mAll = NULL,
                          colWithScore = "score",
                          useNormScore = T,
                          verbose = F){
  if (is.null(mAll)){
    mAll <- usefulLDfunctionsGR::findAllOverlap(myGRsToOverlap,
                                                verbose = verbose)
  }
  n <- ncol(mAll)
  mAll$score <- 0
  for (i in 1:n){
    nameI <- colnames(mAll)[i]
    # Getting the indices where the value is not NA
    nonNAI <- !is.na(mAll[, nameI])
    # Getting the corresponding scores
    scoresToAdd <- mcols(myGRsToOverlap[[nameI]])[mAll[nonNAI, nameI],
                                                  colWithScore]
    # If useNormScore is TRUE, the scores are divided by
    # the sum of all scores of the GRange
    if (useNormScore){
      scoresToAdd <- scoresToAdd /
        sum(mcols(myGRsToOverlap[[nameI]])[, colWithScore])
    }
    mAll$score[nonNAI] <- mAll$score[nonNAI] + scoresToAdd
  }
  # All overlaps are sorted by decreasing scores
  mAllS <- mAll[order(mAll$score, decreasing = T), ]
  # All duplicated values for each of the GRanges are removed.
  # This step is biased as the order may have a high influence.
  mAllSF <- mAllS
  for (i in 1:n){
    nameI <- colnames(mAll)[i]
    mAllSF <- mAllSF[is.na(mAllSF[, nameI]) | !duplicated(mAllSF[, nameI]), ]
  }
  # To compensate this bias, we will try to "rescue"
  # the overlaps which have been deleted but exists.
  mAllSRescue <- mAllS
  # We put in mAllSRescue all overlaps which are not using
  # indices already present in the mAllSF
  for (i in 1:n){
    nameI <- colnames(mAll)[i]
    mAllSRescue <- mAllSRescue[is.na(mAllSRescue[, nameI]) |
                      !mAllSRescue[, nameI] %in% mAllSF[, nameI], ]
  }
  # We remove the duplicated indices in the mAllSRescue
  for (i in 1:n){
    nameI <- colnames(mAll)[i]
    mAllSRescue <- mAllSRescue[is.na(mAllSRescue[, nameI]) |
                        !duplicated(mAllSRescue[, nameI]), ]
  }
  # While there are still overlaps which could be added.
  # We redo the process.
  while (nrow(mAllSRescue) > 0){
    mAllSF <- rbind(mAllSF, mAllSRescue)
    mAllSRescue <- mAllS
    # We put in mAllSRescue all overlaps which are not using
    # indices already present in the mAllSF
    for (i in 1:n){
      nameI <- colnames(mAll)[i]
      mAllSRescue <- mAllSRescue[is.na(mAllSRescue[, nameI]) |
                        !mAllSRescue[, nameI] %in% mAllSF[, nameI], ]
    }
    # We remove the duplicated indices in the mAllSRescue
    for (i in 1:n){
      nameI <- colnames(mAll)[i]
      mAllSRescue <- mAllSRescue[is.na(mAllSRescue[, nameI]) |
                          !duplicated(mAllSRescue[, nameI]), ]
    }
  }
  # We now need to add the indices that are not present.
  mAllF <- mAllSF[, 1:n]
  for (i in 1:n){
    nameI <- colnames(mAll)[i]
    temp.df <- data.frame(1:length(myGRsToOverlap[[nameI]]))
    colnames(temp.df) <- nameI
    mAllF <- merge(mAllF, temp.df, all = T)
  }
  return(mAllF[, setdiff(colnames(mAll), "score")])
}

#' Get a dataframe summarizing the clusters created by overlaps
#' from a GRangesList where each element is only in one cluster
#'
#' @param myGRsToOverlap a GRangesList with items to overlap
#' @return A dataframe with one column per item in the list with colnames identicals to names of the input
#' each row is a cluster of overlapping items
#' each value is a vector with the indices of the items in the cluster.
#' @importFrom GenomicRanges findOverlaps
#' @importFrom igraph graph_from_edgelist vcount add_vertices components
#' @export
#' @examples
#'gr1 <- GenomicRanges::GRanges(seqnames = "chr1",
#'ranges = IRanges::IRanges(start = c(1, 11, 199),
#'                          end = c(10, 12, 200)),
#'score = c(20, 30, 100))
#'gr2 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(5, 101),
#'                                                        end = c(15, 102)),
#'                              score = c(1, 2))
#'gr3 <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(1, 100),
#'                                                        end = c(4, 120)),
#'                              score = c(10, 20))
#'grL <- list(first = gr1, second = gr2, third = gr3)
#'findOverlapAsClusters(grL)
findOverlapAsClusters <- function(myGRsToOverlap){
  # From Seurat:
  # determine pairwise combinations
  combinations <- expand.grid(1:(length(x = myGRsToOverlap)),
                              1:(length(x = myGRsToOverlap)))
  combinations <- combinations[combinations$Var1 < combinations$Var2,
                               , drop = FALSE]
  # determine the proper offsets
  nitems <- sapply(X = myGRsToOverlap, FUN = length)
  breaksItems <- as.vector(x = cumsum(x = c(0, nitems)))
  offsets <- breaksItems[1:length(x = nitems)]
  # We are now building a network
  # First define edges (overlaps)
  all.overlaps <- lapply(1:nrow(combinations), function(row){
    i <- combinations[row, 1]
    j <- combinations[row, 2]
    ov <- as.data.frame(GenomicRanges::findOverlaps(myGRsToOverlap[[i]],
                                                    myGRsToOverlap[[j]]))
    ov[, 1] <- ov[, 1] + offsets[i]
    ov[, 2] <- ov[, 2] + offsets[j]
    return(ov)
  })
  all.overlaps <- do.call(rbind, all.overlaps)
  # Building the network
  network <- igraph::graph_from_edgelist(as.matrix(all.overlaps),
                                         directed = FALSE)
  nvertices.to.add <- sum(nitems) - igraph::vcount(network)
  network <- igraph::add_vertices(network,
                                  nv = nvertices.to.add)
  # Cut the network in connected clusters
  cluster <- igraph::components(network)
  # Make a list where each element is a vector
  # with all items of the cluster
  clusterList <- aggregate(list(id = 1:sum(nitems)),
                           by = list(clusterid = cluster$membership),
                           FUN = c)$id
  # Split the members of the cluster per gr origin
  clusterListByGR <- sapply(clusterList, function(v){
    split(x = v,
          f = cut(x = v,
                  breaks = breaksItems,
                  labels = names(myGRsToOverlap)))
  })

  clusterListWOoffsets <- lapply(1:length(offsets), function(i){
    return(sapply(clusterListByGR[i, ], "-", offsets[i]))
  })

  # Convert it to data.frame
  clusterListAsDF <- as.data.frame(do.call(cbind,
                                           clusterListWOoffsets))
  colnames(clusterListAsDF) <- names(myGRsToOverlap)

  return(clusterListAsDF)
}
