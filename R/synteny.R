#' Find Synteny Blocks
#'
#' Creates synteny blocks from homeolog positions.
#' See the vignette for more information.
#'
#' @param df Data frame with positions.
#' @param eps Numeric, eps parameter for dbscan
#' @param minPts Integer, minPTS parameter for dbscan.
#' @param minLength Integer, minimal block length.
#' Blocks below this threshold will be filtered out.
#' @param maxLength Integer, maximal block length.
#' Blocks above this threshold will be split up.
#' @param gap Numeric, gap size between chromosomes.
#' @export
find_blocks <-
  function(df, eps = 2e+06L, minPts = 100L, minLength = 1e+06L, maxLength = 1e+07L,
           gap = 0L) {
    df$glo1 <- df$pos1
    df$glo2 <- df$pos2
    chrs <- sort(unique(df$chr1))
    start <- end <- numeric(length(chrs))
    genome <- rep(1, length(chrs))
    offset <- 0
    for (i in 1:length(chrs)) {
      if (i > 1)
        offset <- gap + end[i - 1]
      chr1match <- df$chr1 == chrs[i]
      start[i] <- min(df$pos1[chr1match]) + offset
      end[i] <- max(df$pos1[chr1match]) + offset
      df$glo1[chr1match] <- df$glo1[chr1match] + offset
    }
    out <-
      synteny_info$new(chrs = data.frame(
        chr = chrs, start = start, end = end, genome = genome
      ))
    chrs <- sort(unique(df$chr2))
    genome <- rep(2, length(chrs))
    start <- end <- numeric(length(chrs))
    offset <- 0
    for (i in 1:length(chrs)) {
      if (i > 1)
        offset <- gap + end[i - 1]
      chr2match <- df$chr2 == chrs[i]
      start[i] <- min(df$pos2[chr2match]) + offset
      end[i] <- max(df$pos2[chr2match]) + offset
      df$glo2[chr2match] <- df$glo2[chr2match] + offset
    }
    out$chrs <-
      rbind(out$chrs, data.frame(
        chr = chrs, start = start, end = end, genome = genome
      ))
    
    require_package("dbscan")
    dbs <-
      dbscan::dbscan(cbind(df$glo1, df$glo2), eps = eps, minPts = minPts)
    chr1 <- chr2 <- start1 <- start2 <- end1 <- end2 <- c()
    for (i in 1:max(dbs$cluster)) {
      df2 <- df[dbs$cluster == i,]
      chrc1 <- df2$chr1
      chrc2 <- df2$chr2
      for (j in unique(chrc1)) {
        chrl <- which(chrc1 == j)
        chr1 <- c(chr1, j)
        chr2 <- c(chr2, chrc2[stats::median(chrl)])
        mic <- min(chrl)
        mac <- max(chrl)
        start1 <- c(start1, df2$glo1[mic])
        start2 <- c(start2, df2$glo2[mic])
        end1 <- c(end1, df2$glo1[mac])
        end2 <- c(end2, df2$glo2[mac])
      }
      for (j in unique(chrc2)) {
        chrl <- which(chrc2 == j)
        chr1 <- c(chr1, chrc1[stats::median(chrl)])
        chr2 <- c(chr2, j)
        mic <- min(chrl)
        mac <- max(chrl)
        start1 <- c(start1, df2$glo1[mic])
        start2 <- c(start2, df2$glo2[mic])
        end1 <- c(end1, df2$glo1[mac])
        end2 <- c(end2, df2$glo2[mac])
      }
    }
    blocks <- data.frame(chr1, start1, end1, chr2, start2, end2)
    blocks <- unique(blocks)
    if (nrow(blocks) < 1)
      stop("No blocks detected.")
    blocks <- blocks[order(blocks$start1, blocks$start2),]
    for (i in 1:nrow(blocks)) {
      cblock <- blocks[i,]
      m <- max(df$glo1[df$chr1 == cblock$chr1])
      if (cblock$end1 > m){
        blocks$end1[i] <- m
      }
      m <- max(df$glo2[df$chr2 == cblock$chr2])
      if (cblock$end2 > m){
        blocks$end2[i] <- m
      }
      m <- min(df$glo1[df$chr1 == cblock$chr1])
      if (cblock$start1 < m){
        blocks$start1[i] <- m
      }
      m <- min(df$glo1[df$chr2 == cblock$chr2])
      if (cblock$start2 < m){
        blocks$start2[i] <- m
      }
    }

    for (i in unique(blocks$end1)) {
      wb <- which(blocks$end1 == i)
      if (length(wb) > 1) {
        if (blocks$start1[wb][1] > blocks$start1[wb][2]) {
          blocks <- blocks[-wb[2],]
        }else{
          blocks <- blocks[-wb[1],]
        }
      }
    }
    for (i in unique(blocks$start1)) {
      wb <- which(blocks$start1 == i)
      if (length(wb) > 1) {
        if (blocks$end1[wb][1] > blocks$end1[wb][2]) {
          blocks <- blocks[-wb[2],]
        }else{
          blocks <- blocks[-wb[1],]
        }
      }
    }
    del <- c()
    for (i in 1:nrow(blocks)) {
      if (abs(blocks$end1[i] - blocks$start1[i]) < minLength |
          abs(blocks$end2[i] - blocks$start2[i]) < minLength) {
        del <- c(del, i)
      }
    }
    del <- unique(del)
    if (!is.null(del))
      blocks <- blocks[-del,]
    bnames <- names(blocks)
    tmp <-
      apply(blocks, 1, function(x)
        split_sb(x, df[dbs$cluster > 0,], max = maxLength, bnames))
    blocks <- do.call(rbind, tmp)
    names(blocks) <- bnames
    rownames(blocks) <- 1:nrow(blocks)
    out$blocks <- blocks
    out
  }

#' Split synteny blocks
#'
#' @param block Vector,
#' @param df Data.frame,
#' @param max Integer,
#' @param bnames Character,
#' @param smoothing Logical, if block positions should be smoothed
#' @keywords internal
split_sb <- function(block, df, max, bnames, smoothing = TRUE) {
  block <-
    as.data.frame(matrix(block, ncol = 6), stringsAsFactors = FALSE)
  names(block) <- bnames
  block[, c(2, 3, 5, 6)] <- as.numeric(block[, c(2, 3, 5, 6)])
  if ((block$end1 - block$start1) > max |
      block$end2 - block$start2 > max) {
    dfbckp <- df
    df <- df[df$chr1 == block$chr1 & df$chr2 == block$chr2 &
               df$glo1 > block$start1 & df$glo1 < block$end1 &
               df$glo2 > block$start2 & df$glo2 < block$end2,]
    if (nrow(df) < 2)
      return(block)
    df <- df[order(df$glo1),]
    if (smoothing) {
      #dif <- c(0, diff(df$glo2))
      #if(!mean(dif) > 0) dif <- - dif
      #df <- df[dif >= 0, ]
      #if(nrow(df) < 2) return(block)
      df$glo2 <- stats::smooth(df$glo2)
    }
    dif <- mean(c(max(df$glo1), min(df$glo1)))
    center <-
      which(abs(df$glo1 - dif) == min(abs(df$glo1 - dif)))[1]
    block1 <- block
    block1$end1 <- df$glo1[center]
    block1$end2 <- df$glo2[center]
    if (nrow(df) > 1)
      block1 <- split_sb(unlist(block1), dfbckp, max, bnames)
    block2 <- block
    block2$start1 <- df$glo1[center + 1]
    block2$start2 <- df$glo2[center + 1]
    if (nrow(df) > 1)
      block2 <- split_sb(unlist(block2), dfbckp, max, bnames)
    return(rbind(block1, block2))
  }
  block
}

#' Class for synteny information
#' @import R6
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}}.
#' @format \code{\link{R6Class}} object.
#' @docType class
synteny_info <- R6::R6Class(
  classname = "synteny_info",
  public = list(
    blocks = NA,
    chrs = NA,
    initialize = function(blocks, chrs) {
      if (!missing(blocks))
        self$blocks <- blocks
      if (!missing(chrs))
        self$chrs <- chrs
    }
  )
)
