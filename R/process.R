#' Genotype with kmeans using three clusters
#'
#' @param dat List object, containing at least two matrices "intensity" and "theta".
#' Or matrix with raw data.
#' @param delthresh Numeric between 0 and 1. Intensity threshold for deletions.
#' @param drop Logical, if TRUE theta and intensity values are removed to save memory.
#' @param corr Logical, if TRUE all sample medians are corrected to 0.
#'
#' @return Genotypes
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' }
#' @export
geno_baf_rratio <-
  function(dat, delthresh = 0.7, drop = FALSE, corr = FALSE) {
    clusterm <-
      find_peak(as.vector(dat$theta), npeaks = 3, breaks = 1000)
    if (!length(clusterm) == 3)
      clusterm <- c(0.2, 0.5, 0.8)
    require_package("Ckmeans.1d.dp")
    if (.Platform$OS.type == "unix") {
      sink('/dev/null')
    } else {
      sink("NUL")
    }
    res <-
      lapply(1:length(dat$snps), function(x)
        call_geno(dat$theta[x,], dat$intensity[x,], clusterm, delthresh))
    sink()
    dat$baf <- t(sapply(res, function(x)
      x$baf))
    dat$geno <- t(sapply(res, function(x)
      x$geno))
    dat$rratio <- t(sapply(res, function(x)
      x$rratio))
    if (drop) {
      dat$intensity <- NULL
      dat$theta <- NULL
    }
    if (corr) {
      med <- apply(dat$rratio, 2, stats::median, na.rm = TRUE)
      dat$rratio <- sweep(dat$rratio, 2, med, "-")
    }
    dat
  }



#' Genotype with kmeans using three clusters
#'
#' @param theta Vector with theta values.
#' @param intensity Vector with intensity values.
#' @param clusterm Vector containing the cluster means.
#' @param delthresh Threshold for intensity based deletion calling.
#'
#' @return List containing genotypes, b-allele frequencies and R ratios.
#'
#' @keywords internal
call_geno <-
  function(theta, intensity, clusterm = c(0.2, 0.5, 0.8), delthresh = NULL) {
    outlen <- length(theta)
    res <-
      suppressWarnings(Ckmeans.1d.dp::Ckmeans.1d.dp(theta, k = c(1, 3)))
    k <- max(res$cluster)
    lower <- mean(clusterm[1:2])
    upper <- mean(clusterm[2:3])
    dis <- mean(diff(clusterm)) * 0.5
    if (k == 1) {
      out <- rep(NA, outlen)
      return(list(
        geno = out, baf = out, rratio = out
      ))
    }else if (k == 2) {
      if (res$centers[1] < lower &
          res$centers[2] > lower & res$centers[2] < upper &
          res$centers[2] - res$centers[1] > dis) {
        #AA AB
        out <- res$cluster - 1
        centers <- c(res$centers, res$centers[2] + diff(res$centers))
      }else if (res$centers[1] < lower &
                res$centers[2] > upper) {
        #AA BB
        out <- res$cluster
        out[out == 1] <- 0
        centers <-
          c(res$centers[1], mean(res$centers, na.rm = TRUE), res$centers[2])
      }else if (res$centers[1] > lower &
                res$centers[1] < upper & res$centers[2] > upper &
                res$centers[2] - res$centers[1] > dis) {
        #AB BB
        out <- res$cluster
        centers <-
          c(max(res$centers[1] - diff(res$centers), 0, na.rm = TRUE), res$centers)
      }else{
        out <- rep(NA, outlen)
        return(list(
          geno = out, baf = out, rratio = out
        ))
      }
    }else{
      if (res$centers[1] < lower &
          res$centers[2] > lower &
          res$centers[2] < upper & res$centers[3] > upper &
          res$centers[2] - res$centers[1] > dis &
          res$centers[3] - res$centers[2] > dis) {
        out <- res$cluster - 1
        centers <- res$centers
      }else{
        out <- rep(NA, outlen)
        return(list(
          geno = out, baf = out, rratio = out
        ))
      }
    }
    dc <- diff(centers)
    left <- 0.5 * (theta - centers[1]) / dc[1]
    right <- 0.5 + 0.5 * (theta - centers[2]) / dc[2]
    baf <- left
    baf[theta >= centers[2]] <- right[theta >= centers[2]]
    baf[theta <= centers[1]] <- 0
    baf[theta >= centers[3]] <- 1
    out[intensity < delthresh] <- -1
    vcenters <- c(
      mean(intensity[out == 0], na.rm = TRUE),
      mean(intensity[out == 1], na.rm = TRUE),
      mean(intensity[out == 2], na.rm = TRUE)
    )
    vcenters[is.na(vcenters)] <- mean(vcenters, na.rm = TRUE)
    rratio <- log2(intensity / interpol(theta, centers, vcenters))
    rratio[is.infinite(rratio)] <- min(rratio[!is.infinite(rratio)])
    list(geno = out, baf = baf, rratio = rratio)
  }




#' Segmentation of R ratio values
#'
#' Wrapper for \code{DNAcopy::segment}.
#' Other methods (e.g. sliding window approaches) might be added in a future version.
#'
#' @param dat List object, containing at least two matrices "baf" and "rratio" and two vectors "chr" and "pos".
#' @return norm_data object including a CNA object.
#' @examples
#' \dontrun{
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' }
#' }
#' @export
segm <- function(dat) {
  require_package("DNAcopy")
  cna <-
    DNAcopy::CNA(
      genomdat = dat$rratio, chrom = dat$chr, maploc = dat$pos, data.type = "logratio", sampleid = dat$samples
    )
  cna <- DNAcopy::smooth.CNA(cna)
  cna <-
    DNAcopy::segment(
      cna, alpha = 0.05, nperm = 1000, min.width = 5, undo.splits = "none", verbose = 0
    )$output
  dat$cna <- cna
  dat
}

#' Copy Number Variation
#'
#' Assign copy number variations (duplications and deletions) based on the a threshold.
#' The CNVs are assigned to all SNPs in a CNV segment.
#' The CNV segments can be calculated using \code{segm}.
#'
#' @param dat List object, containing at least two matrices "baf" and "rratio" and two vectors "chr" and "pos".
#' @param del Lower threshold, everything below is called deletion.
#' @param dup Upper threshold, everything above is called duplication.
#' @examples
#' \dontrun{
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' dat <- cnv(dat, dup = 0.03, del = -0.06)
#' }
#' }
#' @export
cnv <- function(dat, del = -0.1, dup = 0.1) {
  if (!is.list(dat$cna)) {
    dat <- segm(dat)
  }
  dat$cnv <-
    matrix(0, ncol = length(dat$samples), nrow = length(dat$snps))
  dels <- dat$cna[dat$cna$seg.mean < del,]
  dups <- dat$cna[dat$cna$seg.mean > dup,]
  for (i in unique(dels$ID)) {
    samp <- which(dat$samples == i)
    for (j in unique(dels$chrom)) {
      delchr <- dels[dels$chrom == j & dels$ID == i,]
      if (nrow(delchr) == 0)
        next()
      if (ncol(delchr) == 0)
        stop("Sample name mismatch.")
      chrmatch <- dat$chr == j
      for (k in 1:nrow(delchr)) {
        dat$cnv[chrmatch, samp][dat$pos[chrmatch] > delchr$loc.start[k] &
                                  dat$pos[chrmatch] < delchr$loc.end[k]] <-
          -1
      }
    }
  }
  for (i in unique(dups$ID)) {
    samp <- which(dat$samples == i)
    for (j in unique(dups$chrom)) {
      dupchr <- dups[dups$chrom == j & dups$ID == i,]
      if (nrow(dupchr) == 0)
        next()
      if (ncol(dupchr) == 0)
        stop("Sample name mismatch.")
      chrmatch <- dat$chr == j
      for (k in 1:nrow(dupchr)) {
        dat$cnv[chrmatch, samp][dat$pos[chrmatch] > dupchr$loc.start[k] &
                                  dat$pos[chrmatch] < dupchr$loc.end[k]] <-
          1
      }
    }
  }
  dat
}

#' Find translocations
#'
#' Synteny blocks and CNVs are combined to detect matching translocations.
#'
#' @param dat List object, containing at least two matrices "baf" and "rratio"
#' and two vectors "chr" and "pos".
#' @param sb synteny_info object
#' @param min Integer, minimal number of markers below / above threshold.
#' @return norm_data object including a CNA object.
#' @examples
#' \dontrun{
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' dat <- cnv(dat, dup = 0.03, del = -0.06)
#' data("synteny_blocks", package = "brassicaData", envir = environment())
#' dat <- trans_location(dat, synteny_blocks, min = 5)
#' }
#' }
#' @export
trans_location <- function(dat, sb, min = 1L) {
  sbb <- sb$blocks
  sbc <- sb$chrs
  out <- matrix(0, ncol = length(dat$samples), nrow = nrow(sbb))
  for (i in 1:nrow(sbb)) {
    off1 <- sbc$start[sbc$chr == sbb$chr1[i]]
    off2 <- sbc$start[sbc$chr == sbb$chr2[i]]
    cnv1 <- dat$cnv[dat$chr == sbb$chr1[i] &
                      dat$pos > (sbb$start1[i] - off1) &
                      dat$pos < (sbb$end1[i] - off1),]
    cnv2 <- dat$cnv[dat$chr == sbb$chr2[i] &
                      dat$pos > (sbb$start2[i] - off2) &
                      dat$pos < (sbb$end2[i] - off2),]
    if (nrow(cnv1) < 1 || nrow(cnv2) < 1)
      next()
    #Duplication
    l1 <-
      apply(cnv1, 2, function(x)
        sum(x == 1, na.rm = TRUE) >= min)
    if (any(l1)) {
      l2 <- apply(cnv2, 2, function(x)
        sum(x == -1, na.rm = TRUE) >= min)
      out[i, l1 & l2] <- 1
    }
    l2 <-
      apply(cnv2, 2, function(x)
        sum(x == 1, na.rm = TRUE) >= min)
    if (any(l2)) {
      l1 <- apply(cnv1, 2, function(x)
        sum(x == -1, na.rm = TRUE) >= min)
      out[i, l1 & l2] <- 1
    }
  }
  dat$tl <- out
  dat
}
