#' Filter Samples
#'
#' Failed or negative control samples (e.g. H2O) should be filtered out before the normalization.
#'
#' @param dat raw_data object.
#' @param samples Character vector with sample names, prefixes of samples to
#' filter or sample indices.
#' If prefixes fit many samples, indices are highly recommended.
#' @return raw_data object
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' raw_napus <- filt_samp(raw_napus, check_raw(raw = raw_napus, plot = FALSE, thresh = 28000))
#' }
#' @export
filt_samp <- function(dat, samples){
  if(is.character(samples)){
    cnames <- dat$samples
    filt <- c()
    for(i in samples){
      filt <- c(filt, grep(i, cnames))
    }
  }else{
    filt <- samples
  }
  if(length(filt) > 0){
    dat$samples <- dat$samples[-filt]
    if(is.matrix(dat$raw)) dat$raw <- dat$raw[, -filt]
    if(is.matrix(dat$sd) & dim(dat$sd)[2] == length(samples)) dat$sd <- dat$sd[, -filt]
    if(is.matrix(dat$beads) & dim(dat$beads)[2] == length(samples)) dat$beads <- dat$beads[, -filt]
    if(is.matrix(dat$intensity)) dat$intensity <- dat$intensity[, -filt]
    if(is.matrix(dat$theta)) dat$theta <- dat$theta[, -filt]
    if(is.matrix(dat$baf)) dat$baf <- dat$baf[, -filt]
    if(is.matrix(dat$rratio)) dat$rratio <- dat$rratio[, -filt]
    if(is.matrix(dat$geno)) dat$geno <- dat$geno[, -filt]
    if(is.matrix(dat$cnv)) dat$cnv <- dat$cnv[, -filt]
    if(is.matrix(dat$tl)) dat$tl <- dat$tl[, -filt]
    if(is.data.frame(dat$cna)) dat$cna <- dat$cna[,dat$cna$ID %in% dat$samples]
  }
  dat
}

#' Filter SNPs
#'
#' Some SNPs do not work as well as others and might be filtered out.
#'
#' @param dat Matrix with raw data or list with intensity, theta, position and chromosome objects.
#' @param filt Character or numberic vector with SNP names or rownumbers to filter out.
#'
#' @return Filtered matrix.
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:100)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- filt_snps(dat, dat$snps[is.na(rowMeans(dat$baf, na.rm = TRUE))])
#' }
#' @export filt_snps
#' @rdname filt_snps
filt_snps <- function(dat, filt){
  UseMethod("filt_snps")
}

#' @rdname filt_snps
#' @export
filt_snps.raw_data <- function(dat, filt){
  if(is.character(filt)){
    filt <- which(dat$snps %in% filt)
  }
  if(length(filt) > 0){
    if(is.matrix(dat$raw)) dat$raw <- dat$raw[-filt, ]
    if(is.vector(dat$snps)) dat$snps <- dat$snps[-filt]
    if(is.matrix(dat$sd)) dat$sd <- dat$sd[-filt, ]
    if(is.matrix(dat$beads)) dat$beads <- dat$beads[-filt, ]
    if(is.vector(dat$pos)) dat$pos <- dat$pos[-filt]
    if(is.vector(dat$chr)) dat$chr <- dat$chr[-filt]
  }
  dat
}


#' @rdname filt_snps
#' @export
filt_snps.norm_data <- function(dat, filt){
  if(is.character(filt)){
    filt <- which(dat$snps %in% filt)
  }
  if(length(filt) > 0){
    if(is.matrix(dat$intensity)) dat$intensity <- dat$intensity[-filt, ]
    if(is.matrix(dat$theta)) dat$theta <- dat$theta[-filt, ]
    if(is.matrix(dat$baf)) dat$baf <- dat$baf[-filt, ]
    if(is.matrix(dat$rratio)) dat$rratio <- dat$rratio[-filt, ]
    if(is.matrix(dat$geno)) dat$geno <- dat$geno[-filt, ]
    if(is.matrix(dat$cnv)) dat$cnv <- dat$cnv[-filt, ]
    if(is.vector(dat$pos)) dat$pos <- dat$pos[-filt]
    if(is.vector(dat$chr)) dat$chr <- dat$chr[-filt]
    if(is.vector(dat$snps)) dat$snps <- dat$snps[-filt]
    if(is.matrix(dat$cnv)) dat$cnv <- dat$cnv[-filt, ]
  }
  dat
}

