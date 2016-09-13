#' Preprocess raw values and calculate Intensity and Theta values
#'
#' @param raw Raw_data object.
#' @param norm Method for the normalization. One of "none", "quantile", "median" or "both".
#' @param scaling Logical, if each SNP should be scaled or not.
#' @param transf Method for transformation of the raw values.
#' "none", "log" and "fourth-root" are implemented.
#' @param pn Numeric, p-norm for the intensity calculation.
#' @return List with two matrices "intensity" (signal intensities) and "theta" (genotype value).
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData")
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
#' }
#' dat <- intens_theta(raw_napus)
#' }
#' @export
intens_theta <- function(raw, norm = "quantile",  scaling = "mean",
                         transf = "log", pn = 2) {
  norm <-
    match.arg(arg = norm, choices = c("none", "quantile", "median", "both"))
  transf <-
    match.arg(arg = transf, choices = c("none", "log", "fourthroot"))
  scaling <-
    match.arg(arg = scaling, choices = c("none", "ztrans", "mean"))
  sec <- seq(1,ncol(raw$raw),2)
  out <- list()
  out$samples <- raw$samples[sec]
  out$snps <- raw$snps
  out$chr <- raw$chr
  out$pos <- raw$pos
  class(out) <- "norm_data"
  if (norm == "quantile" | norm == "both") {
    require_package("preprocessCore")
  }
  if (norm == "median" | norm == "both") {
    require_package("limma")
  }
  normalize <- function(raw) {
    for (i in seq(2, ncol(raw), 2)) {
      raw[, c(i - 1, i)] <-
        preprocessCore::normalize.quantiles(raw[, c(i - 1, i)])
    }
    limma::normalizeMedianAbsValues(raw)
  }
  raw$raw <- switch(
    norm,
    quantile = preprocessCore::normalize.quantiles(raw$raw),
    median = limma::normalizeMedianAbsValues(raw$raw),
    both = normalize(raw$raw),
    none = raw$raw
  )
  raw$raw <- switch(
    transf,
    none = raw$raw,
    log = log(raw$raw),
    fourthroot = raw$raw ^ (1 / 4)
  )
  raw$raw[is.infinite(raw$raw)] <- NA
  raw$raw <- switch(
    scaling,
    none = raw$raw,
    ztrans = t(scale(t(raw$raw))),
    mean = sweep(
      raw$raw, 1, rowMeans(raw$raw, na.rm = TRUE) - mean(raw$raw, na.rm = TRUE)
    )
  )
  raw$raw[is.na(raw$raw)] <- 0
  storage.mode(raw$raw) <- "numeric"
  out$intensity <- matrix(unlist(lapply(sec,
                                        function(x)
                                          (raw$raw[, x] ^ pn + raw$raw[, x + 1] ^ pn) ^ (1 / pn))),
                          ncol = ncol(raw$raw) / 2, byrow = FALSE)
  out$theta <- matrix(unlist(lapply(sec,
                                    function(x)
                                      atan2(raw$raw[, x], raw$raw[, x + 1]) / (pi / 2))),
                      ncol = ncol(raw$raw) / 2, byrow = FALSE)
  
  if(is.matrix(raw$beads) && dim(raw$beads) != c(1,1)){
    out$beads <- matrix(unlist(lapply(sec,
                                      function(x) 
                                        (raw$beads[, x] +  raw$beads[, x + 1]) / 2)),
                        ncol = ncol(raw$raw) / 2, byrow = FALSE)
  }
  if(is.matrix(raw$sd) && dim(raw$sd) != c(1,1)){
    out$sd <- matrix(unlist(lapply(sec,
                                      function(x) 
                                        (raw$sd[, x] + raw$sd[, x + 1]) / 2)),
                        ncol = ncol(raw$raw) / 2, byrow = FALSE)
  }
  out
}
