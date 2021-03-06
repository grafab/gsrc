% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filt_snps}
\alias{filt_snps}
\alias{filt_snps.norm_data}
\alias{filt_snps.raw_data}
\title{Filter SNPs}
\usage{
filt_snps(dat, filt)

\method{filt_snps}{raw_data}(dat, filt)

\method{filt_snps}{norm_data}(dat, filt)
}
\arguments{
\item{dat}{Matrix with raw data or list with intensity, theta, position and chromosome objects.}

\item{filt}{Character or numberic vector with SNP names or rownumbers to filter out.}
}
\value{
Filtered matrix.
}
\description{
Some SNPs do not work as well as others and might be filtered out.
}
\examples{
if(require(brassicaData)){
data("raw_napus", package = "brassicaData", envir = environment())
\dontshow{
raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:100)])
raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:10)])
}
dat <- intens_theta(raw_napus)
dat <- remove_suffix(dat, "_Grn")
dat <- geno_baf_rratio(dat, delthresh = 11)
dat <- filt_snps(dat, dat$snps[is.na(rowMeans(dat$baf, na.rm = TRUE))])
}
}

