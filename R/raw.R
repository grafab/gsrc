#' Read multiple idat files.
#'
#' This function is a wrapper to for illuminaio's \code{readIDAT} function.
#' A dictionary will translate the standard IDs from the IDAT file into SNP names.
#' A position file, will filter out all SNPs where the position is unknown.
#' This might save a lot of memory and allow to read in more samples at once.
#' This function is provided to simplify the use of the example data.
#' For data produced with other arrays use respective functions available in CRAN or Bioconductor.
#'
#' @param files A vector of filenames.
#' Two consecutive lines per sample.
#' The order (e.g. green,red) must be consistent for all samples.
#' @param dict A dataframe with columns "idatID" and "name".
#' @param cnames Vector of characters as column names.
#' @param pos Position object with at least the column "name". Is used to filter out SNPs where the position is not known.
#' @param beads Logical, if beads should be read in.
#' @param sd Logical, if standard deviation should be read in.
#' @return An object containing the identities from the idat files.
#' @examples
#' library(brassicaData)
#' files <- list.files(system.file("extdata", package = "brassicaData"),
#' full.names = TRUE, pattern = "idat")
#' samples <- read_sample_sheets(files = list.files(system.file("extdata",
#' package = "brassicaData"), full.names = TRUE, pattern = "csv"))
#' column_names <- sapply(strsplit(files,split="/"), FUN=function(x) x[length(x)])
#' data("dictionary", package = "brassicaData", envir = environment())
#' data("chrPos", package = "brassicaData", envir = environment())
#' raw_data <- read_intensities(files = files, dict = dictionary,
#' cnames = column_names, pos = chrPos)
#' @export
read_intensities <- function(files, dict = NULL, cnames = NULL, pos = NULL,
                             beads = FALSE, sd = FALSE) {
  require_package("illuminaio")
  out <- list(chr = c(), pos = c(), raw = matrix(), snps = c(), beads = matrix(), sd = matrix())
  cols <- 1
  if(sd) cols <- c(cols, 2)
  if(beads) cols <- c(cols, 3)
  if(is.null(dict)){
    out$raw <- matrix(unlist(lapply(files,FUN = function(x) illuminaio::readIDAT(x)$Quants[, cols])),
                      ncol=length(files),byrow = FALSE)
  }else{
    if(!is.null(pos)){
      dict <- dict[dict$name %in% pos$name, ]
      if(is.vector(pos$chromosome) & is.vector(pos$position)){
        tmp<-pos$name %in% dict$name
        out$chr <- pos$chromosome[tmp]
        out$pos <- pos$position[tmp]
        rm(tmp)
      }
    }
    out$raw <- matrix(unlist(lapply(files,FUN = function(x) illuminaio::readIDAT(x)$Quants[as.character(dict$idatID), cols])),
                      ncol=length(files)*length(cols),byrow = FALSE)
    out$snps <- dict$name
  }
  if(beads & sd){
    out$sd <- out$raw[, seq(2, ncol(out$raw), 3)]
    out$raw <-  out$raw[,- seq(2, ncol(out$raw), 3)]
    out$beads <- out$raw[, seq(2, ncol(out$raw), 2)]
    out$raw <-  out$raw[,- seq(2, ncol(out$raw), 2)]
  }else if(sd){
    out$sd <- out$raw[, seq(2, ncol(out$raw), 2)]
    out$raw <-  out$raw[,- seq(2, ncol(out$raw), 2)]
  }else if(beads){
    out$beads <- out$raw[, seq(2, ncol(out$raw), 2)]
    out$raw <-  out$raw[,- seq(2, ncol(out$raw), 2)]
  }

  if(!is.null(cnames)){
    out$samples <- cnames
  }
  class(out) <- "raw_data"
  out
}


#' Read sample sheet(s)
#'
#' Read one or multiple sample sheets in csv format.
#' The information is used to translate between cryptic names
#' (e.g. provided by service provider) and names that can be interpreted.
#'
#' @param files Path to sample sheet.
#' @param skip Integer, lines to skip, before data is read.
#' If not provided, the program looks for entry __[Data]__
#' @param cols Character vector, column names to use.
#' @return A data.frame containing the idat names and the meaningful names of the samples.
#' @examples
#' if(require(brassicaData)){
#' samples <- read_sample_sheets(files = list.files(system.file("extdata",
#' package = "brassicaData"), full.names = TRUE, pattern = "csv"))
#' }
#' @export
read_sample_sheets <- function(files, skip = NULL, cols =
                                 c("SampleID", "SentrixBarcode_A", "SentrixPosition_A")) {
  tab <- c()
  for(i in files){
    if(missing(skip)){
      lines1 <- readLines(i)
      skip_loc <- grep(lines1, pattern = "\\[Data\\]")
    }else{
      skip_loc <- skip
    }
    tab2 <- utils::read.table(i, skip = skip_loc, sep = ";", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    if(ncol(tab2) < length(cols)){
      tab2 <- utils::read.table(i, skip = skip_loc, sep = ",", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    }
    if(!is.null(tab) & is.null(cols)){
      both <- intersect(colnames(tab), colnames(tab2))
      tab <- tab[, both]
      tab2 <- tab2[, both]
    }else if(!is.null(tab)){
      tab <- tab[, 1:length(cols)]
      tab2 <- tab2[, 1:length(cols)]
      colnames(tab) <- colnames(tab2) <- cols
    }else if(is.null(tab)){
      tab2 <- tab2[, 1:length(cols)]
      colnames(tab2) <- cols
    }
    tab <- rbind(tab, tab2)
  }
  return(data.frame(Names = tab$SampleID, ID = paste(tab$SentrixBarcode_A, tab$SentrixPosition_A, sep = "_"),stringsAsFactors = FALSE))
}


#' Check raw data
#'
#' Before preprocessing it is advised to check the data for failed samples.
#' This function allows to visualize the signals and return sample indices of failed samples.
#'
#' @param raw Matrix with raw data.
#' Two consecutive columns per sample.
#' The order (e.g. green,red) must be consistent for all samples.
#' @param plot Logical. If TRUE data will be plotted.
#' @param thresh Threshold for filtering.
#' @param ... arguments are forwarded to \code{hist()}.
#' @return Vector containing Samples below threshold.
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' to_filter <- check_raw(raw_napus, thresh = 28000, breaks = 20)
#' }
#' @export
check_raw <- function(raw, plot = TRUE, thresh = 0, ...){
  cmeans <- colMeans(raw$raw)
  cmeans <- cmeans[seq(1,ncol(raw$raw),2)] + cmeans[seq(2,ncol(raw$raw),2)]
  if(plot){
    graphics::hist(cmeans,  main = "Histogram", xlab = "Mean signal per sample", ...)
    graphics::abline(v=thresh, col = "red")
  }
  out<-which(cmeans<thresh)*2
  sort(c(out-1,out))
}
