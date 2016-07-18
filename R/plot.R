#' Plot raw SNP
#'
#' Plots raw data for specified SNP.
#'
#' @param raw Raw data object.
#' @param n Integer, which SNP should be plotted.
#' @param theta Method for calculation of theta. Currently only "atan2" is implemented.
#' @param transf Method for transformation of the raw values. "none", "log" and "fourth-root" are implemented.
#' @param pn Numeric, p-norm for the intensity calculation.
#' @param nan Numveric to replace NaN for ballele computation (division by zero).
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' plot_raw_snp(raw_napus)
#' }
#' @export
plot_raw_snp <-
  function(raw, n = 1L, transf = "log", pn = 2, theta = "atan2",
           nan = 0.5) {
    #, means = c(0.6, 0.8, 0.9)
    theta <- match.arg(arg = theta, choices = c("atan2", "ballele"))
    row <- raw$raw[n,]
    row <- switch(
      transf,
      none = row,
      log = log(row),
      fourthroot = row ^ (1 / 4)
    )
    storage.mode(row) <- "numeric"
    lr <- length(row)
    sec <- seq(1, lr, 2)
    intensity <- (row[sec] ^ pn + row[sec + 1] ^ pn) ^ (1 / pn)
    if (theta == "atan2") {
      theta <- atan2(row[sec], row[sec + 1]) / (pi / 2)
      #dists <- sapply(means, function(x) abs(theta - x))
      #call <- apply(dists, 1, function(x) which(x == min(x)))
      #call2 <- apply(dists, 1, function(x) which(x == median(x)))
      #cdist<- abs(means[call] - means[call2])
      #theta <- dists[cbind(1:ncol(dists), call2)] / cdist * means[call]
      #theta[intensity < 11] <- NA
    }else{
      row <- row - min(row)
      theta <- row[sec] / rowSums(cbind(row[sec], row[sec + 1]))
      theta[is.na(theta)] <- nan
    }
    require_package("KernSmooth")
    graphics::smoothScatter(theta, intensity)
  }

#' Plot all chromosomes of a sample
#'
#'
#'
#' @param dat List object, containing at least two matrices "baf" and "rratio"
#' and two vectors "chr" and "pos".
#' @param samp Integer, which sample should be plotted.
#' @param ncol Integer, number of colors.
#' @param delthresh Numeric, lower threshold for intensities.
#' @param dupthresh Numeric, upper threshold for intensities.
#' @param ... arguments are forwarded to \code{plot()}.
#' @import graphics
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:100)][-(30000:30100)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' plot_samp(dat, 1)
#' }
#' @export
plot_samp <-
  function(dat, samp, ncol = NULL, delthresh = NULL, dupthresh = NULL, ...) {
    global_pos <- dat$pos
    unichr <- sort(unique(dat$chr))
    if (is.null(ncol)) {
      cols <- grDevices::grey.colors(length(unichr), end = 0.5)
    }else{
      cols <- grDevices::grey.colors(ncol, end = 0.5)
    }
    allcol <- cols[match(dat$chr, unichr) %% ncol + 1]
    if (!is.null(delthresh))
      allcol[dat$rratio[, samp] < delthresh] <- 2
    if (!is.null(dupthresh))
      allcol[dat$rratio[, samp] > dupthresh] <- 3
    
    if (length(unichr) > 1) {
      maxs <- sapply(unichr, function(x)
        max(dat$pos[dat$chr == x]))
      for (i in 2:length(unichr)) {
        global_pos[dat$chr == unichr[i]] <-
          global_pos[dat$chr == unichr[i]] + sum(maxs[1:(i - 1)])
      }
    }else{
      maxs <- max(dat$pos)
    }
    
    csmax <- cumsum(maxs)
    xlim <- c(0, max(global_pos))
    op <- par(
      mfrow = c(2, 1),
      oma = c(5, 4, 6, 2) + 0.1,
      mar = c(0, 0, 1, 1) + 0.1
    )
    plot(
      global_pos, dat$baf[, samp], col = allcol, pch = 19, cex = 0.2, xaxt = "n",
      main = "", xlab = "", ylab = "B Allele Frequency", xlim = xlim, ...
    )
    graphics::axis(
      3, at = csmax - min(maxs) / 2, labels = unichr, tick = FALSE,
      cex.axis = 0.8, las = 2
    )
    abline(v = c(0, csmax))
    title(main = dat$samples[samp], outer = TRUE)
    plot(
      global_pos, dat$rratio[, samp], col = allcol, pch = 19, cex = 0.2, xaxt = "n",
      main = "", xlab = "Position", ylab = "Log2 R Ratio", xlim = xlim, ...
    )
    axis(
      1, at = csmax - min(maxs) / 2, labels = unichr, tick = FALSE,
      cex.axis = 0.8, las = 2
    )
    abline(v = c(0, csmax))
    
    par(op)
  }
#' Plot all chromosomes of a sample
#'
#' Plot all chromosomes of a sample.
#' One subgenome is plotted on top, the other at the bottom.
#' If available, the synteny plots can be plotted between them.
#' Genome structure rearrangements (gsr) (e.g.g deletions or duplications) are
#' highlighted by different colors.
#'
#'
#' @param dat List object, containing at least two matrices "baf"
#' and "rratio" and two vectors "chr" and "pos".
#' @param samp Integer, which sample should be plotted.
#' @param sb Synteny blocks. Data frame with columns start1, start2,
#' end1, end2, chr1 and chr2.
#' @param ncol Number of columns.
#' @param baf Logical, if B-Allele frequency should be plotted.
#' @param tl Logical, if translocations should be highlighted.
#' @param tlcoord Numeric vector of length four. The y-coordinates used to
#' highlight translocations.
#' Top and bottom of top translocation and top and bottom of bottom translocation.
#' @param ... arguments are forwarded to \code{plot()}.
#' @import graphics
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:100)][-(30000:30100)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' dat <- cnv(dat, dup = 0.03, del = -0.06)
#' data("synteny_blocks", package = "brassicaData", envir = environment())
#' plot_gsr(dat, samp = 1, sb = synteny_blocks)
#' }
#' @export
plot_gsr <- function(dat, samp, sb = NULL, ncol = NULL, baf = FALSE,
                     tl = FALSE, tlcoord = c(1.1, 0.68, 0.32,-0.1), ...) {
  if (missing(sb)) {
    require_package("brassicaData")
    synteny_blocks <- NULL
    utils::data("synteny_blocks", package = "brassicaData", envir = environment())
    sb <- synteny_blocks$blocks
  }else if ("synteny_info" %in% class(sb)) {
    sb <- sb$blocks
  }
  unichr <- sort(unique(dat$chr))
  agenome <- grepl("A", dat$chr)
  cgenome <- grepl("C", dat$chr)
  aglobal <- dat$pos[agenome]
  cglobal <- dat$pos[cgenome]
  achr <- grep("A", unichr, value = TRUE)
  cchr <- grep("C", unichr, value = TRUE)
  amaxs <- sapply(achr, function(x)
    max(dat$pos[dat$chr == x]))
  cmaxs <- sapply(cchr, function(x)
    max(dat$pos[dat$chr == x]))
  if (length(achr) > 1) {
    for (i in 2:length(achr)) {
      subs <- dat$chr[agenome] == achr[i]
      aglobal[subs] <- aglobal[subs] + sum(amaxs[1:(i - 1)])
    }
  }
  if (length(cchr) > 1) {
    for (i in 2:length(cchr)) {
      subs <- dat$chr[cgenome] == cchr[i]
      cglobal[subs] <- cglobal[subs] + sum(cmaxs[1:(i - 1)])
    }
  }
  acsmax <- cumsum(amaxs)
  ccsmax <- cumsum(cmaxs)
  axlim <- c(-0.1, max(aglobal) + 0.1)
  cxlim <- c(-0.1, max(cglobal) + 0.1)
  # Define colors
  if (is.null(ncol)) {
    cols <- grDevices::grey.colors(length(achr), end = 0.5)
  }else{
    cols <- grDevices::grey.colors(ncol, end = 0.5)
  }
  allcol <- cols[match(dat$chr, unichr) %% length(cols) + 1]
  if (is.matrix(dat$cnv)) {
    allcol[dat$cnv[, samp] == 1] <- 3
    allcol[dat$cnv[, samp] == -1] <- 2
  }
  if (baf) {
    op <- par(
      mfrow = c(5,1),
      oma = c(5, 4, 6, 2) + 0.1,
      mar = c(0, 0, 0, 0) + 0.1,
      xaxs = "i", yaxs = "r"
    )
  }else{
    op <- par(
      mfrow = c(3,1),
      oma = c(5, 4, 6, 2) + 0.1,
      mar = c(0, 0, 0, 0) + 0.1,
      xaxs = "i", yaxs = "r"
    )
  }
  if (baf) {
    plot(
      aglobal, dat$baf[agenome, samp], pch = 19, col = allcol[agenome],
      cex = 0.2, xaxt = "n",
      main = "", xlab = "", ylab = "B Allele Frequency", xlim = axlim, ylim = 0:1
    )
    axis(
      3, at = acsmax - min(amaxs) / 2, labels = achr, tick = FALSE,
      cex.axis = 0.8, las = 2
    )
    abline(v = c(0, acsmax))
    title(main = dat$samples[samp], outer = TRUE)
  }
  
  plot(
    aglobal, dat$rratio[agenome, samp], pch = 19, col = allcol[agenome],
    cex = 0.2, xaxt = "n",
    main = "", xlab = "", ylab = "Log R Ratio", xlim = axlim, ...
  )
  if (!baf)
    axis(
      3, at = acsmax - min(amaxs) / 2, labels = achr,
      tick = FALSE, cex.axis = 0.8, las = 2
    )
  abline(v = c(0, acsmax))
  if (!baf)
    title(main = dat$samples[samp], outer = TRUE)
  
  # plot synteny blocks
  plot.new()
  cols <-
    grDevices::rainbow(
      n = 10, start = 0, end = 1, alpha = 0.2
    )
  max1 <- max(sb$end1)
  max2 <- max(sb$end2)
  y <- c(1, 1, 0, 0)
  for (i in 1:nrow(sb)) {
    x <- c(sb$start1[i],
           sb$end1[i],
           sb$end2[i],
           sb$start2[i])
    x <- c(x[1:2] / max1, x[3:4] / max2)
    polygon(x, y, col = cols[unique(sb$chr1) %in% sb$chr1[i]])
  }
  
  
  plot(
    cglobal, dat$rratio[cgenome, samp], pch = 19, col = allcol[cgenome],
    cex = 0.2, xaxt = "n",
    main = "", xlab = "", ylab = "Log R Ratio", xlim = cxlim, ...
  )
  if (!baf)
    axis(
      1, at = ccsmax - min(cmaxs) / 2, labels = cchr,
      tick = FALSE, cex.axis = 0.8, las = 2
    )
  abline(v = c(0, ccsmax))
  if (baf) {
    plot(
      cglobal, dat$baf[cgenome, samp], pch = 19, col = allcol[cgenome],
      cex = 0.2, xaxt = "n",
      main = "", xlab = "", ylab = "B Allele Frequency", xlim = cxlim, ylim = 0:1
    )
    axis(
      1, at = ccsmax - min(cmaxs) / 2, labels = cchr,
      tick = FALSE, cex.axis = 0.8, las = 2
    )
    abline(v = c(0, ccsmax))
  }
  # Translocations
  if (tl && is.matrix(dat$tl) && any(as.logical(dat$tl[, samp]))) {
    if (baf) {
      op2 <- par(fig = c(0, 1, 0.2, 0.8), new = TRUE)
    }else{
      op2 <- par(fig = c(0, 1, 0, 1), new = TRUE)
    }
    plot.new()
    for (i in which(as.logical(dat$tl[, samp]))) {
      topx1 <- sb$start1[i] / max1
      topx2 <- sb$end1[i] / max1
      topy1 <- tlcoord[1]
      topy2 <- tlcoord[2]
      botx1 <- sb$start2[i] / max2
      botx2 <- sb$end2[i] / max2
      boty1 <- tlcoord[3]
      boty2 <- tlcoord[4]
      rect(topx1, topy2, topx2, topy1, border = "red")
      rect(botx1, boty2, botx2, boty1, border = "red")
      segments(
        x0 = mean(c(topx1, topx2)), y0 = topy2,
        x1 = mean(c(botx1, botx2)), y1 = boty1, col = "red"
      )
    }
  }
  
  par(op)
}

#' Plot all chromosomes of a population
#'
#' Plot all chromosomes of a sample.
#' The one subgenome is plotted on top, the other at the bottom.
#' If available, the synteny plots can be plotted between them.
#' Deletions and duplications are indicated by different colors.
#'
#'
#' @param dat List object, containing at least two matrices "baf"
#' and "rratio" and two vectors "chr" and "pos".
#' @param sb Synteny blocks. Data frame with columns start1, start2,
#' end1, end2, chr1 and chr2.
#' @param baf Logical, if B-Allele frequency should be plotted.
#' @param cex Size of dots.
#' @param ... arguments are forwarded to \code{plot()}.
#' @import graphics
#' @examples
#' if(require(brassicaData)){
#' data("raw_napus", package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:100)][-(30000:30100)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' dat <- cnv(dat, dup = 0.03, del = -0.06)
#' data("synteny_blocks", package = "brassicaData", envir = environment())
#' plot_global(dat, sb = synteny_blocks)
#' }
#' @export
plot_global <-
  function(dat, sb = NULL, baf = FALSE, cex = 0.2, ...) {
    if (missing(sb)) {
      require_package("brassicaData")
      synteny_blocks <- NULL
      utils::data("synteny_blocks", package = "brassicaData", envir = environment())
      sb <- synteny_blocks$blocks
    }else if ("synteny_info" %in% class(sb)) {
      sb <- sb$blocks
    }
    dotcol <- grDevices::rgb(0,0,0, 0.5)
    unichr <- sort(unique(dat$chr))
    agenome <- grepl("A", dat$chr)
    cgenome <- grepl("C", dat$chr)
    aglobal <- dat$pos[agenome]
    cglobal <- dat$pos[cgenome]
    achr <- grep("A", unichr, value = TRUE)
    cchr <- grep("C", unichr, value = TRUE)
    amaxs <- sapply(achr, function(x)
      max(dat$pos[dat$chr == x]))
    cmaxs <- sapply(cchr, function(x)
      max(dat$pos[dat$chr == x]))
    if (length(achr) > 1) {
      for (i in 2:length(achr)) {
        subs <- dat$chr[agenome] == achr[i]
        aglobal[subs] <- aglobal[subs] + sum(amaxs[1:(i - 1)])
      }
    }
    if (length(cchr) > 1) {
      for (i in 2:length(cchr)) {
        subs <- dat$chr[cgenome] == cchr[i]
        cglobal[subs] <- cglobal[subs] + sum(cmaxs[1:(i - 1)])
      }
    }
    acsmax <- cumsum(amaxs)
    ccsmax <- cumsum(cmaxs)
    axlim <- c(0, max(aglobal))
    cxlim <- c(0, max(cglobal))
    if (baf) {
      op <- par(
        mfrow = c(5,1),
        oma = c(5,4,6,2) + 0.1,
        mar = c(0,0,1,1) + 0.1
      )
    }else{
      op <- par(
        mfrow = c(3,1),
        oma = c(5,4,6,2) + 0.1,
        mar = c(0,0,1,1) + 0.1
      )
    }
    if (baf) {
      plot(
        aglobal, rowMeans(0.5 - abs(dat$baf[agenome,] - 0.5)),
        pch = 19, cex = cex, xaxt = "n",
        main = "", xlab = "", ylab = "B Allele Frequency",
        xlim = axlim, ylim = c(0, 0.5), ...
      )
      axis(
        3, at = acsmax - min(amaxs) / 2, labels = achr,
        tick = FALSE, cex.axis = 0.8, las = 2
      )
      abline(v = c(0, acsmax))
      title(main = "Population mean", outer = TRUE)
    }
    
    plot(
      aglobal, rowMeans(dat$rratio[agenome,]),
      pch = 19, cex = cex, xaxt = "n",
      main = "", xlab = "", ylab = "Log R Ratio",
      xlim = axlim, ...
    )
    if (!baf)
      axis(
        3, at = acsmax - min(amaxs) / 2, labels = achr,
        tick = FALSE, cex.axis = 0.8, las = 2
      )
    abline(v = c(0, acsmax))
    if (!baf)
      title(main = "Population mean", outer = TRUE)
    # plot synteny blocks
    plot.new()
    cols <-
      grDevices::rainbow(
        n = 10, start = 0, end = 1, alpha = 0.2
      )
    max1 <- max(sb$end1)
    max2 <- max(sb$end2)
    y <- c(1, 1, 0, 0)
    for (i in 1:nrow(sb)) {
      x <- c(sb$start1[i],
             sb$end1[i],
             sb$end2[i],
             sb$start2[i])
      x <- c(x[1:2] / max1, x[3:4] / max2)
      polygon(x, y, col = cols[unique(sb$chr1) %in% sb$chr1[i]])
    }
    plot(
      cglobal, rowMeans(dat$rratio[cgenome,]), pch = 19,
      cex = cex, xaxt = "n",
      main = "", xlab = "", ylab = "Log R Ratio", xlim = cxlim, ...
    )
    if (!baf)
      axis(
        1, at = ccsmax - min(cmaxs) / 2, labels = cchr,
        tick = FALSE, cex.axis = 0.8, las = 2
      )
    abline(v = c(0, ccsmax))
    if (baf) {
      plot(
        cglobal, rowMeans(0.5 - abs(dat$baf[cgenome,] - 0.5)),
        pch = 19, cex = cex, xaxt = "n",
        main = "", xlab = "", ylab = "B Allele Frequency",
        xlim = cxlim, ylim = c(0, 0.5), ...
      )
      axis(
        1, at = ccsmax - min(cmaxs) / 2, labels = cchr,
        tick = FALSE, cex.axis = 0.8, las = 2
      )
      abline(v = c(0, ccsmax))
    }
    par(op)
  }

#' Plot translocations
#'
#' Plots translocations for whole population.
#'
#' @param dat List object, containing at least two matrices "intensity"
#' and "theta" and two vectors "chr" and "pos".
#' @param sb synteny block object
#' @import graphics
#' @examples
#' if(require(brassicaData)){
#' data(raw_napus, package = "brassicaData", envir = environment())
#' \dontshow{
#' raw_napus <- filt_samp(raw_napus, raw_napus$samples[-(1:10)])
#' raw_napus <- filt_snps(raw_napus, raw_napus$snps[-(1:100)][-(30000:30100)])
#' }
#' dat <- intens_theta(raw_napus)
#' dat <- remove_suffix(dat, "_Grn")
#' dat <- geno_baf_rratio(dat, delthresh = 11)
#' dat <- segm(dat)
#' dat <- cnv(dat, dup = 0.03, del = -0.06)
#' data("synteny_blocks", package = "brassicaData", envir = environment())
#' dat <- trans_location(dat, synteny_blocks, min = 5)
#' plot_trans_locations(dat, sb = synteny_blocks$blocks)
#' }
#' @export
plot_trans_locations <- function(dat, sb = NULL) {
  if (missing(sb)) {
    require_package("brassicaData")
    synteny_blocks <- NULL
    utils::data("synteny_blocks", package = "brassicaData", envir = environment())
    sb <- synteny_blocks$blocks
  }else if ("synteny_info" %in% class(sb)) {
    sb <- sb$blocks
  }
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  image(dat$tl, xaxt = "n", yaxt = "n", col = 0:1)
  axis(
    1, at = seq(0, 1, length.out = nrow(sb)), labels =
      paste(sb$chr1, " / ", sb$chr2, sep = ""), las = 2
  )
  axis(
    2, at = seq(0, 1, length.out = length(dat$samples)),
    labels = dat$samples, las = 1, cex.axis = 0.3
  )
  par(op)
}
