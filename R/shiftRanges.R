#' Shift ranges in certain distance
#' 
#' @param ranges An \code{\link{GRanges}} object to be shifted.
#' @param shiftrange The distance, given \code{\link{GRanges}} object will be shifted.
#' 
#' @author David Porubsky
#' @export

shiftRanges <- function(ranges=NULL, shiftrange=10000000) {
  isCircular(ranges) <- rep(T, length(seqlevels(ranges)))
  
  shifted.grl <- GRangesList()
  for (i in seqlevels(ranges)) {
    chr.gr <- ranges[seqnames(ranges) == i]
    shift.gr <- shift(chr.gr, shiftrange)
    shift.gr[start(shift.gr) > seqlengths(chr.gr)[i]] <- shift(shift.gr[start(shift.gr) > seqlengths(chr.gr)[i]], -seqlengths(chr.gr)[i])
    shifted.grl[[i]] <- shift.gr
  }
  
  shifted.ranges <- unlist(shifted.grl)
  suppressWarnings( isCircular(shifted.ranges) <- rep(F, length(seqlevels(shifted.ranges))) )
  return(shifted.ranges)
}