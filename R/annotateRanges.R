#' Function searching for overlaps with various annotation ranges.
#' 
#' @param searchGRanges Ranges to be searched for overlaps with annotation ranges.
#' @param maxDist Maximal distance of a search range from the annotation range to be considered as hit.
#' @param blacklist List of regions to be excluded from the analysis (format: chromosome start end).
#' @param ... Files containing various annotation ranges (format: chromosome start end).
#' 
#' @author David Porubsky
#' @export

annotateRanges <- function(searchGRanges=NULL, maxDist=1000, blacklist=NULL, ...) {
  
  annotation <- list(...)
  if (length(annotation)==0) {
    stop()
    return(NULL)
  }
  
  if (!is.null(blacklist)) {
    filt.regions <- read.table(blacklist, header=F)
    filt.regions <- GenomicRanges::GRanges(seqnames=filt.regions$V1, ranges=IRanges(start=filt.regions$V2, end=filt.regions$V3))
    suppressWarnings( filt.hits <- findOverlaps(searchGRanges, filt.regions) )
    #suppressWarnings( filt.hits <- findOverlaps(searchGRanges, filt.regions, maxgap = 1000000) )
    if (length(filt.hits)>0) {
      searchGRanges <- searchGRanges[-queryHits(filt.hits)]
    }  
  }
  
  for (i in 1:length(annotation)) {
    ID <- names(annotation[i])
    bedfile <- annotation[[i]]
    bed <- read.table(bedfile, stringsAsFactors = F)
    bed.gr <- GRanges(seqnames=bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3))
    bed.gr <- reduce(bed.gr)
    
    mcols(searchGRanges)[,eval(ID)] <- rep(0, length(searchGRanges))
    hits <- findOverlaps(bed.gr, searchGRanges, maxgap = maxDist)
    hits <- unique(subjectHits(hits))
    mcols(searchGRanges)[,eval(ID)][hits] <- 1 
  }
  return(searchGRanges)
}
