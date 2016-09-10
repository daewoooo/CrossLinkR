#' Search for hotspots of crosslinked reads
#' 
#' @param crosslinks File containing analyzed links between read-pairs
#' @param maxDist Maximal distance to connect neighbouring reads into hotspots.
#' @param minCov Minimal coverage to filter hotspot locations.
#' @param pval P value cut-off to filter hotspot locations.
#' @param chrom.lengths Vector containing lenghts of all analyzed chromosomes (Info in Bam header).
#' 
#' @author David Porubsky
#' @export

hotLinks <- function(crosslinks=NULL, maxDist=1000, minCov=4, pval=NULL, ntrials=100, chrom.lengths) {
  
  #Helper functions
  connectRanges <- function(gr) {
    hits <- GenomicRanges::findOverlaps(gr, gr, maxgap = maxDist) #connect ranges in distance 1kb
    join.gr <- GenomicRanges::split(gr[subjectHits(hits)], queryHits(hits))
    join.gr <- endoapply(join.gr, mergeGR)
    join.gr <- GenomicRanges::reduce(unlist(join.gr))
    return(join.gr)
  }
  
  mergeGR <- function(gr) {
    new.gr <- GenomicRanges::GRanges(seqnames=as.character(seqnames(gr))[1], ranges=IRanges(start=start(gr[1]), end=end(gr[length(gr)])))
    return(new.gr)
  }
  
  #read in discordant read links
  links <- read.table(crosslinks, stringsAsFactors = F, header=T)
  mate1 <- do.call(rbind, strsplit(links$mate1.coord, ","))
  mate2 <- do.call(rbind, strsplit(links$mate2.coord, ","))
  mate1 <- GenomicRanges::GRanges(seqnames=mate1[,1], ranges=IRanges(start=as.numeric(mate1[,2]), end=as.numeric(mate1[,3])))
  mate2 <- GenomicRanges::GRanges(seqnames=mate2[,1], ranges=IRanges(start=as.numeric(mate2[,2]), end=as.numeric(mate2[,3])))
  
  #merge overlapping read pairs
  hits <- GenomicRanges::findOverlaps(mate1, mate2)
  mask <- queryHits(hits)[queryHits(hits) == subjectHits(hits)]
  mates <- IRanges::append(mate1[mask], mate2[mask])
  mates[seq(from=1, to=length(mates), by=2)] <- mate1[mask]
  mates[seq(from=2, to=length(mates), by=2)] <- mate2[mask]
  mateID <- rep(1:(length(mates)/2), each=2)
  mates.grl <- GenomicRanges::split(mates, mateID)
  merged.mates <- unlist(GenomicRanges::reduce(mates.grl))
  mate1 <- mate1[-mask] 
  mate2 <- mate2[-mask]
  
  #connecting ranges that are maxDist apart from each other
  reduced.gr <- c(mate1, mate2, merged.mates)
  suppressWarnings( joined.gr <- connectRanges(reduced.gr) )
  
  #merge pairs overlapping with single hotspot
  hit.mate1 <- GenomicRanges::findOverlaps(mate1, joined.gr)
  hit.mate2 <- GenomicRanges::findOverlaps(mate2, joined.gr)
  to.merge <- queryHits(hit.mate1) == queryHits(hit.mate2) & subjectHits(hit.mate1) == subjectHits(hit.mate2)
  joined.pairs <- mate1[to.merge] #take only mate1 to represent pairs which overlaps with single hotspot
  mate1 <- mate1[!to.merge] 
  mate2 <- mate2[!to.merge]
  reduced.joined.gr <- c(mate1, mate2, merged.mates, joined.pairs)
  
  counts <- GenomicRanges::countOverlaps(joined.gr, reduced.joined.gr)
  joined.gr$counts <- counts
  
  #normalize counts byt total number of ranges
  norm.counts <- (joined.gr$counts/length(reduced.joined.gr))*10000
  joined.gr$norm.counts <- norm.counts
  
  #sort
  chrom.lengths <- chrom.lengths[names(chrom.lengths) %in% seqlevels(joined.gr)]
  seqlevels(joined.gr) <- names(chrom.lengths)
  seqlengths(joined.gr) <- chrom.lengths
  seqlevels(reduced.joined.gr) <- names(chrom.lengths)
  seqlengths(reduced.joined.gr) <- chrom.lengths
  joined.gr <- sort(joined.gr)
  reduced.joined.gr <- sort(reduced.joined.gr)
  
  if (!is.null(minCov)) {
    hotspots <- joined.gr[joined.gr$counts >= minCov]
  } else if (!is.null(pval)) {
    
    trial.counts <- matrix(nrow=length(joined.gr), ncol=ntrials)
    for (i in 1:ntrials) {
      #message("Working on trial ",i)
      #random.shift <- runif(1, 1000000, 100000000)
      #shifted.hotspots <- shiftRanges(ranges = joined.gr, shiftrange = random.shift)
      #counts <- GenomicRanges::countOverlaps(shifted.hotspots, reduced.joined.gr)
      n <- runLength(seqnames(reduced.joined.gr))
      chrom.lengths.adjust <- chrom.lengths - max(width(reduced.joined.gr))
      random.pos <- round(stats::runif(length(reduced.joined.gr), min=1, max=rep(chrom.lengths.adjust, n)))
      random.reduced.joined.gr <- GRanges(seqnames=seqnames(reduced.joined.gr), ranges=IRanges(start=random.pos, width = width(reduced.joined.gr)))
      counts <- GenomicRanges::countOverlaps(joined.gr, random.reduced.joined.gr)
      trial.counts[,i] <- counts
    }
    
    trial.counts <- as(trial.counts, "vector")
    pVals <- 1-ecdf(trial.counts)(joined.gr$counts)
    #pVals.adjust <- p.adjust(pVals)
    #joined.gr$pVals <- pVals.adjust
    joined.gr$pVals <- pVals
    hotspots <- joined.gr[joined.gr$pVals <= pval]
    
  } else {
    hotspots <- joined.gr
  }   
  
  return(hotspots)
} 