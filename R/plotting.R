#' Plot size distribution between links
#' 
#' @param linkFile Filename containing link coordinates
#' 
#' @author David Porubsky
#' @export

plotlinkSizeDist <- function(linkFile) {
  links <- read.table(linkFile, stringsAsFactors = F, header=T)
  mate1 <- do.call(rbind, strsplit(links$mate1.coord, ","))
  mate2 <- do.call(rbind, strsplit(links$mate2.coord, ","))
  coord <- cbind(as.numeric(mate1[,2]), as.numeric(mate2[,3]))
  coord <- t(apply(coord,1,sort))
  merged.gr <- GenomicRanges::GRanges(seqnames=mate1[,1], ranges=IRanges(start=as.numeric(coord[,1]), end=as.numeric(coord[,2])))
  size.dist <- data.frame(distance = width(merged.gr))
  max <- max(size.dist$distance)
  
  ggplot(size.dist, aes(x=distance)) + geom_histogram(bins = 20) + scale_x_continuous(breaks = round(seq(0, max(size.dist$distance), by = signif(max(size.dist$distance)/5, digits = 1)),1)) + theme_bw()
}


#' Plot inter-chromosomal crosslinks per chromosome
#' 
#' @param linkFile Filename containing link coordinates
#' @param genome Reference genome assembly to use to construct ideogram
#' @importFrom biovizBase getIdeogram
#' @importFrom gridExtra grid.arrange
#' 
#' @author David Porubsky

plotTransLinkCircos <- function(linkFile=NULL, genome='hg38') {
  links <- read.table(linkFile, stringsAsFactors = F, header=T)
  mate1 <- do.call(rbind, strsplit(links$mate1.coord, ","))
  mate2 <- do.call(rbind, strsplit(links$mate2.coord, ","))
  mate1 <- GenomicRanges::GRanges(seqnames=mate1[,1], ranges=IRanges(start=as.numeric(mate1[,2]), end=as.numeric(mate1[,3])))
  mate2 <- GenomicRanges::GRanges(seqnames=mate2[,1], ranges=IRanges(start=as.numeric(mate2[,2]), end=as.numeric(mate2[,3])))
  
  #get ideogram ranges
  #ideo <- GRanges(seqnames=names(chrom.lengths), ranges=IRanges(start=1, end=chrom.lengths))
  hg38Ideogram <- biovizBase::getIdeogram(genome, cytoband = FALSE)
  hg38Ideo <- keepSeqlevels(hg38Ideogram, seqlevels(mate1))
  seqlevels(hg38Ideo) <- gtools::mixedsort(seqlevels(hg38Ideo))
  seqlevels(mate1) <- seqlevels(hg38Ideo)
  seqlengths(mate1) <- seqlengths(hg38Ideo)

  #compare allele present in the read with the reference
  mate1.comp <- rep("none", length(mate1))
  mate1.comp[links$read1.hap == links$read1.hap1] <- "hap1"
  mate1.comp[links$read1.hap == links$read1.hap2] <- "hap2"
  mate2.comp <- rep("none", length(mate2))
  mate2.comp[links$read2.hap == links$read2.hap1] <- "hap1"
  mate2.comp[links$read2.hap == links$read2.hap2] <- "hap2"
  
  mcols(mate1)$mate1.comp <- mate1.comp
  mcols(mate1)$mate2.comp <- mate2.comp
  mcols(mate2)$mate1.comp <- mate1.comp
  mcols(mate2)$mate2.comp <- mate2.comp
  
  #get consensus comparison
  mate2.comp.rev <- mate2.comp
  mate2.comp.rev[mate2.comp=='hap1'] <- 'hap2' 
  mate2.comp.rev[mate2.comp=='hap2'] <- 'hap1'
  comp.consensus <- mate1.comp
  comp.consensus[comp.consensus=='none'] <- mate2.comp.rev[comp.consensus=='none']
  mate1$comp.consensus <- comp.consensus
  
  values(mate1)$to.gr <-  mate2
  mate1 <- mate1[mate1$comp.consensus != 'none']
  
  circos <- list()
  for (chr in seqlevels(hg38Ideo)) {
    message("Plotting chromosome: ", chr)
    chr.gr <- mate1[seqnames(mate1)==chr]
  
    p <- ggplot() + layout_circle(chr.gr, geom = "link", linked.to = "to.gr", radius = 10 ,trackWidth = 0.5, space.skip=0.005, size=0.5, alpha=0.5, aes(color=comp.consensus)) + scale_color_manual(values = c("darkgoldenrod1","deepskyblue4"))
    p <- p + layout_circle(hg38Ideo, geom = "ideo", fill = "gray70", radius = 10.5, trackWidth = 1, space.skip=0.005)
    p <- p + layout_circle(hg38Ideo, geom = "text", aes(label = seqnames), vjust = 0, radius = 11, trackWidth = 1, size=5, space.skip=0.005) + theme(legend.position="none")
    circos[[length(circos)+1]] <- p
  }
  
  n <- length(circos)
  nCol <- floor(sqrt(n))
  pdf("circos.pdf", width=50, height=70)
  do.call("grid.arrange", c(circos, ncol=nCol))
  dev.off()
}  


#' Plot inter-chromosomal crosslinks per chromosomal paur
#' 
#' @param linkFile Filename containing link coordinates
#' @param chromosome.pair Vector containing identifiers of two chromosomes to plot Example: c("chr1", "chr2)
#' 
#' @author David Porubsky

plotPairwiseTransLinks <- function(linkFile=NULL, chromosome.pair=NULL) {

  links <- read.table(linkFile, stringsAsFactors = F, header=T)
  mate1 <- do.call(rbind, strsplit(links$mate1.coord, ","))
  mate2 <- do.call(rbind, strsplit(links$mate2.coord, ","))
  mate1 <- GenomicRanges::GRanges(seqnames=mate1[,1], ranges=IRanges(start=as.numeric(mate1[,2]), end=as.numeric(mate1[,3])))
  mate2 <- GenomicRanges::GRanges(seqnames=mate2[,1], ranges=IRanges(start=as.numeric(mate2[,2]), end=as.numeric(mate2[,3])))
  
  mate1.comp <- rep("none", length(mate1))
  mate1.comp[links$read1.hap == links$read1.hap1] <- "hap1"
  mate1.comp[links$read1.hap == links$read1.hap2] <- "hap2"
  mate2.comp <- rep("none", length(mate2))
  mate2.comp[links$read2.hap == links$read2.hap1] <- "hap1"
  mate2.comp[links$read2.hap == links$read2.hap2] <- "hap2"
  
  #get consensus comparison
  mate2.comp.rev <- mate2.comp
  mate2.comp.rev[mate2.comp=='hap1'] <- 'hap2' 
  mate2.comp.rev[mate2.comp=='hap2'] <- 'hap1'
  comp.consensus <- mate1.comp
  comp.consensus[comp.consensus=='none'] <- mate2.comp.rev[comp.consensus=='none']
  mate1$comp.consensus <- comp.consensus
  values(mate1)$to.gr <-  mate2
  mate1 <- mate1[mate1$comp.consensus != 'none']
  
  plt.gr <- mate1[seqnames(mate1)==chromosome.pair[1] & seqnames(mate1$to.gr) == chromosome.pair[2]]

  df1 <- data.frame(pos=start(plt.gr), level=rep(1, length(plt.gr)), grp=c(1:length(plt.gr)), hap=plt.gr$comp.consensus)
  df2 <- data.frame(pos=start(plt.gr$to.gr), level=rep(5, length(plt.gr$to.gr)), grp=c(1:length(plt.gr$to.gr)), hap=plt.gr$comp.consensus)
  plt.df <- rbind(df1, df2)

  plt <- ggplot(plt.df, aes(pos, level, group = grp, color=hap)) + geom_point(size = 3, shape=1) + geom_line() + scale_color_manual(values = c("darkgoldenrod1","deepskyblue4")) + theme_bw() + scale_y_continuous(breaks = c(1,5), labels = chromosome.pair) + xlab("Genomic position") + theme(axis.title.y=element_blank())

  return(plt)
}
