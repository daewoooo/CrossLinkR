#' This function concatenates strings stored in the list
#' 
#' @param x List of character vectors
#' @param sep Separator to concatenate character vectors 
#' 
#' @author David Porubsky
#' @export

collapse.str <- function(x, sep="") {
  if (length(x) > 1) {
    paste(x, collapse=sep)
  } else if (length(x) == 1) {
    as.character(x)
  } else {
    x <- ""
  }
}


#' Search for links in Hi-C read pairs
#'
#' @param reads An \code{\link{GRanges}} object containing read-pairs sorted by name.
#' @param hap.gr An \code{\link{GRanges}} object containing phased SNVs used as a reference
#' @param readLen Length of the read in bp.
#' @param min.mapq Minimum mapping quality of any given read.
#' @param min.baseq Minimum base quality to consider SNV for crosslink.
#' @param ID Unique identifier for the analysis. (Default: Bamfile name)
#' @importFrom Biostrings writeXStringSet
#' @importFrom GenomicAlignments pileLettersAt
#' @author David Porubsky
#' @export

getLinks <- function(reads, hap.gr, outputfolder, tmpDIR, readLen=100, min.mapq=10, min.baseq=20, ID="") {
  
  makeRandomName <- function(lenght=15) {
    randomString <- paste(sample(c(0:9, letters, LETTERS),lenght, replace=TRUE), collapse="")
    return(randomString)
  }
  
  #create filenames
  concordant.first <- file.path(outputfolder, paste0(ID, "_concordant_first.fasta"))
  concordant.last <- file.path(outputfolder, paste0(ID,"_concordant_last.fasta"))
  discordant.first <- file.path(outputfolder, paste0(ID,"_discordant_first.fasta"))
  discordant.last <- file.path(outputfolder, paste0(ID,"_discordant_last.fasta"))
  read.links.trans <- file.path(outputfolder, paste0(ID, "_readLinks_transChrom.txt"))
  #read.links <- file.path(outputfolder, paste0(ID,"_readLinks.txt"))
  read.links <- file.path(tmpDIR, makeRandomName())
  
  reads <- as(reads, 'GRanges')
  
  #split paired-end reads sorted by name to first and last mate
  data.first <- reads[seq(from=1, to=length(reads), by=2)]
  data.last <- reads[seq(from=2, to=length(reads), by=2)]
  
  #filter reads longer than the used read lenght and by mapping quality
  data.first.filt <- mcols(data.first)$mapq >= min.mapq & width(data.first) <= readLen
  data.last.filt <- mcols(data.last)$mapq >= min.mapq & width(data.last) <= readLen
  mask <- data.first.filt & data.last.filt
  data.first <- data.first[mask]
  data.last <- data.last[mask]
  
  #filter reads that do not overlap with reference SNVs
  hits.first <- findOverlaps(hap.gr, data.first)
  hits.last <- findOverlaps(hap.gr, data.last)
  mask <- intersect(subjectHits(hits.first), subjectHits(hits.last))
  data.first <- data.first[mask]
  data.last <- data.last[mask]
  
  #filter reads mapping to different chromosomes
  #mask <- as(seqnames(data.first), "vector") == as(seqnames(data.last), "vector")
  #data.first.trans <- data.first[!mask]
  #data.last.trans <- data.last[!mask]
  #data.first <- data.first[mask]
  #data.last <- data.last[mask]
  
  #print reads mapping to different chroms
  #data.first.trans <- as(data.first.trans[,0], "data.frame")
  #data.last.trans <- as(data.last.trans[,0], "data.frame")
  #trans <- cbind(data.first.trans, data.last.trans)
  #write.table(trans, file = read.links.trans, row.names = F, col.names = F, quote = F, append = T)
  
  #make sure bot GRanges object uses the same chromosome names
  seqlevels(data.first) <- seqlevels(hap.gr)
  seqlevels(data.last) <- seqlevels(hap.gr)
  
  while (length(data.first)>0) {
    
    #filter out duplicated reads
    dedup.first <- which(duplicated(start(data.first))==TRUE)
    dedup.last <- which(duplicated(start(data.last))==TRUE)
    
    #get indices of all non-overlapping reads for first and second mate
    disjoin.first.idx <- which(disjointBins(data.first)==1)
    disjoin.last.idx <- which(disjointBins(data.last)==1)
    #select all corresponding non-overlapping reads for first and last mate
    mask <- intersect(disjoin.first.idx,disjoin.last.idx)
    
    #remove indices of duplicates from mask
    mask <- setdiff(mask, dedup.first)
    mask <- setdiff(mask, dedup.last)
    
    if (length(mask)==0) {
      disjoin.first <- data.first[1]
      disjoin.last <- data.last[1]
      #subtract selected non-overlapping reads from the original objects
      data.first <- data.first[-1]
      data.last <- data.last[-1]
    } else {
      disjoin.first <- data.first[mask]
      disjoin.last <- data.last[mask]
      #subtract selected non-overlapping reads from the original objects
      data.first <- data.first[-mask]
      data.last <- data.last[-mask]
    }
    
    #Select overlapping SNVs with selected reads for first mate
    hit.first <- GenomicRanges::findOverlaps(hap.gr, disjoin.first)
    first.snvs <- hap.gr[queryHits(hit.first)]
    
    #Pile bases and base qualities at selected SNVs for first mate
    piles1 <- GenomicAlignments::pileLettersAt(disjoin.first$seq, seqnames(disjoin.first), start(disjoin.first), disjoin.first$cigar, first.snvs)
    quals1 <- GenomicAlignments::pileLettersAt(disjoin.first$qual, seqnames(disjoin.first), start(disjoin.first), disjoin.first$cigar, first.snvs)
    #df.piles1 <- as(piles1, "data.frame")
    df.quals1 <- as(quals1, "data.frame") 
    
    #filter bases based on base quality
    quals1 <- sapply(df.quals1$x, function(x) as.numeric(charToRaw(x))-33)
    quals1 <- unlist(quals1)
    piles1 <- as(piles1, "vector")
    #filt.piles1 <- mapply(function(X,Y) { X[Y >= min.baseq] }, X=df.piles1$x, Y=quals1)
    #piles1 <- unlist(lapply(filt.piles1, collapse.str))
    piles1[!quals1 >= min.baseq] <- ""
    
    #subset SNVs covered in each single read and join them in a single string
    piles1.alleles <- split(piles1, subjectHits(hit.first))
    piles1.alleles <- lapply(piles1.alleles, collapse.str)
    disjoin.first <- disjoin.first[,1]
    mcols(disjoin.first)$snvs <- unlist(piles1.alleles)
    
    #filter reference bases that were previously filterted by baseq
    first.snvs.hap1 <- first.snvs$hap1.base
    first.snvs.hap2 <- first.snvs$hap2.base
    first.snvs.hap1[!quals1>=min.baseq] <- ""
    first.snvs.hap2[!quals1>=min.baseq] <- ""
    
    #get corresponding reference alleles to the covered SNVs for each read
    #first.hap1.base <- split(first.snvs$hap1.base, subjectHits(hit.first))
    #first.hap2.base <- split(first.snvs$hap2.base, subjectHits(hit.first))
    first.hap1.base <- split(first.snvs.hap1, subjectHits(hit.first))
    first.hap2.base <- split(first.snvs.hap2, subjectHits(hit.first))
    first.hap1.base <- lapply(first.hap1.base, collapse.str)
    first.hap2.base <- lapply(first.hap2.base, collapse.str)
    mcols(disjoin.first)$hap1.ref <- unlist(first.hap1.base)
    mcols(disjoin.first)$hap2.ref <- unlist(first.hap2.base)
    
    #get positions of covered SNVs
    first.ref.snvs <- split(start(first.snvs), subjectHits(hit.first))
    read1.snv.cov <- sapply(first.ref.snvs, length) #get number of covered SNVs for given read
    first.ref.snvs <- lapply(first.ref.snvs, function(x) collapse.str(x, sep=","))
    mcols(disjoin.first)$ref.pos <- unlist(first.ref.snvs)
    
    #Processing steps taken for first mate will be repeated for last mate
    
    #Select overlapping SNVs with selected reads for last mate
    hit.last <- GenomicRanges::findOverlaps(hap.gr, disjoin.last)
    last.snvs <- hap.gr[queryHits(hit.last)]
    
    #Pile bases and base qualities at selected SNVs for last mate
    piles2 <- GenomicAlignments::pileLettersAt(disjoin.last$seq, seqnames(disjoin.last), start(disjoin.last), disjoin.last$cigar, last.snvs)
    quals2 <- GenomicAlignments::pileLettersAt(disjoin.last$qual, seqnames(disjoin.last), start(disjoin.last), disjoin.last$cigar, last.snvs)
    #df.piles2 <- as(piles2, "data.frame")
    df.quals2 <- as(quals2, "data.frame") 
    
    #filter bases based on base quality
    quals2 <- sapply(df.quals2$x, function(x) as.numeric(charToRaw(x))-33)
    quals2 <- unlist(quals2)
    piles2 <- as(piles2, "vector")
    #filt.piles2 <- mapply(function(X,Y) { X[Y >= min.baseq] }, X=df.piles2$x, Y=quals2)
    #piles2 <- unlist(lapply(filt.piles2, collapse.str))
    piles2[!quals2 >= min.baseq] <- ""
    
    #subset SNVs covered in each single read and join them in a single string
    piles2.alleles <- split(piles2, subjectHits(hit.last))
    piles2.alleles <- lapply(piles2.alleles, collapse.str)
    disjoin.last <- disjoin.last[,1]
    mcols(disjoin.last)$snvs <- unlist(piles2.alleles)
    
    #filter reference bases that were previously filterted by baseq
    last.snvs.hap1 <- last.snvs$hap1.base
    last.snvs.hap2 <- last.snvs$hap2.base
    last.snvs.hap1[!quals2>=min.baseq] <- ""
    last.snvs.hap2[!quals2>=min.baseq] <- ""
    
    #get corresponding reference alleles to the covered SNVs for each read
    #last.hap1.base <- split(last.snvs$hap1.base, subjectHits(hit.last))
    #last.hap2.base <- split(last.snvs$hap2.base, subjectHits(hit.last))
    last.hap1.base <- split(last.snvs.hap1, subjectHits(hit.last))
    last.hap2.base <- split(last.snvs.hap2, subjectHits(hit.last))
    last.hap1.base <- lapply(last.hap1.base, collapse.str)
    last.hap2.base <- lapply(last.hap2.base, collapse.str)
    mcols(disjoin.last)$hap1.ref <- unlist(last.hap1.base)
    mcols(disjoin.last)$hap2.ref <- unlist(last.hap2.base)
    
    #get positions of covered SNVs
    last.ref.snvs <- split(start(last.snvs), subjectHits(hit.last))
    read2.snv.cov <- sapply(last.ref.snvs, length) #get number of covered SNVs for given read
    last.ref.snvs <- lapply(last.ref.snvs, function(x) collapse.str(x, sep=","))
    mcols(disjoin.last)$ref.pos <- unlist(last.ref.snvs)
    
    #remove read pairs with SNVs filtered by baseq
    mask <- which(disjoin.first$snvs != "" & disjoin.last$snvs != "")
    disjoin.first <- disjoin.first[mask]
    disjoin.last <- disjoin.last[mask]
    #get number of covered SNVs for each read pair
    snv.cov <- (read1.snv.cov + read2.snv.cov)[mask]		
    
    #compare covered SNVs with the reference SNVs
    concord.hap1.idx <- which(disjoin.first$snvs == disjoin.first$hap1.ref & disjoin.last$snvs == disjoin.last$hap1.ref)
    concord.hap2.idx <- which(disjoin.first$snvs == disjoin.first$hap2.ref & disjoin.last$snvs == disjoin.last$hap2.ref)
    mcols(disjoin.first)$comparison <- rep('Discordant', length(disjoin.first))
    mcols(disjoin.last)$comparison <- rep('Discordant', length(disjoin.last))
    if (length(concord.hap1.idx)>0) {
      disjoin.first[concord.hap1.idx]$comparison <- 'Concordant.hap1'
      disjoin.last[concord.hap1.idx]$comparison <- 'Concordant.hap1'
    }
    if (length(concord.hap2.idx)>0) {
      disjoin.first[concord.hap2.idx]$comparison <- 'Concordant.hap2'
      disjoin.last[concord.hap2.idx]$comparison <- 'Concordant.hap2'
    }
    
    #print results
    first.coord <-  paste(seqnames(disjoin.first), start(disjoin.first), end(disjoin.first), sep="," ) #select only chrName startPos and endPos
    last.coord <-  paste(seqnames(disjoin.last), start(disjoin.last), end(disjoin.last), sep="," )
    
    #print concordant reads
    reads.first.concord <- disjoin.first$seq[disjoin.first$comparison != 'Discordant']
    names(reads.first.concord) <- first.coord[disjoin.first$comparison != 'Discordant']
    Biostrings::writeXStringSet(reads.first.concord, filepath = concordant.first, format = "fasta", append = T)
    reads.last.concord <- disjoin.last$seq[disjoin.last$comparison != 'Discordant']
    names(reads.last.concord) <- last.coord[disjoin.last$comparison != 'Discordant']
    Biostrings::writeXStringSet(reads.last.concord, filepath = concordant.last, format = "fasta", append = T)
    
    #print discordant reads
    reads.first.discord <- disjoin.first$seq[disjoin.first$comparison == 'Discordant']
    names(reads.first.discord) <- first.coord[disjoin.first$comparison == 'Discordant']
    Biostrings::writeXStringSet(reads.first.discord, filepath = discordant.first, format = "fasta", append = T)
    reads.last.discord <- disjoin.last$seq[disjoin.last$comparison == 'Discordant']
    names(reads.last.discord) <- last.coord[disjoin.last$comparison == 'Discordant']
    Biostrings::writeXStringSet(reads.last.discord, filepath = discordant.last, format = "fasta", append = T)
    
    #print alleles
    read1.snvs <- disjoin.first$ref.pos
    read2.snvs <- disjoin.last$ref.pos
    read1.hap <- disjoin.first$snvs
    read2.hap <- disjoin.last$snvs
    read1.hap1 <- disjoin.first$hap1.ref
    read2.hap1 <- disjoin.last$hap1.ref
    read1.hap2 <- disjoin.first$hap2.ref
    read2.hap2 <- disjoin.last$hap2.ref
    
    print.df <- data.frame(mate1.coord=first.coord, mate2.coord=last.coord, snv.cov=snv.cov, read1.snvs=read1.snvs, read2.snvs=read2.snvs, read1.hap=read1.hap, read2.hap=read2.hap, read1.hap1=read1.hap1, read2.hap1=read2.hap1, read1.hap2=read1.hap2, read2.hap2=read2.hap2, comparison=disjoin.first$comparison)
    write.table(print.df, file = read.links, row.names = F, col.names = F, quote = F, append = T)
  } 
}
