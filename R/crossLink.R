#' Wrapper function for CrossLinkR package
#'
#' @param bamFile Bamfile with paired-end reads sorted by readName.
#' @param vcfFile VCF file containing phased alleles
#' @param outputfolder Output directory. If non-existent it will be created.
#' @param yieldSize Number of reads to be loaded as a chunk to save memory (Default: 1000000)
#' @param readLen Length of the read in bp.
#' @param min.mapq Minimum mapping quality of any given read.
#' @param min.baseq Minimum base quality to consider SNV for crosslink.
#' @param nCPU The number of CPUs to use. Should not be more than available on your machine.
#' @param blacklist List of regions to be excluded from the analysis (format: chromosome start end)
#' @importFrom Rsamtools BamFile ScanBamParam scanBamFlag isOpen
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom BiocParallel MulticoreParam bpiterate
#' @importFrom data.table rbindlist fread
#' @author David Porubsky
#' @export

crossLinkR <- function(bamFile, vcfFile, outputfolder='./CrossLinkR_analysis', yieldSize = 100000, readLen=100, min.mapq=10, min.baseq=20, nCPU=3, blacklist=NULL) {

bamIterator <- function(bf) {
  done <- FALSE
  if (!Rsamtools::isOpen(bf))
    open(bf)
  function() {
    if (done)
      return(NULL)
    yld <- GenomicAlignments::readGAlignments(bf, param=Rsamtools::ScanBamParam(what=c('seq', 'qual','mapq','cigar'), flag=Rsamtools::scanBamFlag(isDuplicate=F)))
    if (length(yld) == 0L) {
      close(bf)
      done <<- TRUE
      NULL
    } else yld
  }
}

message("Importing VCF file... ")
ptm <- proc.time()

#suppressWarnings( hap.gr <- vcf2rangesWH(vcfFile = vcfFile, genotypeField = 1, filterUnphased = T, minimalBlock = NULL) )
suppressWarnings( hap.gr <- vcf2rangesSS(vcfFile = vcfFile, genotypeField = 1) )

#filter out missing alleles and homozygous positions
mask <- hap.gr$hap1.base == "N" | hap.gr$hap2.base == "N" | hap.gr$hap1.base == hap.gr$hap2.base
hap.gr <- hap.gr[!mask]

#filter low quality (score < 1; phredQ > 0.99)
#mask <- hap.gr$score1 >= 1 & hap.gr$score2 >= 1 & hap.gr$phredQ1 > 0.99 & hap.gr$phredQ2 > 0.99
mask <- hap.gr$phredQ1 >= 0.99 & hap.gr$phredQ2 >= 0.99
hap.gr <- hap.gr[mask]

#blacklist regions with mixed directional reads (inversion) as well as segmental duplications
if (!is.null(blacklist)) {
  filt.regions <- read.table(blacklist, header=T)
  filt.regions <- GenomicRanges::GRanges(seqnames=filt.regions$seqnames, ranges=IRanges(start=filt.regions$start, end=filt.regions$end))
  suppressWarnings( filt.hits <- GenomicRanges::findOverlaps(hap.gr, filt.regions) )
  if (length(filt.hits)>0) {
    hap.gr <- hap.gr[-queryHits(filt.hits)]
  }  
}

time <- proc.time() - ptm
message("Time spent: ",round(time[3],2),"s")

#Create outputdirectory if it doesn't exist
if (!file.exists(outputfolder)) {
  dir.create(outputfolder)
}

#Create TMP directory if it doesn't exist
tmpDIR <- file.path(outputfolder, "TMP")
if (!file.exists(tmpDIR)) {
  dir.create(tmpDIR)
}

#remove files for export if the file exists
prefix <- gsub("\\.bam", "", basename(bamFile))
concordant.first <- file.path(outputfolder, paste0(prefix, "_concordant_first.fasta"))
concordant.last <- file.path(outputfolder, paste0(prefix,"_concordant_last.fasta"))
discordant.first <- file.path(outputfolder, paste0(prefix,"_discordant_first.fasta"))
discordant.last <- file.path(outputfolder, paste0(prefix,"_discordant_last.fasta"))
read.links.trans <- file.path(outputfolder, paste0(prefix, "_readLinks_transChrom.txt"))
read.links <- file.path(outputfolder, paste0(prefix,"_readLinks.txt"))

if (file.exists(concordant.first)) file.remove(concordant.first)
if (file.exists(concordant.last)) file.remove(concordant.last)
if (file.exists(discordant.first)) file.remove(discordant.first)
if (file.exists(discordant.last)) file.remove(discordant.last)
if (file.exists(read.links.trans)) file.remove(read.links.trans)
if (file.exists(read.links)) file.remove(read.links)

message("Searching for crosslinks... ")
ptm <- proc.time()

file.header <- Rsamtools::scanBamHeader(bamFile)[[1]]
chrom.lengths <- file.header$targets

bf <- Rsamtools::BamFile(bamFile, yieldSize = yieldSize, obeyQname=TRUE) #open bamfile instance
ITER <- bamIterator(bf) #initialize bamfile iterator
bpparam <- BiocParallel::MulticoreParam(workers = nCPU) #set number of parallel operations
#tC <- tryCatch({
  BiocParallel::bpiterate(ITER, getLinks, hap.gr = hap.gr, outputfolder=outputfolder, tmpDIR=tmpDIR, readLen=readLen, min.mapq=min.mapq, min.baseq=min.baseq, ID=prefix, BPPARAM = bpparam)
#}, error = function(err) {
#  stop(file,'\n',err)
#})

time <- proc.time() - ptm
message("Time spent: ",round(time[3],2),"s") 

#concatenate results from tmpDIR
file.list <- list.files(tmpDIR)
links <- data.table::rbindlist(lapply(file.list, function(x) data.table::fread(file.path(tmpDIR,x), sep = " ", colClasses = c("character","character","numeric","character","character","character","character","character","character","character","character","character"))))
colnames(links) <- c('mate1.coord', 'mate2.coord', 'snv.cov', 'read1.snvs', 'read2.snvs', 'read1.hap', 'read2.hap', 'read1.hap1', 'read2.hap1', 'read1.hap2', 'read2.hap2', 'comparison')
write.table(links, file = read.links, row.names = F, col.names = T, quote = F)
unlink(tmpDIR, recursive = T)

#export discordant pairs into separate files
#links <- read.table(read.links, stringsAsFactors = F, header=T)
links <- links[links$comparison == 'Discordant',]
mate1 <- do.call(rbind, strsplit(links$mate1.coord, ","))
mate2 <- do.call(rbind, strsplit(links$mate2.coord, ","))
mask <- mate1[,1] == mate2[,1]
#print discordant links on the same chromosome
cis.links.discordant <- file.path(outputfolder, paste0(prefix,"_CisLinks_discordant.txt"))
write.table(x = links[mask,], file = cis.links.discordant, quote = F, row.names = F)
#print discordant links on different chromosomes
trans.links.discordant <- file.path(outputfolder, paste0(prefix,"_TransLinks_discordant.txt"))
write.table(x = links[!mask,], file = trans.links.discordant, quote = F, row.names = F)

#search for crosslink hotspots
message("Searching for crosslink hotspots... ")
ptm <- proc.time()

hotspots <- hotLinks(crosslinks=cis.links.discordant, maxDist=1000, minCov=4, pval=NULL, chrom.lengths=chrom.lengths)

crossHotspots <- file.path(outputfolder, paste0(prefix,"_crossHotspots.txt"))
hot.df <- as(hotspots, "data.frame")
write.table(hot.df, file = crossHotspots, quote = F, row.names = F, col.names = F)

time <- proc.time() - ptm
message("Time spent: ",round(time[3],2),"s")

#Annotate hotspots for known genomic features
#message("Annotating hotspots... ")
#ptm <- proc.time()

#annot.hotspots <- annotateRanges(searchGRanges = hotspots, maxDist = 1000, blacklist = NULL, DNAseI="DNAse_clustered_hg38", RegulE="OregAnno.bed", G4motif="g4motif_GRCH38.bed")
#venn.df <- as.data.frame(mcols(annot.hotspots)[,2:4])
#library(limma)
#venn.count <- vennCounts(venn.df)
#vennDiagram(venn.count)

#time <- proc.time() - ptm
#message("Time spent: ",round(time[3],2),"s")

}


