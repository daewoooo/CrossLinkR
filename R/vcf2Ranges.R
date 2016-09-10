#' Read VCF file into a \code{\link{GRanges}} object
#' 
#' @param vcfFile VCF file containing phased alleles
#' @param genotypeField In case of multiple samples phased in a single file (Default: 1)
#' 
#' @author David Porubsky
#' @export

vcf2rangesSS <- function(vcfFile=NULL, genotypeField=1) {
  vcf <- read.table(vcfFile, stringsAsFactors = F)
  
  genotypeField <- genotypeField + 9
  vcf <- vcf[,c(1:9, genotypeField)]
  
  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\|")
  alleles <- do.call(rbind, alleles)
  score1 <- as.numeric(gen.block[,2])
  score2 <- as.numeric(gen.block[,3])
  phredQ1 <- as.numeric(gen.block[,4])
  phredQ2 <- as.numeric(gen.block[,5])
  
  newCols <- cbind(alleles, score1, phredQ1, score2, phredQ2)
  colnames(newCols) <- c('allele1', 'allele2', 'score1', 'phredQ1', 'score2', 'phredQ2')
  vcf <- cbind(vcf[1:9], newCols)
  
  allele1 <- rep(".", nrow(vcf))
  allele2 <- rep(".", nrow(vcf))   
  allele1 <- ifelse(test = vcf$allele1 == 0, yes = vcf$V4, no = vcf$V5)
  allele2 <- ifelse(test = vcf$allele2 == 0, yes = vcf$V4, no = vcf$V5)
  
  hap.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2), hap1.base=allele1, hap2.base=allele2, score1=score1, phredQ1=phredQ1, score2=score2, phredQ2=phredQ2)
  
  return(hap.gr)
}


#' Read VCF file into a \code{\link{GRanges}} object
#' 
#' @param vcfFile VCF file containing phased alleles
#' @param genotypeField In case of multiple samples phased in a single file (Default: 1)
#' @param filterUnphased Filter out unphased alleles (Default: TRUE)
#' @param minimalBlock In case of multiple phased blocks, size of the smallest block (Default: 10000)
#' 
#' @author David Porubsky
#' @export

vcf2rangesWH <- function(vcfFile=NULL, genotypeField=1, filterUnphased=T, minimalBlock=10000) {
  vcf <- read.table(vcfFile, stringsAsFactors = F)
  
  genotypeField <- genotypeField + 9
  vcf <- vcf[,c(1:9, genotypeField)]
  
  if (filterUnphased) {
    mask <- grepl(pattern = "\\|", vcf[,10])
    vcf <- vcf[mask,]
  }
  
  gen.block <- strsplit(as.character(vcf[,10]),':')
  gen.block <- do.call(rbind, gen.block)
  alleles <- strsplit(gen.block[,1],"\\|")
  alleles <- do.call(rbind, alleles)
  
  if (!is.null(minimalBlock)) {
    newCols <- cbind(alleles, gen.block[,2])
    colnames(newCols) <- c('allele1', 'allele2', 'hapBlock')
    vcf <- cbind(vcf[1:9], newCols)
    mask <- which(table(vcf$hapBlock) < minimalBlock)
    vcf <- vcf[vcf$hapBlock != c(names(mask)),]
  } else {
    newCols <- cbind(alleles)
    colnames(newCols) <- c('allele1', 'allele2')
    vcf <- cbind(vcf[1:9], newCols)
  }
  
  allele1 <- rep(".", nrow(vcf))
  allele2 <- rep(".", nrow(vcf))   
  allele1 <- ifelse(test = vcf$allele1 == 0, yes = vcf$V4, no = vcf$V5)
  allele2 <- ifelse(test = vcf$allele2 == 0, yes = vcf$V4, no = vcf$V5)
  
  if (!is.null(minimalBlock)) {
    hap.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2), hap1.base=allele1, hap2.base=allele2, blockID=vcf$hapBlock)
  } else {
    hap.gr <- GenomicRanges::GRanges(seqnames=vcf$V1, ranges=IRanges(start=vcf$V2, end=vcf$V2), hap1.base=allele1, hap2.base=allele2)
  }  
  return(hap.gr)
} 
