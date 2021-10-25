##################
#
# Frederick National Laboratory for Cancer Research
# NIAID Collaborative Bioinformatics Resource (NCBR), 
# National Institute of Allergy and Infectious Diseases, NIH
# Tovah Markowitz, tovah.markowitz@nih.gov
# August 27, 2021
#
##################
### This script is designed to take data from two bed files and a file of chromosome lengths.
# It then uses regioneR to calculate the number of overlaps between the the two bed files
# and randomly redistributes the regions to calculate a p-value.

## PARAMETERS
# bedFile1: the name of the first bed file, this is the file whose data will be manipulated for bootstrapping
# bedFile2: the name of the bed file to be compared to
# chrFile: a two column file of chromosome and length that matches the bed files without a header 
# addFlanks: [optional] number of bases to add to each side of the regions in bedFile1 (default: 0)

## EXAMPLE
# RegioneR_overlap_analysis("enhancer.bed", "HPV_integrants_cervical_ALL_2.bed", "hg19.len", addFlanks=1e5/2)

RegioneR_overlap_analysis <- function(bedFile1, bedFile2, chrFile, addFlanks=0) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    BiocManager::install("GenomicRanges")
  }

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    BiocManager::install("rtracklayer")
  }

  if (!requireNamespace("regioneR", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    BiocManager::install("regioneR")
  }

  library(regioneR)

  chrs <- read.table(chrFile)
  #chrs <- chrs[1:23,]

  bed1 <- sort(rtracklayer::import(bedFile1))
  bed1 <- extendRegions(bed1, extend.start=addFlanks, extend.end=addFlanks)

  bed2 <- sort(rtracklayer::import(bedFile2))

  pt <- permTest(A=bed1, B=bed2, ntimes=10000, genome=chrs,
                 randomize.function=randomizeRegions, allow.overlaps=TRUE,
                 evaluate.function=numOverlaps,count.once=TRUE,
                 alternative="auto")

  return(pt)
}
