##################
#
# Frederick National Laboratory for Cancer Research
# NIAID Collaborative Bioinformatics Resource (NCBR), 
# National Institute of Allergy and Infectious Diseases, NIH
# Tovah Markowitz, tovah.markowitz@nih.gov
# December 9, 2020
#
##################
# Purpose: This script is designed to bootstrap integration breakpoints 
#   to determine if they significantly overlap enhancers or FANCD2 
#   binding sites.
#
# Inputs: 1) Bed files of FANCD2 or enhancer sites. 
#         2) A two column file of chromosome lengths
#         3) An Excel file containing Supplemental Tables S2 and S3.
# Output: a table containing the information in Tables 1 and 2.
#
# Notes: 1) For analysis purposes, breakpoints were treated as 2bp
#   positions with start as the actual breakpoint. 
#        2) Chromosome Y was not included in this analysis.

##################
## PACKAGES
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  	install.packages("BiocManager") }
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  	install.packages("BiocManager") }
  BiocManager::install("rtracklayer")
}

if (!requireNamespace("regioneR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  	install.packages("BiocManager") }
  BiocManager::install("regioneR")
}

if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages('readxl')
}

library(regioneR)
library(GenomicRanges)

##################
## SUBFUNCTIONS

prep_file <- function(inFile, inSheet, addFlanks=0) {
  # Steps include importing the data as a GenomicRanges object, 
  # sorting it, and extending the regions.
  tmp <- data.frame( readxl::read_xlsx(inFile, sheet=inSheet) )
  tmp2 <- data.frame(Patient.ID=tmp$Patient.ID, 
  					 chr=tmp$Chromosome,
                     start=tmp$Integration.Breakpoint, 
                     end=tmp$Integration.Breakpoint +1,
                     Integration.ID=tmp$Integration.ID..,
                     N.breakpoints=as.numeric( gsub("Breakpoints ",  "", 
                                        tmp$X..Breakpoints.per.integration.site) )
                     )
  fileData <- makeGRangesFromDataFrame(tmp2, keep.extra.columns=T)
  fileData <- sort(fileData)
  if (addFlanks != 0) {
    fileData <- extendRegions(fileData, extend.start=addFlanks, 
    						  extend.end=addFlanks)
  }
  return(fileData)
}

collapse_cluster <- function(inGRobj) {
  # Takes a GR object and collapses breakpoints that were previously
  # defined as being together
  clusters <- unique(inGRobj$Integration.ID)
  inDF <- data.frame(inGRobj)
  outDF <- data.frame(seqnames=character( length(clusters) ), 
  					  start=numeric( length(clusters) ),
                      end=numeric( length(clusters) ), 
                      Integration.ID=clusters, 
                      N.breakpoints=numeric( length(clusters) ), 
                      stringsAsFactors=F)
  
  for( i in 1:length(clusters) ) {
    tmp <- inDF[which(inDF$Integration.ID == clusters[i]), ]
    if ( nrow(tmp) == tmp$N.breakpoints[1] ) {
      if ( length( unique(tmp$seqnames) ) == 1 ) {
        outDF$N.breakpoints[i] <- tmp$N.breakpoints[1]
        outDF$seqnames[i] <- as.character(tmp$seqnames[1])
        outDF$start[i] <- min(tmp$start)
        outDF$end[i] <- max(tmp$end)
      } 
    } 
  }
  
  outGR <- makeGRangesFromDataFrame(outDF, keep.extra.columns=T)
  outGR <- sort(outGR)
  return(outGR)
}

##################
## MAIN

inChr <- "hg19.len"
inEnhancerBed <- "Intersect_Merged_Brd4_BroadPeaks_All-W12-subclones-H3K27ac-consensusPeaks.bed"
inFANCD2bed <- "Merged C33A-HeLa ChIP-seq_ChIP-ChIP_COMBINED.bed"
inBreakpoints <- "Updated CESC and HNSCC integration datasets_Dec 2020.xlsx"
CESCsheet <- 1
HNSCCsheet <- 2

# read in chromosome length file, exclude chromosome Y
chrs <- read.table(inChr)
chrs <- chrs[1:23,]

# load overlap feature
enhancers <- sort( rtracklayer::import(inEnhancerBed) )
FANCD2 <- sort( rtracklayer::import(inFANCD2bed) )

# prepare breakpoint or loci datasets, using a 50 kb flank size
CESC <- prep_file(inBreakpoints, CESCsheet, addFlanks=5e4)
CESC2 <- CESC[which(CESC$N.breakpoints != 1)]
CESC3 <- collapse_cluster(CESC2)
CESC2B <- CESC[which(CESC$N.breakpoints == 1)]
CESC4 <- sort( c(CESC2B, CESC3) )

HNSC <- prep_file(inBreakpoints, HNSCCsheet, addFlanks=5e4)
HNSC <- HNSC[which(seqnames(HNSC) != "chrY")]
HNSC2 <- HNSC[which(HNSC$N.breakpoints != 1)]
HNSC3 <- collapse_cluster(HNSC2)
HNSC2B <- HNSC[which(HNSC$N.breakpoints == 1)]
HNSC4 <- sort( c(HNSC2B, HNSC3) )

# set up output object
alldata <- data.frame(Overlap.feature=c(rep("Enhancers",10),rep("FANCD2 sites",10)),
			Dataset=rep(c(rep("CESC",5),rep("HNSCC",5)),2),
            Integration.subgroups=rep(c("All breakpoints",
             						    "Single breakpoints",
                   					    "Clustered breakpoints",
                   						"Condensed clustered breakpoints",
                   					    "Combined single and condensed breakpoints"),2),
            Num.Loci=numeric(10),
            Num.Overlaps=numeric(10), 
            Pct.Overlap=numeric(10), 
            P.value=numeric(10) )

# run regioneR for each combination chosen for the output object
# Row 1: Enhancers versus CESC all breakpoints
pt <- permTest(A=CESC, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[1] <- length(CESC)
alldata$Num.Overlaps[1] <- pt$numOverlaps$observed
alldata$Pct.Overlap[1] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[1]
alldata$P.value[1] <- pt$numOverlaps$pval

# Row 2: Enhancers versus CESC single breakpoints
pt <- permTest(A=CESC2B, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[2] <- length(CESC2B)
alldata$Num.Overlaps[2] <- pt$numOverlaps$observed
alldata$Pct.Overlap[2] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[2]
alldata$P.value[2] <- pt$numOverlaps$pval

# Row 3: Enhancers versus CESC clustered breakpoints
pt <- permTest(A=CESC2, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[3] <- length(CESC2)
alldata$Num.Overlaps[3] <- pt$numOverlaps$observed
alldata$Pct.Overlap[3] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[3]
alldata$P.value[3] <- pt$numOverlaps$pval

# Row 4: Enhancers versus CESC condensed clustered breakpoints
pt <- permTest(A=CESC4, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[4] <- length(CESC4)
alldata$Num.Overlaps[4] <- pt$numOverlaps$observed
alldata$Pct.Overlap[4] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[4]
alldata$P.value[4] <- pt$numOverlaps$pval

# Row 5: Enhancers versus CESC combined single and condensed breakpoints
pt <- permTest(A=CESC3, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[5] <- length(CESC3)
alldata$Num.Overlaps[5] <- pt$numOverlaps$observed
alldata$Pct.Overlap[5] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[5]
alldata$P.value[5] <- pt$numOverlaps$pval

# Row 6: Enhancers versus HNSCC all breakpoints
pt <- permTest(A=HNSC, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[6] <- length(HNSC)
alldata$Num.Overlaps[6] <- pt$numOverlaps$observed
alldata$Pct.Overlap[6] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[6]
alldata$P.value[6] <- pt$numOverlaps$pval

# Row 7: Enhancers versus HNSCC single breakpoints
pt <- permTest(A=HNSC2B, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[7] <- length(HNSC2B)
alldata$Num.Overlaps[7] <- pt$numOverlaps$observed
alldata$Pct.Overlap[7] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[7]
alldata$P.value[7] <- pt$numOverlaps$pval

# Row 8: Enhancers versus HNSCC clustered breakpoints
pt <- permTest(A=HNSC2, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[8] <- length(HNSC2)
alldata$Num.Overlaps[8] <- pt$numOverlaps$observed
alldata$Pct.Overlap[8] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[8]
alldata$P.value[8] <- pt$numOverlaps$pval

# Row 9: Enhancers versus HNSCC condensed clustered breakpoints
pt <- permTest(A=HNSC4, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[9] <- length(HNSC4)
alldata$Num.Overlaps[9] <- pt$numOverlaps$observed
alldata$Pct.Overlap[9] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[9]
alldata$P.value[9] <- pt$numOverlaps$pval

# Row 10: Enhancers versus HNSCC combined single and condensed breakpoints
pt <- permTest(A=HNSC3, B=enhancers, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[10] <- length(HNSC3)
alldata$Num.Overlaps[10] <- pt$numOverlaps$observed
alldata$Pct.Overlap[10] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[10]
alldata$P.value[10] <- pt$numOverlaps$pval

# Row 11: FANCD2 versus CESC all breakpoints
pt <- permTest(A=CESC, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[11] <- length(CESC)
alldata$Num.Overlaps[11] <- pt$numOverlaps$observed
alldata$Pct.Overlap[11] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[11]
alldata$P.value[11] <- pt$numOverlaps$pval

# Row 12: FANCD2 versus CESC single breakpoints
pt <- permTest(A=CESC2B, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[12] <- length(CESC2B)
alldata$Num.Overlaps[12] <- pt$numOverlaps$observed
alldata$Pct.Overlap[12] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[12]
alldata$P.value[12] <- pt$numOverlaps$pval

# Row 13: FANCD2 versus CESC clustered breakpoints
pt <- permTest(A=CESC2, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[13] <- length(CESC2)
alldata$Num.Overlaps[13] <- pt$numOverlaps$observed
alldata$Pct.Overlap[13] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[13]
alldata$P.value[13] <- pt$numOverlaps$pval

# Row 14: FANCD2 versus CESC condensed clustered breakpoints
pt <- permTest(A=CESC4, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[14] <- length(CESC4)
alldata$Num.Overlaps[14] <- pt$numOverlaps$observed
alldata$Pct.Overlap[14] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[14]
alldata$P.value[14] <- pt$numOverlaps$pval

# Row 15: FANCD2 versus CESC combined single and condensed breakpoints
pt <- permTest(A=CESC3, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[15] <- length(CESC3)
alldata$Num.Overlaps[15] <- pt$numOverlaps$observed
alldata$Pct.Overlap[15] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[15]
alldata$P.value[15] <- pt$numOverlaps$pval

# Row 16: FANCD2 versus HNSCC all breakpoints
pt <- permTest(A=HNSC, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[16] <- length(HNSC)
alldata$Num.Overlaps[16] <- pt$numOverlaps$observed
alldata$Pct.Overlap[16] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[16]
alldata$P.value[16] <- pt$numOverlaps$pval

# Row 17: FANCD2 versus HNSCC single breakpoints
pt <- permTest(A=HNSC2B, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[17] <- length(HNSC2B)
alldata$Num.Overlaps[17] <- pt$numOverlaps$observed
alldata$Pct.Overlap[17] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[17]
alldata$P.value[17] <- pt$numOverlaps$pval

# Row 18: FANCD2 versus HNSCC clustered breakpoints
pt <- permTest(A=HNSC2, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[18] <- length(HNSC2)
alldata$Num.Overlaps[18] <- pt$numOverlaps$observed
alldata$Pct.Overlap[18] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[18]
alldata$P.value[18] <- pt$numOverlaps$pval

# Row 19: FANCD2 versus HNSCC condensed clustered breakpoints
pt <- permTest(A=HNSC4, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[19] <- length(HNSC4)
alldata$Num.Overlaps[19] <- pt$numOverlaps$observed
alldata$Pct.Overlap[19] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[19]
alldata$P.value[19] <- pt$numOverlaps$pval

# Row 20: FANCD2 versus HNSCC combined single and condensed breakpoints
pt <- permTest(A=HNSC3, B=FANCD2, ntimes=10000, genome=chrs, 
               randomize.function=randomizeRegions, allow.overlaps=TRUE,
               evaluate.function=numOverlaps, count.once=T,
               alternative="auto")

alldata$Num.loci[20] <- length(HNSC3)
alldata$Num.Overlaps[20] <- pt$numOverlaps$observed
alldata$Pct.Overlap[20] <- 100 * pt$numOverlaps$observed / alldata$Num.loci[10]
alldata$P.value[20] <- pt$numOverlaps$pval

print(alldata)