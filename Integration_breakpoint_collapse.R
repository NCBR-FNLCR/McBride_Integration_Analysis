##################
#
# Frederick National Laboratory (NCBR)
# Tovah Markowitz, tovah.markowitz@nih.gov
# April 19, 2021
#
##################
### This script is designed to take data from an Excel file 
# whose first three columns are chromosome, start, and end (or a bed file). 
# It will then collapse rows together based upon the distance between the End value
# of one row and the Start value of the previous row. It will produce a bed-like
# txt file as output with the collapsed region and the number of rows that were 
# combined.

## PARAMETERS
# inFile: the name of the input file
# maxDist: the minimum distance at which two features will no longer be collapsed
# outFile: the name of the output file
# sheetNum: [optional] the position of the sheet in the Excel file to process

## EXAMPLE
# collapse_integrations_into_hotspots(inFile="Hotspots from CESC_single and 5'-3' breakpoints.xlsx", 
#                                     maxDist=3e6, 
#                                     outFile="hotspots.txt", 
#                                     sheetNum=1)

collapse_integrations_into_hotspots <- function(inFile, maxDist, outFile, sheetNum) {
  
  if (grepl(".xlsx",inFile)) {
    # needs readxl to load the Excel file
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("R package 'readxl' needed for this function to work. Please install it.\n",
           "install.packages('readxl')", call. = FALSE)
    }
    inData <- readxl::read_xlsx(inFile, sheet=sheetNum)
  } else {
    inData <- read.delim(inFile,header=F)
  }
  
  names(inData)[1:3] <- c("Chromosome","Start","End")
  
  # pre-define features for loop based upon the first row
  chrA <- inData$Chromosome[1]
  startA <- inData$Start[1]
  endA <- inData$End[1]
  counter <- 1
  
  outData <- c("Chromosome","Start","End","N_integrations")
  
  # for each row: 
  # if chromosome number matches the one in chrA:
    # if the start of the row is within maxDist of endA:
       # replace endA with end value of current row
       # iterate counter
  # if either of the above statements are False, 
    # save the data from the pre-defined features 
    # and refill them with information from the current row
  for (i in 2:(nrow(inData))) {
    if (inData$Chromosome[i] == chrA) {
      if ((inData$Start[i] - endA) < maxDist) {
        endA <- inData$End[i]
        counter <- counter + 1
      } else {
        outData <- rbind(outData,c(chrA,startA,endA,counter))
        chrA <- inData$Chromosome[i]
        startA <- inData$Start[i]
        endA <- inData$End[i]
        counter <- 1
      }
    } else {
      outData <- rbind(outData,c(chrA,startA,endA,counter))
      chrA <- inData$Chromosome[i]
      startA <- inData$Start[i]
      endA <- inData$End[i]
      counter <- 1
    }
  }
  
  # grab final hotspot
  outData <- rbind(outData,c(chrA,startA,endA,counter))
  
  # export the results to a txt file
  outData2 <- data.frame(outData[2:nrow(outData),])
  names(outData2) <- outData[1,]
  write.table(outData2,outFile,quote=F,sep="\t",row.names=F)
}
