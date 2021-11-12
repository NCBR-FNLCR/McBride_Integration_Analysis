# McBride_Integration_Analysis
### Recurrent Integration of Human Papillomavirus Genomes at Transcriptional Regulatory Hubs

#### Overlap_analysis_genomic_features_collapsing_and_flanks.R 
Purpose: This script is designed to bootstrap integration breakpoints to determine if they significantly overlap enhancers or FANCD2 binding sites.  
  
Inputs:  
* Bed files of FANCD2 or enhancer sites. 
* A two column file of chromosome lengths
* An Excel file containing Supplemental Tables S2 and S3.  

Output: a table containing the information in Tables 1 and 2.  

#### RegioneR_overlap_analysis.R
This script is designed to take data from two bed files and a file of chromosome lengths.
It then uses regioneR to calculate the number of overlaps between the the two bed files
and randomly redistributes the regions to calculate a p-value.  
  
PARAMETERS  
* bedFile1: the name of the first bed file, this is the file whose data will be manipulated for bootstrapping
* bedFile2: the name of the bed file to be compared to
* chrFile: a two column file of chromosome and length that matches the bed files without a header 
* addFlanks: [optional] number of bases to add to each side of the regions in bedFile1 (default: 0)

EXAMPLE: `RegioneR_overlap_analysis("enhancer.bed", "HPV_integrants_cervical_ALL_2.bed", "hg19.len", addFlanks=1e5/2)`

#### Integration_breakpoint_collapse.R
This script is designed to take data from an Excel file whose first three columns are chromosome, start, and end (or a bed file). 
It will then collapse rows together based upon the distance between the End value of one row and the Start value of the previous row. It will produce a bed-like
txt file as output with the collapsed region and the number of rows that were combined.  
  
PARAMETERS  
* inFile: the name of the input file
* maxDist: the minimum distance at which two features will no longer be collapsed
* outFile: the name of the output file
* sheetNum: [optional] the position of the sheet in the Excel file to process  

EXAMPLE: `collapse_integrations_into_hotspots(inFile="Hotspots from CESC_single and 5'-3' breakpoints.xlsx", maxDist=3e6, outFile="hotspots.txt", 
sheetNum=1)`

