# POCP-VS-AAI
Evaluation of POCP VS AAI: This is a R script to parse all the results from EzAAI and Bio-PI POCP marix into summary tables, analyze them and create some plots. 

## Dependencies
The script has been tested in R v.4.2.2. It has the following dependendencies: 
  * dplyr
  * tidyr
  * vioplot
  * ggplot2
  * stringr
  * car
  * xlsx

## How to run it
Download the repository. Open the script "AAI VS POCP.R" with R or Rstudio and change the "path_POCP_AAI" in the second line to the folder with the contents of the repository (it should contain a folder called pocp and aai with all the output files from the comparisons) 
Once the "path_POCP_AAI" is changed run the script. It should create a folder in the path_POCP_AAI path called results, with the summary files in XLSX and the plots in PDF. 

## Additional files in the repository not necessary for the script to run
 * Accession table.txt: contains all the NCBI accession numbers of the genomes used in the analysis
 * results: this folder contains all the outputs from the script. 
