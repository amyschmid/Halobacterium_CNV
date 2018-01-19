# 
# Calculating arCOG enrichment

Calculating enrichment in arCOG categories for Dulmage et al. " Growth phase specific modulation of cell morphology and gene expression by an archaeal histone protein."

- This repository contains the code necessary to calculate statistically significant enrichment in arCOG functional categories. 
- The code was used to generate tables 2 and 3 in the manuscript.
- Raw cog data file for Halobacterium salinarum and other archaea can be downloaded from ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/ (see also Wolf et al., 2012, Biol Direct 7:46).

## **Dependencies**

- R (version 2.15.2 or higher) base package.

## **Steps to recreate tables 2 and 3 (or your own arCOG enrichments):**

- Load R data file containing cog file for H. salinarum (NRC-1\_arCOG\_2012\_final.Rdata) or download from Wolf et al., FTP site.
- Run the command 'cogtest2' to generate arCOG enrichment table
- Run the command 'cogset' to determine which genes from your list of interest are contained within each enriched arCOG category.
