# QC_affymetrix
R script for affymetrix quality check

This script evaluate RLE,NUSE and RNA degradation for a given list of cel Files
The RLE and NUSE median values and the slopes of the RNA degradation curves are evaluated for each samples
Outliers are identified by computing the statistics on th boxplot of their pooled values
The variable out will contain the name of all the samples that are outliers for at least two of the methods.

```{r}
fileFolder = "cel_files/" #where to read the cel files
RDATAfileFolder = "RDATA/" #folder in which save all the intermediate and final results in RData file format
RLE_fig_folder =  "figures/RLE/"
NUSE_fig_folder = "figures/NUSE/"
RNADeg_fig_folder = "figures/RNADeg/"
maxFiles = 100 #maximum number of files to be preprocessed at the same time. It depends on the machine RAM capability on which you run this code. 
isParallel = TRUE
nCores = 40

source("src/evaluate_goodness.R")

```
