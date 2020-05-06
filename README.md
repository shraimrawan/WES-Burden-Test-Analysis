# WES-Burden-Test-Analysis
This folder contains code to run a burdentest on data in vcf format 

## QCVariants 
This folder contains optional filtering options to remove low quality reads 

## EthnicityAnlaysis
To get a homogenous population, we wanted to filter to only caucasians. This folder has scripts to filter the vcf to only caucasian race. 

## SNPEff Analysis
Annotates and predicts the effects of genetic variants on genes and proteins 

## BurdentTest.sh
Script to run burden test on vcf data. Will result in a txt file that can be fed into the R script below 

## BurdentestVariantCount.R
Get counts of each variant in each group (in this case responder vs. nonresponder) 
