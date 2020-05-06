# Run these steps to filter out low quality variants 
Steps are numbered based on the order it needs to be run
## Decompose data 

## Normalize data 

## GATK
DP < 15 -G-filter-name "LowDP" 
GQ <= 20" -G-filter-name "LowGQ"
