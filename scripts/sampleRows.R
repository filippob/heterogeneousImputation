## Script to randomly sample rows within breed/group (1st column from the ped file)
## Reading the .ped file in from the shared server folder 
##Ensure that R knows that the SNP calls are characters, otherwise will read T as TRUE etc
## run as: Rscript --vanilla sampleRows.R <pedfile> <sampleSize>
#ped <- read.table("/storage/share/jody/data/cowSubset.ped", colClasses = c("character"), header = FALSE)

library("tidyverse")
library("data.table")

args = commandArgs(trailingOnly = TRUE)

print(paste("arg1: ",args[1],sep=" "))
print(paste("arg2: ", args[2],sep=" "))
print(paste("arg3: ",args[3], sep=" "))

pedFile = args[1]
sampleSize = as.numeric(args[2])
pathMain = args[3]

# load functions to inject missing
source(paste(pathMain,"heterogeneousImputation/scripts/functions.R",sep="/"))

#read first two columns from ped file
print("Reading in the first two columns of the ped file ...")# get n. of samples
tfam <- fread("transposed.tfam", header = FALSE, select = c(1,2))

n <- nrow(tfam)
print(paste(n,"rows read from transposed", pedFile,sep=" "))

## generate a keep file for Plink: two columns, Family ID and Sample ID 
print("Sampling rows ...")
vec <- sample(n,sampleSize)
keepID <- tfam[vec,]

print("Writing out sampled rows ...")
write.table(keepID,file="keepIDs.txt",quote = FALSE, col.names = FALSE, row.names = FALSE)

print(paste("File keepIDs.txt written out in",getwd(),sep=" "))
print("DONE!")
