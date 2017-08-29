## Measure computation performance of functions to inject missing data
library("ggplot2")
library("microbenchmark")


# load functions to inject missing
source("functions.R")

## Reading the .ped file in from the shared server folder 
##Ensure that R knows that the SNP calls are characters, otherwise will read T as TRUE etc
## run as: Rscript --vanilla measurePerformance.R <pedfile> <proportionMissing>
#ped <- read.table("/storage/share/jody/data/cowSubset.ped", colClasses = c("character"), header = FALSE)

args = commandArgs(trailingOnly = TRUE)

print(paste("arg1: ",args[1],sep=" "))
print(paste("arg2: ", args[2],sep=" "))

pedFile = args[1]
proportionMissing = as.numeric(args[2])

## A check for the progress of the script
print("Reading in the ped file ...")
ped <- read.table(pedFile, header = FALSE, colClasses = c("character"), na.strings = c("0"))
print("Ped file read in")

# ped <- read.table("/storage/share/jody/data/cowChr12.ped", header=FALSE, stringsAsFactors = FALSE, colClasses = c("character"))

ped <- ped[,-c(1:6)]
n <- nrow(ped)
m <- ncol(ped)/2

print(paste("N. of samples:",n,"Number of markers:",m,sep=" "))

print("Running microbenchmarking ...")
mbm <- microbenchmark("missIndex" = { idx <- getIndexVector(n = n, m = m, p = proportionMissing); M <- setToMissing_colInds(ped,idx) },
                      "missDouble" = {
                        idx <- getIndexVector(n = n, m = m, p = proportionMissing);
                        M <- setToMissing_doubleM(ped,idx);
                      })

print("Saving results from microbenchmarking ...")
save(mbm,file = "performance.RData")

pdf("missingInjection_performance.pdf")
autoplot(mbm)
dev.off()

print("DONE!")