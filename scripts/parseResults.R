## Reading the .ped file in from the shared server folder 
##Ensure that R knows that the SNP calls are characters, otherwise will read T as TRUE etc
## run as: Rscript --vanilla parseResults.R <pedfile> <proportionMissing>
#ped <- read.table("/storage/share/jody/data/cowSubset.ped", colClasses = c("character"), header = FALSE)

## load libraries
library("plyr")

# load functions to inject missing
source("functions.R")

## read input from command line
args = commandArgs(trailingOnly = TRUE)

print(paste("arg1: ",args[1],sep=" "))
print(paste("arg2: ", args[2],sep=" "))

originalRaw_file = args[1]
impRaw_file = args[2]
idx_file = args[3]

# originalRaw_file = "/storage/share/jody/data/cowChr12.raw"
# impRaw_file =  "/storage/share/jody/filippo/heterogeneousImputation/scripts/devR/imputed.raw"
# idx_file = "/storage/share/jody/filippo/heterogeneousImputation/scripts/devR/imputed.raw"

###########################################################################################################

## read in the original data in .raw format
originalRaw <- read.table(originalRaw_file, header = TRUE)
originalRaw <- originalRaw[,-c(1:6)]

n <- nrow(originalRaw)
m <- ncol(originalRaw) 

print(paste(n,"samples and",m,"markers read from",originalRaw_file,sep=" "))

## read in the indexes
idx <- read.table(idx_file, header = FALSE, colClasses = c("numeric"))
idx <- idx$V1

print(paste("Proportion of missing genotypes:", round(length(idx)/(n*m),3), sep= " "))

## read in the imputed geotypes in .raw format
impRaw <- read.table(impRaw_file,header = TRUE)
impRaw <- impRaw[,-c(1:6)]

## retrieve original and imputed genotypes (0/1/2)
originalGenotypes <-retrieveGenotypes(originalRaw,idx)
imputedGenotypes <- retrieveGenotypes(impRaw,idx)

# remove large .raw files
rm(originalRaw,impRaw)

print(paste(length(imputedGenotypes),"imputed genotypes retrieved from the original and imputed ped files", sep=" "))

## dataframe with results
res <- cbind.data.frame(originalGenotypes,imputedGenotypes)

## measure total accuracy of imputation
totalAccuracy <- sum(res$originalGenotypes==res$imputedGenotypes)/nrow(res)

D <- ddply(res,"originalGenotypes",function(x) {
  
  accuracy=sum(x$originalGenotypes==x$imputedGenotypes)/nrow(x)
  return("accuracy"=accuracy)
})

##preparing results
print("Preparing dataframe with results")

ergebnisse <- data.frame(
  "file_name"=NULL,
  "sample_size"=NULL,
  "avgMAF"=NULL,
  "injectedMissing"=NULL,
  "totalAccuracy"=NULL,
  "accuracyAA"=NULL,
  "accuracyAB"=NULL,
  "accuracyBB"=NULL
)

filename <- gsub("\\.raw$","",basename(originalRaw_file))

ergebnisse <- data.frame(
  "file_name"=filename,
  "sample_size"=n,
  "avgMAF"=NA,
  "injectedMissing"=length(idx),
  "totalAccuracy"=totalAccuracy,
  "accuracyAA"=D[D$originalGenotypes==0,"V1"],
  "accuracyAB"=D[D$originalGenotypes==1,"V1"],
  "accuracyBB"=D[D$originalGenotypes==2,"V1"]
)

## writing out results to 'results.csv'
print("Writing out results to results.csv")

if(!file.exists("results.csv")){
  write.table(ergebnisse,"results.csv",col.names = TRUE,row.names = FALSE, sep=",")
} else write.table(res,"results.csv",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)

## to be done: measure accuracy for AA (0), AB (1) and BB (2) genotypes. Beware of possible NA's (not with the cow data, but with the sheep data)

