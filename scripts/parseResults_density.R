## Reading the .ped file in from the shared server folder 
##Ensure that R knows that the SNP calls are characters, otherwise will read T as TRUE etc
## run as: Rscript --vanilla parseResults.R <original_ped.raw> <imputed_ped.raw> <indexes> <experiment_name>
#ped <- read.table("/storage/share/jody/data/cowSubset.ped", colClasses = c("character"), header = FALSE)

end.time <- Sys.time()
anfangZeit <- scan("anfangZeit")

## load libraries
library("plyr")
library("data.table")


## read input from command line
args = commandArgs(trailingOnly = TRUE)

print(paste("Original raw file: ",args[1],sep=" "))
print(paste("Combined (missing) raw file: ", args[2],sep=" "))
print(paste("Imputed raw file: ", args[3],sep=" "))
#print(paste("experiment name: ", args[4],sep=" "))
print(paste("N. of LD samples: ", args[5],sep=" "))
print(paste("Name of low-density SNP file: ",args[6] ,sep=" "))
print(paste("Main path is: ", args[6], sep=" "))

originalRaw_file = args[1]
combinedRaw_file = args[2]
impRaw_file = args[3]
ldSize = args[4]
ldFile = args[5]
pathMain = args[6]

experiment = paste(ldFile,ldSize,sep="_")

# load functions to inject missing
source(paste(pathMain,"heterogeneousImputation/scripts/functions.R",sep="/"))
#source("/storage/share/jody/software/scripts/functions.R")

# originalRaw_file = "/storage/share/jody/filippo/density/prova/DENSITYIMP_chr22.100_20_283948240.04-09-2017/originalRaw.raw"
# combinedRaw_file = "/storage/share/jody/filippo/density/prova/DENSITYIMP_chr22.100_20_283948240.04-09-2017/combinedRaw.raw"
# impRaw_file =  "/storage/share/jody/filippo/density/prova/DENSITYIMP_chr22.100_20_283948240.04-09-2017/imputedRaw.raw"
# idx_file = "/storage/share/jody/filippo/heterogeneousImputation/scripts/devR/imputed.raw"

###########################################################################################################

## read in the original data in .raw format
originalRaw <- fread(originalRaw_file, header = TRUE)
originalRaw <- originalRaw[,7:ncol(originalRaw), with = FALSE]

n <- nrow(originalRaw)
m <- ncol(originalRaw) 

print(paste(n,"samples and",m,"markers read from",originalRaw_file,sep=" "))

## read in the combined missing data in .raw format
combinedRaw <- fread(combinedRaw_file, header = TRUE)
combinedRaw <- combinedRaw[,7:ncol(combinedRaw), with = FALSE]

n <- nrow(combinedRaw)
m <- ncol(combinedRaw) 

## read in the imputed geotypes in .raw format
impRaw <- fread(impRaw_file,header = TRUE)
impRaw <- impRaw[,7:ncol(impRaw), with=FALSE]

originalRaw <- as.matrix(originalRaw)
impRaw <- as.matrix(impRaw)

## retrieve original and imputed genotypes (0/1/2)
## corresponding to missing values in the combinedRaw file (LD + HD datasets)
originalGenotypes <- originalRaw[is.na(combinedRaw)]
imputedGenotypes <- impRaw[is.na(combinedRaw)]

proportionMissing <- round(length(imputedGenotypes)/(n*m),3)
print(paste("Proportion of missing genotypes:", proportionMissing, sep= " "))

# remove large .raw files
rm(originalRaw,impRaw,combinedRaw)

print(paste(length(imputedGenotypes),"imputed genotypes retrieved from the original and imputed ped files", sep=" "))

## dataframe with results
res <- cbind.data.frame(originalGenotypes,imputedGenotypes)
res <- na.omit(res)

## measure total accuracy of imputation
totalAccuracy <- sum(res$originalGenotypes==res$imputedGenotypes)/nrow(res)

D <- ddply(res,"originalGenotypes",function(x) {
  
  accuracy=sum(x$originalGenotypes==x$imputedGenotypes)/nrow(x)
  return("accuracy"=accuracy)
})

dd <- ddply(res,"originalGenotypes",function(x) {
  
  n= nrow(x)
  to0= sum(x$imputedGenotypes==0)
  to1= sum(x$imputedGenotypes==1)
  to2= sum(x$imputedGenotypes==2)
  
  return(c("n"=n,"to0"=to0,"to1"=to1,"to2"=to2))
})

##preparing results
print("Preparing dataframe with results")

# elapsed_time <- scan("time_results")
elapsed_time <- (end.time-as.POSIXct(anfangZeit, origin="1970-01-01"))
# filename <- gsub("\\.raw$","",basename(originalRaw_file))

##maf
freq <- read.table("freq.frq",header = TRUE)
print("MAF read from freq.frq (Plink)")

ergebnisse <- data.frame(
  "experiment_name"=experiment,
  "sample_size"=n,
  "scaling_up"=paste(nrow(ldFile),"to",m,sep="_"),
  "proportion_missing"=proportionMissing,
  "avgMAF"=mean(freq$MAF,na.rm = TRUE),
  "nLD"=ldSize,
  "totalAccuracy"=totalAccuracy,
  "accuracyAA"=D[D$originalGenotypes==0,"V1"],
  "accuracyAB"=D[D$originalGenotypes==1,"V1"],
  "accuracyBB"=D[D$originalGenotypes==2,"V1"],
  "elapsed_time"=elapsed_time,
  "nAA"=dd[dd$originalGenotypes==0,"n"],
  "nAB"=dd[dd$originalGenotypes==1,"n"],
  "nBB"=dd[dd$originalGenotypes==2,"n"],
  "AAtoAB"=dd[dd$originalGenotypes==0,"to1"],
  "AAtoBB"=dd[dd$originalGenotypes==0,"to2"],
  "ABtoAA"=dd[dd$originalGenotypes==1,"to0"],
  "ABtoBB"=dd[dd$originalGenotypes==1,"to2"],
  "BBtoAA"=dd[dd$originalGenotypes==2,"to0"],
  "BBtoAB"=dd[dd$originalGenotypes==2,"to1"],
  "currentDate"= as.numeric(Sys.time())
)

## writing out results to 'results.csv'
print("Writing out results to results.csv")

if(!file.exists("results.csv")){
  write.table(ergebnisse,"results.csv",col.names = TRUE,row.names = FALSE, sep=",")
} else write.table(ergebnisse,"results.csv",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)

## to be done: measure accuracy for AA (0), AB (1) and BB (2) genotypes. Beware of possible NA's (not with the cow data, but with the sheep data)

