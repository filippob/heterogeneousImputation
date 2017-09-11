## Script to randomly sample rows within breed/group (1st column from the ped file)
## Reading the .ped file in from the shared server folder 
##Ensure that R knows that the SNP calls are characters, otherwise will read T as TRUE etc
## run as: Rscript --vanilla sampleRows.R <pedfile> <sampleSize>
#ped <- read.table("/storage/share/jody/data/cowSubset.ped", colClasses = c("character"), header = FALSE)

library("plyr")

args = commandArgs(trailingOnly = TRUE)

print(paste("arg1: ",args[1],sep=" "))
print(paste("arg2: ", args[2],sep=" "))
print(paste("arg3: ",args[3], sep=" "))

pedFile = args[1]
sampleSize = as.numeric(args[2])
pathMain = args[3];

# load functions to inject missing
source(paste(pathMain,"heterogeneousImputation/scripts/functions.R",sep="/"))

#get n. of columns
ncols <-as.numeric(unlist(strsplit(trim(system2("head",paste("-1", pedFile,  "| wc", sep=" "), stdout = TRUE)),split = "\\s+"))[2])

#read first two columns from ped file
print("Reading in the first two columns of the ped file ...")# get n. of samples
ped <- read.table(pedFile, header = FALSE, colClasses = c(rep("character",2),rep("NULL",(ncols-2))))

n <- nrow(ped)
print(paste(n,"rows read from",pedFile,sep=" "))

# dd <- ddply(ped,"V1",nrow)/n
# round(dd*sampleSize,0)

print("Sampling rows ...")
keepID <- ddply(ped,"V1",function(x,n=nrow(ped),s=sampleSize) {
  
  ns <- round((nrow(x)/n)*s) #n. of individuals to be sampled in each group
  return(x[sample(nrow(x),ns),])
})

print("Writing out sampled rows ...")
write.table(keepID,file="keepIDs.txt",quote = FALSE, col.names = FALSE, row.names = FALSE)

print(paste("File keepIDs.txt written out in",getwd(),sep=" "))
print("DONE!")
