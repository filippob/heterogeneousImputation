
## function to extract genotypes from a .raw Plink file, corresponding to indexes in a vector
retrieveGenotypes <- function(rawPed,idx) {
  
  if(all(c("FID","IID","PAT","MAT","SEX","PHENOTYPE") %in% names(rawPed))) rawPed <- rawPed[,-c(1:6)]
  #get n. of samples and markers
  n <- nrow(rawPed)
  m <- ncol(rawPed)
  
  vec <- rep(FALSE,n*m)
  vec[idx] <- TRUE
  
  M <- matrix(vec,nrow=n,byrow = TRUE)
  inds <- which(M, arr.ind = TRUE)
  
  genotypes <- rawPed[inds]
  return(genotypes)
}

###########################################################################################################

## read in the original data in .raw format
originalRaw <- read.table("/storage/share/jody/data/cowChr12.raw", header = TRUE)

## read in the indexes
idx <- read.table("indexes.txt", header = FALSE, colClasses = c("numeric"))
idx <- idx$V1

print(paste("Proportion of missing genotypes:", round(length(idx)/(n*m),3), sep= " "))

## read in the imputed geotypes in .raw format
impRaw <- read.table("imputed.raw",header = TRUE)

## retrieve original and imputed genotypes (0/1/2)
originalGenotypes <-retrieveGenotypes(originalRaw,idx)
imputedGenotypes <- retrieveGenotypes(impRaw,idx)

## measure total accuracy of imputation
totalAccuracy <- sum(originalGenotypes==imputedGenotypes)/length(idx)

## to be done: measure accuracy for AA (0), AB (1) and BB (2) genotypes. Beware of possible NA's (not with the cow data, but with the sheep data)

