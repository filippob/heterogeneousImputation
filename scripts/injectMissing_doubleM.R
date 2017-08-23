####################################
# SCRIPT TO INJECT ARTIFICIAL
# MISSING GENOTYPES INTO A PED FILE
# DUPLICATE MISSING MATRIX
####################################

## EXAMPLE! (real data to follow)

# synthetic dataset
n <- 4
m <- 5
proportionMissing <- 0.20

# we create a (n*2m) ped file (without the first 6 descriptive columns)
ped <- matrix(rnorm(n*2*m),nrow=n, byrow = TRUE)

###############################################################################################
# FUNCTIONS
###############################################################################################
#function to generate the vector of indexes for artificial missing genotypes (to be save out as text file)
getIndexVector <- function(n,m,p=proportionMissing) {
  
  vec <- rep(FALSE,n*m)
  #now we randomly sample some artificial missing "genotypes" 
  idx <- sample(length(vec),p*(n*m))
  vec[idx] <- TRUE #use the randomly sampled indexes to set some of the FALSEs to TRUE
  return(which(vec))
}

#function to set some genotypes to missing as specified by getIndexVector()
#method of doubling the missing matrix
setToMissing_doubleM <- function(ped,inds) {
  
  n <- nrow(ped)
  m <- ncol(ped)/2
  vec <- rep(FALSE,n*m)
  vec[inds] <- TRUE
  
  #matrix of missing values
  missing <- matrix(vec,nrow = n, byrow = TRUE)
  
  # (n*m) --> (n*2m) TO let it known that the m is actually 2 columns per SNP 
  missing <- missing[,rep(1:ncol(missing),rep(2,ncol(missing)))]
  
  #assign missing to NAs
  ped[missing] <- NA
  return(ped)
}


###############################################################################################

idx <- getIndexVector(n = nrow(ped), m = ncol(ped)/2, p = proportionMissing)
M <- setToMissing_doubleM(ped,idx)


