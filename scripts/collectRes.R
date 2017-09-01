## Script to collect results from experiments
## Rscript --vanilla collectRes.R (no arguments)

# parameters to modify
res.folder = '/storage/share/jody/pilots'
outfile = 'collected_res.csv'

#first we move there
print(paste("Looking for results in",res.folder,sep=" "))
setwd(res.folder)

#then we list all available result folders
res.list = list.files(pattern='POPIMP_*')

data = NULL
for (res in res.list){
  #opening the result file
  res.file = file.path(res, 'results.csv')
  
  #if file does not exist, we silently skip the folder
  if (!file.exists(res.file)) next
  
  data = rbind(data, read.csv(res.file, stringsAsFactors = FALSE))
}

write.csv(data, outfile)
print(paste("Results collected and written to", outfile, sep=" "))