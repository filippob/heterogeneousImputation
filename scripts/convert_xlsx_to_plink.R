## R script to convert from Excel ffile format to Plink ped/map files
## kinship matrix used to account for population structure in the data
## input: Plink .raw and .map files + phenotype file
# run as Rscript --vanilla convert_xlsx_to_plink.R excel_file=path_to_xlsx sheet_index=sheet_index_number label=population_label


# library("snpStats")
# browseVignettes("snpStats")


library("GenABEL")
library("tidyverse")
library("data.table")

print("Converting from Excel file to Plink ped/map format")

###################################
## read arguments from command line
###################################
allowed_parameters = c(
  'excel_file',
  'sheet_index',
  'label'
)

args <- commandArgs(trailingOnly = TRUE)

print(args)
for (p in args){
  pieces = strsplit(p, '=')[[1]]
  #sanity check for something=somethingElse
  if (length(pieces) != 2){
    stop(paste('badly formatted parameter:', p))
  }
  if (pieces[1] %in% allowed_parameters)  {
    assign(pieces[1], pieces[2])
    next
  }
  
  #if we get here, is an unknown parameter
  stop(paste('bad parameter:', pieces[1]))
}

# excel_file = "data/pesco/SNP_array/MASPES_18K_genomeV2.xlsx"
# sheet_index = 1
# label = "pop001"

## convert sheet index to numeric/integer
sheet_index = as.integer(sheet_index)

print(paste("excel file name:", excel_file))
print(paste("sheet index number:", sheet_index))
print(paste("dataset label:", label))


#### read data #######
print("reading data")
temp <- readxl::read_xlsx(excel_file, sheet = sheet_index)
temp <- temp %>% rename(chr = chrom, nsame = rs) %>% mutate(strand = "+")
temp <- temp %>% as_tibble() %>% relocate(strand, .after = pos)
temp[temp == "--"] <- "00"

path = dirname(excel_file)
fname = paste(path, "/", "geno", ".illu", sep="")
fwrite(x = temp, file = fname, col.names = TRUE, sep = "\t")

print(paste("genotypes written to", fname))

gtname = paste(path, "/", "geno", ".raw", sep="")
convert.snp.illumina(inf=fname,
                       out=gtname,
                       strand="file")

## generating file with pseudo-phenotypes
print("generating pseudophenotypes")
df = data.frame("id"=names(temp)[-c(1:4)],"sex"=0)
df$phenotype = rnorm(nrow(df))
phname = paste(path, "/", "mock_phenotypes.txt", sep="")
fwrite(x = df, file = phname, sep = "\t")

#### load into GenABEL #######
print("loading data into GenABEL")
gendf <- load.gwaa.data(phe=phname,
                       gen=gtname,
                       force=TRUE)

### convert to Plink #####
print("exporting to Plink")
GenABEL::export.plink(data = gendf, filebasename = paste(path,label,sep="/"), 
                      phenotypes = NULL, transpose = FALSE,
                      export012na = FALSE, extendedmap = FALSE)


print("DONE!")
