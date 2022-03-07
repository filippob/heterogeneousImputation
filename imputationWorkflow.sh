#!/bin/bash

# run as: bash imputationWorkflow.sh -f <plink_filename> -s <species> -p <proportion_missing> -n <sample_size> -m <maf> -o <outdir>
# <plink_filename>: Plink name without .ped/.map extension
# <proportion_missing>: proportion of markers
# <species> species identifier to be used in Plink
# <sample_size> n of individuals to be sampled randomly from the ped file
# <maf> minimum MAF to filter data (imputation from ped/map files does not work with monomorphic sites (no info on the alternative allele) 
# <outdir> root output directory [must be an existing directory]
# Use -gt 1 to consume two arguments per pass in
# the loop (e.g. each
# argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).

#######################
## PARAMETERS
#######################
GENO=0.50 ## safety threshold to remove loci with excess "natural" missing values
MIND=0.50 ## safety threshold to remove samples with excess "natural" missing values


Help()
{
   # Display Help
   echo "This script runs the pipeline to measure the within-dataset imputation accuracy (residual missing values).\n"
   echo "Paths to required software packages (e.g. R, Plink, Java, Beagle) are to be set in the file pathNames.txt"
   echo
   echo "Syntax: imputationWorkflow.sh [-h|f|s|p|n|m|o]"
   echo "options:"
   echo "h     print this help"
   echo "f     Plink filename (path to) [required]"
   echo "s     species [required]"
   echo "p     proportion of missing values to inject [required]"
   echo "n     sample size (to be sampled randomly from the dataset) [required]"
   echo "m     MAF threshold to filter data [required]"
   echo "o     output directory [required]"
   echo
}

# Get the options
while getopts ":h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
   esac
done


while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -f|--input_file)
    INPUTFILE="$2"
    shift # past argument
    ;;
    -s|--species)
    SPECIES="$2"
    shift
    ;;
    -p|--proportion_missing)
    MISSING="$2"
    shift # past argument
    ;;
    -n|--sample_size)
    SAMPLESIZE="$2"
    shift # past argument
    ;;
    -m|--maf)
    MAF="$2"
    shift
    ;;
    -o|--outdir)
    OUTDIR="$2"
    shift # past argument
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


echo "########################################################"
echo "## WRAPPER BASH SCRIPT TO RUN THE IMPUTATION EXPERIMENTS"
echo "########################################################"

echo INPUT FILE  = "${INPUTFILE}"
echo SPECIES  = "${SPECIES}"
echo PROPORTION OF MISSING     = "${MISSING}"
echo SAMPLE SIZE     = "${SAMPLESIZE}"
echo MAF THRESHOLD     = "${MAF}"
echo OUT FOLDER     = "${OUTDIR}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

### Hard-coded main paths
echo "########################################################"
echo "## MAIN PATHS TO SOFTWARE - FROM pathNames.txt         #"
echo "########################################################"
cwd=`pwd`
echo "current directory is $cwd"
source heterogeneousImputation/pathNames.txt

echo "Main path to software is ${MAINPATH}"
echo "Path to Rscript is ${RPATH}"
echo "Path to Plink is ${PLINKPATH}"
echo "Path to Bagle is ${BEAGLEPATH}"

echo "#######################################"
echo "## STEP -1 "
echo "## create unique folders for each run"
echo "#######################################"
PREFIX="GAPIMP"

#tmstmp=$(date +%N)
## fix for MAC OS
#if [[ "$OSTYPE" == "darwin"* ]]; then
#  tmstmp=$(date +%s)
#fi

tmstmp=$(date +%s)

currDate=$(date +%d-%m-%Y)
folderName=${PREFIX}_$( basename $INPUTFILE).${SAMPLESIZE}_${MISSING}_${tmstmp}.${currDate}

echo "Folder name is:$folderName"
cd $OUTDIR

if [ -d "$folderName" ]; then
>&2 echo "!! ERROR: folder $folderName already exists !!"
exit
fi

mkdir $folderName
cd $folderName

echo "#######################################"
echo "## STEP 0"
echo "## sample individuals from the ped file"
echo "## NB: if sample size is small, it is  "
echo "## safer to use a higher MAF filter    "
echo "## (otherwise all minor alleles could  "
echo "## be set to missing by chance no ALT) "
echo "#######################################"
$PLINKPATH --$SPECIES --allow-extra-chr --bfile ${MAINPATH}/${INPUTFILE} --recode transpose --out transposed
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE $MAINPATH
$PLINKPATH "--$SPECIES" --allow-extra-chr --bfile ${INPUTFILE} --keep keepIDs.txt --maf $MAF --bp-space 1 --snps-only 'just-acgt' --not-chr 0 --geno $GENO --mind $MIND --recode --out subset

echo "#######################################"
echo "## STEP 0.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
$PLINKPATH --$SPECIES --file subset --freq --out subset
$PLINKPATH --$SPECIES --file subset --recode A --out originalRaw

echo "#######################################"
echo "## STEP 1"
echo "## injecting artificial missing"
echo "#######################################"
## STEP 1
## injecting artificial missing genotypes
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/injectMissing.R subset.ped $MISSING $MAINPATH
cp subset.map artificialMissing.map # copy the map file to where the injected ped is created

echo "#######################################"
echo "## STEP 2"
echo "## imputation of missing genotypes"
echo "#######################################"
## STEP 2
## Imputation of missing genotypes
$PLINKPATH --$SPECIES --file artificialMissing --recode vcf --out artificialMissing
java -Xss5m -Xmx4g -jar $BEAGLEPATH gt=artificialMissing.vcf out=imputed
$PLINKPATH --$SPECIES --vcf imputed.vcf.gz --recode --out imputed
$PLINKPATH --$SPECIES --file imputed --recode A --out imputed

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
$PLINKPATH --$SPECIES --file imputed --freq --out freq 

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/parseResults.R originalRaw.raw imputed.raw indexes.txt $( basename $INPUTFILE) $MAINPATH
rm artificialMissing.ped artificialMissing.map originalRaw.raw imputed.raw *.vcf imputed.ped imputed.map *.gz subset.ped subset.log subset.map *.frq *.nosex originalRaw.log freq.log transposed.*

