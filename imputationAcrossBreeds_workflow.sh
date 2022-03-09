#!/bin/bash

# run as: bash imputationAcrossBreeds_workflow.sh -f <plink_filename> -d <low_density_array> -n <sample_size> -b <breeds> -s <species> -o <outdir> -c <config>
# <plink_filename>: Plink name without .ped/.map extension
# <low_density_array>: list of SNP names from the desired low density SNP array 
# <sample_size>: n. of individuals that will be sampled initially from the original file 
# <breeds>: breed(s) that will be assigned the low density SNP array 
# <species>: Plink species identifier (e.g. cow, sheep etc.)
# <outdir> root output directory
# <config> path to config file with parameters
# Use -gt 1 to consume two arguments per pass in
# the loop (e.g. each
# argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).

Help()
{
   # Display Help
   echo "This script runs the pipeline to measure the imputation accuracy from one breed to the other.\n"
   echo "Samples from one or more breeds are assigned to the LD SNP array and imputed based on complete\n"
   echo "records from the other breed(s)"
   echo "Paths to required software packages (e.g. R, Plink, Java, Beagle) are to be set in the file pathNames.txt"
   echo
   echo "Syntax: imputationAcrossBreeds_Workflow.sh [-h|f|s|d|n|b|o]"
   echo "options:"
   echo "h     print this help"
   echo "f     Plink filename (path to) [required]"
   echo "s     species [required]"
   echo "d     low density array (list of SNP names) [required]"
   echo "n     sample size (to be sampled randomly from the dataset) [required]"
   echo "b     breed(s) assigned to the low density SNP array [required]"
   echo "o     output directory [required]"
   echo "c     path to config file [if not provided, default is used]"
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
    -d|--low_density_array)
    LOWDENSITY="$2"
    shift # past argument
    ;;
    -n|--sample_size)
    SAMPLESIZE="$2"
    shift # past argument
    ;;
    -b|--breeds)
    BREEDS="$2"
    shift # past argument
    ;;
    -s|--species)
    SPECIES="$2"
    shift
    ;; 
    -o|--outdir)
    OUTDIR="$2"
    shift # past argument
    ;;
    -c|--config)
    CONFIG="$2"
    shift
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
echo LOW DENSITY FILE     = "${LOWDENSITY}"
echo SAMPLE SIZE     = "${SAMPLESIZE}"
echo BREEDS     = "${BREEDS}"
echo OUT FOLDER     = "${OUTDIR}"
configFile="${CONFIG:-pathNames.txt}"
echo CONFIG FILE = "${configFile}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

cwd=`pwd`
echo "current directory is $cwd"
source $configFile

### Parameters from the config file
echo "#######################"
echo "## EXPERIMENT         #"
echo "#######################"

echo "Experiment type: ${PREFIX}"

echo "########################################################"
echo "## MAIN PATHS TO SOFTWARE - FROM pathNames.txt         #"
echo "########################################################"

echo "Main path to software is ${MAINPATH}"
echo "Path to Rscript is ${RPATH}"
echo "Path to Plink is ${PLINKPATH}"
echo "Path to Bagle is ${BEAGLEPATH}"

echo "#################################################"
echo "## GENOTYPE FILTERING THRESHOLDS               ##"
echo "#################################################"

echo "MAF threshold (to remove unimputable monomorphic loci): $MAF"
echo "MIND threshold (to remove samples with excess missing data): $MIND"
echo "GENO threshold (to remove loci with excess missing data): $GENO"


echo "#######################################"
echo "## STEP -1"
echo "## create unique folders for each run"
echo "#######################################"
tmstmp=$(date +%N)
currDate=$(date +%d-%m-%Y)
folderName=${PREFIX}_$( basename $INPUTFILE).${SAMPLESIZE}_$(echo $BREEDS | sed 's/,/-/')_${tmstmp}.${currDate}

echo "Folder name is: $folderName"
cd $OUTDIR

if [ -d "$folderName" ]; then
>&2 echo "!! ERROR: folder $folderName already exists !!"
exit
fi

mkdir $folderName
folderpath=`realpath $folderName`
echo "The following folder has just been created: $folderpath" 

cd $folderName

echo "#######################################"
echo "## STEP 0"
echo "## sample individuals from the ped file"
echo "#######################################"

if [ ! -f "${MAINPATH}/${INPUTFILE}.fam" ]; then
    echo "${MAINPATH}/${INPUTFILE}.fam does not exist"
    exit
fi

$PLINKPATH --$SPECIES --bfile ${MAINPATH}/${INPUTFILE} --recode transpose --out transposed
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE $MAINPATH
$PLINKPATH --${SPECIES} --bfile ${MAINPATH}/${INPUTFILE} --keep keepIDs.txt --maf $MAF --bp-space 1 --recode --out subset


echo "#############################################"
echo "## STEP 1                                    "
echo "## sample individuals for the HD and LD array"
echo "#############################################"
## STEP 1
## Prepare file with HD and LD (missing values) SNP genotypes

if [ ! -f "subset.ped" ]; then
	echo "file subset.ped/map does not exist"
	exit
fi

# get breed IDs to be assigned the LD array
RASSEN=$(echo $BREEDS | sed 's/,/\\|/') #transforms comma-separated breed names into a regular expression for grep
echo "RASSEN     = $RASSEN"
cut -f1-2 -d' ' subset.ped | grep ${RASSEN} > keep.ids
# create low-density and high-density subsets
$PLINKPATH --${SPECIES} --file subset --keep keep.ids --extract ${LOWDENSITY} --recode --out subsetLD
$PLINKPATH --${SPECIES} --file subset --remove keep.ids --recode --out subsetHD
# put HD and LD subsets together into a combined file
$PLINKPATH --${SPECIES} --file subsetHD --merge subsetLD --recode --out combined
rm subsetHD* subsetLD.ped

echo "#######################################"
echo "## STEP 1.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
$PLINKPATH --${SPECIES} --file subset --freq --out subset
$PLINKPATH --${SPECIES} --file subset --recode A --out originalRaw
$PLINKPATH --${SPECIES} --file combined --recode A --out combinedRaw
rm combinedRaw.nosex combinedRaw.log combined.bed combined.bim originalRaw.nosex originalRaw.log

echo "#######################################"
echo "## STEP 2"
echo "## imputation of missing genotypes"
echo "#######################################"
## STEP 2
## Imputation of missing genotypes
date +%s > anfangZeit
$PLINKPATH --$SPECIES --file combined --recode vcf --out combined
java -Xss5m -Xmx4g -jar $BEAGLEPATH gt=combined.vcf out=imputed
$PLINKPATH --$SPECIES --vcf imputed.vcf.gz --recode --out imputed
$PLINKPATH --$SPECIES --file imputed --recode A --out imputed
rm imputed.nosex imputed.log

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
$PLINKPATH --${SPECIES} --file imputed --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/parseResults_across.R originalRaw.raw combinedRaw.raw imputed.raw $( basename $INPUTFILE) ${BREEDS} $LOWDENSITY $MAINPATH
rm originalRaw.raw imputed.raw combinedRaw.raw subset.ped subset.log subset.map freq.frq combined.* subset.nosex


