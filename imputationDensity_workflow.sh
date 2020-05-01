#!/bin/bash

# run as: bash imputationDensity_Workflow.sh -f <plink_filename> -s <species> -d <low_density_array> -n <sample_size> -m <maf> -l <ld_size> -o <outdir>
# <plink_filename>: Plink name without .ped/.map extension
# <species> species identifier to be used in Plink
# <low_density_array>: list of SNP names from the desired low density SNP array 
# <sample_size>: n. of individuals that will be sampled initially from the original file 
# <ld_size>: n. of individuals that will be assigned the low density SNP array 
# <maf> minimum MAF to filter data (imputation from ped/map files does not work with monomorphic sites (no info on the alternative allele) 
# <outdir> root output directory
# Use -gt 1 to consume two arguments per pass in
# the loop (e.g. each
# argument has a corresponding value to go with it).
# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).

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
    -d|--low_density_array)
    LOWDENSITY="$2"
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
    -l|--ld_size)
    LDSIZE="$2"
    shift # past argument
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
echo LOW DENSITY FILE     = "${LOWDENSITY}"
echo SAMPLE SIZE     = "${SAMPLESIZE}"
echo LD SAMPLE SIZE     = "${LDSIZE}"
echo MAF THRESHOLD     = "${MAF}"
echo OUT FOLDER     = "${OUTDIR}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

echo "################################################"
echo "## MAIN PATHS TO SOFTWARE - from pathNames.txt #"
echo "################################################"
source pathNames.txt

echo "Main path to software is $MAINPATH"
echo "Path to Rscript is $RPATH"
echo "Path to Plink is $PLINKPATH"
echo "Path to Bagle is ${BEAGLEPATH}"

echo "#######################################"
echo "## STEP -1"
echo "## create unique folders for each run"
echo "#######################################"
PREFIX="DENSITYIMP"
tmstmp=$(date +%N)
currDate=$(date +%d-%m-%Y)
folderName=${PREFIX}_$( basename $INPUTFILE).${SAMPLESIZE}_$( basename $LDSIZE)_${tmstmp}.${currDate}

echo "Folder name is: $folderName"
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
echo "#######################################"
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE $MAINPATH
$PLINKPATH --$SPECIES --file ${INPUTFILE} --keep keepIDs.txt --maf $MAF --bp-space 1 --recode --out subset
$PLINKPATH --$SPECIES --file subset --freq --out subset

echo "#######################################"
echo "## STEP 1"
echo "## sample individuals from the ped file"
echo "#######################################"
## STEP 1
## Prepare file with HD and LD (missing values) SNP genotypes
# randomly sample individuals to assign to the LD array
cut -f1-2 -d' ' subset.ped | shuf -n $LDSIZE > keep.ids
# create low-density and high-density subsets
$PLINKPATH --$SPECIES --file subset --keep keep.ids --extract ${LOWDENSITY} --recode --out subsetLD
$PLINKPATH --$SPECIES --file subset --remove keep.ids --recode --out subsetHD
# put HD and LD subsets together into a combined file
$PLINKPATH --$SPECIES --file subsetHD --merge subsetLD --allow-no-sex --recode --out combined
rm subsetHD* subsetLD.ped

echo "#######################################"
echo "## STEP 1.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
$PLINKPATH --$SPECIES --file subset --recode A --out originalRaw
$PLINKPATH --$SPECIES --file combined --recode A --out combinedRaw
rm combinedRaw.nosex combinedRaw.log combined.bed combined.bim originalRaw.nosex originalRaw.log

echo "#######################################"
echo "## STEP 2"
echo "## imputation of missing genotypes"
echo "#######################################"
## STEP 2
## Imputation of missing genotypes
date +%s > anfangZeit
## Imputation of missing genotypes
$PLINKPATH --$SPECIES --file combined --recode vcf --out combined
java -Xss5m -Xmx4g -jar $BEAGLEPATH gt=combined.vcf out=imputed
$PLINKPATH --$SPECIES --vcf imputed.vcf.gz --recode --out imputed
$PLINKPATH --$SPECIES --file imputed --recode A --out imputed
rm *.vcf *.gz

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
$PLINKPATH --$SPECIES --file imputed --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/parseResults_density.R originalRaw.raw combinedRaw.raw imputed.raw ${LDSIZE} $( basename $INPUTFILE) ${MAINPATH}
rm originalRaw.raw imputed.raw imputed.map imputed.ped combinedRaw.raw freq.frq combined.* subset.* subsetLD.*

