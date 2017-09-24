#!/bin/bash

# run as: bash imputationAcrossBreeds_Workflow.sh -f <plink_filename> -d <low_density_array> -n <sample_size> -b <breeds> -s <species> -o <outdir>
# <plink_filename>: Plink name without .ped/.map extension
# <low_density_array>: list of SNP names from the desired low density SNP array 
# <sample_size>: n. of individuals that will be sampled initially from the original file 
# <breeds>: breed(s) that will be assigned the low density SNP array 
# <species>: Plink species identifier (e.g. cow, sheep etc.)
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

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi


### Hard-coded main paths
echo "########################################################"
echo "## MAIN PATHS TO SOFTWARE - FROM pathNames.txt         #"
echo "########################################################"
source pathNames.txt ##

echo "Main path to software is ${MAINPATH}"
echo "Path to Rscript is ${RPATH}"
echo "Path to Plink is ${PLINKPATH}"

echo "#######################################"
echo "## STEP -1"
echo "## create unique folders for each run"
echo "#######################################"
PREFIX="ACROSSBREEDIMP"
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
cd $folderName

echo "#######################################"
echo "## STEP 0"
echo "## sample individuals from the ped file"
echo "#######################################"
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE $MAINPATH
$PLINKPATH --${SPECIES} --file ${INPUTFILE} --keep keepIDs.txt --maf 0.01 --bp-space 1 --recode --out subset

echo "#############################################"
echo "## STEP 1                                    "
echo "## sample individuals for the HD and LD array"
echo "#############################################"
## STEP 1
## Prepare file with HD and LD (missing values) SNP genotypes
# get breed IDs to be assigned the LD array
RASSEN=$(echo $BREEDS | sed 's/,/\\|/') #transforms comma-separated breed names into a regular expression for grep
echo "RASSEN     = $RASSEN"
cut -f1-2 -d' ' subset.ped | grep ${RASSEN} > keep.ids
# create low-density and high-density subsets
$PLINKPATH --${SPECIES} --file subset --keep keep.ids --extract ${LOWDENSITY} --recode --out subsetLD
$PLINKPATH --${SPECIES} --file subset --remove keep.ids --recode --out subsetHD
# put HD and LD subsets together into a combined file
$PLINKPATH --${SPECIES} --file subsetHD --merge subsetLD --recode --out combined
rm subsetHD* subsetLD*

echo "#######################################"
echo "## STEP 1.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
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
cp ${MAINPATH}/Zanardi/PARAMFILE_DENSITY.txt PARAMFILE.txt
#(/usr/bin/time --format "%e" python /storage/share/jody/software/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4) > imputation_step.log 2> time_results
python ${MAINPATH}/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4 > imputation_step.log
$PLINKPATH --${SPECIES} --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --recode A --out imputedRaw
rm imputedRaw.nosex imputedRaw.log

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
$PLINKPATH --${SPECIES} --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/parseResults_density.R originalRaw.raw combinedRaw.raw imputedRaw.raw $( basename $INPUTFILE) ${BREEDS} $MAINPATH
rm originalRaw.raw imputedRaw.raw combinedRaw.raw subset.ped subset.log subset.map freq.frq combined.* subset.nosex
rm -r OUTPUT/


