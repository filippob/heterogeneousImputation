#!/bin/bash

# run as: bash imputationWorkflow.sh -f <plink_filename> -p <proportion_missing> -n <sample_size> -o <outdir>
# <plink_filename>: Plink name without .ped/.map extension
# <proportion_missing>: proportion of markers 
# <outdir> root output directory [must be an existing directory]
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
    -p|--proportion_missing)
    MISSING="$2"
    shift # past argument
    ;;
    -n|--sample_size)
    SAMPLESIZE="$2"
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
echo PROPORTION OF MISSING     = "${MISSING}"
echo SAMPLE SIZE     = "${SAMPLESIZE}"
echo OUT FOLDER     = "${OUTDIR}"

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

### Hard-coded main paths
echo "########################################################"
echo "## MAIN PATHS TO SOFTWARE - FROM pathNames.txt         #"
echo "########################################################"
source pathNames.txt

echo "Main path to software is ${MAINPATH}"
echo "Path to Rscript is ${RPATH}"
echo "Path to Plink is ${PLINKPATH}"


echo "#######################################"
echo "## STEP -1"
echo "## create unique folders for each run"
echo "#######################################"
PREFIX="GAPIMP"
tmstmp=$(date +%N)
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
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE $MAINPATH
$PLINKPATH --cow --file ${INPUTFILE} --keep keepIDs.txt --maf 0.05 --bp-space 1 --recode --out subset

echo "#######################################"
echo "## STEP 0.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
$PLINKPATH --cow --file subset --recode A --out originalRaw
rm originalRaw.nosex originalRaw.log

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
cp ${MAINPATH}/Zanardi/PARAMFILE.txt .
#(/usr/bin/time --format "%e" python /storage/share/jody/software/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4) > imputation_step.log 2> time_results
python ${MAINPATH}/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4 > imputation_step.log
$PLINKPATH --cow --file OUTPUT/BEAGLE_OUT_beagle4_IMPUTED --recode A --out imputedRaw
rm imputedRaw.nosex imputedRaw.log

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
$PLINKPATH --cow --file OUTPUT/BEAGLE_OUT_beagle4_IMPUTED --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
$RPATH --vanilla ${MAINPATH}/heterogeneousImputation/scripts/parseResults.R originalRaw.raw imputedRaw.raw indexes.txt $( basename $INPUTFILE) $MAINPATH
rm artificialMissing.ped artificialMissing.map originalRaw.raw imputedRaw.raw subset.ped subset.log subset.map freq.frq
rm -r OUTPUT/


