#!/bin/bash

# run as: bash imputationWorkflow.sh -f <plink_filename> -p <proportion_missing> -n <sample_size> -o <outdir>
# <plink_filename>: Plink name without .ped/.map extension
# <proportion_missing>: proportion of markers 
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

echo "#######################################"
echo "## STEP -1"
echo "## create unique folders for each run"
echo "#######################################"
tmstmp=$(date +%N)
currDate=$(date +%d-%m-%Y)
folderName=$( basename $INPUTFILE).${SAMPLESIZE}_${MISSING}_${tmstmp}.${currDate}

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
echo "#######################################"
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/heterogeneousImputation/scripts/sampleRows.R ${INPUTFILE}.ped $SAMPLESIZE
/storage/software/plink --cow --file ${INPUTFILE} --keep keepIDs.txt --recode --out subset

echo "#######################################"
echo "## STEP 0.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
/storage/software/plink --cow --file subset --recode A --out originalRaw
rm originalRaw.nosex originalRaw.log

echo "#######################################"
echo "## STEP 1"
echo "## injecting artificial missing"
echo "#######################################"
## STEP 1
## injecting artificial missing genotypes
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/heterogeneousImputation/scripts/injectMissing.R subset.ped $MISSING
cp subset.map artificialMissing.map # copy the map file to where the injected ped is created

echo "#######################################"
echo "## STEP 2"
echo "## imputation of missing genotypes"
echo "#######################################"
## STEP 2
## Imputation of missing genotypes
cp /storage/share/jody/software/Zanardi/PARAMFILE.txt .
(/usr/bin/time --format "%e" python /storage/share/jody/software/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4) > imputation_step.log 2> time_results
/storage/software/plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --recode A --out imputedRaw
rm imputedRaw.nosex imputedRaw.log

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
/storage/software/plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/heterogeneousImputation/scripts/parseResults.R originalRaw.raw imputedRaw.raw indexes.txt $( basename $INPUTFILE)


