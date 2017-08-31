#!/bin/bash

# Use -gt 1 to consume two arguments per pass in the loop (e.g. each
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
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
echo INPUT FILE  = "${INPUTFILE}"
echo PROPORTION OF MISSING     = "${MISSING}"
if [[-n $1]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 $1
fi

echo "#######################################"
echo "## STEP 0"
echo "## sample individuals from the ped file"
echo "#######################################"
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/scripts/sampleRows.R ${INPUTFILE}.ped 100
plink --cow --file /storage/share/jody/data/experimentalData/${INPUTFILE} --keep keepIDs.txt --recode --out subset

echo "#######################################"
echo "## STEP 0.5"
echo "## recode the ped file into a .raw file"
echo "#######################################"
plink --cow --file subset --recode A --out originalRaw
rm originalRaw.nosex originalRaw.log

echo "#######################################"
echo "## STEP 1"
echo "## injecting artificial missing"
echo "#######################################"
## STEP 1
## injecting artificial missing genotypes
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/scripts/injectMissing.R subset.ped $MISSING
cp subset.map artificialMissing.map # copy the map file to where the injected ped is created

echo "#######################################"
echo "## STEP 2"
echo "## imputation of missing genotypes"
echo "#######################################"
## STEP 2
## Imputation of missing genotypes
cp /storage/share/jody/software/Zanardi/PARAMFILE.txt .
(/usr/bin/time --format "%e" python /storage/share/jody/software/Zanardi/Zanardi.py --param=PARAMFILE.txt --beagle4) > imputation_step.log 2> time_results
plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --recode A --out imputedRaw
rm imputedRaw.nosex imputedRaw.log

echo "#######################################"
echo "## STEP 3"
echo "## Caclulate MAF"
echo "#######################################"
## STEP 3
## MAF calculation
plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --freq --out freq 
rm freq.log freq.nosex

echo "#######################################"
echo "## STEP 4"
echo "## parsing results"
echo "#######################################"
## STEP 4
## parsing results
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/scripts/parseResults.R originalRaw.raw imputedRaw.raw indexes.txt


