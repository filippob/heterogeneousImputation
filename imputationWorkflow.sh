## STEP 0
## recode the ped file into a .raw file
plink --cow --file /storage/share/jody/data/cowChr25 --recode A --out originalRaw
rm originalRaw.nosex originalRaw.log

## STEP 1
## injecting artificial missing genotypes
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/scripts/injectMissing.R /storage/share/jody/data/cowChr25.ped 0.025
cp /storage/share/jody/data/cowChr25.map /storage/share/jody/TEMP/artificialMissing.map # copy the map file to where the injected ped is created

## STEP 2
## Imputation of missing genotypes
(/usr/bin/time --format "%e" python /storage/share/jody/software/Zanardi/Zanardi.py --param=/storage/share/jody/software/Zanardi/PARAMFILE.txt --beagle4) > imputation_step.log 2> time_results
plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --recode A --out imputedRaw
rm imputedRaw.nosex imputedRaw.log

## STEP 3
## MAF calculation
plink --cow --file OUTPUT/BEAGLE_OUT_stsm_IMPUTED --freq --out freq 
rm freq.log freq.nosex

## STEP 4
## parsing results
/storage/biscarinif/R-3.1.1/bin/Rscript --vanilla /storage/share/jody/software/scripts/parseResults.R originalRaw.raw imputedRaw.raw indexes.txt


