#!/bin/bash

prjfolder="$HOME/Documents/chiara/imputation/heterogeneousImputation/simulate_data"
expfolder="r_simdata"
mapfile="lm_mrk_001.txt"
pedfile="Line 1_mrk_001.txt"
label="line1"

cd $prjfolder

## map file
bsname="${mapfile%.*}"
tail -n +2 ${expfolder}/${mapfile} | awk '{print $2, $1, 0, $3}' > ${expfolder}/${bsname}.map
mv ${expfolder}/${bsname}.map ${expfolder}/${label}.map

## ped file
fname="Line 1_mrk_001.txt"
echo $fname
tail -n +2 ${expfolder}/"${fname}" | awk '{printf "POP " $1; printf " 0 0 0 -9 " ; for(i=2; i<=NF; i++) printf $i" "; print""}' > "${expfolder}/${label}.ped"

echo "DONE!"
