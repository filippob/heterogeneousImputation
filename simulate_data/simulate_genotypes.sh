#!/bin/bash

prjfolder="$HOME/Documents/chiara/imputation/heterogeneousImputation/simulated_data"

cd $prjfolder

echo " - running simulations ... "
./QMSim simdata.prm -o 

echo "DONE!"
