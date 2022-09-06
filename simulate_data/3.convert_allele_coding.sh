#!/bin/sh

## script that converts alleles simulated by QMsim
## from 1/2 coding to arbitrary nucleotides (e.g. A/C)
## would this work? (operationally yes, but what about possible biases?)

prefix="line1"
targetpath="simulate_data/r_simdata"

if [ ! -d $targetpath ]; then
	mkdir -p $targetpath
fi

echo "processing file line1.ped"
inpf=${prefix}.ped
cut -f7- -d' ' ${targetpath}/$inpf > ${targetpath}/temp
cut -f1-6 -d' ' ${targetpath}/$inpf > ${targetpath}/metadata

sed -i 's/1/A/g' ${targetpath}/temp
sed -i 's/2/C/g' ${targetpath}/temp
	
outf=${prefix}_recoded.ped
paste -d' ' ${targetpath}/metadata ${targetpath}/temp > ${targetpath}/${outf}

## rename map file to make pair of equally named files
cp ${targetpath}/${prefix}.map ${targetpath}/${prefix}_recoded.map
rm ${targetpath}/temp ${targetpath}/metadata

echo "DONE!"
