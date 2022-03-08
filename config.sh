#!/bin/sh

## parameters of the experiment
PREFIX="GAPIMP"

## genotype filtering threshold
MIND=0.50 #safety threshold to remove samples with excess missing
GENO=0.50 #safety threshold to remove loci with excess missing
MAF=0.01 #safety threshold to make sure that there are no monomorphic markers (impossible to impute)

## main paths
MAINPATH="/home/biscarinif/imputation/" #main project folder
RPATH="/home/biscarinif/.conda/envs/imputation/bin/Rscript" #path to Rscript
PLINKPATH="/home/biscarinif/software/plink/plink" #path to plink1.9
BEAGLEPATH="/home/biscarinif/software/beagle/beagle.5.1.jar" #path to Beagle

