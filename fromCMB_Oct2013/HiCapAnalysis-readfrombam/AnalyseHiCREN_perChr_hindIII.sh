#!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
cd /mnt/crick/pelina/bin/HiCapAnalysis-readfrombam

export LD_LIBRARY_PATH=/mnt/crick/pelina/bin/bamtools/lib/:${LD_LIBRARY_PATH}

echo ./AnalyseHiCap_perChr_hindIII Experiments_HiC.REN.txt 1 1000 mE_HiCREN_thr1_hindIII_$chrname\_ $chrname >mE_HiCREN_thr1_hindIII_$chrname.out 
./AnalyseHiCap_perChr_hindIII Experiments_HiC.REN.txt 1 1000 mE_HiCREN_thr1_hindIII_$chrname\_ $chrname >mE_HiCREN_thr1_hindIII_$chrname.out &



