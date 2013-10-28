#!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
cd /mnt/crick/pelina/bin/HiCapAnalysis-readfrombam

export LD_LIBRARY_PATH=/mnt/crick/pelina/bin/bamtools/lib/:${LD_LIBRARY_PATH}

echo ./AnalyseHiCap_perChr Experiments_HiC.txt 3 1000 mE_HiC_$chrname\_ $chrname
./AnalyseHiCap_perChr Experiments_HiC.txt 3 1000 mE_HiC_$chrname\_ $chrname >mE_HiC_$chrname\.out &

