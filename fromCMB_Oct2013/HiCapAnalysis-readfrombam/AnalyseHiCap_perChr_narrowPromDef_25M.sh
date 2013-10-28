#!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
cd /mnt/crick/pelina/bin/HiCapAnalysis-readfrombam

export LD_LIBRARY_PATH=/mnt/crick/pelina/bin/bamtools/lib/:${LD_LIBRARY_PATH}

./AnalyseHiCap_perChr_narrowPromDef Experiments_HiCap.25Mreads.txt 1 1000 mE_HiCap_narrowPromDef_25M_$chrname\_ $chrname >mE_HiCap_narrowPromDef_25M_$chrname\.out &

