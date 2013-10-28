#!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
cd /mnt/crick/pelina/bin/HiCapAnalysis-readfrombam

export LD_LIBRARY_PATH=/mnt/crick/pelina/bin/bamtools/lib/:${LD_LIBRARY_PATH}

./AnalyseHiCap_perChr_narrowPromDef Experiments_HiCap.merged.txt 3 1000 mE_HiCap_narrowPromDef_$chrname\_ $chrname >mE_HiCap_narrowPromDef_$chrname\.out &

