#!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr_HG19 chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
cd /mnt/crick/pelina/bin/HiCapAnalysis-readfrombam

export LD_LIBRARY_PATH=/mnt/crick/pelina/bin/bamtools/lib/:${LD_LIBRARY_PATH}


echo ./AnalyseHiCap_perChr_HG19 Experiments_ChIAPET_25Mreads.txt 1 1000 HG19_ChiAPET_25M_thr1_$chrname\_ $chrname >HG19_ChiAPET_25M_thr1_$chrname.out
./AnalyseHiCap_perChr_HG19 Experiments_ChIAPET_25Mreads.txt 1 1000 HG19_ChiAPET_25M_thr1_$chrname\_ $chrname >HG19_ChiAPET_25M_thr1_$chrname.out &



