#!/bin/bash -l

#SBATCH -A b2011029
#SBATCH -p devel
#SBATCH -t 1:00:00
#SBATCH -J EnhInt
#SBATCH -o EnhInt_HiCap.out
#SBATCH -e EnhInt_HiCap.err
#SBATCH --mail-user pelin.akan@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

# Arguments:  DetectEnhInts_perChr chrname

chrname=$1

#HiCap Analysis 
cd /bubo/home/h20/pelin/3Cproj/bin/HiCapAnalysis/

echo ./DetectEnhInts_perChr Experiments_3C.txt 3 1000 mE_HiCap_merged_EnhInt_$chrname\_ $chrname
./DetectEnhInts_perChr Experiments_3C.txt 3 1000 mE_HiCap_merged_EnhInt_$chrname\_ $chrname
 
