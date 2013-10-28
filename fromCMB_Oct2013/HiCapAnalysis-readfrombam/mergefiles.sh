#! /bin/bash


space=" ";

for file in `ls mE_HiCREN_25M_thr1*_SignificantInteractions_Promoters.txt`; do filelist+=($file$space); done;
echo ${filelist[@]}

asize=${#filelist[@]}
echo $asize

for ((i=0;i<asize;i++))
do awk <${filelist[i]} 'NR!=1 {print}' >${filelist[i]}\.noheader;
done;

for file2 in `ls mE_HiCREN_25M_thr1*_SignificantInteractions_Promoters.txt.noheader`; do filelist2+=($file2$space); done;
echo ${filelist2[@]}

cat ${filelist2[@]} >merged.file

awk <${filelist[1]} 'NR==1' >header.file

cat  header.file merged.file >merged.file.wheader

mv merged.file.wheader mE_HiCREN_25M_thr1_SignificantInteractions_Promoters.txt
rm merged.file
rm header.file
rm *.noheader
unset file
unset filelist
unset asize
unset file2
unset filelist2


########################
for file in `ls mE_HiCREN_25M_thr1_*_Interactions_PromProm.txt`; do filelist+=($file$space); done;
echo ${filelist[@]}

asize=${#filelist[@]}
echo $asize

for ((i=0;i<asize;i++))
do awk <${filelist[i]} 'NR!=1 {print}' >${filelist[i]}\.noheader;
done;

for file2 in `ls mE_HiCREN_25M_thr1_*_Interactions_PromProm.txt.noheader`; do filelist2+=($file2$space); done;
echo ${filelist2[@]}

cat ${filelist2[@]} >merged.file

awk <${filelist[1]} 'NR==1' >header.file

cat  header.file merged.file >merged.file.wheader

mv merged.file.wheader mE_HiCREN_25M_thr1_Interactions_PromProm.txt

rm merged.file
rm header.file
rm *.noheader
unset file
unset filelist
unset asize
unset file2
unset filelist2

###################################################
for file in `ls mE_HiCREN_25M_thr1_*_SignificantInteractions_NegCtrls.txt`; do filelist+=($file$space); done;
echo ${filelist[@]}

asize=${#filelist[@]}
echo $asize

for ((i=0;i<asize;i++))
do awk <${filelist[i]} 'NR!=1 {print}' >${filelist[i]}\.noheader;
done;

for file2 in `ls mE_HiCREN_25M_thr1_*_SignificantInteractions_NegCtrls.txt.noheader`; do filelist2+=($file2$space); done;
echo ${filelist2[@]}

cat ${filelist2[@]} >merged.file

awk <${filelist[1]} 'NR==1' >header.file

cat  header.file merged.file >merged.file.wheader

mv merged.file.wheader mE_HiCREN_25M_thr1_SignificantInteractions_NegCtrls.txt
rm merged.file
rm header.file
rm *.noheader
unset file
unset filelist
unset asize
unset file2
unset filelist2

#########################################
for file in `ls mE_HiCREN_25M_thr1_*_SignificantInteractions_NegCtrlstoNegCtrls.txt`; do filelist+=($file$space); done;
echo ${filelist[@]}

asize=${#filelist[@]}
echo $asize

for ((i=0;i<asize;i++))
do awk <${filelist[i]} 'NR!=1 {print}' >${filelist[i]}\.noheader;
done;

for file2 in `ls mE_HiCREN_25M_thr1_*_SignificantInteractions_NegCtrlstoNegCtrls.txt.noheader`; do filelist2+=($file2$space); done;
echo ${filelist2[@]}

cat ${filelist2[@]} >merged.file

awk <${filelist[1]} 'NR==1' >header.file

cat  header.file merged.file >merged.file.wheader

mv merged.file.wheader mE_HiCREN_25M_thr1_Interactions_NegCtrlstoNegCtrls.txt

rm merged.file
rm header.file
rm *.noheader
unset file
unset filelist
unset asize
unset file2
unset filelist2




