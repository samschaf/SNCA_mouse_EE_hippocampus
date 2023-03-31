#!/bin/bash

code="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methratio_RRBS/methratio_fixed.py"
cd "/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_mapq_filtered/PX0876_BS_pool4/"
outdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/filtered_methratio/PX0876_BS_pool4/"
adaptdir="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/"
adaptfiles=($(ls /home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/))
adapters=(ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC)
files=($(ls *))
refgenome="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/mm10_hSNCA.fa"

for n in 0 1 2 3 4 5 6 7; do
eval seq=${files[$n]};
echo $seq;
eval outfile="${outdir}${adapters[n]}_mapqfiltered_methratio.txt";
python $code -d $refgenome -o $outfile $seq;
let n=n+1;
done

#Paired-end files mapped to mm10 genome
#specify .txt outfile

