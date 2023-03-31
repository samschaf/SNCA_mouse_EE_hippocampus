#!/bin/bash

#sorting
cd "/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/2-aligned_RRBS/PX1175_oxBS_pool2_redo/"
outdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_sorted/PX1175_oxBS_pool2_redo/"
adaptdir="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/adapters_1-8/"
adaptfiles=($(ls /home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/adapters_1-8/))
adapters=(ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC)
files=($(ls *_aligned_bsmap.bam))

for n in 1; do
eval seq=${files[$n]};
eval outfile=$outdir${adapters[n]}_sorted.bam;
echo $seq; 
samtools sort $seq -o $outfile;
done

#indexing
#cd "/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_sorted/"
#outdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_indexed/"
#adaptdir="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/"
#adaptfiles=($(ls /home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/))
#adapters=(ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC)
#files=($(ls *))

#for n in 0; do
#eval seq=${files[$n]};
#eval outfile=$outdir${adapters[n]}_indexed.bai;
#echo $seq; 
#samtools index $seq $outfile;
#done

#filtering for MAPQ > 10
cd "/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_sorted/PX1175_oxBS_pool2_redo/"
outdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/082022/samtools_mapq_filtered/PX1175_oxBS_pool2_redo/"
adaptdir="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/adapters_1-8/"
adaptfiles=($(ls /home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/adapters_1-8/))
adapters=(ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC)
files=($(ls *))

for n in 1; do
eval seq=${files[$n]};
eval outfile=$outdir${adapters[n]}_filtered.bam;
echo $seq; 
samtools view -q 10 -b $seq -o $outfile;
done

#returns files with minimum MAPQ 10
#prints aligned samples

