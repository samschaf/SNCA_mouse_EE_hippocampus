#!/bin/bash

seqdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/Data/Hippocampus/RRBS/raw_RRBS_data/hipp_oxBS_pool2_redo/PX1175/CD5MEANXX_6/"
cd $seqdir


#cd $seqdir/forward/
forwardfiles=($(ls CD5MEANXX_6_1_*))
#cd $seqdir/reverse/
revfiles=($(ls CD5MEANXX_6_2_*))

outdir="/home/BCRICWH.LAN/sschaffner/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/2-aligned_RRBS/PX1175_oxBS_pool2_redo/"

refgenome="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/mm10.fa"

adaptdir="/home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/"
adaptfiles=($(ls /home/BCRICWH.LAN/sschaffner/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/align_scripts_final/DecipherPD_adapters/))
adapters=(ACAGTG ACTTGA ATCACG CAGATC CGATGT GCCAAT TGACCA TTAGGC)

num_mismatch=5

for n in 0 1 2 3 4 5 6 7; do
eval outfile=$outdir${adapters[n]}_aligned_bsmap.bam;
eval adapter=($(cat $adaptdir${adaptfiles[n]}));
bsmap -a $seqdir/${forwardfiles[$n]} -b $seqdir/${revfiles[$n]} -d $refgenome -o $outfile -v $num_mismatch -g 0 -x 500 -m 20 -A $adapter -D C-CGG -r 0 -p 2;   
done


#Paired-end files mapped to mm10_hSNCA genome
#specify BAM outfile
#5 mismatches allowed
#no gaps 
#20-500bp
#adapter trimming
#MspI-digested RRBS mode
#report only uniquely mapped reads
