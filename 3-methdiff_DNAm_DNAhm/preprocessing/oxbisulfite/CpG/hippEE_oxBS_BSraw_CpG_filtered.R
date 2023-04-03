#!/usr/bin/env Rscript

#DecipherPD Mouse Enriched Environment Hippocampus Data Pre-Processing
#================================

##### Analyst: Samantha Schaffner
##### Date: Sept 15, 2022

#This project investigates the effects of genotype and lifestyle factors on the DNA methylome of C57BL/6 mice. Specifically, wild type and human SNCA #transgenic mice (with a BAC insertion containing full-length human SNCA and all its regulatory elements) were reared to the age of either 6 or 12 #months before sacrificing. A subset of 6 month old mice from both genotypes were exposed to a chronic stress paradigm, while a subset of 12 month old #mice from both genotypes were exposed to an enriched environment paradigm. From all groups, hippocampus and striatum tissues were collected.

#DNA was extracted from mouse brain tissues and prepped for reduced representation bisulfite sequenicng (RRBS) using either a bisulfite conversion (mC + hmC) or oxidatitve bisulfite conversion (mC only) protocol. Samples were pooled in groups of 8 and sequenced on an Illumina HiSeq2500 (75bp, paired-end) at the BC Genome Sciences Centre. Alignment and methylation ratio calculations were performed using the BSMAP pipeline.

#In this script:

#Data from all samples was read into R and combined into BSraw/BSrel objects containing methylation information for the common sites covered across all samples. The data was additionally filtered for read coverage (>=10X and <=99.9th percentile).

## Working directory and libraries

.libPaths(c(.libPaths(), "/mnt/scratch/KoborLab/R_Libs/4.2", "/mnt/scratch/KoborLab/Personal_Folders/sschaffner/Rstudio/R/x86_64-pc-linux-gnu-library/4.0/"))
library(BiSeq)
library(GenomicRanges)
library(gplots)
library(RColorBrewer)
library(lumi)

## Quality check: bisulfite conversion efficiency
#Before continuing with analysis, I will look at the average DNA methylation across lambda phage spike-in DNA to determine bisulfite conversion efficiency. Phage DNA should all be unmethylated, and high quality samples should show an average methylation of 2% or less.

### WT 12m std
#WT_12m_std_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTSE/ACAGTG_methratio_lambda_only.txt")
#mean(WT_12m_std_1$ratio) #0.006708479
#WT_12m_std_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTSE/ACTTGA_methratio_lambda_only.txt")
#mean(WT_12m_std_2$ratio) #0.007476594
#WT_12m_std_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTSE/ATCACG_methratio_lambda_only.txt")
#mean(WT_12m_std_3$ratio) #0.01060752
#WT_12m_std_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTSE/CAGATC_methratio_lambda_only.txt")
#mean(WT_12m_std_4$ratio) #0.009279661

###WT 12m EE
#WT_12m_EE_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTEE/CGATGT_methratio_lambda_only.txt")
#mean(WT_12m_EE_1$ratio) #0.007564729
#WT_12m_EE_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTEE/GCCAAT_methratio_lambda_only.txt")
#mean(WT_12m_EE_2$ratio) #0.008773356
#WT_12m_EE_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTEE/TGACCA_methratio_lambda_only.txt")
#mean(WT_12m_EE_3$ratio) #0.0079125
#WT_12m_EE_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/WTEE/TTAGGC_methratio_lambda_only.txt")
#ean(WT_12m_EE_4$ratio) #0.008823061

###TG 12m std
#TG_12m_std_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGSE/ACAGTG_methratio_lambda_only.txt")
#mean(TG_12m_std_1$ratio) #0.007886425
#TG_12m_std_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGSE/ACTTGA_methratio_lambda_only.txt")
#ean(TG_12m_std_2$ratio) #0.008868263
#G_12m_std_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGSE/ATCACG_methratio_lambda_only.txt")
#mean(TG_12m_std_3$ratio) #0.008858986
#TG_12m_std_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGSE/CAGATC_methratio_lambda_only.txt")
#mean(TG_12m_std_4$ratio) #0.008204137

###TG 12m EE
#TG_12m_EE_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGEE/CGATGT_methratio_lambda_only.txt")
#mean(TG_12m_EE_1$ratio) #0.006677606
#TG_12m_EE_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGEE/GCCAAT_methratio_lambda_only.txt")
#mean(TG_12m_EE_2$ratio) #0.00653367
#TG_12m_EE_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGEE/TGACCA_methratio_lambda_only.txt")
#mean(TG_12m_EE_3$ratio) #0.007410984
#TG_12m_EE_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGEE/TTAGGC_methratio_lambda_only.txt")
#mean(TG_12m_EE_4$ratio) #0.006608163


#All the oxBS samples had successful bisulfite conversion (average DNAm in lambda spike-in was no greater than 1%).

## Load in data and subset to common sites
#I'll be looking at the methylation-only data first, contained within PX0869_oxBS_pool1 and PX0870_oxBS_pool2.
#Note: each file is very large (>1GB, ~20 million reads), so only a couple can be loaded into the environment at a time without affecting RStudio's performance.

#I'll obtain all the cytosine data and subset it to common sites.

#Enriched environment, WT and TG
WT_12m_EE_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/ACAGTG_mapqfiltered_methratio.txt")
#subset to CpG sites
#WT_12m_EE_1 <- WT_12m_EE_1[WT_12m_EE_1$context=="site",]
#create a unique identifier for each site (chromosome and position)
WT_12m_EE_1$site <- paste(WT_12m_EE_1$chr, WT_12m_EE_1$pos, sep=": ")
#filter for coverage >= 10
WT_12m_EE_1 <- WT_12m_EE_1[WT_12m_EE_1$CT_count>=10 & WT_12m_EE_1$CT_count<=quantile(WT_12m_EE_1$CT_count, 0.999) & WT_12m_EE_1$context=="CG",]
summary(WT_12m_EE_1$CT_count)
#rm(WT_12m_EE_1)

WT_12m_EE_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/ACTTGA_mapqfiltered_methratio.txt")
#WT_12m_EE_2 <- WT_12m_EE_2[WT_12m_EE_2$context=="site",]
WT_12m_EE_2$site <- paste(WT_12m_EE_2$chr, WT_12m_EE_2$pos, sep=": ")
WT_12m_EE_2 <- WT_12m_EE_2[WT_12m_EE_2$CT_count>=10 & WT_12m_EE_2$CT_count<=quantile(WT_12m_EE_2$CT_count, 0.999) & WT_12m_EE_2$context=="CG",]
#check coverage >= 10
summary(WT_12m_EE_2$CT_count)
#rm(WT_12m_EE_2)

WT_12m_EE_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/CAGATC_mapqfiltered_methratio.txt")
#WT_12m_EE_3 <- WT_12m_EE_3[WT_12m_EE_3$context=="site",]
WT_12m_EE_3$site <- paste(WT_12m_EE_3$chr, WT_12m_EE_3$pos, sep=": ")
WT_12m_EE_3 <- WT_12m_EE_3[WT_12m_EE_3$CT_count>=10 & WT_12m_EE_3$CT_count<=quantile(WT_12m_EE_3$CT_count, 0.999) & WT_12m_EE_3$context=="CG",]
#check coverage >= 10
summary(WT_12m_EE_3$CT_count)
#rm(WT_12m_EE_3)

WT_12m_EE_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/GCCAAT_mapqfiltered_methratio.txt")
#WT_12m_EE_4 <- WT_12m_EE_4[WT_12m_EE_4$context=="site",]
WT_12m_EE_4$site <- paste(WT_12m_EE_4$chr, WT_12m_EE_4$pos, sep=": ")
WT_12m_EE_4 <- WT_12m_EE_4[WT_12m_EE_4$CT_count>=10 & WT_12m_EE_4$CT_count<=quantile(WT_12m_EE_4$CT_count, 0.999) & WT_12m_EE_4$context=="CG",]
#check coverage >= 10
summary(WT_12m_EE_4$CT_count)
#rm(WT_12m_EE_4)

TG_12m_EE_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/ACAGTG_mapqfiltered_methratio.txt")
#TG_12m_EE_1 <- TG_12m_EE_1[TG_12m_EE_1$context=="site",]
TG_12m_EE_1$chr <- paste("chr", TG_12m_EE_1$chr, sep="")
TG_12m_EE_1$site <- paste(TG_12m_EE_1$chr, TG_12m_EE_1$pos, sep=": ")
TG_12m_EE_1 <- TG_12m_EE_1[TG_12m_EE_1$CT_count>=10 & TG_12m_EE_1$CT_count<=quantile(TG_12m_EE_1$CT_count, 0.999) & TG_12m_EE_1$context=="CG",]
#check coverage >= 10
summary(TG_12m_EE_1$CT_count)
#rm(TG_12m_EE_1)

TG_12m_EE_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/ACTTGA_mapqfiltered_methratio.txt")
#TG_12m_EE_2 <- TG_12m_EE_2[TG_12m_EE_2$context=="site",]
TG_12m_EE_2$chr <- paste("chr", TG_12m_EE_2$chr, sep="")
TG_12m_EE_2$site <- paste(TG_12m_EE_2$chr, TG_12m_EE_2$pos, sep=": ")
TG_12m_EE_2 <- TG_12m_EE_2[TG_12m_EE_2$CT_count>=10 & TG_12m_EE_2$CT_count<=quantile(TG_12m_EE_2$CT_count, 0.999) & TG_12m_EE_2$context=="CG",]
#check coverage >= 10
summary(TG_12m_EE_2$CT_count)
#rm(TG_12m_EE_2)

TG_12m_EE_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/CAGATC_mapqfiltered_methratio.txt")
#TG_12m_EE_3 <- TG_12m_EE_3[TG_12m_EE_3$context=="site",]
TG_12m_EE_3$chr <- paste("chr", TG_12m_EE_3$chr, sep="")
TG_12m_EE_3$site <- paste(TG_12m_EE_3$chr, TG_12m_EE_3$pos, sep=": ")
TG_12m_EE_3 <- TG_12m_EE_3[TG_12m_EE_3$CT_count>=10 & TG_12m_EE_3$CT_count<=quantile(TG_12m_EE_3$CT_count, 0.999) & TG_12m_EE_3$context=="CG",]
#check coverage >= 10
summary(TG_12m_EE_3$CT_count)
#rm(TG_12m_EE_3)

TG_12m_EE_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/GCCAAT_mapqfiltered_methratio.txt")
#TG_12m_EE_4 <- TG_12m_EE_4[TG_12m_EE_4$context=="site",]
TG_12m_EE_4$chr <- paste("chr", TG_12m_EE_4$chr, sep="")
TG_12m_EE_4$site <- paste(TG_12m_EE_4$chr, TG_12m_EE_4$pos, sep=": ")
TG_12m_EE_4 <- TG_12m_EE_4[TG_12m_EE_4$CT_count>=10 & TG_12m_EE_4$CT_count<=quantile(TG_12m_EE_4$CT_count, 0.999) & TG_12m_EE_4$context=="CG",]
#check coverage >= 10
summary(TG_12m_EE_4$CT_count)
#rm(TG_12m_EE_4)

#get common sites
common_sites <- TG_12m_EE_1$site[TG_12m_EE_1$site %in% TG_12m_EE_2$site]
common_sites <- common_sites[common_sites %in% TG_12m_EE_3$site]
common_sites <- common_sites[common_sites %in% TG_12m_EE_4$site]
common_sites <- common_sites[common_sites %in% WT_12m_EE_1$site]
common_sites <- common_sites[common_sites %in% WT_12m_EE_2$site]
common_sites <- common_sites[common_sites %in% WT_12m_EE_3$site]
common_sites <- common_sites[common_sites %in% WT_12m_EE_4$site]
length(common_sites) #5,379,072

#subset to common sites
TG_12m_EE_1 <- TG_12m_EE_1[TG_12m_EE_1$site %in% common_sites,]
TG_12m_EE_2 <- TG_12m_EE_2[TG_12m_EE_2$site %in% common_sites,]
TG_12m_EE_3 <- TG_12m_EE_3[TG_12m_EE_3$site %in% common_sites,]
TG_12m_EE_4 <- TG_12m_EE_4[TG_12m_EE_4$site %in% common_sites,]
WT_12m_EE_1 <- WT_12m_EE_1[WT_12m_EE_1$site %in% common_sites,]
WT_12m_EE_2 <- WT_12m_EE_2[WT_12m_EE_2$site %in% common_sites,]
WT_12m_EE_3 <- WT_12m_EE_3[WT_12m_EE_3$site %in% common_sites,]
WT_12m_EE_4 <- WT_12m_EE_4[WT_12m_EE_4$site %in% common_sites,]

#Standard environment, WT and TG
WT_12m_std_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/ATCACG_mapqfiltered_methratio.txt")
#WT_12m_std_1 <- WT_12m_std_1[WT_12m_std_1$context=="site",]
WT_12m_std_1$site <- paste(WT_12m_std_1$chr, WT_12m_std_1$pos, sep=": ")
WT_12m_std_1 <- WT_12m_std_1[WT_12m_std_1$CT_count>=10 & WT_12m_std_1$CT_count<=quantile(WT_12m_std_1$CT_count, 0.999) & WT_12m_std_1$context=="CG",]
#check coverage >= 10
summary(WT_12m_std_1$CT_count)
#rm(WT_12m_std_1)

WT_12m_std_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/CGATGT_mapqfiltered_methratio.txt")
#WT_12m_std_2 <- WT_12m_std_2[WT_12m_std_2$context=="site",]
WT_12m_std_2$site <- paste(WT_12m_std_2$chr, WT_12m_std_2$pos, sep=": ")
WT_12m_std_2 <- WT_12m_std_2[WT_12m_std_2$CT_count>=10 & WT_12m_std_2$CT_count<=quantile(WT_12m_std_2$CT_count, 0.999) & WT_12m_std_2$context=="CG",]
#check coverage >= 10
summary(WT_12m_std_2$CT_count)
#rm(WT_12m_std_2)

WT_12m_std_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/TGACCA_mapqfiltered_methratio.txt")
#WT_12m_std_3 <- WT_12m_std_3[WT_12m_std_3$context=="site",]
WT_12m_std_3$site <- paste(WT_12m_std_3$chr, WT_12m_std_3$pos, sep=": ")
WT_12m_std_3 <- WT_12m_std_3[WT_12m_std_3$CT_count>=10 & WT_12m_std_3$CT_count<=quantile(WT_12m_std_3$CT_count, 0.999) & WT_12m_std_3$context=="CG",]
#check coverage >= 10
summary(WT_12m_std_3$CT_count)
#rm(WT_12m_std_3)

WT_12m_std_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0869_oxBS_pool1/TTAGGC_mapqfiltered_methratio.txt")
#WT_12m_std_4 <- WT_12m_std_4[WT_12m_std_4$context=="site",]
WT_12m_std_4$site <- paste(WT_12m_std_4$chr, WT_12m_std_4$pos, sep=": ")
WT_12m_std_4 <- WT_12m_std_4[WT_12m_std_4$CT_count>=10 & WT_12m_std_4$CT_count<=quantile(WT_12m_std_4$CT_count, 0.999) & WT_12m_std_4$context=="CG",]
#check coverage >= 10
summary(WT_12m_std_4$CT_count)
#rm(WT_12m_std_4)

TG_12m_std_1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/ATCACG_mapqfiltered_methratio.txt")
#TG_12m_std_1 <- TG_12m_std_1[TG_12m_std_1$context=="site",]
TG_12m_std_1$chr <- paste("chr", TG_12m_std_1$chr, sep="")
TG_12m_std_1$site <- paste(TG_12m_std_1$chr, TG_12m_std_1$pos, sep=": ")
TG_12m_std_1 <- TG_12m_std_1[TG_12m_std_1$CT_count>=10 & TG_12m_std_1$CT_count<=quantile(TG_12m_std_1$CT_count, 0.999) & TG_12m_std_1$context=="CG",]
#check coverage >= 10
summary(TG_12m_std_1$CT_count)
#rm(TG_12m_std_1)

TG_12m_std_2 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/CGATGT_mapqfiltered_methratio.txt")
#TG_12m_std_2 <- TG_12m_std_2[TG_12m_std_2$context=="site",]
TG_12m_std_2$chr <- paste("chr", TG_12m_std_2$chr, sep="")
TG_12m_std_2$site <- paste(TG_12m_std_2$chr, TG_12m_std_2$pos, sep=": ")
TG_12m_std_2 <- TG_12m_std_2[TG_12m_std_2$CT_count>=10 & TG_12m_std_2$CT_count<=quantile(TG_12m_std_2$CT_count, 0.999) & TG_12m_std_2$context=="CG",]
#check coverage >= 10
summary(TG_12m_std_2$CT_count)
#rm(TG_12m_std_2)

TG_12m_std_3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TGACCA_mapqfiltered_methratio.txt")
#TG_12m_std_3 <- TG_12m_std_3[TG_12m_std_3$context=="site",]
TG_12m_std_3$chr <- paste("chr", TG_12m_std_3$chr, sep="")
TG_12m_std_3$site <- paste(TG_12m_std_3$chr, TG_12m_std_3$pos, sep=": ")
TG_12m_std_3 <- TG_12m_std_3[TG_12m_std_3$CT_count>=10 & TG_12m_std_3$CT_count<=quantile(TG_12m_std_3$CT_count, 0.999) & TG_12m_std_3$context=="CG",]
#check coverage >= 10
summary(TG_12m_std_3$CT_count)
#rm(TG_12m_std_3)

TG_12m_std_4 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/bsmap_methratio_mapq_filtered/PX0870_oxBS_pool2/TTAGGC_mapqfiltered_methratio.txt")
#TG_12m_std_4 <- TG_12m_std_4[TG_12m_std_4$context=="site",]
TG_12m_std_4$chr <- paste("chr", TG_12m_std_4$chr, sep="")
TG_12m_std_4$site <- paste(TG_12m_std_4$chr, TG_12m_std_4$pos, sep=": ")
TG_12m_std_4 <- TG_12m_std_4[TG_12m_std_4$CT_count>=10 & TG_12m_std_4$CT_count<=quantile(TG_12m_std_4$CT_count, 0.999) & TG_12m_std_4$context=="CG",]
#check coverage >= 10
summary(TG_12m_std_4$CT_count)
#rm(TG_12m_std_4)


#get common sites between all groups
all_12m_common <- common_sites[common_sites %in% TG_12m_std_1$site]
all_12m_common <- all_12m_common[all_12m_common %in% TG_12m_std_2$site]
all_12m_common <- all_12m_common[all_12m_common %in% TG_12m_std_3$site]
all_12m_common <- all_12m_common[all_12m_common %in% TG_12m_std_4$site]
all_12m_common <- all_12m_common[all_12m_common %in% WT_12m_std_1$site]
all_12m_common <- all_12m_common[all_12m_common %in% WT_12m_std_2$site]
all_12m_common <- all_12m_common[all_12m_common %in% WT_12m_std_3$site]
all_12m_common <- all_12m_common[all_12m_common %in% WT_12m_std_4$site]
length(all_12m_common) #914,215

#subset to common sites
TG_12m_EE_1 <- TG_12m_EE_1[TG_12m_EE_1$site %in% all_12m_common,]
TG_12m_EE_2 <- TG_12m_EE_2[TG_12m_EE_2$site %in% all_12m_common,]
TG_12m_EE_3 <- TG_12m_EE_3[TG_12m_EE_3$site %in% all_12m_common,]
TG_12m_EE_4 <- TG_12m_EE_4[TG_12m_EE_4$site %in% all_12m_common,]
WT_12m_EE_1 <- WT_12m_EE_1[WT_12m_EE_1$site %in% all_12m_common,]
WT_12m_EE_2 <- WT_12m_EE_2[WT_12m_EE_2$site %in% all_12m_common,]
WT_12m_EE_3 <- WT_12m_EE_3[WT_12m_EE_3$site %in% all_12m_common,]
WT_12m_EE_4 <- WT_12m_EE_4[WT_12m_EE_4$site %in% all_12m_common,]
TG_12m_std_1 <- TG_12m_std_1[TG_12m_std_1$site %in% all_12m_common,]
TG_12m_std_2 <- TG_12m_std_2[TG_12m_std_2$site %in% all_12m_common,]
TG_12m_std_3 <- TG_12m_std_3[TG_12m_std_3$site %in% all_12m_common,]
TG_12m_std_4 <- TG_12m_std_4[TG_12m_std_4$site %in% all_12m_common,]
WT_12m_std_1 <- WT_12m_std_1[WT_12m_std_1$site %in% all_12m_common,]
WT_12m_std_2 <- WT_12m_std_2[WT_12m_std_2$site %in% all_12m_common,]
WT_12m_std_3 <- WT_12m_std_3[WT_12m_std_3$site %in% all_12m_common,]
WT_12m_std_4 <- WT_12m_std_4[WT_12m_std_4$site %in% all_12m_common,]

# create BSraw object
#The BSraw class allows storage of position, group, and read coverage information per site.

all_EE <- data.frame(WT_12m_std_1$ratio, WT_12m_std_2$ratio, WT_12m_std_3$ratio, WT_12m_std_4$ratio, WT_12m_EE_1$ratio, WT_12m_EE_2$ratio, WT_12m_EE_3$ratio, WT_12m_EE_4$ratio, TG_12m_std_1$ratio, TG_12m_std_2$ratio, TG_12m_std_3$ratio, TG_12m_std_4$ratio, TG_12m_EE_1$ratio, TG_12m_EE_2$ratio, TG_12m_EE_3$ratio, TG_12m_EE_4$ratio)
colnames(all_EE) <- c("WTSE1", "WTSE2", "WTSE3", "WTSE4", "WTEE1", "WTEE2", "WTEE3", "WTEE4", "TGSE1", "TGSE2", "TGSE3", "TGSE4", "TGEE1", "TGEE2", "TGEE3", "TGEE4")
rownames(all_EE) <- paste(WT_12m_EE_1$chr, WT_12m_EE_1$pos, sep=".")
save(all_EE, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/all_EE_mC_betas_CpG_10X.RData")

total_reads <- data.frame(WT_12m_std_1$CT_count, WT_12m_std_2$CT_count, WT_12m_std_3$CT_count, WT_12m_std_4$CT_count, WT_12m_EE_1$CT_count, WT_12m_EE_2$CT_count, WT_12m_EE_3$CT_count, WT_12m_EE_4$CT_count, TG_12m_std_1$CT_count, TG_12m_std_2$CT_count, TG_12m_std_3$CT_count, TG_12m_std_4$CT_count, TG_12m_EE_1$CT_count, TG_12m_EE_2$CT_count, TG_12m_EE_3$CT_count, TG_12m_EE_4$CT_count)
colnames(total_reads) <- c("WTSE1","WTSE2","WTSE3","WTSE4","WTEE1","WTEE2","WTEE3","WTEE4","TGSE1","TGSE2","TGSE3","TGSE4","TGEE1","TGEE2","TGEE3","TGEE4")
meth_reads <- data.frame(WT_12m_std_1$C_count, WT_12m_std_2$C_count, WT_12m_std_3$C_count, WT_12m_std_4$C_count, WT_12m_EE_1$C_count, WT_12m_EE_2$C_count, WT_12m_EE_3$C_count, WT_12m_EE_4$C_count, TG_12m_std_1$C_count, TG_12m_std_2$C_count, TG_12m_std_3$C_count, TG_12m_std_4$C_count, TG_12m_EE_1$C_count, TG_12m_EE_2$C_count, TG_12m_EE_3$C_count, TG_12m_EE_4$C_count)
colnames(meth_reads) <- c("WTSE1","WTSE2","WTSE3","WTSE4","WTEE1","WTEE2","WTEE3","WTEE4","TGSE1","TGSE2","TGSE3","TGSE4","TGEE1","TGEE2","TGEE3","TGEE4")

metadata <- list(Sequencer="HiSeq2500", Year="2018")
rowRanges <- GRanges(seqnames=WT_12m_std_1$chr, ranges=IRanges(start=WT_12m_std_1$pos, end=WT_12m_std_1$pos))
colData <- DataFrame(group=c(rep("WTSE",4), rep("WTEE",4), rep("TGSE",4), rep("TGEE",4)), row.names=c("WTSE1", "WTSE2", "WTSE3", "WTSE4", "WTEE1", "WTEE2", "WTEE3", "WTEE4", "TGSE1", "TGSE2", "TGSE3", "TGSE4", "TGEE1", "TGEE2", "TGEE3", "TGEE4"))

oxBS_BSraw_CpG <- BSraw(metadata=metadata, rowRanges=rowRanges, colData=colData, totalReads=as.matrix(total_reads), methReads=as.matrix(meth_reads))

save(oxBS_BSraw_CpG, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/hippEE_oxBS_BSraw_CpG_10X.RData")
