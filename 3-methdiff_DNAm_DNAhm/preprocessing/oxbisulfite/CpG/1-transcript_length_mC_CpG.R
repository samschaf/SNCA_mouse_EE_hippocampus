#!/usr/bin/env Rscript

###Creating longest transcript annotation for protein-coding genes in a mouse RRBS dataset
#Author: Samantha Schaffner
#Date: Dec 9, 2022

#Most cytosine sites are annotated to multiple features. In this script, I reference the #AnnotationHub database to retrieve the length of transcripts associated with each site, and #determine which transcript is the longest.

#Note that each comparison had the same list of CpG sites to start with, so this annotation can be #used for all comparisons.

## Load in your methylation object and annotate all possible genomic contexts

.libPaths(c(.libPaths(), "/mnt/scratch/KoborLab/R_Libs/4.2", "/mnt/scratch/KoborLab/Personal_Folders/sschaffner/Rstudio/R/x86_64-pc-linux-gnu-library/4.0/"))
library(lifecycle, lib.loc="/mnt/scratch/KoborLab/R_Libs/4.2")
library(scales, lib.loc="/mnt/scratch/KoborLab/R_Libs/4.2")
library(vctrs, lib.loc="/mnt/scratch/KoborLab/R_Libs/4.2")

library(annotatr)
library(AnnotationHub)
library(BiSeq)
library(GenomicRanges)
library(dplyr)

#this worked in Rstudio but not in the script submitted to GPCC for some reason
#will save and load in Rscript
#Annotations <- builtin_annotations()
#(MM10_Annos <- Annotations[grep("mm10", Annotations)]) #19 total
#MM10_Annos
#MM10_Annos <- MM10_Annos[-which(MM10_Annos %in% c("mm10_genes_cds", "mm10_enhancers_fantom","mm10_cpgs","mm10_basicgenes"))]
#Mouse_Annos <- build_annotations(genome='mm10', annotations = MM10_Annos)

#get the site list from the BSraw object used as input in each regression
#load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/oxBS_BSraw_CpG_10X_noout_auto_var.RData")

#BSraw_anno <- data.frame(annotate_regions(regions=GRanges(ranges=ranges(BSraw_auto_var), seqnames=seqnames(BSraw_auto_var)), annotations=Mouse_Annos))
#BSraw_anno$site <- paste(BSraw_anno$seqnames, BSraw_anno$start, sep=": ")
#save(BSraw_anno, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/mC_BSraw_anno_CpG.RData")

load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/mC_BSraw_anno_CpG.RData")
unique_cpgs <- unique(BSraw_anno$site) #329,507


## Determine the longest transcript for each site
#In this step we will (1) subset the data to sites that have multiple annotations, and (2) for #those sites, determine which is the longest open reading frame.

#subset to sites that are annotated to transcripts (tx_id column does not contain NAs)
BSraw_anno_multi_genes <- BSraw_anno[complete.cases(BSraw_anno$annot.tx_id),]

#subset further to only the columns we are interested in
BSraw_anno_multi_genes_sub <- BSraw_anno_multi_genes[,c(1:3,12:16)]
#remove any duplicate rows
BSraw_anno_multi_genes_sub <- distinct(BSraw_anno_multi_genes_sub)

#access the mm10 gene database from AnnotationHub
ah <- AnnotationHub() #snapshotDate(): 2022-04-25
query(ah, c("Mus musculus", "TxDb"))

#"knownGene" has transcript IDs in matching format
knownGenedb <- ah[["AH52263"]] #TxDb.Mmusculus.UCSC.mm10.ensGene.sqlite
knownGene_txlengths <- transcriptLengths(knownGenedb)

#remove "." extension in annotated BSraw object
BSraw_anno_multi_genes_sub$annot.tx_id <- sapply(1:nrow(BSraw_anno_multi_genes_sub), function(x) { unlist(strsplit(BSraw_anno_multi_genes_sub$annot.tx_id[x], split="\\."))[1]})

nrow(knownGene_txlengths_sub <- knownGene_txlengths[knownGene_txlengths$tx_name %in% BSraw_anno_multi_genes_sub$annot.tx_id,]) #60,294 entries (transcripts only contained in your dataset)

#column "annot.id" in results_anno_multi_genes is unique for every entry. Each gene has a few #different transcripts associated with it, and several features per transcript that might overlap #- e.g. "exon" and "first exon."
#first, determine which is the longest transcript, and subset to those.
#then, within each longest transcript annotation, determine what features would over-ride each #other...

longest_transcript <- sapply(1:10000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript2 <- sapply(10001:20000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript2)

longest_transcript3 <- sapply(20001:30000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript3)

longest_transcript4 <- sapply(30001:40000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript4)

longest_transcript5 <- sapply(40001:50000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript5)

longest_transcript6 <- sapply(50001:60000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript6)

longest_transcript7 <- sapply(60001:70000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript7)

longest_transcript8 <- sapply(70001:80000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript8)

longest_transcript9 <- sapply(80001:90000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript9)

longest_transcript10 <- sapply(90001:100000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript10)

longest_transcript11 <- sapply(100001:110000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript11)

longest_transcript12 <- sapply(110001:120000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript12)

longest_transcript13 <- sapply(120001:130000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript13)

longest_transcript14 <- sapply(130001:140000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript14)

longest_transcript15 <- sapply(140001:150000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript15)

longest_transcript16 <- sapply(150001:160000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript16)

longest_transcript17 <- sapply(160001:170000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript17)

longest_transcript18 <- sapply(170001:180000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript18)

longest_transcript19 <- sapply(180001:190000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript19)

longest_transcript20 <- sapply(190001:200000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript20)

longest_transcript21 <- sapply(200001:210000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript21)

longest_transcript22 <- sapply(210001:220000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript22)

longest_transcript23 <- sapply(220001:230000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript23)

longest_transcript24 <- sapply(230001:240000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript24)

longest_transcript25 <- sapply(240001:250000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript25)

longest_transcript26 <- sapply(250001:260000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript26)

longest_transcript27 <- sapply(260001:270000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript27)

longest_transcript28 <- sapply(270001:280000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript28)

longest_transcript29 <- sapply(280001:290000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript29)

longest_transcript30 <- sapply(290001:300000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript30)

longest_transcript31 <- sapply(300001:310000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript31)

longest_transcript32 <- sapply(310001:320000, function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript32)

longest_transcript33 <- sapply(320001:length(unique_cpgs), function(i) {
  
  #print(i)
  cpg_set <- BSraw_anno_multi_genes_sub[BSraw_anno_multi_genes_sub$site==unique_cpgs[i],]
  cpg_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% cpg_set$annot.tx_id,]
  #cpg_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=cpg_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(cpg_set_genes)>1){ max_length_gene = cpg_set_genes[cpg_set_genes$tx_len==max(cpg_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(cpg_set_genes)==1){ return(cpg_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

longest_transcript <- c(longest_transcript, longest_transcript33)

names(longest_transcript) <- unique_cpgs

save(longest_transcript, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/oxbisulfite/CpG/longest_transcript_mC_CpG.RData")

