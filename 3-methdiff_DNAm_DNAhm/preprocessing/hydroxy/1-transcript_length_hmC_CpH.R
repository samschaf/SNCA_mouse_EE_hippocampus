#!/usr/bin/env Rscript

###Creating longest transcript annotation for protein-coding genes in a mouse RRBS dataset
#Author: Samantha Schaffner
#Date: Jan 5, 2023

#Most cytosine sites are annotated to multiple features. In this script, I reference the #AnnotationHub database to retrieve the length of transcripts associated with each site, and #determine which transcript is the longest.

#Note that each comparison had the same list of CpH sites to start with, so this annotation can be #used for all comparisons.

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
Annotations <- builtin_annotations()
(MM10_Annos <- Annotations[grep("mm10", Annotations)]) #19 total
MM10_Annos
MM10_Annos <- MM10_Annos[-which(MM10_Annos %in% c("mm10_genes_cds", "mm10_enhancers_fantom","mm10_CpHs","mm10_basicgenes"))]
Mouse_Annos <- build_annotations(genome='mm10', annotations = MM10_Annos)

#get the site list from the hmC object used as input in each regression
load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/hydroxy/CpH/all_EE_hmC_betas_CpH_10X_auto_noout.RData")
chr <- sapply(1:nrow(betas_variable), function(x) unlist(strsplit(rownames(betas_variable)[x], split="\\."))[[1]])
pos <- sapply(1:nrow(betas_variable), function(x) as.integer(unlist(strsplit(rownames(betas_variable)[x], split="\\."))[[2]]))

hmC_anno <- data.frame(annotate_regions(regions=GRanges(ranges=IRanges(start=pos,end=pos), seqnames=chr), annotations=Mouse_Annos))
hmC_anno$site <- paste(hmC_anno$seqnames, hmC_anno$start, sep=": ")
save(hmC_anno, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/hydroxy/CpH/hmC_anno_CpH.RData")

load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/hydroxy/CpH/hmC_anno_CpH.RData")
unique_CpHs <- unique(hmC_anno$site)


## Determine the longest transcript for each site
#In this step we will (1) subset the data to sites that have multiple annotations, and (2) for #those sites, determine which is the longest open reading frame.

#subset to sites that are annotated to transcripts (tx_id column does not contain NAs)
hmC_anno_multi_genes <- hmC_anno[complete.cases(hmC_anno$annot.tx_id),]

#subset further to only the columns we are interested in
hmC_anno_multi_genes_sub <- hmC_anno_multi_genes[,c(1:3,12:16)]
#remove any duplicate rows
hmC_anno_multi_genes_sub <- distinct(hmC_anno_multi_genes_sub)

#access the mm10 gene database from AnnotationHub
ah <- AnnotationHub() #snapshotDate(): 2022-04-25
query(ah, c("Mus musculus", "TxDb"))

#"knownGene" has transcript IDs in matching format
knownGenedb <- ah[["AH52263"]] #TxDb.Mmusculus.UCSC.mm10.ensGene.sqlite
knownGene_txlengths <- transcriptLengths(knownGenedb)

#remove "." extension in annotated hmC object
hmC_anno_multi_genes_sub$annot.tx_id <- sapply(1:nrow(hmC_anno_multi_genes_sub), function(x) { unlist(strsplit(hmC_anno_multi_genes_sub$annot.tx_id[x], split="\\."))[1]})

nrow(knownGene_txlengths_sub <- knownGene_txlengths[knownGene_txlengths$tx_name %in% hmC_anno_multi_genes_sub$annot.tx_id,]) #274 entries (transcripts only contained in your dataset)

#column "annot.id" in results_anno_multi_genes is unique for every entry. Each gene has a few #different transcripts associated with it, and several features per transcript that might overlap #- e.g. "exon" and "first exon."
#first, determine which is the longest transcript, and subset to those.
#then, within each longest transcript annotation, determine what features would over-ride each #other...

longest_transcript <- sapply(1:length(unique_CpHs), function(i) {
  
  #print(i)
  CpH_set <- hmC_anno_multi_genes_sub[hmC_anno_multi_genes_sub$site==unique_CpHs[i],]
  CpH_set_genes <- knownGene_txlengths_sub[knownGene_txlengths_sub$tx_name %in% CpH_set$annot.tx_id,]
  #CpH_set_length <- getBM(mart=mouse_mart, filters="external_gene_name", values=CpH_set$annot.symbol, attributes=c("external_gene_name", "transcript_length"))
  
  if (nrow(CpH_set_genes)>1){ max_length_gene = CpH_set_genes[CpH_set_genes$tx_len==max(CpH_set_genes$tx_len),"tx_name"]
  if (length(max_length_gene)==1) { return(max_length_gene) }
  else { return(max_length_gene[1]) }}
  
  else if (nrow(CpH_set_genes)==1){ return(CpH_set_genes$tx_name) }
  else { return("none") }
  #else { stop() }
})

names(longest_transcript) <- unique_CpHs

save(longest_transcript, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/preprocessing/hydroxy/CpH/longest_transcript_hmC_CpH.RData")

