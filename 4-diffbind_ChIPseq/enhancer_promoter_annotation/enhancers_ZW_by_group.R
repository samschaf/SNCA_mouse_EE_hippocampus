#####Defining active vs poised enhancers in DecipherPD EE mice
#Samantha Schaffner
#Oct 21 2021

#Reference: Creyghton et al., PNAS 2010 (https://www.pnas.org/content/107/50/21931)

#Enhancer definition (Creyghton et al.):
#- Peak of H3K4me1 enriched region
#- At least 1kb away from H3K4me3 peak or TSS (from UCSC table browser)
#- Merged enhancer peaks within 500bp into one enhancer
#- “Active” = H3K27ac peak within 3kb of an H3K4me1 peak (consolidate the two regions into one enhancer)

#Enhancer definiton (Zinah Wassouf reccomendation):
#- H3K4me1 consensus peaks (>=2 replicates)
#- 2kb away from H3K4me3 peaks
#- (2kb? 1kb?) away from TSS
#- "active" = within 4kb of H3K27ac peak
#- "poised" = H3K4me1 only

#For now, trying with enhancers defined in all groups. Can do them separately after.

K4me1 <- read.csv("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/H3K4me1_consPeaks_2.csv")
head(K4me1) #first 3 columns are chr, start, and stop
#K4me1$cons_all <- (K4me1$WTSE>0 & K4me1$TGSE>0 & K4me1$WTEE>0 & K4me1$TGEE>0)
#summary(K4me1$cons_all) #46,123 cons peaks

K4me3 <- read.csv("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/H3K4me3_consPeaks_2.csv")
head(K4me3) #first 3 columns are chr, start, and stop
#K4me3$cons_all <- (K4me3$WTSE>0 & K4me3$TGSE>0 & K4me3$WTEE>0 & K4me3$TGEE>0)
#summary(K4me3$cons_all) #20,635 cons peaks

K27ac <- read.csv("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/H3K27ac_consPeaks_2.csv")
head(K27ac) #first 3 columns are chr, start, and stop
#K27ac$cons_all <- (K27ac$WTSE>0 & K27ac$TGSE>0 & K27ac$WTEE>0 & K27ac$TGEE>0)
#summary(K27ac$cons_all) #38,579 cons peaks

#Annotating K4me1 regions within 2kb of a TSS
#mm10 bed file (Ensembl UCSC Table Browser) 
load("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/multi_omic/mm10_bed.RData")
K4me1$TSS_2kb <- sapply(1:nrow(K4me1), function(x){
  if (nrow(mm10_bed[mm10_bed$V3==K4me1$seqnames[x] & mm10_bed$V5>=(K4me1$start[x]-2000) & mm10_bed$V5<=(K4me1$end[x]+2000),])>0) { return(TRUE) }
  else { return(FALSE) }
})
#K4me1_all <- K4me1[K4me1$cons_all==TRUE,]

#Annotating K4me1 regions within 4kb of a K27ac peak
#K27ac_all <- K27ac[K27ac$cons_all==TRUE,]
K27ac$up_4kb <- K27ac$start-4000
K27ac$down_4kb <- K27ac$end+4000

library(GenomicRanges)
K27ac_4kb_gr <- GRanges(seqnames=K27ac$seqnames, ranges=IRanges(start=K27ac$up_4kb, end=K27ac$down_4kb), mcols=K27ac[,2:10])
K4me1_all_gr <- GRanges(seqnames=K4me1$seqnames, ranges=IRanges(start=K4me1$start, end=K4me1$end), mcols=K4me1[,4:10])

K4me1_K27ac <- as.data.frame(subsetByOverlaps(K4me1_all_gr, K27ac_4kb_gr))
K4me1_K27ac <- K4me1_K27ac[,c(1:5,8:12)]
colnames(K4me1_K27ac)[6:10] <- sapply(colnames(K4me1_K27ac)[6:10], function(x) unlist(strsplit(x, split="mcols."))[2])

K4me1_K27ac$coords <- paste(K4me1_K27ac$seqnames, paste(K4me1_K27ac$start, K4me1_K27ac$end, sep="-"), sep=": ")
K4me1$coords <- paste(K4me1$seqnames, paste(K4me1$start, K4me1$end, sep="-"), sep=": ")
K4me1$K27ac_4kb <- (K4me1$coords %in% K4me1_K27ac$coords)

#Subsetting to enhancers
#Poised vs active
nrow(enhancers_poised <- K4me1[K4me1$TSS_2kb==FALSE & K4me1$K27ac_4kb==FALSE,]) #22,270
nrow(enhancers_active <- K4me1_K27ac[K4me1_K27ac$TSS_2kb==FALSE,]) #72,899

write.table(enhancers_poised, file="~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/enhancers_poised_by_group.txt", sep="\t")
write.table(enhancers_active, file="~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/enhancers_active_by_group.txt", sep="\t")
