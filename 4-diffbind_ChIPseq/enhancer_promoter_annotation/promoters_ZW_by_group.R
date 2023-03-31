#####Defining active promoters in DecipherPD EE mice
#Samantha Schaffner
#Oct 21, 2021

#Reference: Creyghton et al., PNAS 2010 (https://www.pnas.org/content/107/50/21931)

#Promoter definiton (Zinah Wassouf reccomendation):
#- H3K27ac consensus peaks (>=2 replicates)
#- H3K4me3 peaks within 2kb OR
#- TSS within 2kb

#For now, trying with promoters defined in all groups. Can do them separately after.

K4me3 <- read.csv("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/H3K4me3_consPeaks_2.csv")
head(K4me3) #first 3 columns are chr, start, and stop
#K4me3$cons_all <- (K4me3$WTSE>0 & K4me3$TGSE>0 & K4me3$WTEE>0 & K4me3$TGEE>0)
#summary(K4me3$cons_all) #20,635 cons peaks

K27ac <- read.csv("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/H3K27ac_consPeaks_2.csv")
head(K27ac) #first 3 columns are chr, start, and stop
#K27ac$cons_all <- (K27ac$WTSE>0 & K27ac$TGSE>0 & K27ac$WTEE>0 & K27ac$TGEE>0)
#summary(K27ac$cons_all) #38,579 cons peaks

#Annotating K27ac regions within 2kb of a TSS
#mm10 bed file (Ensembl UCSC Table Browser) 
load("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/multi_omic/mm10_bed.RData")
K27ac$TSS_2kb <- sapply(1:nrow(K27ac), function(x){
  if (nrow(mm10_bed[mm10_bed$V3==K27ac$seqnames[x] & mm10_bed$V5>=(K27ac$start[x]-2000) & mm10_bed$V5<=(K27ac$end[x]+2000),])>0) { return(TRUE) }
  else { return(FALSE) }
})
#K27ac_all <- K27ac[K27ac$cons_all==TRUE,]

#Subsetting to active promoters
nrow(promoters <- K27ac[K27ac$TSS_2kb==TRUE,]) #22,165

write.table(promoters, file="~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/ChIPseq/Zinah_Wassouf_data/consensus_peaks/promoters_by_group.txt", sep="\t")
