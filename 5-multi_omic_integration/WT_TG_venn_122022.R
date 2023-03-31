##Samantha Schaffner
# Dec 9, 2022
.libPaths(c(.libPaths(), "/mnt/scratch/KoborLab/R_Libs/4.2", "/mnt/scratch/KoborLab/Personal_Folders/sschaffner/Rstudio/R/x86_64-pc-linux-gnu-library/4.0/"))
library(VennDiagram)

####DNAm
lm_WTSE_TGSE <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/methylation/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_mC <- lm_WTSE_TGSE[lm_WTSE_TGSE$threshold==TRUE,] #5 CpGs

lm_WTEE_TGEE <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/methylation/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_mC <- lm_WTEE_TGEE[lm_WTEE_TGEE$threshold==TRUE,] #16 CpGs

length(WTEE_TGEE_mC$site[WTEE_TGEE_mC$site %in% WTSE_TGSE_mC$site]) #0

x=list(WTSE_TGSE_mC$site, WTEE_TGEE_mC$site)

venn.diagram(x=list(WTEE_TGEE_mC$site, WTSE_TGSE_mC$site), 
             category.names=c("TGEE vs WTEE", "TGSE vs WTSE"), 
             col="transparent",
             fill=c("palegreen3","grey48"),
             alpha=0.5,
             print.mode=c("raw","percent"),
             imagetype="png",
             output=FALSE,
             ext.text=FALSE, 
             margin=0.4, cat.just=list(c(-0.4,18), c(1,11)), 
             cat.fontfamily="arial", fontfamily="arial",
             width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="DNA methylation", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/figures/Fig2/venn_DNAm.png")

#permuting lack of overlap - expected or depleted?
source("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/overlap_permutation_2way.R")
overlap_permutation_2way(WTSE_TGSE_mC$site, WTEE_TGEE_mC$site, background.probes=lm_WTEE_TGEE$site, 1000, Group1="SE", Group2="EE")
#[1] "Permutation P values for enrichment and depletion"
#[1] "Enrichment SE-EE: 1; Depletion SE-EE: 1"

###DNAhm
lm_WTSE_TGSE <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/hydroxy/WTSE_TGSE/linear_model/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_hmC <- lm_WTSE_TGSE[lm_WTSE_TGSE$threshold==TRUE,] #1583 DHM CpGs

lm_WTEE_TGEE <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/hydroxy/WTEE_TGEE/linear_model/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_hmC <- lm_WTEE_TGEE[lm_WTEE_TGEE$threshold==TRUE,] #1414 CpGs

length(WTEE_TGEE_hmC$site[WTEE_TGEE_hmC$site %in% WTSE_TGSE_hmC$site]) #823

x=list(WTSE_TGSE_hmC$site, WTEE_TGEE_hmC$site)

venn.diagram(x=list(WTEE_TGEE_hmC$site, WTSE_TGSE_hmC$site), 
             category.names=c("TGEE vs WTEE", "TGSE vs WTSE"), 
             col="transparent",
             fill=c("palegreen3","grey48"),
             alpha=0.5,
             print.mode=c("raw","percent"),
             imagetype="png",
             output=FALSE,
             ext.text=FALSE, 
             margin=0.4, cat.just=list(c(1,5), c(0,5)), 
             cat.fontfamily="arial", fontfamily="arial",
             width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="DNA hydroxymethylation", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/venn_DNAm.png")

#permuting overlaps
source("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/overlap_permutation_2way.R")
overlap_permutation_2way(WTSE_TGSE_hmC$site, WTEE_TGEE_hmC$site, background.probes=lm_WTEE_TGEE$site, 1000, Group1="SE", Group2="EE")
#[1] "Permutation P values for enrichment and depletion"
#[1] "Enrichment SE-EE: 0; Depletion SE-EE: 1"

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)

####H3K4me1

WTSE_TGSE_H3K4me1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K4me1 <- WTSE_TGSE_H3K4me1[-which(WTSE_TGSE_H3K4me1$FDR==min(WTSE_TGSE_H3K4me1$FDR)),]
WTSE_TGSE_H3K4me1$coord <- paste(WTSE_TGSE_H3K4me1$chr, paste(WTSE_TGSE_H3K4me1$start, WTSE_TGSE_H3K4me1$end, sep="-"),sep=": ")

WTEE_TGEE_H3K4me1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTEE_TGEE_H3K4me1 <- WTEE_TGEE_H3K4me1[-which(WTEE_TGEE_H3K4me1$FDR==min(WTEE_TGEE_H3K4me1$FDR)),]
WTEE_TGEE_H3K4me1$coord <- paste(WTEE_TGEE_H3K4me1$chr, paste(WTEE_TGEE_H3K4me1$start, WTEE_TGEE_H3K4me1$end, sep="-"),sep=": ")

venn.diagram(x=list(WTSE_TGSE_H3K4me1$coord, WTEE_TGEE_H3K4me1$coord), 
             category.names=c("TGSE vs WTSE", "TGEE vs WTEE"), 
             col="transparent",
             fill=c("grey48","palegreen3"),
             alpha=0.5,
             print.mode=c("raw","percent"),
             imagetype="png",
             output=FALSE,
             ext.text=FALSE, 
             margin=c(0.4,0.4), cat.just=list(c(1,5), c(0,5)), 
             cat.fontfamily="arial", fontfamily="arial",
             width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K4me1", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/venn_H3K4me1_test.png")


####H3K27ac

WTSE_TGSE_H3K27ac <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K27ac <- WTSE_TGSE_H3K27ac[-which(WTSE_TGSE_H3K27ac$FDR==min(WTSE_TGSE_H3K27ac$FDR)),]
WTSE_TGSE_H3K27ac$coord <- paste(WTSE_TGSE_H3K27ac$chr, paste(WTSE_TGSE_H3K27ac$start, WTSE_TGSE_H3K27ac$end, sep="-"),sep=": ")

WTEE_TGEE_H3K27ac <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTEE_TGEE_H3K27ac <- WTEE_TGEE_H3K27ac[-which(WTEE_TGEE_H3K27ac$FDR==min(WTEE_TGEE_H3K27ac$FDR)),]
WTEE_TGEE_H3K27ac$coord <- paste(WTEE_TGEE_H3K27ac$chr, paste(WTEE_TGEE_H3K27ac$start, WTEE_TGEE_H3K27ac$end, sep="-"),sep=": ")

venn.diagram(x=list(WTSE_TGSE_H3K27ac$coord, WTEE_TGEE_H3K27ac$coord), 
             category.names=c("TGSE vs WTSE", "TGEE vs WTEE"), 
             col="transparent",
             fill=c("grey48","palegreen3"),
             alpha=0.5,
             print.mode=c("raw","percent"),
             imagetype="png",
             output=FALSE,
             ext.text=FALSE, 
             margin=c(0.4,0.4), cat.just=list(c(1,5), c(0,5)), 
             cat.fontfamily="arial", fontfamily="arial",
             width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K27ac", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/8-multi_omic_integration/overlap/venn_H3K27ac.png")


####H3K4me3

WTSE_TGSE_H3K4me3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K4me3 <- WTSE_TGSE_H3K4me3[-which(WTSE_TGSE_H3K4me3$FDR==min(WTSE_TGSE_H3K4me3$FDR)),]
WTSE_TGSE_H3K4me3$coord <- paste(WTSE_TGSE_H3K4me3$chr, paste(WTSE_TGSE_H3K4me3$start, WTSE_TGSE_H3K4me3$end, sep="-"),sep=": ")

WTEE_TGEE_H3K4me3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTEE_TGEE_H3K4me3 <- WTEE_TGEE_H3K4me3[-which(WTEE_TGEE_H3K4me3$FDR==min(WTEE_TGEE_H3K4me3$FDR)),]
WTEE_TGEE_H3K4me3$coord <- paste(WTEE_TGEE_H3K4me3$chr, paste(WTEE_TGEE_H3K4me3$start, WTEE_TGEE_H3K4me3$end, sep="-"),sep=": ")

venn.diagram(x=list(WTSE_TGSE_H3K4me3$coord, WTEE_TGEE_H3K4me3$coord), 
             category.names=c("TGSE vs WTSE", "TGEE vs WTEE"), 
             col="transparent",
             fill=c("grey48","palegreen3"),
             alpha=0.5,
             print.mode=c("raw","percent"),
             imagetype="png",
             output=FALSE,
             ext.text=TRUE, ext.length=c(0.9,0.5),
             margin=c(0.4,0.4), cat.just=list(c(1,5), c(0,5)), 
             cat.fontfamily="arial", fontfamily="arial",
             width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K4me3", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/8-multi_omic_integration/overlap/venn_H3K4me3.png")

####RNA-seq data
WTSE_TGSE_RNA <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/6-processed_RNAseq_data/mmus_hippo12m_TG_STD_vs_WT_STD.csv")
WTSE_TGSE_RNA <- WTSE_TGSE_RNA[abs(WTSE_TGSE_RNA$log2FoldChange)>=0.3 & WTSE_TGSE_RNA$padj<=0.15,]

WTEE_TGEE_RNA <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/6-processed_RNAseq_data/mmus_hippo12m_TG_EE_vs_WT_EE.csv")
WTEE_TGEE_RNA <- WTEE_TGEE_RNA[abs(WTEE_TGEE_RNA$log2FoldChange)>=0.3 & WTEE_TGEE_RNA$padj<=0.15,]

venn.diagram(x=list(WTSE_TGSE_RNA$ensembl, WTEE_TGEE_RNA$ensembl), 
                           category.names=c("TGSE vs WTSE", "TGEE vs WTEE"), 
                           col="transparent",
                          fill=c("grey48","palegreen3"),
                           alpha=0.5,
                           print.mode=c("raw","percent"),
                           imagetype="png",
                           output=FALSE,
                           ext.text=TRUE, ext.length=c(0.4,0.9),
                           margin=c(0.4,0.4), cat.just=list(c(0,25), c(1,10)), 
                          cat.fontfamily="arial", fontfamily="arial",
                          width=500, height=500, cex=0.2, cat.cex=0.2, ext.line.lwd=0.2,
                          resolution=600,
                          main="mRNA", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
                          filename="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/8-multi_omic_integration/overlap/venn_mRNA.png")
