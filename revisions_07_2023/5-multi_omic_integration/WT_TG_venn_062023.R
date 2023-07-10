##Samantha Schaffner
# June 13, 2023: padj < 0.05 for DNA and RNA
#.libPaths(c(.libPaths(), "/mnt/scratch/KoborLab/R_Libs/4.2", "/mnt/scratch/KoborLab/Personal_Folders/sschaffner/Rstudio/R/x86_64-pc-linux-gnu-library/4.0/"))
library(VennDiagram)
source("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/overlap_permutation_2way.R")

####DNAm (CpH)
lm_WTSE_TGSE <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTSE_TGSE/linear_model/auto/CpH/lm_WTSE_TGSE_CpH.csv")
nrow(WTSE_TGSE_mC <- lm_WTSE_TGSE[lm_WTSE_TGSE$threshold==TRUE & lm_WTSE_TGSE$adjP<=0.05,]) #14 CpHs

lm_WTEE_TGEE <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTEE_TGEE/linear_model/auto/CpH/lm_WTEE_TGEE_CpH.csv")
nrow(WTEE_TGEE_mC <- lm_WTEE_TGEE[lm_WTEE_TGEE$threshold==TRUE & lm_WTEE_TGEE$adjP<=0.05,]) #28 CpHs

length(WTEE_TGEE_mC$site[WTEE_TGEE_mC$site %in% WTSE_TGSE_mC$site]) #6

WTSE_TGSE_mC[WTSE_TGSE_mC$site %in% WTEE_TGEE_mC$site,c("site","gene","genomic_anno","DB","adjP")]
#site gene          genomic_anno        DB       adjP
#4422  chr1.177421430 <NA> mm10_genes_intergenic 0.1214298 0.03266175
#6794  chr10.65938345 <NA> mm10_genes_intergenic 0.3286373 0.02071479
#20415 chr13.99788771 <NA> mm10_genes_intergenic 0.3118476 0.01734535
#22242 chr14.53929892 <NA> mm10_genes_intergenic 0.3350568 0.04176936
#34160 chr18.51073430 <NA> mm10_genes_intergenic 0.1232450 0.02973066
#58730  chr6.15971331 <NA> mm10_genes_intergenic 0.3151683 0.04176936

WTEE_TGEE_mC[WTEE_TGEE_mC$site %in% WTSE_TGSE_mC$site,c("site","gene","genomic_anno","DB","adjP")]
#site gene          genomic_anno        DB        adjP
#4422  chr1.177421430 <NA> mm10_genes_intergenic 0.1196461 0.027010572
#6794  chr10.65938345 <NA> mm10_genes_intergenic 0.2916247 0.023709136
#20415 chr13.99788771 <NA> mm10_genes_intergenic 0.2390866 0.028353712
#22242 chr14.53929892 <NA> mm10_genes_intergenic 0.3758910 0.021828734
#34160 chr18.51073430 <NA> mm10_genes_intergenic 0.1991364 0.001518654
#58730  chr6.15971331 <NA> mm10_genes_intergenic 0.4335454 0.007575420

#All 6 CpHs have increased DNAm in TG mice in both environments

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
             width=500, height=500, cex=0.3, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="CpH methylation", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_DNAm_CpH_p0.05_v2.png")

#permuting lack of overlap - expected or depleted?
overlap_permutation_2way(WTSE_TGSE_mC$site, WTEE_TGEE_mC$site, background.probes=lm_WTEE_TGEE$site, 1000, Group1="SE", Group2="EE")
#[1] "Permutation P values for enrichment and depletion"
#[1] "Enrichment EE-SE: 0; Depletion EE-SE: 1"

###DNAhm (CpG)
lm_WTSE_TGSE <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
nrow(WTSE_TGSE_hmC <- lm_WTSE_TGSE[lm_WTSE_TGSE$threshold==TRUE & lm_WTSE_TGSE$adjP<=0.05,]) #707 DHM CpGs

lm_WTEE_TGEE <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
nrow(WTEE_TGEE_hmC <- lm_WTEE_TGEE[lm_WTEE_TGEE$threshold==TRUE & lm_WTEE_TGEE$adjP<=0.05,]) #370 CpGs

length(WTEE_TGEE_hmC$site[WTEE_TGEE_hmC$site %in% WTSE_TGSE_hmC$site]) #248

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
             width=500, height=500, cex=0.3, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="CpG hydroxymethylation", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_DNAhm_CpG_p0.05_v2.png")

#permuting overlaps
overlap_permutation_2way(WTSE_TGSE_hmC$site, WTEE_TGEE_hmC$site, background.probes=lm_WTEE_TGEE$site, 1000, Group1="SE", Group2="EE")
#[1] "Permutation P values for enrichment and depletion"
#[1] "Enrichment SE-EE: 0; Depletion SE-EE: 1"

####RNA-seq data
WTSE_TGSE_RNA_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/4-processed_RNAseq_data/mmus_hippo12m_TG_STD_vs_WT_STD.csv")
nrow(WTSE_TGSE_RNA <- WTSE_TGSE_RNA_bg[abs(WTSE_TGSE_RNA_bg$log2FoldChange)>=0.3 & WTSE_TGSE_RNA_bg$padj<=0.15,]) #329
nrow(WTSE_TGSE_RNA <- WTSE_TGSE_RNA[abs(WTSE_TGSE_RNA$log2FoldChange)>=0.3 & WTSE_TGSE_RNA$padj<=0.05,]) #215

WTEE_TGEE_RNA_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/4-processed_RNAseq_data/mmus_hippo12m_TG_EE_vs_WT_EE.csv")
nrow(WTEE_TGEE_RNA <- WTEE_TGEE_RNA_bg[abs(WTEE_TGEE_RNA_bg$log2FoldChange)>=0.3 & WTEE_TGEE_RNA_bg$padj<=0.15,]) #14
nrow(WTEE_TGEE_RNA <- WTEE_TGEE_RNA[abs(WTEE_TGEE_RNA$log2FoldChange)>=0.3 & WTEE_TGEE_RNA$padj<=0.05,]) #3

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
                          main="mRNA expression", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
                          filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_mRNA_p0.05.png")

#permuting overlaps
overlap_permutation_2way(WTSE_TGSE_RNA$ensembl, WTEE_TGEE_RNA$ensembl, background.probes=WTEE_TGEE_RNA_bg$ensembl, 1000, Group1="SE", Group2="EE")
#[1] "Permutation P values for enrichment and depletion"
#[1] "Enrichment SE-EE: 0; Depletion SE-EE: 1"

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)

####H3K4me1

WTSE_TGSE_H3K4me1 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K4me1 <- WTSE_TGSE_H3K4me1[-which(WTSE_TGSE_H3K4me1$FDR==min(WTSE_TGSE_H3K4me1$FDR)),]
WTSE_TGSE_H3K4me1$coord <- paste(WTSE_TGSE_H3K4me1$chr, paste(WTSE_TGSE_H3K4me1$start, WTSE_TGSE_H3K4me1$end, sep="-"),sep=": ")

WTEE_TGEE_H3K4me1 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
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
             width=500, height=500, cex=0.3, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K4me1", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_H3K4me1_v2.png")


####H3K27ac

WTSE_TGSE_H3K27ac <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K27ac <- WTSE_TGSE_H3K27ac[-which(WTSE_TGSE_H3K27ac$FDR==min(WTSE_TGSE_H3K27ac$FDR)),]
WTSE_TGSE_H3K27ac$coord <- paste(WTSE_TGSE_H3K27ac$chr, paste(WTSE_TGSE_H3K27ac$start, WTSE_TGSE_H3K27ac$end, sep="-"),sep=": ")

WTEE_TGEE_H3K27ac <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
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
             width=500, height=500, cex=0.3, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K27ac", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_H3K27ac_v2.png")


####H3K4me3

WTSE_TGSE_H3K4me3 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
#removing Gulp1
WTSE_TGSE_H3K4me3 <- WTSE_TGSE_H3K4me3[-which(WTSE_TGSE_H3K4me3$FDR==min(WTSE_TGSE_H3K4me3$FDR)),]
WTSE_TGSE_H3K4me3$coord <- paste(WTSE_TGSE_H3K4me3$chr, paste(WTSE_TGSE_H3K4me3$start, WTSE_TGSE_H3K4me3$end, sep="-"),sep=": ")

WTEE_TGEE_H3K4me3 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
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
             width=500, height=500, cex=0.3, cat.cex=0.2, ext.line.lwd=0.2,
             resolution=600,
             main="H3K4me3", main.cex=0.3, main.fontfamily="arial", main.pos=c(0.5,0.7),
             filename="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/venn_H3K4me3_v2.png")
