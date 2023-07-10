##making a nested list of all genes

####DNAm/DNAhm data
lm_WTEE_TGEE_mC_CpG <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_mC_CpG <- lm_WTEE_TGEE_mC_CpG[lm_WTEE_TGEE_mC_CpG$threshold==TRUE & lm_WTEE_TGEE_mC_CpG$adjP<=0.05,]

lm_WTEE_TGEE_hmC_CpG <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_hmC_CpG <- lm_WTEE_TGEE_hmC_CpG[lm_WTEE_TGEE_hmC_CpG$threshold==TRUE & lm_WTEE_TGEE_hmC_CpG$adjP<=0.05,]

lm_WTEE_TGEE_mC_CpH <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTEE_TGEE/linear_model/auto/CpH/lm_WTEE_TGEE_CpH.csv")
WTEE_TGEE_mC_CpH <- lm_WTEE_TGEE_mC_CpH[lm_WTEE_TGEE_mC_CpH$threshold==TRUE & lm_WTEE_TGEE_mC_CpH$adjP<=0.05,]

lm_WTEE_TGEE_hmC_CpH <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTEE_TGEE/linear_model/auto/CpH/lm_WTEE_TGEE_CpH_auto.csv")
WTEE_TGEE_hmC_CpH <- lm_WTEE_TGEE_hmC_CpH[lm_WTEE_TGEE_hmC_CpH$threshold==TRUE & lm_WTEE_TGEE_hmC_CpH$adjP<=0.05,]

#retrieving DM/DHM genes
length(CpG_mC <- unique(WTEE_TGEE_mC_CpG$gene[complete.cases(WTEE_TGEE_mC_CpG$gene)])) #0 genes
length(CpG_hmC <- unique(WTEE_TGEE_hmC_CpG$gene[complete.cases(WTEE_TGEE_hmC_CpG$gene)])) #87 genes
length(CpH_mC <- unique(WTEE_TGEE_mC_CpH$gene[complete.cases(WTEE_TGEE_mC_CpH$gene)])) #1 gene
length(CpH_hmC <- unique(WTEE_TGEE_hmC_CpH$gene[complete.cases(WTEE_TGEE_hmC_CpH$gene)])) #0 genes

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)
WTEE_TGEE_H3K4me1 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTEE_TGEE_H3K27ac <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTEE_TGEE_H3K4me3 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")

#subsetting to genes
H3K4me1 <- unique(as.character(WTEE_TGEE_H3K4me1$SYMBOL[complete.cases(WTEE_TGEE_H3K4me1$SYMBOL)])) #3,633 genes
H3K4me3 <- unique(as.character(WTEE_TGEE_H3K4me3$SYMBOL[complete.cases(WTEE_TGEE_H3K4me3$SYMBOL)])) #197 genes
H3K27ac <- unique(as.character(WTEE_TGEE_H3K27ac$SYMBOL[complete.cases(WTEE_TGEE_H3K27ac$SYMBOL)])) #1,062 genes

####RNA-seq data
WTEE_TGEE_RNA_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/4-processed_RNAseq_data/mmus_hippo12m_TG_EE_vs_WT_EE.csv")
WTEE_TGEE_RNA <- WTEE_TGEE_RNA_bg[abs(WTEE_TGEE_RNA_bg$log2FoldChange)>=0.3 & WTEE_TGEE_RNA_bg$padj<=0.05,]
mRNA <- unique(as.character(WTEE_TGEE_RNA$mgi_symbol[complete.cases(WTEE_TGEE_RNA$mgi_symbol)])) #3 genes

##putting it all together
omics <- list(CpH_mC=CpH_mC, CpG_hmC=CpG_hmC, H3K4me3=H3K4me3, H3K4me1=H3K4me1, H3K27ac=H3K27ac, mRNA=mRNA)

##########Upset plot
library(UpSetR)
#upset(fromList(omics), nsets=6, order.by="freq")

omics_meta <- data.frame(sets=names(omics))
omics_meta$type <- c("DNA_mod", "DNA_mod", rep("Histone_PTM",3), "mRNA")

upset <- upset(fromList(omics), nsets=6, order.by="freq", set.metadata=list(data=omics_meta, plots=list(list(type="matrix_rows", column="type", colors=c(DNA_mod="blue", Histone_PTM="red", mRNA="yellow")))))

##########permuting overlaps
source('~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/overlap_permutation_6way.R')
genes.RNA=unique(as.character(WTEE_TGEE_RNA_bg$mgi_symbol[complete.cases(WTEE_TGEE_RNA_bg$mgi_symbol)])) 
genes.DNAm=unique(lm_WTEE_TGEE_mC_CpH$gene[complete.cases(lm_WTEE_TGEE_mC_CpH$gene)])
genes.DNAhm=unique(lm_WTEE_TGEE_hmC_CpG$gene[complete.cases(lm_WTEE_TGEE_hmC_CpG$gene)])
H3K4me1_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K4me1/H3K4me1_WTEE-TGEE_anno.csv")
genes.K4me1 <- unique(H3K4me1_bg$SYMBOL[complete.cases(H3K4me1_bg$SYMBOL)])
H3K27ac_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K27ac/H3K27ac_WTEE-TGEE_anno.csv")
genes.K27ac <- unique(H3K27ac_bg$SYMBOL[complete.cases(H3K27ac_bg$SYMBOL)])
H3K4me3_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K4me3/H3K4me3_WTEE-TGEE_anno.csv")
genes.K4me3 <- unique(H3K4me3_bg$SYMBOL[complete.cases(H3K4me3_bg$SYMBOL)])

perm <- overlap_permutation_6way(mRNA, genes.RNA, 
  CpH_mC, genes.DNAm,
  CpG_hmC, genes.DNAhm,
  H3K4me1, genes.K4me1,
  H3K4me3, genes.K4me3,
  H3K27ac, genes.K27ac, 
  "mRNA","CpH_mC","CpG_hmC",
  "H3K4me1","H3K4me3","H3K27ac",
  permutation_number=1000)

perm[perm$real_overlap>0,]
#                      Category real_overlap enrich_padj depletion_padj
#38 CpG_hmC - H3K4me1 - H3K27ac            4       0.000          0.985
#40 H3K4me1 - H3K4me3 - H3K27ac            9       0.000          1.000
#50           CpG_hmC - H3K4me1           19       0.072          0.884
#52           CpG_hmC - H3K27ac            5       0.140          0.742
#53           H3K4me1 - H3K4me3           34       0.312          0.619
#54           H3K4me1 - H3K27ac          473       0.000          1.000
#55           H3K4me3 - H3K27ac            5       0.769          0.167

