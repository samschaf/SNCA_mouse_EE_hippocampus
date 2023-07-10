##making a nested list of all genes

####DNAm/DNAhm data
lm_WTSE_TGSE_mC_CpG <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_mC_CpG <- lm_WTSE_TGSE_mC_CpG[lm_WTSE_TGSE_mC_CpG$threshold==TRUE & lm_WTSE_TGSE_mC_CpG$adjP<=0.05,]

lm_WTSE_TGSE_hmC_CpG <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_hmC_CpG <- lm_WTSE_TGSE_hmC_CpG[lm_WTSE_TGSE_hmC_CpG$threshold==TRUE & lm_WTSE_TGSE_hmC_CpG$adjP<=0.05,]

lm_WTSE_TGSE_mC_CpH <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/methylation/WTSE_TGSE/linear_model/auto/CpH/lm_WTSE_TGSE_CpH.csv")
WTSE_TGSE_mC_CpH <- lm_WTSE_TGSE_mC_CpH[lm_WTSE_TGSE_mC_CpH$threshold==TRUE & lm_WTSE_TGSE_mC_CpH$adjP<=0.05,]

lm_WTSE_TGSE_hmC_CpH <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/3-methdiff_DNAm_DNAhm/hydroxy/WTSE_TGSE/linear_model/auto/CpH/lm_WTSE_TGSE_CpH_auto.csv")
WTSE_TGSE_hmC_CpH <- lm_WTSE_TGSE_hmC_CpH[lm_WTSE_TGSE_hmC_CpH$threshold==TRUE & lm_WTSE_TGSE_hmC_CpH$adjP<=0.05,]

#retrieving DM/DHM genes
length(CpG_mC <- unique(WTSE_TGSE_mC_CpG$gene[complete.cases(WTSE_TGSE_mC_CpG$gene)])) #0 genes
length(CpG_hmC <- unique(WTSE_TGSE_hmC_CpG$gene[complete.cases(WTSE_TGSE_hmC_CpG$gene)])) #147 genes
length(CpH_mC <- unique(WTSE_TGSE_mC_CpH$gene[complete.cases(WTSE_TGSE_mC_CpH$gene)])) #2 genes
length(CpH_hmC <- unique(WTSE_TGSE_hmC_CpH$gene[complete.cases(WTSE_TGSE_hmC_CpH$gene)])) #0 genes

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)
WTSE_TGSE_H3K4me1 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTSE_TGSE_H3K27ac <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTSE_TGSE_H3K4me3 <- read.delim("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/archive/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")

#subsetting to genes
H3K4me1 <- unique(as.character(WTSE_TGSE_H3K4me1$SYMBOL[complete.cases(WTSE_TGSE_H3K4me1$SYMBOL)])) #3,437 genes
H3K4me3 <- unique(as.character(WTSE_TGSE_H3K4me3$SYMBOL[complete.cases(WTSE_TGSE_H3K4me3$SYMBOL)])) #10 genes
H3K27ac <- unique(as.character(WTSE_TGSE_H3K27ac$SYMBOL[complete.cases(WTSE_TGSE_H3K27ac$SYMBOL)])) #1,024 genes

####RNA-seq data
WTSE_TGSE_RNA_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/4-processed_RNAseq_data/mmus_hippo12m_TG_STD_vs_WT_STD.csv")
WTSE_TGSE_RNA <- WTSE_TGSE_RNA_bg[abs(WTSE_TGSE_RNA_bg$log2FoldChange)>=0.3 & WTSE_TGSE_RNA_bg$padj<=0.05,]
mRNA <- unique(as.character(WTSE_TGSE_RNA$mgi_symbol[complete.cases(WTSE_TGSE_RNA$mgi_symbol)])) #166 genes

##putting it all together
omics <- list(CpH_mC=CpH_mC, CpG_hmC=CpG_hmC, H3K4me3=H3K4me3, H3K4me1=H3K4me1, H3K27ac=H3K27ac, mRNA=mRNA)

##########Upset plot
library(UpSetR)
#upset(fromList(omics), nsets=6, order.by="freq")

omics_meta <- data.frame(sets=names(omics))
omics_meta$type <- c("DNA_mod", "DNA_mod", rep("Histone_PTM",3), "mRNA")

upset <- upset(fromList(omics), nsets=6, order.by="freq", set.metadata=list(data=omics_meta, plots=list(list(type="matrix_rows", column="type", colors=c(DNA_mod="blue", Histone_PTM="red", mRNA="yellow")))))

##########permuting overlaps
source('~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/overlap/overlap_permutation_6way.R')
genes.RNA=unique(as.character(WTSE_TGSE_RNA_bg$mgi_symbol[complete.cases(WTSE_TGSE_RNA_bg$mgi_symbol)])) 
genes.DNAm=unique(lm_WTSE_TGSE_mC_CpH$gene[complete.cases(lm_WTSE_TGSE_mC_CpH$gene)])
genes.DNAhm=unique(lm_WTSE_TGSE_hmC_CpG$gene[complete.cases(lm_WTSE_TGSE_hmC_CpG$gene)])
H3K4me1_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K4me1/H3K4me1_WTSE-TGSE_anno.csv")
genes.K4me1 <- unique(H3K4me1_bg$SYMBOL[complete.cases(H3K4me1_bg$SYMBOL)])
H3K27ac_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K27ac/H3K27ac_WTSE-TGSE_anno.csv")
genes.K27ac <- unique(H3K27ac_bg$SYMBOL[complete.cases(H3K27ac_bg$SYMBOL)])
H3K4me3_bg <- read.csv("~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/5-ChIPseq/diffbind_2_reps_all_conspeaks/H3K4me3/H3K4me3_WTSE-TGSE_anno.csv")
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
#29    mRNA - H3K4me1 - H3K27ac            2       0.211          0.534
#38 CpG_hmC - H3K4me1 - H3K27ac            4       0.014          0.945
#40 H3K4me1 - H3K4me3 - H3K27ac            1       0.007          0.911
#43              mRNA - H3K4me1           34       0.045          0.939
#45              mRNA - H3K27ac            6       0.488          0.380
#50           CpG_hmC - H3K4me1           30       0.084          0.889
#52           CpG_hmC - H3K27ac            9       0.095          0.830
#53           H3K4me1 - H3K4me3            1       0.062          0.899
#54           H3K4me1 - H3K27ac          455       0.000          1.000
#55           H3K4me3 - H3K27ac            1       0.006          0.960

### GO enrichment for H3K4me1/mRNA overlap
library(gprofiler2)
length(background_genes <- unique(genes.K4me1[genes.K4me1 %in% genes.RNA])) #8,059
length(K4me1_mRNA_genes <- H3K4me1[H3K4me1 %in% mRNA]) #36
length(K4me1_mRNA_genes <- K4me1_mRNA_genes[-which(K4me1_mRNA_genes %in% H3K27ac)]) #34

gr_H3K4me1_mRNA <- gost(K4me1_mRNA_genes, organism="mmusculus", user_threshold = 0.05, correction_method="fdr",
                        domain_scope="custom", custom_bg=background_genes, evcodes=TRUE)
gr_H3K4me1_mRNA$result[,c("term_name","source","p_value","intersection_size","intersection")]
#term_name source    p_value intersection_size intersection
#1 ActRIIA-ActRIB-Smad3-Arip1 complex  CORUM 0.01527761                 1       Acvr2a
#2            Ng2-Grip1-Glur2 complex  CORUM 0.01527761                 1        Cspg4
#3            Heg1-Krit1-Ccm2 complex  CORUM 0.01527761                 1         Heg1
#4          Ctbp1-Glis2-Hdac3 complex  CORUM 0.01527761                 1        Glis2
#5                   mmu-miR-1843a-5p  MIRNA 0.04272093                 2 Lrch3,Rnf169

#with FDR 0.2 and intersection size > 2
gr_H3K4me1_mRNA$result[gr_H3K4me1_mRNA$result$intersection_size>2,c("term_name","source","p_value","intersection_size","intersection")]
#term_name source   p_value intersection_size
#35 developmental maturation  GO:BP 0.1959418                 4
#70        metal ion binding  GO:MF 0.1258488                14
#73           cation binding  GO:MF 0.1352242                14
#83      calcium ion binding  GO:MF 0.1791895                 4
#intersection
#35                                                                 Rps6ka2,Cspg4,Snx10,Fat4
#70 Neurl1b,Hpcal1,Otof,Glis2,Rps6ka2,Mgat5b,Tns1,Hivep1,Heg1,Rnf169,Fat4,Acvr2a,Msmo1,Zfhx2
#73 Neurl1b,Hpcal1,Otof,Glis2,Rps6ka2,Mgat5b,Tns1,Hivep1,Heg1,Rnf169,Fat4,Acvr2a,Msmo1,Zfhx2
#83                                                                    Hpcal1,Otof,Heg1,Fat4


# Table of genes with H3K4me1/mRNA change in SE
nrow(tbl.K4me1 <- WTSE_TGSE_H3K4me1[WTSE_TGSE_H3K4me1$SYMBOL %in% K4me1_mRNA_genes,]) #45 peaks
tbl.K4me1$coord <- paste(tbl.K4me1$chr, paste(tbl.K4me1$start, tbl.K4me1$end, sep="-"), sep=": ")
write.csv(tbl.K4me1, file="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/H3K4me1_peaks_mRNA_overlap.csv", row.names=F)

nrow(tbl.mRNA <- WTSE_TGSE_RNA[WTSE_TGSE_RNA$mgi_symbol %in% K4me1_mRNA_genes,]) #34 genes
write.csv(tbl.mRNA, file="~/KoborLab/Projects/DecipherPD_mouse/EE_asyn/sschaffner/6-multi_omic_integration/upset_plots_venn_diagrams/DEGs_H3K4me1_overlap.csv", row.names=F)
       