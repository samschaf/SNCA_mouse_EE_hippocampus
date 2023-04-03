##making a nested list of all genes

####DNAm/DNAhm data
lm_WTEE_TGEE_mC <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/methylation/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_mC <- lm_WTEE_TGEE_mC[lm_WTEE_TGEE_mC$threshold==TRUE,]

lm_WTEE_TGEE_hmC <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/hydroxy/WTEE_TGEE/linear_model/auto/CpG/lm_WTEE_TGEE.csv")
WTEE_TGEE_hmC <- lm_WTEE_TGEE_hmC[lm_WTEE_TGEE_hmC$threshold==TRUE,]

#retrieving DM/DHM genes
DNAm <- unique(WTEE_TGEE_mC$gene[complete.cases(WTEE_TGEE_mC$gene)]) #1 gene
DNAhm <- unique(WTEE_TGEE_hmC$gene[complete.cases(WTEE_TGEE_hmC$gene)]) #300 genes

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)
WTEE_TGEE_H3K4me1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTEE_TGEE_H3K27ac <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTEE_TGEE_H3K4me3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTEE-TGEE_cons-2_FDR0_05_FC0_5_annot.tsv")

#subsetting to genes
H3K4me1 <- unique(as.character(WTEE_TGEE_H3K4me1$SYMBOL[complete.cases(WTEE_TGEE_H3K4me1$SYMBOL)])) #3,633 genes
H3K4me3 <- unique(as.character(WTEE_TGEE_H3K4me3$SYMBOL[complete.cases(WTEE_TGEE_H3K4me3$SYMBOL)])) #197 genes
H3K27ac <- unique(as.character(WTEE_TGEE_H3K27ac$SYMBOL[complete.cases(WTEE_TGEE_H3K27ac$SYMBOL)])) #1,062 genes

####RNA-seq data
WTEE_TGEE_RNA_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/6-processed_RNAseq_data/mmus_hippo12m_TG_EE_vs_WT_EE.csv")
WTEE_TGEE_RNA <- WTEE_TGEE_RNA_bg[abs(WTEE_TGEE_RNA_bg$log2FoldChange)>=0.3 & WTEE_TGEE_RNA_bg$padj<=0.15,]
mRNA <- unique(as.character(WTEE_TGEE_RNA$mgi_symbol[complete.cases(WTEE_TGEE_RNA$mgi_symbol)])) #13 genes

##putting it all together
omics <- list(DNAm=DNAm, DNAhm=DNAhm, H3K4me3=H3K4me3, H3K4me1=H3K4me1, H3K27ac=H3K27ac, mRNA=mRNA)

##########Upset plot
library(UpSetR)

omics_meta <- data.frame(sets=names(omics))
omics_meta$type <- c("DNA_mod", "DNA_mod", rep("Histone_PTM",3), "mRNA")

upset(fromList(omics), nsets=6, order.by="freq", set.metadata=list(data=omics_meta, plots=list(list(type="matrix_rows", column="type", colors=c(DNA_mod="blue", Histone_PTM="red", mRNA="yellow")))))

##########permuting overlaps
source('/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/8-multi_omic_integration/overlap/overlap_permutation_6way.R')
genes.RNA=unique(as.character(WTEE_TGEE_RNA_bg$mgi_symbol[complete.cases(WTEE_TGEE_RNA_bg$mgi_symbol)])) 
genes.DNAm=unique(lm_WTEE_TGEE_mC$gene[complete.cases(lm_WTEE_TGEE_mC$gene)])
genes.DNAhm=unique(lm_WTEE_TGEE_hmC$gene[complete.cases(lm_WTEE_TGEE_hmC$gene)])
H3K4me1_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K4me1/H3K4me1_WTEE-TGEE_anno.csv")
genes.K4me1 <- unique(H3K4me1_bg$SYMBOL[complete.cases(H3K4me1_bg$SYMBOL)])
H3K27ac_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K27ac/H3K27ac_WTEE-TGEE_anno.csv")
genes.K27ac <- unique(H3K27ac_bg$SYMBOL[complete.cases(H3K27ac_bg$SYMBOL)])
H3K4me3_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K4me3/H3K4me3_WTEE-TGEE_anno.csv")
genes.K4me3 <- unique(H3K4me3_bg$SYMBOL[complete.cases(H3K4me3_bg$SYMBOL)])

perm2 <- overlap_permutation_6way(mRNA, genes.RNA, 
                                 DNAm, genes.DNAm,
                                 DNAhm, genes.DNAhm,
                                 H3K4me1, genes.K4me1,
                                 H3K4me3, genes.K4me3,
                                 H3K27ac, genes.K27ac, 
                                 permutation_number=1000)
save(perm2, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/WTEE_TGEE_overlap_perm.RData")

#enrichments: *** K27ac/K4me1, *** K27ac/K4me1/K4me3, *** DNAhm/K27ac/K4me1, *** DNAhm/K4me1/K4me3

#lists of genes

################## 1 gene per group
#####4 'omics
#H3K4me3, H3K27ac, DNAhm, H3K4me1
set1 <- H3K4me3[H3K4me3 %in% DNAhm & H3K4me3 %in% H3K4me1 & H3K4me3 %in% H3K27ac] #"Spred2"

#####3 'omics
#H3K4me3, DNAm, H3K4me1
set3 <- H3K4me1[H3K4me1 %in% H3K4me3 & H3K4me1 %in% DNAm] #"Sptbn4"
intersect <- c("Sptbn4", "Spred2")
#H3K4me3, DNAm, H3K27ac
set7 <- H3K4me3[H3K4me3 %in% H3K27ac & H3K4me3 %in% DNAm] #"Hlcs"
#set7[-which(set7 %in% intersect)] #"Hpcal1"
intersect <- unique(c(intersect, set7))

#####2 'omics
#mRNA, H3K4me1
set4 <- H3K4me1[H3K4me1 %in% mRNA] #"Tmsb10"
#set4[-which(set4 %in% intersect)] #"Eef1a1"
intersect <- unique(c(intersect, set4))

################## 2-10 genes per group
#####4 'omics
#DNAm, H3K27ac, DNAhm, H3K4me1
set5 <- DNAm[DNAm %in% DNAhm & DNAm %in% H3K4me1 & DNAm %in% H3K27ac]
# [1] "Pkp4"    "Casz1"   "Scarb1"  "Slco3a1" "Ntm"     "Ccdc6"   "Rgs6"    "Mast4"   "Mb21d2" 
#[10] "Slc26a8"
#set5[-which(set5 %in% intersect)] #"Eef1a1"
intersect <- unique(c(intersect, set5))

#####3 'omics
#H3K4me3, DNAhm, H3K4me1
set6 <- H3K4me1[H3K4me1 %in% H3K4me3 & H3K4me1 %in% DNAhm]
save(set6, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_H3K4me1_DNAhm_H3K4me3_EE.RData")
set6[-which(set6 %in% intersect)] #"Mogat2"  "Zfp423"  "Col13a1" "Sox5"
intersect <- unique(c(intersect, set6))
#DNAm, H3K27ac, DNAhm
set7 <- H3K27ac[H3K27ac %in% DNAm & H3K27ac %in% DNAhm]
set7[-which(set7 %in% intersect)] #"Csmd1"  "Dclk3"  "Synm"   "Itpr3"  "Prex1"  "Katnb1"
intersect <- unique(c(intersect, set7))
#H3K4me3, H3K27ac, H3K4me1
set8 <- H3K27ac[H3K27ac %in% H3K4me1 & H3K27ac %in% H3K4me3]
save(set8, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_H3K4me1_H3K27ac_H3K4me3_EE.RData")
set8[-which(set8 %in% intersect)]
#"Gulp1"    "Pim3"     "Mir100hg" "Arid1b"   "Smim23"   "Nfib"     "Chd7"     "Foxk1" 
intersect <- unique(c(intersect, set8))

#####2 'omics
#H3K4me3, DNAm
set9 <- H3K4me3[H3K4me3 %in% DNAm]
set9[-which(set9 %in% intersect)] #"Gm5934" "Wnt3a" 
intersect <- unique(c(intersect, set9))
#H3K4me3, H3K27ac
set10 <- H3K4me3[H3K4me3 %in% H3K27ac]
set10[-which(set10 %in% intersect)] #"Lmo4"   "Gstm1"  "Grrp1"  "Ube2v1"
intersect <- unique(c(intersect, set10))
#H3K4me3, DNAhm
set11 <- H3K4me3[H3K4me3 %in% DNAhm]
set11[-which(set11 %in% intersect)]
#[1] "Nol4l"         "Cyth1"         "Plec"          "Ppp1r16b"      "Inf2"          "Vwa7"         
#[7] "9530036O11Rik" "Abcb9"         "Cnot3"  
intersect <- unique(c(intersect, set11))

################## 12-65 genes per group
#####3 'omics
#DNAm, H3K27ac, H3K4me1
set12 <- H3K4me1[H3K4me1 %in% H3K27ac & H3K4me1 %in% DNAm]
length(set12 <- set12[-which(set12 %in% intersect)]) #16
intersect <- unique(c(intersect, set12))
#DNAm, DNAhm, H3K4me1
set13 <- H3K4me1[H3K4me1 %in% DNAhm & H3K4me1 %in% DNAm]
length(set13 <- set13[-which(set13 %in% intersect)]) #25
intersect <- unique(c(intersect, set13))
#DNAhm, H3K27ac, H3K4me1
set14 <- H3K4me1[H3K4me1 %in% DNAhm & H3K4me1 %in% H3K27ac] #11 genes
save(set14, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_DNAhm_H3K4me1_H3K27ac_EE.RData")

#how many overlap with the SE group?
genes_EE <- set16
load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_DNAhm_H3K4me1_H3K27ac_SE.RData")
genes_SE <- set14
length(genes_SE[genes_SE %in% genes_EE]) #2
#"Cacna1c" "Mgat5" 

length(set14 <- set14[-which(set14 %in% intersect)]) #65
intersect <- unique(c(intersect, set14))

#####2 'omics
#DNAm, H3K27ac
set15 <- DNAm[DNAm %in% H3K27ac]
length(set15 <- set15[-which(set15 %in% intersect)]) #13
intersect <- unique(c(intersect, set15))
#H3K4me3, H3K4me1
set16 <- H3K4me3[H3K4me3 %in% H3K4me1]
length(set16 <- set16[-which(set16 %in% intersect)]) #29
intersect <- unique(c(intersect, set16))
#H3K27ac, DNAhm
set17 <- DNAhm[DNAhm %in% H3K27ac]
length(set17 <- set17[-which(set17 %in% intersect)]) #40
intersect <- unique(c(intersect, set17))
#DNAm, DNAhm
set18 <- DNAhm[DNAhm %in% DNAm]
length(set18 <- set18[-which(set18 %in% intersect)]) #45
intersect <- unique(c(intersect, set18))
#DNAm, H3K4me1
set19 <- DNAm[DNAm %in% H3K4me1]
length(set19 <- set19[-which(set19 %in% intersect)]) #56
intersect <- unique(c(intersect, set19))

#####1 'omic
set20 <- mRNA #13 genes
#length(set20[-which(set20 %in% intersect)]) #56
intersect <- unique(c(intersect, set20))

################## 138-2814 genes per group
#####2 'omics
#DNAhm, H3K4me1
set21 <- DNAhm[DNAhm %in% H3K4me1]
length(set21 <- set21[-which(set21 %in% intersect)]) #217
intersect <- unique(c(intersect, set21))
#H3K27ac, H3K4me1
set22 <- H3K27ac[H3K27ac %in% H3K4me1] #486 genes
save(set22, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_H3K4me1_H3K27ac_EE.RData")

#how many overlap with the SE group?
genes_EE <- set22
load("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_H3K4me1_H3K27ac_SE.RData")
genes_SE <- set24 #462 genes
length(genes_SE[genes_SE %in% genes_EE]) #134

length(set22 <- set22[-which(set22 %in% intersect)]) #386
intersect <- unique(c(intersect, set22))

#####1 'omic
#H3K4me3
set23 <- H3K4me3
length(set23 <- set23[-which(set23 %in% intersect)]) #101
intersect <- unique(c(intersect, set23))
#DNAm
set24 <- DNAm
length(set24 <- set24[-which(set24 %in% intersect)]) #226
intersect <- unique(c(intersect, set24))
#H3K27ac
set25 <- H3K27ac
length(set25 <- set25[-which(set25 %in% intersect)]) #512
intersect <- unique(c(intersect, set25))
#DNAhm
set26 <- DNAhm
length(set26 <- set26[-which(set26 %in% intersect)]) #675
intersect <- unique(c(intersect, set26))
#H3K4me1
set27 <- H3K4me1
length(set27 <- set27[-which(set27 %in% intersect)]) #2814
intersect <- unique(c(intersect, set27))

############GO enrichment
######H3K4me1 & H3K27ac (386 genes)
gr_H3K4me1_H3K27ac <- gost(query = set22, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_H3K27ac$source))
#CORUM GO:BP GO:CC GO:MF  KEGG MIRNA  REAC    TF 
#2   544   112    55    12     2     6   306 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
library(dplyr)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
library(ggplot2)
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC

######H3K4me1 & DNAhm (217 genes)
#How many with gain/loss of DNAhm in TG mice?
load("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/methdiff/hydroxy/BiSeq_site_specific/WTEE_TGEE/BH/betaResults_hydroxy_WTEE_TGEE_BH.RData")
#anno <- read.delim("~/KoborLab/kobor_space/sschaffner/DecipherPD_RRBS/SNCA_OE/Hippocampus/12m_EE/methdiff/hydroxy/BiSeq_site_specific/cpg_anno_chrom.txt")
#anno_sub <- anno[anno$gene %in% set21,]
#length(unique(anno_sub$gene)) #153 genes?
betaResults_sub <- betaResults_WTEE_TGEE_hmC[betaResults_WTEE_TGEE_hmC$gene %in% set21,]
length(unique(betaResults_sub$gene)) #217 genes
betaResults_sub <- betaResults_sub[betaResults_sub$threshold==TRUE,]
length(unique(betaResults_sub$gene)) #217 genes
summary(as.factor(betaResults_sub$DNAm_change))
#Decrease Increase 
#224       37

#How many with gain/loss of H3K4me1 in TG mice?
H3K4me1_sub <- WTEE_TGEE_H3K4me1[WTEE_TGEE_H3K4me1$SYMBOL %in% set21,]
length(unique(H3K4me1_sub$SYMBOL)) #217
summary(as.factor(sign(H3K4me1_sub$Fold)))
#-1   1 
#223 105

gr_H3K4me1_DNAhm <- gost(query = set21, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_DNAhm$source))
#GO:BP GO:CC GO:MF MIRNA  REAC    TF 
#65    29    37     1     1   297 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
#library(dplyr)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
#library(ggplot2)
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC

######H3K4me1, H3K27ac, & DNAhm (65 genes)
gr_H3K4me1_DNAhm <- gost(query = set14, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_DNAhm$source))
#CORUM GO:BP GO:CC GO:MF    HP  KEGG    TF    WP 
#18    46    13   104     1    34   136     1 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
#library(dplyr)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
#library(ggplot2)
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC

######H3K4me1 & DNAm (56 genes)
gr_H3K4me1_DNAm <- gost(query = set19, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_DNAm$source))
#CORUM  KEGG MIRNA    WP 
#3    44     2     2 

#plotting top 10 KEGG terms
KEGG <- gr_H3K4me1_DNAm[gr_H3K4me1_DNAm$source=="KEGG",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
KEGG <- KEGG %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
KEGG$term_name <- factor(KEGG$term_name, levels=KEGG$term_name)
#Num terms on x-axis, coloured by adjP
#library(ggplot2)
KEGG <- ggplot(KEGG[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="palegreen2",high="palegreen4", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("KEGG")
KEGG