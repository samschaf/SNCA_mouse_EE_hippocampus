##making a nested list of all genes

####DNAm/DNAhm data
lm_WTSE_TGSE_mC <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/methylation/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_mC <- lm_WTSE_TGSE_mC[lm_WTSE_TGSE_mC$threshold==TRUE,]

lm_WTSE_TGSE_hmC <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/5-methdiff_DNAm_DNAhm/hydroxy/WTSE_TGSE/linear_model/auto/CpG/lm_WTSE_TGSE.csv")
WTSE_TGSE_hmC <- lm_WTSE_TGSE_hmC[lm_WTSE_TGSE_hmC$threshold==TRUE,]

#retrieving DM/DHM genes
DNAm <- unique(WTSE_TGSE_mC$gene[complete.cases(WTSE_TGSE_mC$gene)]) #0 genes
DNAhm <- unique(WTSE_TGSE_hmC$gene[complete.cases(WTSE_TGSE_hmC$gene)]) #315 genes

####Histone PTM ChIP-seq data (only includes significant differentially bound regions, FDR < 0.05)
WTSE_TGSE_H3K4me1 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me1/H3K4me1_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTSE_TGSE_H3K27ac <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K27ac/H3K27ac_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")
WTSE_TGSE_H3K4me3 <- read.delim("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_sig_only/H3K4me3/H3K4me3_WTSE-TGSE_cons-2_FDR0_05_FC0_5_annot.tsv")

#subsetting to genes
H3K4me1 <- unique(as.character(WTSE_TGSE_H3K4me1$SYMBOL[complete.cases(WTSE_TGSE_H3K4me1$SYMBOL)])) #3,437 genes
H3K4me3 <- unique(as.character(WTSE_TGSE_H3K4me3$SYMBOL[complete.cases(WTSE_TGSE_H3K4me3$SYMBOL)])) #10 genes
H3K27ac <- unique(as.character(WTSE_TGSE_H3K27ac$SYMBOL[complete.cases(WTSE_TGSE_H3K27ac$SYMBOL)])) #1,024 genes

####RNA-seq data
WTSE_TGSE_RNA_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/6-processed_RNAseq_data/mmus_hippo12m_TG_STD_vs_WT_STD.csv")
WTSE_TGSE_RNA <- WTSE_TGSE_RNA_bg[abs(WTSE_TGSE_RNA_bg$log2FoldChange)>=0.3 & WTSE_TGSE_RNA_bg$padj<=0.15,]
mRNA <- unique(as.character(WTSE_TGSE_RNA$mgi_symbol[complete.cases(WTSE_TGSE_RNA$mgi_symbol)])) #244 genes

##putting it all together
omics <- list(DNAm=DNAm, DNAhm=DNAhm, H3K4me3=H3K4me3, H3K4me1=H3K4me1, H3K27ac=H3K27ac, mRNA=mRNA)

##########Upset plot
library(UpSetR)
#upset(fromList(omics), nsets=6, order.by="freq")

omics_meta <- data.frame(sets=names(omics))
omics_meta$type <- c("DNA_mod", "DNA_mod", rep("Histone_PTM",3), "mRNA")

upset <- upset(fromList(omics), nsets=6, order.by="freq", set.metadata=list(data=omics_meta, plots=list(list(type="matrix_rows", column="type", colors=c(DNA_mod="blue", Histone_PTM="red", mRNA="yellow")))))

##########permuting overlaps
source('/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/8-multi_omic_integration/overlap/overlap_permutation_6way.R')
genes.RNA=unique(as.character(WTSE_TGSE_RNA_bg$mgi_symbol[complete.cases(WTSE_TGSE_RNA_bg$mgi_symbol)])) 
genes.DNAm=unique(lm_WTSE_TGSE_mC$gene[complete.cases(lm_WTSE_TGSE_mC$gene)])
genes.DNAhm=unique(lm_WTSE_TGSE_hmC$gene[complete.cases(lm_WTSE_TGSE_hmC$gene)])
H3K4me1_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K4me1/H3K4me1_WTSE-TGSE_anno.csv")
genes.K4me1 <- unique(H3K4me1_bg$SYMBOL[complete.cases(H3K4me1_bg$SYMBOL)])
H3K27ac_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K27ac/H3K27ac_WTSE-TGSE_anno.csv")
genes.K27ac <- unique(H3K27ac_bg$SYMBOL[complete.cases(H3K27ac_bg$SYMBOL)])
H3K4me3_bg <- read.csv("/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/7-ChIPseq/Zinah_Wassouf_data/diffbind_2_reps_all_conspeaks/H3K4me3/H3K4me3_WTSE-TGSE_anno.csv")
genes.K4me3 <- unique(H3K4me3_bg$SYMBOL[complete.cases(H3K4me3_bg$SYMBOL)])

perm <- overlap_permutation_6way(mRNA, genes.RNA, 
  DNAm, genes.DNAm,
  DNAhm, genes.DNAhm,
  H3K4me1, genes.K4me1,
  H3K4me3, genes.K4me3,
  H3K27ac, genes.K27ac, 
  permutation_number=1000)
save(perm, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/WTSE_TGSE_overlap_perm_122022.RData")

##########Plotting examples of co-regulated genes
#Any IEGs?
IEGs <- c("Egr1", "Nr4a2", "Nurr1", "Tyro3", "Arc", "Homer1a")
IEG_ind <- sapply(IEGs, function(x) names(omics[grep(x, omics)]))
#$Nr4a2
#[1] "H3K4me1"

#####H3K4me1 & H3K27ac (167)
#GO enrichment
library(gprofiler2)

#lists of genes

################## 1 gene per group
#####4 'omics
#mRNA, DNAhm, H3K4me1, H3K27ac
set1 <- mRNA[mRNA %in% DNAhm & mRNA %in% H3K4me1 & mRNA %in% H3K27ac] #"Slit3"
#DNAm, mRNA, DNAhm, H3K4me1
set2 <- mRNA[mRNA %in% DNAhm & mRNA %in% H3K4me1 & mRNA %in% DNAm] #"Prdm16"

#####3 'omics
#H3K4me3, H3K4me1, H3K27ac
set3 <- H3K4me1[H3K4me1 %in% H3K4me3 & H3K4me1 %in% H3K27ac] #"Gulp1"
intersect <- c("Slit3", "Prdm16", "Gulp1")
#DNAm, mRNA, H3K4me1
set7 <- H3K4me1[H3K4me1 %in% mRNA & H3K4me1 %in% DNAm]
set7[-which(set7 %in% intersect)] #"Hpcal1"
intersect <- unique(c(intersect, set7))

#####2 'omics
#H3K4me3, H3K4me1
set4 <- H3K4me1[H3K4me1 %in% H3K4me3]
set4[-which(set4 %in% intersect)] #"Eef1a1"
intersect <- unique(c(intersect, set4))
#H3K4me3, H3K27ac
set5 <- H3K27ac[H3K27ac %in% H3K4me3]
set5[-which(set5 %in% intersect)] #"Abcg1"
intersect <- unique(c(intersect, set5))
#H3K4me3, DNAhm
set6 <- H3K4me3[H3K4me3 %in% DNAhm] #"Hnmt"
#set6[-which(set6 %in% intersect)] #"Abcg1"
intersect <- unique(c(intersect, set6))

################## 2-10 genes per group
#####3 'omics
#mRNA, DNAhm, H3K27ac
set8 <- H3K27ac[H3K27ac %in% mRNA & H3K27ac %in% DNAhm]
set8[-which(set8 %in% intersect)] #"Gpr83"  "Als2cl"
intersect <- unique(c(intersect, set8))
#mRNA, H3K27ac, H3K4me1
set9 <- H3K27ac[H3K27ac %in% mRNA & H3K27ac %in% H3K4me1]
set9[-which(set9 %in% intersect)] #"Tmem204" "Syne1"   "Gbp3"    "Mdga1" 
intersect <- unique(c(intersect, set9))
#DNAm, H3K27ac, H3K4me1
set10 <- H3K27ac[H3K27ac %in% DNAm & H3K27ac %in% H3K4me1] 
#"Pced1b"        "Galnt18"       "Xkr6"          "Slc20a2"       "1700057H15Rik"
#set10[-which(set10 %in% intersect)]
intersect <- unique(c(intersect, set10))
#DNAm, DNAhm, H3K4me1
set11 <- H3K4me1[H3K4me1 %in% DNAm & H3K4me1 %in% DNAhm]
set11[-which(set11 %in% intersect)]
#"Wwc1"    "Dscam"   "Cdh23"   "Zfp629"  "Zfp64"   "Mapk4"   "Map7"    "Ccdc88c" "Efhd1" 
intersect <- unique(c(intersect, set11))
#mRNA, DNAhm, H3K4me1
set12 <- H3K4me1[H3K4me1 %in% mRNA & H3K4me1 %in% DNAhm]
set12[-which(set12 %in% intersect)]
#[1] "Otof"      "Arhgef10l" "Cdh4"      "Flnb"      "Mgat5b"    "Atp10a"    "Cspg4"     "Snx10"    
#[9] "Srrm4"     "Cacna1i"
intersect <- unique(c(intersect, set12))
save(set12, file="/mnt/scratch/KoborLab/DecipherPD_mouse/sschaffner/ee_asyn/hipp/8-multi_omic_integration/overlap/genes_mRNA_H3K4me1_DNAhm.RData")


#####2 'omics
#DNAm, H3K27ac
set13 <- H3K27ac[H3K27ac %in% DNAm]
set13[-which(set13 %in% intersect)] #"Hif3a"   "Ptprf"   "Foxk1"   "Adora2b"
intersect <- unique(c(intersect, set13))
#mRNA, H3K27ac
set14 <- H3K27ac[H3K27ac %in% mRNA]
save(set14, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_mRNA_H3K27ac.RData")
set14[-which(set14 %in% intersect)]
# "Rgs4"    "Apod"    "Mmd"     "Pitpnm3" "Epas1"   "Nptx2"   "Kcnk13"  "Megf11"
intersect <- unique(c(intersect, set14))

#####1 'omic
#H3K4me3
set15 <- H3K4me3
set15[-which(set15 %in% intersect)]
#"1700081H04Rik" "Cbx4"          "Zfp985"        "Efna5"         "Slitrk2"       "Cyp3a16"   
intersect <- unique(c(intersect, set15))

################## 13-78 genes per group
#####3 'omics
#DNAhm, H3K4me1, H3K27ac
set16 <- H3K27ac[H3K27ac %in% DNAhm & H3K27ac %in% H3K4me1] #14 genes
save(set16, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_DNAhm_H3K4me1_H3K27ac.RData")

length(set16 <- set16[-which(set16 %in% intersect)]) #62
intersect <- unique(c(intersect, set16))

#####2 'omics
#mRNA, DNAhm
set17 <- DNAhm[DNAhm %in% mRNA]
length(set17 <- set17[-which(set17 %in% intersect)]) #13
save(set17, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_mRNA_DNAhm.RData")
intersect <- unique(c(intersect, set17))                   
#DNAm, DNAhm
set18 <- DNAhm[DNAhm %in% DNAm]
length(set18 <- set18[-which(set18 %in% intersect)]) #13
intersect <- unique(c(intersect, set18))         
#DNAm, H3K4me1
set19 <- DNAm[DNAm %in% H3K4me1]
length(set19 <- set19[-which(set19 %in% intersect)]) #17
intersect <- unique(c(intersect, set19))        
#DNAhm, H3K27ac
set20 <- DNAhm[DNAhm %in% H3K27ac]
length(set20 <- set20[-which(set20 %in% intersect)]) #36
intersect <- unique(c(intersect, set20))  
#mRNA, H3K4me1
set21 <- mRNA[mRNA %in% H3K4me1]
length(set21 <- set21[-which(set21 %in% intersect)]) #47
intersect <- unique(c(intersect, set21))   
save(set21, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_mRNA_H3K4me1_47.RData")


#####1 'omic
#DNAm
set22 <- DNAm
length(set22 <- set22[-which(set22 %in% intersect)]) #78
intersect <- unique(c(intersect, set22))   

################## 167-2732 genes per group
#####2 'omics
#DNAhm, H3K4me1
set23 <- DNAhm[DNAhm %in% H3K4me1]
length(set23 <- set23[-which(set23 %in% intersect)]) #167
intersect <- unique(c(intersect, set23))   
#H3K4me1, H3K27ac
set24 <- H3K27ac[H3K27ac %in% H3K4me1] #462 genes
save(set24, file="/mnt/scratch/KoborLab/DecipherPD_mouse/ee_asyn/hipp/sschaffner/082022/8-multi_omic_integration/overlap/genes_H3K4me1_H3K27ac.RData")

length(set24 <- set24[-which(set24 %in% intersect)]) #389
intersect <- unique(c(intersect, set24))   

#####1 'omic
#mRNA
set25 <- mRNA
length(set25 <- set25[-which(set25 %in% intersect)]) #167
intersect <- unique(c(intersect, set25))   
#H3K27ac
set26 <- H3K27ac
length(set26 <- set26[-which(set26 %in% intersect)]) #511
intersect <- unique(c(intersect, set26))   
#DNAhm
set27 <- DNAhm
length(set27 <- set27[-which(set27 %in% intersect)]) #616
intersect <- unique(c(intersect, set27))   
#H3K4me1
set28 <- H3K4me1
length(set28 <- set28[-which(set28 %in% intersect)]) #2732
intersect <- unique(c(intersect, set28))   

############GO enrichment
######H3K4me1 & H3K27ac (389 genes)
gr_H3K4me1_H3K27ac <- gost(query = set24, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_H3K27ac$source))
#GO:BP GO:CC GO:MF  KEGG MIRNA  REAC    TF 
#572   119    57    34     2     1   329 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
library(dplyr)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
library(ggplot2)
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_H3K27ac[gr_H3K4me1_H3K27ac$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC

######H3K4me1 & DNAhm (167 genes)
gr_H3K4me1_DNAhm <- gost(query = set23, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_DNAhm$source))
#CORUM GO:BP GO:CC GO:MF  KEGG    TF 
#17    28    40     5     1   238 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:5,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC

######H3K4me1, H3K27ac, & DNAhm (62 genes)
gr_H3K4me1_DNAhm <- gost(query = set16, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result
summary(as.factor(gr_H3K4me1_DNAhm$source))
#CORUM GO:BP GO:CC GO:MF  KEGG    TF 
#12   262    26    43     2    38 

#plotting top 10 GO:BP terms
GOBP <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:BP",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOBP <- GOBP %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)
#Num terms on x-axis, coloured by adjP
GOBP <- ggplot(GOBP[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:BP")

#plotting top 10 GO:CC terms
GOCC <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:CC",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOCC <- GOCC %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOCC$term_name <- factor(GOCC$term_name, levels=GOCC$term_name)
#Num terms on x-axis, coloured by adjP
GOCC <- ggplot(GOCC[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:CC")

#plotting top 10 GO:MF terms
GOMF <- gr_H3K4me1_DNAhm[gr_H3K4me1_DNAhm$source=="GO:MF",]
#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOMF <- GOMF %>% arrange(desc(p_value)) #reverse order, so most sig will be at top of bar graph
GOMF$term_name <- factor(GOMF$term_name, levels=GOMF$term_name)
#Num terms on x-axis, coloured by adjP
GOMF <- ggplot(GOMF[1:10,], aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank()) + ggtitle("GO:MF")

GOBP
GOMF
GOCC






#####H3K4me1, H3K27ac, & DNAhm (62)
overlap <- H3K4me1[H3K4me1 %in% H3K27ac & H3K4me1 %in% DNAhm]
#example: Ank3

gr <- gost(query = overlap, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result

summary(as.factor(gr$source))
#CORUM GO:BP GO:CC GO:MF  KEGG    TF 
#12   268    25    43     2    46

#plotting top 10 GO:BP terms
GOBP <- gr[gr$source=="GO:BP",]

#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOBP <- GOBP[1:10,] %>% arrange(desc(p_value))#reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)

#Num terms on x-axis, coloured by adjP
ggplot(GOBP, aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank())

#####H3K4me1 & mRNA (37)
overlap <- H3K4me1[H3K4me1 %in% mRNA]

gr <- gost(query = overlap, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result

summary(as.factor(gr$source))
#CORUM GO:MF    TF 
#6     1    86

#####H3K27ac & DNAhm (36)
overlap <- H3K27ac[H3K27ac %in% DNAhm]

gr <- gost(query = overlap, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result

summary(as.factor(gr$source))
#CORUM GO:BP GO:CC GO:MF  KEGG  REAC    TF 
#19   313    31    48     1     1   121 

#plotting top 10 GO:BP terms
GOBP <- gr[gr$source=="GO:BP",]

#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOBP <- GOBP[1:10,] %>% arrange(desc(p_value))#reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)

#Num terms on x-axis, coloured by adjP
ggplot(GOBP, aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank())

#H3K4me1 & DNAm (17)
overlap <- H3K4me1[H3K4me1 %in% DNAhm]

gr <- gost(query = overlap, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", sources = NULL, evcodes=TRUE)$result

summary(as.factor(gr$source))
#GO:BP GO:CC GO:MF  KEGG    TF    WP 
#264    72    46     2   307     1 

#plotting top 10 GO:BP terms
GOBP <- gr[gr$source=="GO:BP",]

#arrange by p-value and factor GO term description (so that levels will be ordered correctly for plotting)
GOBP <- GOBP[1:10,] %>% arrange(desc(p_value))#reverse order, so most sig will be at top of bar graph
GOBP$term_name <- factor(GOBP$term_name, levels=GOBP$term_name)

#Num terms on x-axis, coloured by adjP
ggplot(GOBP, aes(x=term_name, y=intersection_size, fill=p_value))  + scale_fill_continuous(low="gray48",high="gray87", name="adjP") + geom_bar(stat="identity", position="dodge")  +
  coord_flip() + theme_classic() + scale_y_continuous(name="Number of Genes")  + scale_colour_manual(values="black", guide="none")  + theme(axis.title.y=element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y=element_blank())
