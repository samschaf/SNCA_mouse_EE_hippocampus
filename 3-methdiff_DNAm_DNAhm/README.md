This folder contains R scripts used to pre-process the DNAm and DNAhm data used in the manuscript (subfolder "preprocessing"), and to analyze differential methylation (subfolder "methylation") and hydroxymethylation (subfolder "hydroxy").

lambda_bisulfite_check.Rmd: Used to read in methratio files aligned to the phage lambda genome and check for successful bisulfite conversion in each sample (<1% methylation level).

CGI_gene_permutation_enrichment.R: Used to permute enrichment/depletion of differentially (hydroxy)methylated CpGs for gene and CpG island features, relative to the background of CpGs assayed.

overlap_permutation_2way.R: Used to permute enrichment/depletion of the overlap between differentially (hydroxy)methylated CpGs from two comparisons (TGSE vs WTSE and TGEE vs WTEE).
