CGI_Gene_permutation_enrichment<-function(CpG_list,background.probes, permutation_number, plotmin, plotmax){
  
  #create a data frame with your hit list of probes assigned to genomic contexts, and remove duplicates
  Genes_correlated_CpGs<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$site%in%CpG_list),]
  Genes_correlated_CpGs<-Genes_correlated_CpGs[!duplicated(Genes_correlated_CpGs),]
  
  #calculate the number of probes in your hit list that occur in each genomic context and make a data frame
  #if no probes map to a specific feature, add this as a blank entry in the data frame
  Gene_hits_regionMeans<-tapply(Genes_correlated_CpGs$site, Genes_correlated_CpGs$anno, length)
  Gene_hits_regionMeans[is.na(Gene_hits_regionMeans)]<-0
  if (length(grep("mm10_genes_intergenic", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_intergenic"=0)
  }
  if (length(grep("mm10_enhancers_fantom", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_enhancers_fantom"=0)
  }
  if (length(grep("mm10_genes_1to5kb", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_1to5kb"=0)
  }
  if (length(grep("mm10_genes_3UTRs", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_3UTRs"=0)
  }
  if (length(grep("mm10_genes_5UTRs", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_5UTRs"=0)
  }
  if (length(grep("mm10_genes_introns", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_introns"=0)
  }
  if (length(grep("enhancers_poised", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "enhancers_poised"=0)
  }
  if (length(grep("mm10_genes_exonintronboundaries", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_exonintronboundaries"=0)
  }
  if (length(grep("mm10_genes_exons", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_exons"=0)
  }
  if (length(grep("mm10_genes_intronexonboundaries", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_intronexonboundaries"=0)
  }  
  if (length(grep("enhancers_active", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "enhancers_active"=0)
  }
  if (length(grep("mm10_genes_promoters", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_genes_promoters"=0)
  }
  if (length(grep("promoters_active", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "promoters_active"=0)
  }
  if (length(grep("mm10_lncrna_gencode", names(Gene_hits_regionMeans)))==0){
    Gene_hits_regionMeans <- c(Gene_hits_regionMeans, "mm10_lncrna_gencode"=0)
  }
  Gene_hits_regionMeans<-data.frame(Site_Count=as.numeric(Gene_hits_regionMeans), Region=names(Gene_hits_regionMeans))
  
  ## Boot strapping (to see if hits more in feature than expected)
  bootstrap_genes<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    #take a random sample of all probes, same size as your hit list
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    #get the probe --> genomic context association for the random sample of probes
    Gene_rnd<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$site%in%rnd_CpGs),]
    Gene_rnd<-Gene_rnd[!duplicated(Gene_rnd),]
    
    #caluclate the number of probes in random list that occur in each genomic context
    Gene_rnd_regionMeans<-tapply(Gene_rnd$site, Gene_rnd$anno, length)
    Gene_rnd_regionMeans[is.na(Gene_rnd_regionMeans)]<-0
    #check for missing annotations and add as blank columns
    if (length(grep("mm10_genes_intergenic", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_intergenic"=0)
    }
    if (length(grep("mm10_enhancers_fantom", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_enhancers_fantom"=0)
    }
    if (length(grep("mm10_genes_1to5kb", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_1to5kb"=0)
    }
    if (length(grep("mm10_genes_3UTRs", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_3UTRs"=0)
    }
    if (length(grep("mm10_genes_5UTRs", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_5UTRs"=0)
    }
    if (length(grep("mm10_genes_introns", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_introns"=0)
    }
    if (length(grep("enhancers_poised", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "enhancers_poised"=0)
    }
    if (length(grep("mm10_genes_exonintronboundaries", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_exonintronboundaries"=0)
    }
    if (length(grep("mm10_genes_exons", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_exons"=0)
    }
    if (length(grep("mm10_genes_intronexonboundaries", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_intronexonboundaries"=0)
    }  
    if (length(grep("enhancers_active", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "enhancers_active"=0)
    }
    if (length(grep("mm10_genes_promoters", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_genes_promoters"=0)
    }
    if (length(grep("promoters_active", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "promoters_active"=0)
    }
    if (length(grep("mm10_lncrna_gencode", names(Gene_rnd_regionMeans)))==0){
      Gene_rnd_regionMeans <- c(Gene_rnd_regionMeans, "mm10_lncrna_gencode"=0)
    }
    
    Gene_rnd_regionMeans <- Gene_rnd_regionMeans[match(Gene_hits_regionMeans$Region, names(Gene_rnd_regionMeans))]
    Gene_rnd_regionMeans
  })
  #bind together each permutation
  bootstrap_genes<-do.call(rbind, bootstrap_genes)
  
  print("FDR Adjusted Permutation P values for enrichment and depletion")
  Region_results<-sapply(1:nrow(Gene_hits_regionMeans), function(x){
    #loop over each genomic context
    #get the number of probes in that genomic region which showed up in hit list
    real_CpG_in_region<-Gene_hits_regionMeans$Site_Count[x]
    #how many iterations of bootstrapping found MORE cpgs in that region than the real data? Divide this by the number of permutations to get a p-value, then adjust this for multiple testing across the number of genomic regions tested
    Adjusted_enrich_p<-p.adjust(length(which(bootstrap_genes[,x]>=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Gene_hits_regionMeans))
    #how many iterations of bootstrapping found LESS cpgs in that region than the real data? Divide this by the number of permutations to get a p-value, then adjust this for multiple testing across the number of genomic regions tested
    Adjusted_depletion_p<-p.adjust(length(which(bootstrap_genes[,x]<=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Gene_hits_regionMeans))
    print(paste("Enrichment: ", Adjusted_enrich_p, "; Depletion ", Adjusted_depletion_p, "; Feature: ", Gene_hits_regionMeans$Region[x],sep=""))
  })
  
  
  ## CGI
  CGI<-cpg_anno_cgi
  
  #get the probe - CGI associations for your hit list
  Resort_hits<-CGI[which(CGI$site%in%CpG_list),]
  
  #get the number of probes that match to each CGI context for hit list
  #Resorts_hits_featureMeans should be a data frame with a column for Site_Count
  #There should be a row for every feature, regardless of if it occurs in your data (e.g. Site_Count column would be 0 in this case)
  Resort_hits_featureMeans<-tapply(Resort_hits$site, Resort_hits$anno, length)
  Resort_hits_featureMeans[is.na(Resort_hits_featureMeans)]<-0
  
  #add in the "zero" features at this stage
  if (!("mm10_cpg_islands" %in% names(Resort_hits_featureMeans))){
    feature <- 0
    names(feature) <- "mm10_cpg_islands"
    Resort_hits_featureMeans <- c(Resort_hits_featureMeans, feature)
  }
  if (!("mm10_cpg_inter" %in% names(Resort_hits_featureMeans))){
    feature <- 0
    names(feature) <- "mm10_cpg_inter"
    Resort_hits_featureMeans <- c(Resort_hits_featureMeans, feature)
  }
  if (!("mm10_cpg_shelves" %in% names(Resort_hits_featureMeans))){
    feature <- 0
    names(feature) <- "mm10_cpg_shelves"
    Resort_hits_featureMeans <- c(Resort_hits_featureMeans, feature)
  }  
  if (!("mm10_cpg_shores" %in% names(Resort_hits_featureMeans))){
    feature <- 0
    names(feature) <- "mm10_cpg_shores"
    Resort_hits_featureMeans <- c(Resort_hits_featureMeans, feature)
  }
  
  
  Resort_hits_featureMeans<-data.frame(Site_Count=as.numeric(Resort_hits_featureMeans),
                                       Feature=names(Resort_hits_featureMeans))
  levels(Resort_hits_featureMeans$Feature)<-c("mm10_cpg_inter","mm10_cpg_islands","mm10_cpg_shelves","mm10_cpg_shores")
  
  ## Boot strapping (to see if hits more in feature than expected)
  bootstrap_CGI<-lapply(1:permutation_number, function(x){
    set.seed(x)
    Hit_number<-length(CpG_list)
    rnd_CpGs<-background.probes[sample(1:length(background.probes),Hit_number)]
    Resort_rnd<-CGI[which(CGI$site%in%rnd_CpGs),]
    Resort_rnd_featureMeans<-tapply(Resort_rnd$site, Resort_rnd$anno, length)
    Resort_rnd_featureMeans[is.na(Resort_rnd_featureMeans)]<-0
    #add in the "zero" features at this stage
    if (!("mm10_cpg_islands" %in% names(Resort_rnd_featureMeans))){
      feature <- 0
      names(feature) <- "mm10_cpg_islands"
      Resort_rnd_featureMeans <- c(Resort_rnd_featureMeans, feature)
    }
    if (!("mm10_cpg_inter" %in% names(Resort_rnd_featureMeans))){
      feature <- 0
      names(feature) <- "mm10_cpg_inter"
      Resort_rnd_featureMeans <- c(Resort_rnd_featureMeans, feature)
    }
    if (!("mm10_cpg_shelves" %in% names(Resort_rnd_featureMeans))){
      feature <- 0
      names(feature) <- "mm10_cpg_shelves"
      Resort_rnd_featureMeans <- c(Resort_rnd_featureMeans, feature)
    }  
    if (!("mm10_cpg_shores" %in% names(Resort_rnd_featureMeans))){
      feature <- 0
      names(feature) <- "mm10_cpg_shores"
      Resort_rnd_featureMeans <- c(Resort_rnd_featureMeans, feature)
    }
    Resort_rnd_featureMeans <- Resort_rnd_featureMeans[match(Resort_hits_featureMeans$Feature, names(Resort_rnd_featureMeans))]
    Resort_rnd_featureMeans
  })
  bootstrap_CGI<-do.call(rbind,bootstrap_CGI)
  
  
  CGI_results<-sapply(1:nrow(Resort_hits_featureMeans), function(x){
    real_CpG_in_region<-Resort_hits_featureMeans$Site_Count[x]
    Adjusted_enrich_p<-p.adjust(length(which(bootstrap_CGI[,x]>=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Resort_hits_featureMeans))
    Adjusted_depletion_p<-p.adjust(length(which(bootstrap_CGI[,x]<=real_CpG_in_region))/permutation_number, method="fdr", n=nrow(Resort_hits_featureMeans))
    print(paste("Enrichment: ", Adjusted_enrich_p, "; Depletion: ", Adjusted_depletion_p, "; Feature: ", Resort_hits_featureMeans$Feature[x], sep=""))
  })
  
  #c(Region_results,CGI_results)
  
  
}