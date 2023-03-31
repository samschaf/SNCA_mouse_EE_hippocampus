overlap_permutation_2way <- function(CpG_list_1, CpG_list_2, background.probes, permutation_number, Group1, Group2){
  
  #real overlaps
  real_overlap12 <- length(CpG_list_1[CpG_list_1 %in% CpG_list_2])
  
  #collapsing to numbers
  #real_overlap12 <- length(real_overlap12[-which(real_overlap12 %in% c(real_overlap1234, real_overlap123, real_overlap124))])
  
  ################ Boot strapping for Group1-Group2-Group3-Group4 overlap
  bootstrap_overlap12 <-lapply(1:permutation_number, function(x){
    set.seed(x)
    
    #take 4 random samples of all probes, same size as your hit lists
    hitlist1 <-background.probes[sample(1:length(background.probes),length(CpG_list_1))]
    hitlist2 <-background.probes[sample(1:length(background.probes),length(CpG_list_2))]
  
    #get lists of overlaps for the random sample of probes
    rnd_overlap12 <- length(hitlist1[hitlist1 %in% hitlist2])
    #rnd_overlap13 <- hitlist1[hitlist1 %in% hitlist3]
    #rnd_overlap14 <- hitlist1[hitlist1 %in% hitlist4]
    #rnd_overlap23 <- hitlist2[hitlist2 %in% hitlist3]
    #rnd_overlap24 <- hitlist2[hitlist2 %in% hitlist4]
    #rnd_overlap123 <- rnd_overlap12[rnd_overlap12 %in% hitlist3]
    #rnd_overlap124 <- rnd_overlap12[rnd_overlap12 %in% hitlist4]
    #rnd_overlap134 <- rnd_overlap13[rnd_overlap13 %in% hitlist4]
    #rnd_overlap234 <- rnd_overlap23[rnd_overlap23 %in% hitlist4]
    #rnd_overlap1234 <- rnd_overlap123[rnd_overlap123 %in% hitlist4]
    
    #get number of overlaps for the random sample of probes
    #rnd_overlap12 <- length(rnd_overlap12[-which(rnd_overlap12 %in% c(rnd_overlap1234, rnd_overlap123, rnd_overlap124))])
    #rnd_overlap13 <- length(rnd_overlap13[-which(rnd_overlap13 %in% c(rnd_overlap1234, rnd_overlap123, rnd_overlap134))])
    #rnd_overlap14 <- length(rnd_overlap14[-which(rnd_overlap14 %in% c(rnd_overlap1234, rnd_overlap124, rnd_overlap134))])
    #rnd_overlap23 <- length(rnd_overlap23[-which(rnd_overlap23 %in% c(rnd_overlap1234, rnd_overlap123, rnd_overlap234))])
    #rnd_overlap123 <- length(rnd_overlap123[-which(rnd_overlap123 %in% rnd_overlap1234)])
    #rnd_overlap124 <- length(rnd_overlap124[-which(rnd_overlap124 %in% rnd_overlap1234)])
    #rnd_overlap134 <- length(rnd_overlap134[-which(rnd_overlap134 %in% rnd_overlap1234)])
    #rnd_overlap234 <- length(rnd_overlap234[-which(rnd_overlap234 %in% rnd_overlap1234)])
    #rnd_overlap1234 <- length(rnd_overlap1234)
    
    rnd_overlap12
  })
  #bind together each permutation
  bootstrap_overlap12<-do.call(rbind, bootstrap_overlap12)
  
  print("Permutation P values for enrichment and depletion")
  
  #how many iterations of bootstrapping found MORE overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  enrich_p12 <-length(which(bootstrap_overlap12>=real_overlap12))/permutation_number
  
  #how many iterations of bootstrapping found LESS overlapping probes than the real number? Divide this by the number of permutations to get a p-value
  depletion_p12 <-length(which(bootstrap_overlap12<=real_overlap12))/permutation_number
  
  print(paste("Enrichment ", Group1, "-", Group2, ": ", enrich_p12, "; Depletion ", Group1, "-", Group2, ": ", depletion_p12, sep=""))
  
}

