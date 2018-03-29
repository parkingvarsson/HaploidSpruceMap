#setwd("~/Dropbox/IntraSpruceDiversity/Exome_data/For_Hap_map")

#install.packages("pegas")
library(pegas)

#install.packages("vcfR")
library(vcfR)

###Data
probe_genome <- ape::read.dna("../../Genomes/Diploid_probes_genome/probes_RG_0802-diploid_scaffolds.fa",format="fasta")
vcf_ann <- read.vcfR("../PopGen_extended_SNPs_biallelic_filtered_LowQual_GATK-hard-filtering_minQD5_minMQ50_maxDP16k_minDP3k_remove_replicated_6samples_missing0.25.ann.vcf.gz", verbose = FALSE)
bed <- read.table("../bed_file/overlapping_probes_with_map_new_with_probe_and_scaffold_info.bed", sep="\t", quote = "",head=T)

###Own script for calculating genetic diversity and segregating sites
check_DNAdist<-function(x){
  
  reduced_matrix<-as.data.frame(as.character(x[,apply(x,2,function(y){length(unique(y))})>1]))
  if (dim(reduced_matrix)[2]==0){
    return(c(0,NA))
  }else{
  
    SNPs<-apply(summary(reduced_matrix),2,function(SNP){as.data.frame(SNP)})
    S<-length(SNPs)
    site_frequency<-rep(0,(dim(reduced_matrix)[1]*2)-1)
  
    het<-c("m","s","r","k","w","y")
  
    for (i in 1:length(SNPs)){
      SNP<-as.data.frame(t(as.data.frame(apply(SNPs[[i]],1,function(x){strsplit(x,split=":")}))))
      if (length(which(is.na(SNP[,1])))>0){
        SNP<-as.data.frame(SNP[!(is.na(SNP[,1])),])
      }
      SNP[,2]<-as.integer(as.character(SNP[,2]))
      SNP[,1]<-as.character(SNP[,1])
    
      if (dim(SNP)[1]==2){
        site_frequency[min(SNP[!(SNP[,1] %in% het),2]*2 + SNP[SNP[,1] %in% het,2], SNP[SNP[,1] %in% het,2])]<- site_frequency[min(SNP[!(SNP[,1] %in% het),2]*2 + SNP[SNP[,1] %in% het,2], SNP[SNP[,1] %in% het,2])] + 1
      }else{
        site_frequency[min(SNP[!(SNP[,1] %in% het),2][1]*2 + SNP[SNP[,1] %in% het,2], SNP[!(SNP[,1] %in% het),2][2]*2 + SNP[SNP[,1] %in% het,2])]<- site_frequency[min(SNP[!(SNP[,1] %in% het),2][1]*2 + SNP[SNP[,1] %in% het,2], SNP[!(SNP[,1] %in% het),2][2]*2 + SNP[SNP[,1] %in% het,2])] + 1
      }
    }
  
    dist<-0
    n<-2*dim(reduced_matrix)[1]
    i<-seq(1:(n-1))
    dist<-sum(i*(n-i)*site_frequency[1:(n-1)])
  
    average_dist<-dist/(n*(n-1)/2)
  
    return(c(S,average_dist))
  }
}


## hacked version of TajimaD function from pegas
neutrality.stats<-function (x) 
{
  n <- if (is.list(x)) 
    length(x)
  else dim(x)[1]
  
  div<-check_DNAdist(x)
  
  S <- div[1]
  if (S==0) {
    warning("no segregating sites")
    return(list(S = 0, thetaW = 0, pi = 0, D = NaN))
  }
  
  if (n < 4) {
    warning("Tajima test requires at least 4 sequences")
    return(list(S = 0, thetaW = 0, pi = 0, D = NaN))
  }
  
  khat <- div[2]
  
  
  tmp <- 1:(n - 1)
  a1 <- sum(1/tmp)
  a2 <- sum(1/tmp^2)
  b1 <- (n + 1)/(3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1 <- b1 - 1/a1
  c2 <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1 <- c1/a1
  e2 <- c2/(a1^2 + a2)
  D <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  
  list(S = S,thetaW = (S/a1)/dim(x)[2], pi = khat/dim(x)[2], D = D)
}


### Function to calculate nautrality stats for each probe separately
calculate_probe_stats<-function(vcf,bed,fasta){
  Stats<-NULL
  Stats$probe<-NULL
  Stats$pi<-NULL
  Stats$S<-NULL
  Stats$D<-NULL
  Stats$thetaW<-NULL
  
  index<-0
  
  for(i in 1:length(unique(bed$scaffold))){
    scaffold<-fasta[grep(paste(as.character(unique(bed$scaffold)[i])," ",sep=""), names(fasta))]
    names(scaffold) <- as.character(unique(bed$scaffold)[i])
    scaffold <- as.matrix(scaffold)
    
    scaffold_probes<-bed[bed$scaffold == unique(bed$scaffold)[i],]
    
    for(j in 1:dim(scaffold_probes)[1]){  
      print(c(i,j))
      if(scaffold_probes$stop[j] > length(scaffold)){
        scaffold_probes$stop[j] <-length(scaffold)
      }
      probe_vcf<-vcfR2DNAbin(vcf[which(vcf@fix ==as.character(unique(bed$scaffold)[i])),], consensus = TRUE, extract.haps = FALSE, ref.seq = scaffold[,scaffold_probes[j,2]:scaffold_probes[j,3]], start.pos = scaffold_probes[j,2], verbose = FALSE)
      
      index<-index +1  
      
      Stats$probe[index]<-as.character(scaffold_probes$scaf_probe[j])
      
      ##calculate individuals with missing data
      idx <- apply(probe_vcf,1,function(x){sum(x==240)})>0
      
      if(dim(probe_vcf[!idx,])[1]>0){
        
        D<-neutrality.stats(probe_vcf[!idx,])
        Stats$S[index]<-D$S
        Stats$thetaW[index]<-D$thetaW
        Stats$pi[index]<-D$pi
        Stats$D[index]<-D$D
      }else{
        Stats$S[index]<-0
        Stats$thetaW[index]<-0
        Stats$pi[index]<-0
        Stats$D[index]<-NA
      }
    }
  }
  Stats<-as.data.frame(Stats)
  return(Stats)
}

All_probes_Stats<-calculate_probe_stats(vcf_ann,bed,probe_genome)

save.image("Probe_neutrality.RData")

write.table(All_probes_Stats, "All_probes_Stats.txt", quote=F,row.names = F, col.names = T, sep="\t")




###Just to try script on a specific scaffold and probe
#scaffold<-probe_genome[grep(paste(as.character(unique(bed$scaffold)[1])," ",sep=""), names(probe_genome))]
#names(scaffold) <- as.character(unique(bed$scaffold)[1])
#scaffold <- as.matrix(scaffold)
#scaffold_probes<-bed[bed$scaffold == as.character(unique(bed$scaffold))[1],]
#probe_vcf<-vcfR2DNAbin(vcf_ann[which(vcf_ann@fix ==as.character(unique(scaffold_probes$scaffold))),], consensus = TRUE, extract.haps = FALSE, ref.seq = scaffold[,scaffold_probes[1,2]:scaffold_probes[1,3]], start.pos = scaffold_probes[1,2], verbose = FALSE)
#idx <- apply(probe_vcf,1,function(x){sum(x==240)})>0
#reduced_matrix<-as.data.frame(as.character(probe_vcf[!idx,apply(probe_vcf[!idx,],2,function(y){length(unique(y))})>1]))









