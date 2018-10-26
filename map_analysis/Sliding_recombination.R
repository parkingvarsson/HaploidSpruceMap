###Recombination hotspots/coldspots
##Read in consensus map
read.table("../Consensus_maps.txt",head=T)->maps

##Read in scaffold length
read.table("haploid_scaffolds.txt",head=F)->scaf_length
names(scaf_length)<-c("Scaffold","length")

### read in gene models
read.table("~/Dropbox/work/Spruce_projects/reference_genome/Eugene.gff3")->all_models
read.table("~/Dropbox/work/Spruce_projects/reference_genome/HC_genes.gff3")->HC_models
read.table("~/Dropbox/work/Spruce_projects/reference_genome/MC_genes.gff3")->MC_models
read.table("~/Dropbox/work/Spruce_projects/reference_genome/LC_genes.gff3")->LC_models

HC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(HC_models$V9),";")))[,1])),"=")))[,2])
MC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(MC_models$V9),";")))[,1])),"=")))[,2])
LC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(LC_models$V9),";")))[,1])),"=")))[,2])

library(gdata)
HC_models<-drop.levels(HC_models[HC_models$V1 %in% as.character(maps$Scaffold),])
MC_models<-drop.levels(MC_models[MC_models$V1 %in% as.character(maps$Scaffold),])
LC_models<-drop.levels(LC_models[LC_models$V1 %in% as.character(maps$Scaffold),])
map_models<-rbind(HC_models,MC_models,LC_models)


## Function that counts the number of markers, scaffolds and total lenght of all scaffolds in a predefined window size along the LGs
Recomb_spots<-function(map,LG,window_size){
  Pos<-NULL
  Markers<-NULL
  Scaffolds<-NULL
  Length<-NULL
  Genes<-NULL
  tmp<-map[map$LG==LG,]
  for(i in seq(window_size/2,as.integer(max(tmp$Consensus)),window_size/2)){
    tmp2<-tmp[tmp$Consensus >= i-window_size/2 & tmp$Consensus <= i + window_size/2,]  
    Pos<-c(Pos,i)
    #Markers<-c(Markers,dim(tmp2)[1]/dim(tmp)[1])
    #Scaffolds<-c(Scaffolds,length(unique(tmp2$Scaffold))/length(unique(tmp$Scaffold)))
    #Length<-c(Length,sum(scaf_length$length[scaf_length$Scaffold %in% tmp2$Scaffold])/sum(scaf_length$length[scaf_length$Scaffold %in% tmp$Scaffold]))
    #Genes<-c(Genes,length(unique(map_models$gene[map_models$V1 %in% tmp2$Scaffold]))/length(unique(map_models$gene[map_models$V1 %in% tmp$Scaffold])))
    Markers<-c(Markers,dim(tmp2)[1])
    Scaffolds<-c(Scaffolds,length(unique(tmp2$Scaffold)))
    Length<-c(Length,sum(scaf_length$length[scaf_length$Scaffold %in% tmp2$Scaffold]))
    Genes<-c(Genes,length(unique(map_models$gene[map_models$V1 %in% tmp2$Scaffold])))
    
     }
  Markers_mean<-rep(mean(Markers,na.rm=T),length(Markers))
  Markers_upperSD<-rep(mean(Markers,na.rm=T)+sd(Markers,na.rm=T),length(Markers))
  Markers_lowerSD<-rep(mean(Markers,na.rm=T)-sd(Markers,na.rm=T),length(Markers))
  Length_mean<-rep(mean(Length,na.rm=T),length(Length))
  Length_upperSD<-rep(mean(Length,na.rm=T)+sd(Length,na.rm=T),length(Length))
  Length_lowerSD<-rep(mean(Length,na.rm=T)-sd(Length,na.rm=T),length(Length))
  Scaffolds_mean<-rep(mean(Scaffolds,na.rm=T),length(Scaffolds))
  Scaffolds_upperSD<-rep(mean(Scaffolds,na.rm=T)+sd(Scaffolds,na.rm=T),length(Scaffolds))
  Scaffolds_lowerSD<-rep(mean(Scaffolds,na.rm=T)-sd(Scaffolds,na.rm=T),length(Scaffolds))
  Genes_mean<-rep(mean(Genes,na.rm=T),length(Genes))
  Genes_upperSD<-rep(mean(Genes,na.rm=T)+sd(Genes,na.rm=T),length(Genes))
  Genes_lowerSD<-rep(mean(Genes,na.rm=T)-sd(Genes,na.rm=T),length(Genes))
  
  return(as.data.frame(cbind(Pos,Markers,Markers_mean,Markers_lowerSD,Markers_upperSD,Scaffolds,Scaffolds_mean,Scaffolds_lowerSD,Scaffolds_upperSD,Length,Length_mean,Length_lowerSD,Length_upperSD,Genes,Genes_mean,Genes_lowerSD,Genes_upperSD)))
}  

## Create sliding window of 10 cM with 5 cM incremental steps through all 12 LGs
Recomb_LGI<-Recomb_spots(maps,1,10)
Recomb_LGII<-Recomb_spots(maps,2,10)
Recomb_LGIII<-Recomb_spots(maps,3,10)
Recomb_LGIV<-Recomb_spots(maps,4,10)
Recomb_LGV<-Recomb_spots(maps,5,10)
Recomb_LGVI<-Recomb_spots(maps,6,10)
Recomb_LGVII<-Recomb_spots(maps,7,10)
Recomb_LGVIII<-Recomb_spots(maps,8,10)
Recomb_LGIX<-Recomb_spots(maps,9,10)
Recomb_LGX<-Recomb_spots(maps,10,10)
Recomb_LGXI<-Recomb_spots(maps,11,10)
Recomb_LGXII<-Recomb_spots(maps,12,10)


#confidence_interval <- function(vector, interval) {
#  # Standard deviation of sample
#  vec_sd <- sd(vector)
#  # Sample size
#  n <- length(vector)
#  # Mean of sample
#  vec_mean <- mean(vector)
#  # Error according to t distribution
#  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
#  # Confidence interval as a vector
#  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
#  return(result)
#}

Recomb_maps<-rbind(Recomb_LGI,Recomb_LGII,Recomb_LGIII,Recomb_LGIV,Recomb_LGV,Recomb_LGVI,Recomb_LGVII,Recomb_LGVIII,Recomb_LGIX,Recomb_LGX,Recomb_LGXI,Recomb_LGXII)
Recomb_maps$LG<-c(rep(1,dim(Recomb_LGI)[1]),rep(2,dim(Recomb_LGII)[1]),rep(3,dim(Recomb_LGIII)[1]),rep(4,dim(Recomb_LGIV)[1]),rep(5,dim(Recomb_LGV)[1]),
              rep(6,dim(Recomb_LGVI)[1]),rep(7,dim(Recomb_LGVII)[1]),rep(8,dim(Recomb_LGVIII)[1]),rep(9,dim(Recomb_LGIX)[1]),rep(10,dim(Recomb_LGX)[1]),
              rep(11,dim(Recomb_LGXI)[1]),rep(12,dim(Recomb_LGXII)[1]))

##read in scripts for creating manhattan plots
pdf("sliding-10cM-window_5cM-steps_Markers.pdf")
par(mfrow=c(1,1),mar=c(4,4,2,0)+0.1)
#manhattan(x = Recomb_maps,chr="LG",p="Scaffolds",bp="Pos",mean="Scaffolds_mean",logp = F,ymax = max(Recomb_maps$Scaffolds)+10,point = F,ylab="",xlab="",main="Number of Scaffolds")
manhattan(x = Recomb_maps,chr="LG",p="Markers",bp="Pos",mean="Markers_mean",logp = F,ymax = max(Recomb_maps$Markers)+10,point = F,ylab="",xlab="Linkage group",main="Number of Markers per window")
#manhattan(x = Recomb_maps,chr="LG",p="Length",bp="Pos",mean="Length_mean",logp = F,ymax = max(Recomb_maps$Length)+10,point = F,ylab="",xlab="",main="Length (bp)")
#manhattan(x = Recomb_maps,chr="LG",p="Genes",bp="Pos",mean="Genes_mean",logp = F,ymax = max(Recomb_maps$Genes)+10,point = F,ylab="",xlab="Linkage group",main="Number of Genes per window")
dev.off()


