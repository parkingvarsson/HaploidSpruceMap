C3_LG10_dropped$flipped_Size<-sapply(C3_LG10_dropped$Size,function(x) {max(C3_LG10_dropped$Size)-C3_LG10_dropped$Size})[,1] ##LG I
C2_LG12_dropped$flipped_Size<-sapply(C2_LG12_dropped$Size,function(x) {max(C2_LG12_dropped$Size)-C2_LG12_dropped$Size})[,1] ##LG II
C2_LG1_dropped$flipped_Size<-sapply(C2_LG1_dropped$Size,function(x) {max(C2_LG1_dropped$Size)-C2_LG1_dropped$Size})[,1] ##LG III
C3_LG3_dropped$flipped_Size<-sapply(C3_LG3_dropped$Size,function(x) {max(C3_LG3_dropped$Size)-C3_LG3_dropped$Size})[,1] ##LG V
C1_LG6_dropped$flipped_Size<-sapply(C1_LG6_dropped$Size,function(x) {max(C1_LG6_dropped$Size)-C1_LG6_dropped$Size})[,1] ##LG VI
C3_LG8_dropped$flipped_Size<-sapply(C3_LG8_dropped$Size,function(x) {max(C3_LG8_dropped$Size)-C3_LG8_dropped$Size})[,1] ##LG VII
C2_LG8_dropped$flipped_Size<-sapply(C2_LG8_dropped$Size,function(x) {max(C2_LG8_dropped$Size)-C2_LG8_dropped$Size})[,1] ##LG VIII
C1_LG9_dropped$flipped_Size<-sapply(C1_LG9_dropped$Size,function(x) {max(C1_LG9_dropped$Size)-C1_LG9_dropped$Size})[,1] ## LG IX
C2_LG10_dropped$flipped_Size<-sapply(C2_LG10_dropped$Size,function(x) {max(C2_LG10_dropped$Size)-C2_LG10_dropped$Size})[,1] ##LG X


##compare orders between clusters

Cluster_comp<-function(C1_LG,C2_LG,C3_LG){
  C1_vsC2<-which(C1_LG$Marker %in% C2_LG$Marker)
  C2_vsC1<-NULL
  for(i in 1:length(C1_vsC2)){
    C2_vsC1[i]<-which(as.character(C2_LG$Marker) == as.character(C1_LG$Marker[C1_vsC2[i]]))
  }
  C1_vsC3<-which(C1_LG$Marker %in% C3_LG$Marker)
  C3_vsC1<-NULL
  for(i in 1:length(C1_vsC3)){
    C3_vsC1[i]<-which(as.character(C3_LG$Marker) == as.character(C1_LG$Marker[C1_vsC3[i]]))
  }
  C2_vsC3<-which(C2_LG$Marker %in% C3_LG$Marker)
  C3_vsC2<-NULL
  for(i in 1:length(C2_vsC3)){
    C3_vsC2[i]<-which(as.character(C3_LG$Marker) == as.character(C2_LG$Marker[C2_vsC3[i]]))
  }
  combined<-list(C1_vsC2,C2_vsC1,C1_vsC3,C3_vsC1,C2_vsC3,C3_vsC2)
  return(combined)
}


par(mfrow=c(4,3))
plot(1,type="n",xlim=c(0,510),ylim=c(0,510),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG I")
abline(v=seq(0,510,10),h=seq(0,510,10),col="grey95")
points(C1_LG1_dropped$Size[C1LG1_C2LG11_C3LG10_dropped[1][[1]]],C2_LG11_dropped$Size[C1LG1_C2LG11_C3LG10_dropped[2][[1]]],pch=19,col="black")
points(C1_LG1_dropped$Size[C1LG1_C2LG11_C3LG10_dropped[3][[1]]],C3_LG10_dropped$flipped_Size[C1LG1_C2LG11_C3LG10_dropped[4][[1]]],pch=19,col="red")
points(C2_LG11_dropped$Size[C1LG1_C2LG11_C3LG10_dropped[5][[1]]],C3_LG10_dropped$flipped_Size[C1LG1_C2LG11_C3LG10_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,350),ylim=c(0,350),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG II")
abline(v=seq(0,350,10),h=seq(0,350,10),col="grey95")
points(C1_LG2_dropped$Size[C1LG2_C2LG12_C3LG1_dropped[1][[1]]],C2_LG12_dropped$flipped_Size[C1LG2_C2LG12_C3LG1_dropped[2][[1]]],pch=19,col="black")
points(C1_LG2_dropped$Size[C1LG2_C2LG12_C3LG1_dropped[3][[1]]],C3_LG1_dropped$Size[C1LG2_C2LG12_C3LG1_dropped[4][[1]]],pch=19,col="red")
points(C2_LG12_dropped$flipped_Size[C1LG2_C2LG12_C3LG1_dropped[5][[1]]],C3_LG1_dropped$Size[C1LG2_C2LG12_C3LG1_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,430),ylim=c(0,430),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG III")
abline(v=seq(0,430,10),h=seq(0,430,10),col="grey95")
points(C1_LG3_dropped$Size[C1LG3_C2LG1_C3LG2_dropped[1][[1]]],C2_LG1_dropped$flipped_Size[C1LG3_C2LG1_C3LG2_dropped[2][[1]]],pch=19,col="black")
points(C1_LG3_dropped$Size[C1LG3_C2LG1_C3LG2_dropped[3][[1]]],C3_LG2_dropped$Size[C1LG3_C2LG1_C3LG2_dropped[4][[1]]],pch=19,col="red")
points(C2_LG1_dropped$flipped_Size[C1LG3_C2LG1_C3LG2_dropped[5][[1]]],C3_LG2_dropped$Size[C1LG3_C2LG1_C3LG2_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,420),ylim=c(0,420),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG IV")
abline(v=seq(0,420,10),h=seq(0,420,10),col="grey95")
points(C1_LG4_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[1][[1]]],C2_LG3_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[2][[1]]],pch=19,col="black")
points(C1_LG4_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[3][[1]]],C3_LG5_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[4][[1]]],pch=19,col="red")
points(C2_LG3_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[5][[1]]],C3_LG5_dropped$Size[C1LG4_C2LG3_C3LG5_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,480),ylim=c(0,480),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG V")
abline(v=seq(0,480,10),h=seq(0,480,10),col="grey95")
points(C1_LG5_dropped$Size[C1LG5_C2LG9_C3LG3_dropped[1][[1]]],C2_LG9_dropped$Size[C1LG5_C2LG9_C3LG3_dropped[2][[1]]],pch=19,col="black")
points(C1_LG5_dropped$Size[C1LG5_C2LG9_C3LG3_dropped[3][[1]]],C3_LG3_dropped$flipped_Size[C1LG5_C2LG9_C3LG3_dropped[4][[1]]],pch=19,col="red")
points(C2_LG9_dropped$Size[C1LG5_C2LG9_C3LG3_dropped[5][[1]]],C3_LG3_dropped$flipped_Size[C1LG5_C2LG9_C3LG3_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,330),ylim=c(0,330),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG VI")
abline(v=seq(0,330,10),h=seq(0,330,10),col="grey95")
points(C1_LG6_dropped$flipped_Size[C1LG6_C2LG5_C3LG7_dropped[1][[1]]],C2_LG5_dropped$Size[C1LG6_C2LG5_C3LG7_dropped[2][[1]]],pch=19,col="black")
points(C1_LG6_dropped$flipped_Size[C1LG6_C2LG5_C3LG7_dropped[3][[1]]],C3_LG7_dropped$Size[C1LG6_C2LG5_C3LG7_dropped[4][[1]]],pch=19,col="red")
points(C2_LG5_dropped$Size[C1LG6_C2LG5_C3LG7_dropped[5][[1]]],C3_LG7_dropped$Size[C1LG6_C2LG5_C3LG7_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,450),ylim=c(0,450),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG VII")
abline(v=seq(0,450,10),h=seq(0,450,10),col="grey95")
points(C1_LG7_dropped$Size[C1LG7_C2LG6_C3LG8_dropped[1][[1]]],C2_LG6_dropped$Size[C1LG7_C2LG6_C3LG8_dropped[2][[1]]],pch=19,col="black")
points(C1_LG7_dropped$Size[C1LG7_C2LG6_C3LG8_dropped[3][[1]]],C3_LG8_dropped$flipped_Size[C1LG7_C2LG6_C3LG8_dropped[4][[1]]],pch=19,col="red")
points(C2_LG6_dropped$Size[C1LG7_C2LG6_C3LG8_dropped[5][[1]]],C3_LG8_dropped$flipped_Size[C1LG7_C2LG6_C3LG8_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,420),ylim=c(0,420),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG VIII")
abline(v=seq(0,420,10),h=seq(0,420,10),col="grey95")
points(C1_LG8_dropped$Size[C1LG8_C2LG8_C3LG11_dropped[1][[1]]],C2_LG8_dropped$flipped_Size[C1LG8_C2LG8_C3LG11_dropped[2][[1]]],pch=19,col="black")
points(C1_LG8_dropped$Size[C1LG8_C2LG8_C3LG11_dropped[3][[1]]],C3_LG11_dropped$Size[C1LG8_C2LG8_C3LG11_dropped[4][[1]]],pch=19,col="red")
points(C2_LG8_dropped$flipped_Size[C1LG8_C2LG8_C3LG11_dropped[5][[1]]],C3_LG11_dropped$Size[C1LG8_C2LG8_C3LG11_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,370),ylim=c(0,370),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG IX")
abline(v=seq(0,370,10),h=seq(0,370,10),col="grey95")
points(C1_LG9_dropped$flipped_Size[C1LG9_C2LG4_C3LG6_dropped[1][[1]]],C2_LG4_dropped$Size[C1LG9_C2LG4_C3LG6_dropped[2][[1]]],pch=19,col="black")
points(C1_LG9_dropped$flipped_Size[C1LG9_C2LG4_C3LG6_dropped[3][[1]]],C3_LG6_dropped$Size[C1LG9_C2LG4_C3LG6_dropped[4][[1]]],pch=19,col="red")
points(C2_LG4_dropped$Size[C1LG9_C2LG4_C3LG6_dropped[5][[1]]],C3_LG6_dropped$Size[C1LG9_C2LG4_C3LG6_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,410),ylim=c(0,410),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG X")
abline(v=seq(0,410,10),h=seq(0,410,10),col="grey95")
points(C1_LG10_dropped$Size[C1LG10_C2LG10_C3LG12_dropped[1][[1]]],C2_LG10_dropped$flipped_Size[C1LG10_C2LG10_C3LG12_dropped[2][[1]]],pch=19,col="black")
points(C1_LG10_dropped$Size[C1LG10_C2LG10_C3LG12_dropped[3][[1]]],C3_LG12_dropped$Size[C1LG10_C2LG10_C3LG12_dropped[4][[1]]],pch=19,col="red")
points(C2_LG10_dropped$flipped_Size[C1LG10_C2LG10_C3LG12_dropped[5][[1]]],C3_LG12_dropped$Size[C1LG10_C2LG10_C3LG12_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,350),ylim=c(0,350),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG XI")
abline(v=seq(0,350,10),h=seq(0,350,10),col="grey95")
points(C1_LG11_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[1][[1]]],C2_LG2_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[2][[1]]],pch=19,col="black")
points(C1_LG11_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[3][[1]]],C3_LG4_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[4][[1]]],pch=19,col="red")
points(C2_LG2_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[5][[1]]],C3_LG4_dropped$Size[C1LG11_C2LG2_C3LG4_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)
plot(1,type="n",xlim=c(0,430),ylim=c(0,430),xlab="First Cluster (cM)",ylab="Second Cluster (cM)",main="LG XII")
abline(v=seq(0,430,10),h=seq(0,430,10),col="grey95")
points(C1_LG12_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[1][[1]]],C2_LG7_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[2][[1]]],pch=19,col="black")
points(C1_LG12_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[3][[1]]],C3_LG9_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[4][[1]]],pch=19,col="red")
points(C2_LG7_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[5][[1]]],C3_LG9_dropped$Size[C1LG12_C2LG7_C3LG9_dropped[6][[1]]],pch=19,col="blue")
legend("topleft", c("C1 vs C2","C1 vs C3","C2 vs C3"), text.col = c("black","red","blue"), bty="n", cex=0.9)



### CREATE CONSENSUS MAPS
LGI_ConsensusAll<-LPmerge(list(C1_LG1_dropped[,2:3],C2_LG11_dropped[,2:3],C3_LG10_dropped[,c(2,4)]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGI_ConsensusAll,file="LGI_ConsensusAll.RData")
LGII_ConsensusAll<-LPmerge(list(C1_LG2_dropped[,2:3],C2_LG12_dropped[,c(2,4)],C3_LG1_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGII_ConsensusAll,file="LGII_ConsensusAll.RData")
LGIII_ConsensusAll<-LPmerge(list(C1_LG3_dropped[,2:3],C2_LG1_dropped[,c(2,4)],C3_LG2_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGIII_ConsensusAll,file="LGIII_ConsensusAll.RData")
LGIV_ConsensusAll<-LPmerge(list(C1_LG4_dropped[,2:3],C2_LG3_dropped[,2:3],C3_LG5_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGIV_ConsensusAll,file="LGIV_ConsensusAll.RData")
LGV_ConsensusAll<-LPmerge(list(C1_LG5_dropped[,2:3],C2_LG9_dropped[,2:3],C3_LG3_dropped[,c(2,4)]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGV_ConsensusAll,file="LGV_ConsensusAll.RData")
LGV_ConsensusC1C3<-LPmerge(list(C1_LG5_dropped[,2:3],C3_LG3_dropped[,c(2,4)]),max.interval=1:5,weights=c(0.27,0.73)) ##interval 2
save(LGV_ConsensusC1C3,file="LGV_ConsensusC1C3.RData")
LGVI_ConsensusAll<-LPmerge(list(C1_LG6_dropped[,c(2,4)],C2_LG5_dropped[,2:3],C3_LG7_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGVI_ConsensusAll,file="LGVI_ConsensusAll.RData")
LGVII_ConsensusAll<-LPmerge(list(C1_LG7_dropped[,2:3],C2_LG6_dropped[,2:3],C3_LG8_dropped[,c(2,4)]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGVII_ConsensusAll,file="LGVII_ConsensusAll.RData")
LGVIII_ConsensusAll<-LPmerge(list(C1_LG8_dropped[,2:3],C2_LG8_dropped[,c(2,4)],C3_LG11_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGVIII_ConsensusAll,file="LGVIII_ConsensusAll.RData")
LGIX_ConsensusAll<-LPmerge(list(C1_LG9_dropped[,c(2,4)],C2_LG4_dropped[,2:3],C3_LG6_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGIX_ConsensusAll,file="LGIX_ConsensusAll.RData")
LGX_ConsensusAll<-LPmerge(list(C1_LG10_dropped[,2:3],C2_LG10_dropped[,c(2,4)],C3_LG12_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGX_ConsensusAll,file="LGX_ConsensusAll.RData")
LGXI_ConsensusAll<-LPmerge(list(C1_LG11_dropped[,2:3],C2_LG2_dropped[,2:3],C3_LG4_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGXI_ConsensusAll,file="LGXI_ConsensusAll.RData")
LGXII_ConsensusAll<-LPmerge(list(C1_LG12_dropped[,2:3],C2_LG7_dropped[,2:3],C3_LG9_dropped[,2:3]),max.interval=1:5,weights=c(0.22,0.19,0.59))
save(LGXII_ConsensusAll,file="LGXII_ConsensusAll.RData")



all_splitted_scaffolds$Consensus_LG<-NULL
for (i in 1:dim(all_splitted_scaffolds)[1]){
  if(all_splitted_scaffolds$Cluster[i]=="C1"){
    all_splitted_scaffolds$Consensus_LG[i]<-as.character(Cluster_comp$LG[Cluster_comp$C1 == all_splitted_scaffolds$LG[i]])
  } 
  if(all_splitted_scaffolds$Cluster[i]=="C2"){
    all_splitted_scaffolds$Consensus_LG[i]<-as.character(Cluster_comp$LG[Cluster_comp$C2 == all_splitted_scaffolds$LG[i]])
  }
  if(all_splitted_scaffolds$Cluster[i] == "C3"){
    all_splitted_scaffolds$Consensus_LG[i]<-as.character(Cluster_comp$LG[Cluster_comp$C3 == all_splitted_scaffolds$LG[i]])  
  }
}

output<-data.frame(matrix(ncol = 13, nrow = 157))
names(output)<-c("Scaffold","LG I","LG II","LG III","LG IV","LG V","LG VI","LG VII","LG VIII","LG IX","LG X","LG XI","LG XII")
for(i in 1:length(levels(all_splitted_scaffolds$scaffold))){
  output[i,1]<-levels(all_splitted_scaffolds$scaffold)[i]
  for(j in 1:length(levels(all_splitted_scaffolds$Consensus_LG))){
    if (length(unique(all_splitted_scaffolds$probe[all_splitted_scaffolds$scaffold == levels(all_splitted_scaffolds$scaffold)[i] & all_splitted_scaffolds$Consensus_LG==levels(all_splitted_scaffolds$Consensus_LG)[j]])) >0){
      tmp<-as.character(unique(all_splitted_scaffolds$probe[all_splitted_scaffolds$scaffold == levels(all_splitted_scaffolds$scaffold)[i] & all_splitted_scaffolds$Consensus_LG==levels(all_splitted_scaffolds$Consensus_LG)[j]]))
      output[i,j+1]<-paste(tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],sep=",")    
    }
    else{ 
      output[i,j+1]<-"-"
     }  
  }
}