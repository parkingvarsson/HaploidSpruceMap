load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGI_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGII_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGIII_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGIV_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGV_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGVI_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGVII_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGVIII_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGIX_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGX_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGXI_Consensus.RData")
load("~/Dropbox/work/Spruce_projects/megagametophytes/own_SNP_calls/final_calls/mother_vs_megs/3_clusters_from_pca/Onemap/probe_level_maps/BatchMap_2018/Cluster_comparison/Consensus/LGXII_Consensus.RData")



##Shortest map within the lowest range of mean RMSE

pdf("Comparison-Component_vs_Consensus.pdf")
plotCompMap_new(LGI_Consensus[[8]],"I",450)
plotCompMap_new(LGII_Consensus[[2]],"II",330)
plotCompMap_new(LGIII_Consensus[[10]],"III",400)
plotCompMap_new(LGIV_Consensus[[10]],"IV",360)
plotCompMap_new(LGV_Consensus[[4]],"V",380)
plotCompMap_new(LGVI_Consensus[[7]],"VI",320)
plotCompMap_new(LGVII_Consensus[[10]],"VII",380)
plotCompMap_new(LGVIII_Consensus[[8]],"VIII",380)
plotCompMap_new(LGIX_Consensus[[10]],"IX",350)
plotCompMap_new(LGX_Consensus[[6]],"X",320)
plotCompMap_new(LGXI_Consensus[[10]],"XI",280)
plotCompMap_new(LGXII_Consensus[[1]],"XII",370)
dev.off()

Cluster_comp<-function(Consensus_list){
 Consensus_order<-seq(1:length(Consensus_list$marker))
  C1_order<-NULL
  tmp<-Consensus_list[order(Consensus_list$C1,na.last=NA),]
  for(i in 1:length(Consensus_list$marker)){
    if(Consensus_list$marker[i] %in% tmp$marker){
      C1_order[i]<-which(tmp$marker == Consensus_list$marker[i])
    } else {C1_order[i]<-NA}  
  }
  C2_order<-NULL
  tmp1<-Consensus_list[order(Consensus_list$C2,na.last=NA),]
  for(i in 1:length(Consensus_list$marker)){
    if(Consensus_list$marker[i] %in% tmp1$marker){
      C2_order[i]<-which(tmp1$marker == Consensus_list$marker[i])
    } else {C2_order[i]<-NA}  
  }
  C3_order<-NULL
  tmp2<-Consensus_list[order(Consensus_list$C3,na.last=NA),]
  for(i in 1:length(Consensus_list$marker)){
    if(Consensus_list$marker[i] %in% tmp2$marker){
      C3_order[i]<-which(tmp2$marker == Consensus_list$marker[i])
    } else {C3_order[i]<-NA}  
  }
  
  combined<-as.data.frame(cbind(Consensus_order,C1_order,C2_order,C3_order))
  return(combined)
}

##order comparison for chosen interval in consensus maps
LGI_orderComp<-Cluster_comp(LGI_Consensus[[8]])
LGII_orderComp<-Cluster_comp(LGII_Consensus[[2]])
LGIII_orderComp<-Cluster_comp(LGIII_Consensus[[10]])
LGIV_orderComp<-Cluster_comp(LGIV_Consensus[[10]])
LGV_orderComp<-Cluster_comp(LGV_Consensus[[4]])
LGVI_orderComp<-Cluster_comp(LGVI_Consensus[[7]])
LGVII_orderComp<-Cluster_comp(LGVII_Consensus[[10]])
LGVIII_orderComp<-Cluster_comp(LGVIII_Consensus[[8]])
LGIX_orderComp<-Cluster_comp(LGIX_Consensus[[10]])
LGX_orderComp<-Cluster_comp(LGX_Consensus[[6]])
LGXI_orderComp<-Cluster_comp(LGXI_Consensus[[10]])
LGXII_orderComp<-Cluster_comp(LGXII_Consensus[[1]])


## Order correlation between consensus and frammework maps for chosen interval 
LGI_corr<-NULL
LGII_corr<-NULL
LGIII_corr<-NULL
LGIV_corr<-NULL
LGV_corr<-NULL
LGVI_corr<-NULL
LGVII_corr<-NULL
LGVIII_corr<-NULL
LGIX_corr<-NULL
LGX_corr<-NULL
LGXI_corr<-NULL
LGXII_corr<-NULL

for(i in 1:3){
  for(j in 2:4){
    if(i != j & i < j){
      LGI_corr<-c(LGI_corr,cor.test(LGI_orderComp[,i],LGI_orderComp[,j],method="kendall")$estimate[[1]])
      LGII_corr<-c(LGII_corr,cor.test(LGII_orderComp[,i],LGII_orderComp[,j],method="kendall")$estimate[[1]])
      LGIII_corr<-c(LGIII_corr,cor.test(LGIII_orderComp[,i],LGIII_orderComp[,j],method="kendall")$estimate[[1]])
      LGIV_corr<-c(LGIV_corr,cor.test(LGIV_orderComp[,i],LGIV_orderComp[,j],method="kendall")$estimate[[1]])
      LGV_corr<-c(LGV_corr,cor.test(LGV_orderComp[,i],LGV_orderComp[,j],method="kendall")$estimate[[1]])
      LGVI_corr<-c(LGVI_corr,cor.test(LGVI_orderComp[,i],LGVI_orderComp[,j],method="kendall")$estimate[[1]])
      LGVII_corr<-c(LGVII_corr,cor.test(LGVII_orderComp[,i],LGVII_orderComp[,j],method="kendall")$estimate[[1]])
      LGVIII_corr<-c(LGVIII_corr,cor.test(LGVIII_orderComp[,i],LGVIII_orderComp[,j],method="kendall")$estimate[[1]])
      LGIX_corr<-c(LGIX_corr,cor.test(LGIX_orderComp[,i],LGIX_orderComp[,j],method="kendall")$estimate[[1]])
      LGX_corr<-c(LGX_corr,cor.test(LGX_orderComp[,i],LGX_orderComp[,j],method="kendall")$estimate[[1]])
      LGXI_corr<-c(LGXI_corr,cor.test(LGXI_orderComp[,i],LGXI_orderComp[,j],method="kendall")$estimate[[1]])
      LGXII_corr<-c(LGXII_corr,cor.test(LGXII_orderComp[,i],LGXII_orderComp[,j],method="kendall")$estimate[[1]])
    }
  }  
}  

Order_correlations<-as.data.frame(rbind(LGI_corr,LGII_corr,LGIII_corr,LGIV_corr,LGV_corr,LGVI_corr,LGVII_corr,LGVIII_corr,LGIX_corr,LGX_corr,LGXI_corr,LGXII_corr))
names(Order_correlations)<-c("Consensus_C1","Consensus_C2","Consensus_C3","C1_C2","C1_C3","C2_C3")

plot_correlations<-function(LG_order,order_corr,LG,name){
  par(mfrow=c(4,4),mar=c(0,0,0,0)+0.1,oma=c(4,4,3,2))
  Distorted<-LG_order[which(as.character(Consensus_maps$Marker[Consensus_maps$LG==LG]) %in% Distorted_Scaffolds_in_map$Marker[Distorted_Scaffolds_in_map$Distorted_C3==T]),]
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center","Consensus",col="black",bty="n",cex=2)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center",legend=bquote(.(round(order_corr[1],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center",legend=bquote(.(round(order_corr[2],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center",legend=bquote(.(round(order_corr[3],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(LG_order$Consensus_order,LG_order$C1_order,pch=19,cex=0.5,ylab="",xlab="",xaxt="n",las=1)
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center","Cluster 1",col="black",bty="n",cex=2)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center",legend=bquote(.(round(order_corr[4],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="")
  legend("center",legend=bquote(.(round(order_corr[5],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(LG_order$Consensus_order,LG_order$C2_order,pch=19,cex=0.5,xlab="",ylab="",xaxt="n",las=1)
  plot(LG_order$C1_order,LG_order$C2_order,pch=19,cex=0.5,xlab="",ylab="",yaxt="n",xaxt="n",las=1)
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="")
  legend("center","Cluster 2",col="black",bty="n",cex=2)
  box()
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",las=1)
  legend("center",legend=bquote(.(round(order_corr[6],digit=3))),col="black",bty="n",cex=1.8)
  box()
  plot(LG_order$Consensus_order,LG_order$C3_order,pch=19,cex=0.5,xlab="",ylab="",las=1)
  points(Distorted$Consensus_order,Distorted$C3_order,pch=19,cex=0.5,col=5)
  plot(LG_order$C1_order,LG_order$C3_order,pch=19,cex=0.5,xlab="",ylab="",yaxt="n")
  points(Distorted$C1_order,Distorted$C3_order,pch=19,cex=0.5,col=5)
  plot(LG_order$C2_order,LG_order$C3_order,pch=19,cex=0.5,xlab="",ylab="",yaxt="n")
  points(Distorted$C2_order,Distorted$C3_order,pch=19,cex=0.5,col=5)
  plot(1,type="n",pch=19,axes=F,ylab="",xlab="",yaxt="n")
  legend("center","Cluster 3",col="black",bty="n",cex=2)
  mtext(side=3,text=paste("LG ",name,sep=""),outer=T)
  box()
}


pdf("Framework-framework-consensus_comparision.pdf")
plot_correlations(LGI_orderComp,Order_correlations[1,],1,"I")
plot_correlations(LGII_orderComp,Order_correlations[2,],2,"II")
plot_correlations(LGIII_orderComp,Order_correlations[3,],3,"III")
plot_correlations(LGIV_orderComp,Order_correlations[4,],4,"IV")
plot_correlations(LGV_orderComp,Order_correlations[5,],5,"V")
plot_correlations(LGVI_orderComp,Order_correlations[6,],6,"VI")
plot_correlations(LGVII_orderComp,Order_correlations[7,],7,"VII")
plot_correlations(LGVIII_orderComp,Order_correlations[8,],8,"VIII")
plot_correlations(LGIX_orderComp,Order_correlations[9,],9,"IX")
plot_correlations(LGX_orderComp,Order_correlations[10,],10,"X")
plot_correlations(LGXI_orderComp,Order_correlations[11,],11,"XI")
plot_correlations(LGXII_orderComp,Order_correlations[12,],12,"XII")
dev.off()


###Extract consensus maps for downstream analysis



Consensus_maps<-rbind(LGI_Consensus[[8]],LGII_Consensus[[2]],
      LGIII_Consensus[[10]],
      LGIV_Consensus[[10]],
      LGV_Consensus[[4]],
      LGVI_Consensus[[7]],
      LGVII_Consensus[[10]],
      LGVIII_Consensus[[8]],
      LGIX_Consensus[[10]],
      LGX_Consensus[[6]],
      LGXI_Consensus[[10]],
      LGXII_Consensus[[1]])

Consensus_maps<-cbind(c(rep(1,dim(LGI_Consensus[[8]])[1]),
                        rep(2,dim(LGII_Consensus[[2]])[1]),
                        rep(3,dim(LGIII_Consensus[[10]])[1]),
                        rep(4,dim(LGIV_Consensus[[10]])[1]),
                        rep(5,dim(LGV_Consensus[[4]])[1]),
                        rep(6,dim(LGVI_Consensus[[7]])[1]),
                        rep(7,dim(LGVII_Consensus[[10]])[1]),
                        rep(8,dim(LGVIII_Consensus[[8]])[1]),
                        rep(9,dim(LGIX_Consensus[[10]])[1]),
                        rep(10,dim(LGX_Consensus[[6]])[1]),
                        rep(11,dim(LGXI_Consensus[[10]])[1]),
                        rep(12,dim(LGXII_Consensus[[1]])[1])),Consensus_maps)
names(Consensus_maps)<-c("LG","Marker","Consensus","C3","C2","C1")
Consensus_maps$Scaffold<-as.factor(t(as.data.frame(strsplit(as.character(Consensus_maps$Marker),split=":")))[,1])
Consensus_maps$Probe<-as.factor(t(as.data.frame(strsplit(as.character(Consensus_maps$Marker),split=":")))[,2])

write.table(Consensus_maps,"Consensus_maps.txt",quote=F,row.names = F,sep="\t")
save.image("LPmerge_outputs.RData")


###Where are distorted C3 markers distributed in the consensus map?
read.table("Distorted_C3_markers.txt",head=F)->Distorted_C3
names(Distorted_C3)<-"Marker"
Distorted_C3$Scaffold<-as.factor(t(as.data.frame(strsplit(as.character(Distorted_C3$Marker),split=":")))[,1])
Distorted_C3$Probe<-as.factor(t(as.data.frame(strsplit(as.character(Distorted_C3$Marker),split=":")))[,2])


Distorted_Scaffolds_in_map<-Consensus_maps[Consensus_maps$Scaffold %in% as.character(Distorted_C3$Scaffold),]
Distorted_Scaffolds_in_map$Distorted_C3<-"FALSE"
for(i in 1:dim(Distorted_Scaffolds_in_map)[1]){
  if(Distorted_Scaffolds_in_map$Marker[i] %in% as.character(Distorted_C3$Marker)){
    Distorted_Scaffolds_in_map$Distorted_C3[i]<-TRUE
  }
}

library(gdata)
Distorted_Scaffolds_in_map <- drop.levels(Distorted_Scaffolds_in_map)

for(i in 1:181){
  if (length(which(table(Distorted_Scaffolds_in_map$Scaffold,Distorted_Scaffolds_in_map$LG)[i,] > 0)) >1){
    print(as.data.frame(table(Distorted_Scaffolds_in_map$Scaffold,Distorted_Scaffolds_in_map$LG))[i,])
  }
}  




