##Lowest mean RMSE

pdf("Comparison-Component_vs_Consensus.pdf")
plotCompMap_new(LGI_Consensus[[2]],"I",450)
plotCompMap_new(LGII_Consensus[[3]],"II",350)
plotCompMap_new(LGIII_Consensus[[2]],"III",400)
plotCompMap_new(LGIV_Consensus[[3]],"IV",350)
plotCompMap_new(LGV_Consensus[[9]],"V",400)
plotCompMap_new(LGVI_Consensus[[2]],"VI",300)
plotCompMap_new(LGVII_Consensus[[10]],"VII",400)
plotCompMap_new(LGVIII_Consensus[[6]],"VIII",350)
plotCompMap_new(LGIX_Consensus[[3]],"IX",350)
plotCompMap_new(LGX_Consensus[[9]],"X",350)
plotCompMap_new(LGXI_Consensus[[10]],"XI",300)
plotCompMap_new(LGXII_Consensus[[2]],"XII",400)
dev.off()

Cluster_comp<-function(Consensus_list,interval){
 Consensus_order<-seq(1:length(Consensus_list[[interval]]$marker))
  C1_order<-NULL
  tmp<-Consensus_list[[interval]][order(Consensus_list[[interval]]$X3,na.last=NA),]
  for(i in 1:length(Consensus_list[[interval]]$marker)){
    if(Consensus_list[[interval]]$marker[i] %in% tmp$marker){
      C1_order[i]<-which(tmp$marker == Consensus_list[[interval]]$marker[i])
    } else {C1_order[i]<-NA}  
  }
  C2_order<-NULL
  tmp1<-Consensus_list[[interval]][order(Consensus_list[[interval]]$X2,na.last=NA),]
  for(i in 1:length(Consensus_list[[interval]]$marker)){
    if(Consensus_list[[interval]]$marker[i] %in% tmp1$marker){
      C2_order[i]<-which(tmp1$marker == Consensus_list[[interval]]$marker[i])
    } else {C2_order[i]<-NA}  
  }
  C3_order<-NULL
  tmp2<-Consensus_list[[interval]][order(Consensus_list[[interval]]$X1,na.last=NA),]
  for(i in 1:length(Consensus_list[[interval]]$marker)){
    if(Consensus_list[[interval]]$marker[i] %in% tmp2$marker){
      C3_order[i]<-which(tmp2$marker == Consensus_list[[interval]]$marker[i])
    } else {C3_order[i]<-NA}  
  }
  
  combined<-as.data.frame(cbind(Consensus_order,C1_order,C2_order,C3_order))
  return(combined)
}

##order com of all interval consensus maps
LGI_orderCompAllIntervals<-NULL
LGII_orderCompAllIntervals<-NULL
LGIII_orderCompAllIntervals<-NULL
LGIV_orderCompAllIntervals<-NULL
LGV_orderCompAllIntervals<-NULL
LGVI_orderCompAllIntervals<-NULL
LGVII_orderCompAllIntervals<-NULL
LGVIII_orderCompAllIntervals<-NULL
LGIX_orderCompAllIntervals<-NULL
LGX_orderCompAllIntervals<-NULL
LGXI_orderCompAllIntervals<-NULL
LGXII_orderCompAllIntervals<-NULL
for (i in 1:10){
  LGI_orderCompAllIntervals[[i]]<-Cluster_comp(LGI_Consensus,i)
  LGII_orderCompAllIntervals[[i]]<-Cluster_comp(LGII_Consensus,i)
  LGIII_orderCompAllIntervals[[i]]<-Cluster_comp(LGIII_Consensus,i)
  LGIV_orderCompAllIntervals[[i]]<-Cluster_comp(LGIV_Consensus,i)
  LGV_orderCompAllIntervals[[i]]<-Cluster_comp(LGV_Consensus,i)
  LGVI_orderCompAllIntervals[[i]]<-Cluster_comp(LGVI_Consensus,i)
  LGVII_orderCompAllIntervals[[i]]<-Cluster_comp(LGVII_Consensus,i)
  LGVIII_orderCompAllIntervals[[i]]<-Cluster_comp(LGVIII_Consensus,i)
  LGIX_orderCompAllIntervals[[i]]<-Cluster_comp(LGIX_Consensus,i)
  LGX_orderCompAllIntervals[[i]]<-Cluster_comp(LGX_Consensus,i)
  LGXI_orderCompAllIntervals[[i]]<-Cluster_comp(LGXI_Consensus,i)
  LGXII_orderCompAllIntervals[[i]]<-Cluster_comp(LGXII_Consensus,i)
}

RMSE_corr<-as.data.frame(matrix(NA,nrow=10,ncol=3))
for (i in 1:dim(RMSE_corr)[1]){
  RMSE_corr[i,1]<-cor.test(as.numeric(LGs_orderCompRMSE[[i]][,1]),as.numeric(LGs_orderCompRMSE[[i]][,2]),method="pearson")$estimate[[1]]
  RMSE_corr[i,2]<-cor.test(as.numeric(LGs_orderCompRMSE[[i]][,1]),as.numeric(LGs_orderCompRMSE[[i]][,3]),method="pearson")$estimate[[1]]
  RMSE_corr[i,3]<-cor.test(as.numeric(LGs_orderCompRMSE[[i]][,1]),as.numeric(LGs_orderCompRMSE[[i]][,4]),method="pearson")$estimate[[1]]
}

library(PerformanceAnalytics)
chart.Correlation(LGI_orderCompAllIntervals[[2]],digit=3,method="kendall")
chart.Correlation(LGII_orderCompAllIntervals[[3]],digit=3,method="kendall")
chart.Correlation(LGIII_orderCompAllIntervals[[2]],digit=3,method="kendall")
chart.Correlation(LGIV_orderCompAllIntervals[[3]],digit=3,method="kendall")
chart.Correlation(LGV_orderCompAllIntervals[[9]],digit=3,method="kendall")
chart.Correlation(LGVI_orderCompAllIntervals[[2]],digit=3,method="kendall")
chart.Correlation(LGVII_orderCompAllIntervals[[10]],digit=3,method="kendall")
chart.Correlation(LGVIII_orderCompAllIntervals[[6]],digit=3,method="kendall")
chart.Correlation(LGIX_orderCompAllIntervals[[3]],digit=3,method="kendall")
chart.Correlation(LGX_orderCompAllIntervals[[9]],digit=3,method="kendall")
chart.Correlation(LGXI_orderCompAllIntervals[[10]],digit=3,method="kendall")
chart.Correlation(LGXII_orderCompAllIntervals[[2]],digit=3,method="kendall")


plot_Order.vs.cM<-function(LG,intervall,max.size,main){
  plot(LG[[intervall]]$consensus,pch=21,col="black",type="l",ylim=c(0,max.size),ylab="Position in cM",xlab="Consensus order",main=main) 
  points(LG[[intervall]]$X1,pch=21,col="red",type="l") 
  points(LG[[intervall]]$X2,pch=21,col="green",type="l")
  points(LG[[intervall]]$X3,pch=21,col="blue",type="l")
  legend(0,max.size, c("Consensus","Cluster 1", "Cluster 2","Cluster 3"),pch=21,col=c("black","blue","green","red"),bty="n",cex=0.7)
}

par(mfrow=c(1,1))
#plot_Order.vs.cM(LGI_Consensus,?,400,,"LGI")
plot_Order.vs.cM(LGII_Consensus,3,300,"LGII")
plot_Order.vs.cM(LGIII_Consensus,2,400,"LGIII")
plot_Order.vs.cM(LGIV_Consensus,3,350,"LGIV")
plot_Order.vs.cM(LGV_Consensus,9,400,"LGV")
plot_Order.vs.cM(LGVI_Consensus,2,300,"LGVI")
plot_Order.vs.cM(LGVII_Consensus,10,400,"LGVII")
plot_Order.vs.cM(LGVIII_Consensus,6,350,"LGVIII")
plot_Order.vs.cM(LGIX_Consensus,3,350,"LGIX")
plot_Order.vs.cM(LGX_Consensus,9,350,"LGX")
plot_Order.vs.cM(LGXI_Consensus,10,300,"LGXI")
plot_Order.vs.cM(LGXII_Consensus,2,400,"LGXII")

par(mfrow=c(1,1))
plotConsensusMap(list(LGI_Consensus[[2]][,2],LGII_Consensus[[3]][,2],LGIII_Consensus[[2]][,2],LGVI_Consensus[[3]][,2],LGV_Consensus[[9]][,2],LGVI_Consensus[[2]][,2],LGVII_Consensus[[10]][,2],
     LGVIII_Consensus[[6]][,2],LGIX_Consensus[[3]][,2],LGX_Consensus[[9]][,2],LGXI_Consensus[[10]][,2],LGXII_Consensus[[2]][,2]),360)
