plotCompMap_new<-function(CompMap,LG,max_length){
  par(oma=c(2,2,6,2))
  layout(matrix(c(1,1,1,1,2,3,4,5),4,2,byrow=F))
  par(mar=c(2,2,0,2))
  par(las=1)
  plot(x=c(1,2),y=c(-25,max_length),xlim=c(0.8,3.2),type="n",yaxt="n",xaxt="n",xlab="", ylab="")
  
  lines(x=c(2,2),y=c(0,max(CompMap$consensus,na.rm=T)))
  lines(x=c(1,1),y=c(0,max(CompMap$X1,na.rm=T)/2))
  lines(x=c(1,1),y=c(max(CompMap$X1,na.rm=T)/2+50,max(CompMap$X1,na.rm=T)/2+50 + max(CompMap$X3,na.rm=T)/2))
  lines(x=c(3,3),y=c(0,max(CompMap$X2,na.rm=T)/2))
  text(c(1,1,2,3),c(-10,max(CompMap$X1,na.rm=T)/2+40,-10,-10), labels="0 cM" ,cex=0.8 )
  text(c(1,1,2,3),c(max(CompMap$X1,na.rm=T)/2+10,max(CompMap$X1,na.rm=T)/2 + max(CompMap$X3,na.rm=T)/2 +60,max(CompMap$consensus,na.rm=T)+10,max(CompMap$X2,na.rm=T)/2+10), labels=c(paste(max(CompMap$X1,na.rm=T)," cM"),paste(max(CompMap$X3,na.rm=T)," cM"),paste(max(CompMap$consensus,na.rm=T)," cM"),paste(max(CompMap$X2,na.rm=T)," cM")) ,cex=0.8 )
  text(c(0.8,0.8,3.2), c(max(CompMap$X1,na.rm=T)/4,max(CompMap$X1,na.rm=T)/2+50 + max(CompMap$X3,na.rm=T)/4,max(CompMap$X2,na.rm=T)/4),labels=c("Cluster 3", "Cluster 1", "Cluster 2"),cex=0.9,srt=90)
  
  for(i in 1:length(CompMap$marker)){
    lines(x=c(1.95,2.05),y=c(CompMap$consensus[i],CompMap$consensus[i]),col="black")
    lines(x=c(0.95,1.05),y=c(CompMap$X1[i]/2,CompMap$X1[i]/2),col="gray")
    lines(x=c(2.95,3.05),y=c(CompMap$X2[i]/2,CompMap$X2[i]/2),col="gray")
    lines(x=c(0.95,1.05),y=c(max(CompMap$X1,na.rm=T)/2+50 + CompMap$X3[i]/2,max(CompMap$X1,na.rm=T)/2+50 + CompMap$X3[i]/2),col="gray")
    if(!is.na(CompMap$X1[i])){
      lines(x=c(1.15,1.85),y=c(CompMap$X1[i]/2,CompMap$consensus[i]),lty=1,col="green")
    }
    if(!is.na(CompMap$X2[i])){
      lines(x=c(2.15,2.85),y=c(CompMap$consensus[i],CompMap$X2[i]/2),lty=1,col="red")
    }
    if(!is.na(CompMap$X3[i])){
      lines(x=c(1.15,1.85),y=c(max(CompMap$X1,na.rm=T)/2+50 + CompMap$X3[i]/2,CompMap$consensus[i]),lty=1,col="blue")
    }
  }

  Consensus_order<-seq(1:length(CompMap$marker))
  C1_order<-NULL
  tmp<-CompMap[order(CompMap$X3,na.last=NA),]
  for(i in 1:length(CompMap$marker)){
    if(CompMap$marker[i] %in% tmp$marker){
      C1_order[i]<-which(tmp$marker == CompMap$marker[i])
    } else {C1_order[i]<-NA}  
  }
  C2_order<-NULL
  tmp1<-CompMap[order(CompMap$X2,na.last=NA),]
  for(i in 1:length(CompMap$marker)){
    if(CompMap$marker[i] %in% tmp1$marker){
      C2_order[i]<-which(tmp1$marker == CompMap$marker[i])
    } else {C2_order[i]<-NA}  
  }
  C3_order<-NULL
  tmp2<-CompMap[order(CompMap$X1,na.last=NA),]
  for(i in 1:length(CompMap$marker)){
    if(CompMap$marker[i] %in% tmp2$marker){
      C3_order[i]<-which(tmp2$marker == CompMap$marker[i])
    } else {C3_order[i]<-NA}  
  }
  
  C1<-cor.test(as.numeric(C1_order),as.numeric(Consensus_order),method="kendall")$estimate[[1]]
  C2<-cor.test(as.numeric(C2_order),as.numeric(Consensus_order),method="kendall")$estimate[[1]]
  C3<-cor.test(as.numeric(C3_order),as.numeric(Consensus_order),method="kendall")$estimate[[1]]
  
  par(mar=c(0,4,0,2))
  plot(CompMap$consensus,CompMap$X1,na.rm=T,yaxt="n",xaxt="n",ylab="Cluster 3",xlab="")
  legend("topleft",legend=bquote(rho == .(round(C3,digits=3))),bty="n")
  plot(CompMap$consensus,CompMap$X2,na.rm=T,yaxt="n",xaxt="n",ylab="Cluster 2",xlab="")
  legend("topleft",legend=bquote(rho == .(round(C2,digits=3))),bty="n")
  plot(CompMap$consensus,CompMap$X3,na.rm=T,yaxt="n",xaxt="n",ylab="Cluster 1",xlab="")
  legend("topleft",legend=bquote(rho == .(round(C1,digits=3))),bty="n")
  plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(0,max(CompMap$consensus)))
  lines(c(0,max(CompMap$consensus)),c(1,1),type="l")
  for(i in 1:length(CompMap$consensus)){
    lines(x=c(CompMap$consensus[i],CompMap$consensus[i]),y=c(0.9,1.1))
  }
  title(main=paste("LG",LG),outer=T)
}
