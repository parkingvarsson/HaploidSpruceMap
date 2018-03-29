###Recombination hotspots/coldspots

##Read in scaffold length
read.table("haploid_scaffolds.txt",head=F)->scaf_length
names(scaf_length)<-c("Scaffold","length")

### read in gene models
read.table("Eugene.gff3")->all_models
read.table("HC_genes.gff3")->HC_models
read.table("MC_genes.gff3")->MC_models
read.table("LC_genes.gff3")->LC_models

HC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(HC_models$V9),";")))[,1])),"=")))[,2])
MC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(MC_models$V9),";")))[,1])),"=")))[,2])
LC_models$gene<-as.factor(t(as.data.frame(strsplit(as.character(as.factor(t(as.data.frame(strsplit(as.character(LC_models$V9),";")))[,1])),"=")))[,2])

library(gdata)
HC_models<-drop.levels(HC_models[HC_models$V1 %in% maps$scaf,])
MC_models<-drop.levels(MC_models[MC_models$V1 %in% maps$scaf,])
LC_models<-drop.levels(LC_models[LC_models$V1 %in% maps$scaf,])
map_models<-rbind(HC_models,MC_models,LC_models)


##Read in consensus map
read.table("Consensus-maps.txt",head=T)->maps
maps$scaf<-as.factor(t(as.data.frame(strsplit(as.character(maps$marker),":")))[,1])


## Function that counts the number of markers, scaffolds and total lenght of all scaffolds in a predefined window size along the LGs
Recomb_spots<-function(map,LG,window_size){
  Recomb_markers<-NULL
  Recomb_scaffolds<-NULL
  Recomb_length<-NULL
  Recomb_genes<-NULL
  tmp<-map[map$LG==LG,]
  for(i in window_size:as.integer(max(tmp$consensus))){
    tmp2<-tmp[tmp$consensus >i-window_size & tmp$consensus <i,]  
    Recomb_markers<-c(Recomb_markers,dim(tmp2)[1])
    Recomb_scaffolds<-c(Recomb_scaffolds,length(unique(tmp2$scaf)))
    Recomb_length<-c(Recomb_length,sum(scaf_length$length[scaf_length$Scaffold %in% tmp2$scaf]))
    Recomb_genes<-c(Recomb_genes,length(unique(map_models$gene[map_models$V1 %in% tmp2$scaf])))
  }
  return(as.data.frame(cbind(Recomb_markers,Recomb_scaffolds,Recomb_length,Recomb_genes)))
}  

## Create sliding window of 5 cM through all 12 LGs
Recomb_LGI<-Recomb_spots(maps,1,5)
Recomb_LGII<-Recomb_spots(maps,2,5)
Recomb_LGIII<-Recomb_spots(maps,3,5)
Recomb_LGIV<-Recomb_spots(maps,4,5)
Recomb_LGV<-Recomb_spots(maps,5,5)
Recomb_LGVI<-Recomb_spots(maps,6,5)
Recomb_LGVII<-Recomb_spots(maps,7,5)
Recomb_LGVIII<-Recomb_spots(maps,8,5)
Recomb_LGIX<-Recomb_spots(maps,9,5)
Recomb_LGX<-Recomb_spots(maps,10,5)
Recomb_LGXI<-Recomb_spots(maps,11,5)
Recomb_LGXII<-Recomb_spots(maps,12,5)


par(mfrow=c(4,3))
par(mar=c(5,5,2,3))
plot(Recomb_LGI$Recomb_genes,type="l",ylim=c(0,290),ylab="Number of genes",xlab="",main="LGI, Total physical size 35.3 Mb",las=1)
abline(h=mean(Recomb_LGI$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGI$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGI$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,4,2,4))
plot(Recomb_LGII$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGII, Total physical size 26.0 Mb",las=1)
abline(h=mean(Recomb_LGII$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGII$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGII$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,3,2,5))
plot(Recomb_LGIII$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGIII, Total physical size 32.2 Mb",las=1)
abline(h=mean(Recomb_LGIII$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGIII$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGIII$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
mtext(side=4,line=3,"Number of scaffolds",las=0,cex=0.7)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,5,2,3))
plot(Recomb_LGIV$Recomb_genes,type="l",ylim=c(0,290),ylab="Number of genes",xlab="",main="LGIV, Total physical size 27 Mb",las=1)
abline(h=mean(Recomb_LGIV$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGIV$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGIV$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,4,2,4))
plot(Recomb_LGV$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGV, Total physical size 30.1 Mb",las=1)
abline(h=mean(Recomb_LGV$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGV$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGV$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,3,2,5))
plot(Recomb_LGVI$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGVI, Total physical size 28 Mb",las=1)
abline(h=mean(Recomb_LGVI$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGVI$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGVI$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
mtext(side=4,line=3,"Number of scaffolds",las=0,cex=0.7)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,5,2,3))
plot(Recomb_LGVII$Recomb_genes,type="l",ylim=c(0,290),ylab="Number of genes",xlab="",main="LGVII, Total physical size 29.3 Mb",las=1)
abline(h=mean(Recomb_LGVII$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGVII$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGVII$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,4,2,4))
plot(Recomb_LGVIII$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGVIII, Total physical size 26.1 Mb",las=1)
abline(h=mean(Recomb_LGVIII$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGVIII$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGVIII$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,3,2,5))
plot(Recomb_LGIX$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="",main="LGIX, Total physical size 29.6 Mb",las=1)
abline(h=mean(Recomb_LGIX$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGIX$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGIX$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
mtext(side=4,line=3,"Number of scaffolds",las=0,cex=0.7)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,5,2,3))
plot(Recomb_LGX$Recomb_genes,type="l",ylim=c(0,290),ylab="Number of genes",xlab="LG in cM",main="LGX, Total physical size 28.6 Mb",las=1)
abline(h=mean(Recomb_LGX$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGX$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGX$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,4,2,4))
plot(Recomb_LGXI$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="LG in cM",main="LGXI, Total physical size 28.1 Mb",las=1)
abline(h=mean(Recomb_LGXI$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGXI$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGXI$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
par(mar=c(5,3,2,5))
plot(Recomb_LGXII$Recomb_genes,type="l",ylim=c(0,290),ylab="",xlab="LG in cM",main="LGXII, Total physical size 28.4 Mb",las=1)
abline(h=mean(Recomb_LGXII$Recomb_genes),lty=3)
par(new=T)
plot(Recomb_LGXII$Recomb_scaffolds,type="l",col="red",ylim=c(0,290),axes=F,xlab=NA,ylab=NA)
abline(h=mean(Recomb_LGXII$Recomb_scaffolds),lty=3,col="red")
axis(side=4,las=1)
mtext(side=4,line=3,"Number of scaffolds",las=0,cex=0.7)
legend("topright",legend=c("N genes","N scaffolds"),lty=1,col=c("black","red"),bty="n")
