######################################
##Make a circular plot of linkage groups and connect markers from the same scaffolds
######################################
source("https://bioconductor.org/biocLite.R")
biocLite("OmicCircos")

read.table("Consensus-maps.txt",head=T)->Full_consensus
Full_consensus$scaf<-as.factor(t(as.data.frame(strsplit(as.character(Full_consensus$marker),":")))[,1])

library (OmicCircos)
options (stringsAsFactors=FALSE) ;
set.seed(1234) ;

##Segment frame
cons.name <- paste("LG",Full_consensus$LG,sep=" ")
cons.start <- as.character(round(Full_consensus$consensus),digit=3)
cons.end <- as.character(round(Full_consensus$consensus),digit=3)
cons.f <- data.frame(seg.name=cons.name, seg.start=cons.start,seg.end=cons.end,the.v=runif(length(cons.end)),Note=as.character(Full_consensus$marker))
db <- segAnglePo(cons.f,seg=unique(cons.name));

##Segment mapping
cons.v<-Full_consensus[,c(1,3)]
names(cons.v)<-c("seg.name","seg.pos")
cons.v$seg.pos<-as.character(cons.v$seg.pos)
cons.v$seg.name<-paste("LG",cons.v$seg.name,sep=" ")

###Segment links
cons.links <- NULL
cons.links$chr1 <-NULL
cons.links$pos1 <-NULL
cons.links$marker1 <-NULL
cons.links$chr2 <-NULL
cons.links$pos2 <-NULL
cons.links$marker2 <-NULL
for (i in 1:length(levels(Full_consensus$scaf))){
  tmp <- Full_consensus[Full_consensus$scaf == levels(Full_consensus$scaf)[i],]
  if (dim(tmp)[1] > 1){
    for (j in 2:dim(tmp)[1]){
      cons.links$chr1 <- c(cons.links$chr1,as.character(paste("LG",tmp$LG[j-1],sep=" ")))
      cons.links$chr2 <- c(cons.links$chr2,as.character(paste( "LG", tmp$LG[j],sep=" ")))
      cons.links$pos1 <- c(cons.links$pos1,as.character(tmp$consensus[j-1]))
      cons.links$pos2 <- c(cons.links$pos2,as.character(tmp$consensus[j]))
      cons.links$marker1 <- c(cons.links$marker1,as.character(tmp$marker[j-1]))
      cons.links$marker2 <- c(cons.links$marker2,as.character(tmp$marker[j]))
    }
  }
}

cons.links <- as.data.frame(cons.links)
cons.links <- cons.links[,c(1,3,5,2,4,6)]

intra_cons.links <- cons.links[cons.links$chr1 == cons.links$chr2,]
inter_cons.links <- cons.links[cons.links$chr1 != cons.links$chr2,]

col.intra.links<-NULL
for ( i in 1:dim(intra_cons.links)[1]){
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) >5)){
    col.intra.links[i] <- "grey35"
  }
  else{
    col.intra.links[i] <- "grey90"
  }
}


col.inter.links<-NULL
tmp<-t(as.data.frame(strsplit(inter_cons.links$marker1,":"))[1,])
tmp2<-t(as.data.frame(strsplit(inter_cons.links$marker2,":"))[1,])
tmp<-table(tmp)
tmp2<-table(tmp2)
tmp<-unique(names(c(tmp[tmp>1],tmp2[tmp2>1])))
for (i in 1:dim(inter_cons.links)[1]){
  if (strsplit(inter_cons.links$marker1[i],":")[[1]][1] %in% tmp){
    col.inter.links[i]<-"dark blue"
  } 
  else{
    col.inter.links[i]<-"orange"  
  } 
}

inter_cons.links2<-inter_cons.links[which(col.inter.links == "dark blue"),]
intra_cons.links2<-intra_cons.links[which(col.intra.links== "grey35"),]


###Add recombination hotspots/coldspots
#load("recomb_spots.RData")
Recomb_LGI$LG<-"LG 1"
Recomb_LGI$pos<-row.names(Recomb_LGI)
Recomb_LGII$LG<-"LG 2"
Recomb_LGII$pos<-row.names(Recomb_LGII)
Recomb_LGIII$LG<-"LG 3"
Recomb_LGIII$pos<-row.names(Recomb_LGIII)
Recomb_LGIV$LG<-"LG 4"
Recomb_LGIV$pos<-row.names(Recomb_LGIV)
Recomb_LGV$LG<-"LG 5"
Recomb_LGV$pos<-row.names(Recomb_LGV)
Recomb_LGVI$LG<-"LG 6"
Recomb_LGVI$pos<-row.names(Recomb_LGVI)
Recomb_LGVII$LG<-"LG 7"
Recomb_LGVII$pos<-row.names(Recomb_LGVII)
Recomb_LGVIII$LG<-"LG 8"
Recomb_LGVIII$pos<-row.names(Recomb_LGVIII)
Recomb_LGIX$LG<-"LG 9"
Recomb_LGIX$pos<-row.names(Recomb_LGIX)
Recomb_LGX$LG<-"LG 10"
Recomb_LGX$pos<-row.names(Recomb_LGX)
Recomb_LGXI$LG<-"LG 11"
Recomb_LGXI$pos<-row.names(Recomb_LGXI)
Recomb_LGXII$LG<-"LG 12"
Recomb_LGXII$pos<-row.names(Recomb_LGXII)
Recomb_all_LGs<-rbind(Recomb_LGI,Recomb_LGII,Recomb_LGIII,Recomb_LGIV,Recomb_LGV,Recomb_LGVI,Recomb_LGVII,Recomb_LGVIII,Recomb_LGIX,Recomb_LGX,Recomb_LGXI,Recomb_LGXII)

Recomb_all_LGs<-Recomb_all_LGs[,c(5:6,1:4)]
Recomb_all_LGs$Mb<-Recomb_all_LGs$Recomb_length/1e6
Recomb_all_LGs$genesPerMb<-Recomb_all_LGs$Recomb_genes/Recomb_all_LGs$Mb
Recomb_all_LGs$genesPerMb[Recomb_all_LGs$genesPerMb == "NaN"]<-0
Recomb_all_LGs$genesPerScaffold<-Recomb_all_LGs$Recomb_genes/Recomb_all_LGs$Recomb_scaffolds
Recomb_all_LGs$genesPerScaffold[Recomb_all_LGs$genesPerScaffold == "NaN"]<-1
Recomb_all_LGs$avg_Scaf_size<-Recomb_all_LGs$Mb/Recomb_all_LGs$Recomb_scaffolds
Recomb_all_LGs$avg_Scaf_size[Recomb_all_LGs$avg_Scaf_size == "NaN"]<-0

pdf("Circosplot_with_scaf-number_genes-per-scaffold_kb_genes-per-kb.pdf")
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
circos(R=350, cir=db, type="chr",col="grey50", print.chr.lab=TRUE, W=4, scale=FALSE)
circos(R=310, cir=db, mapping=cons.v, type="b3",col="black", W=40 ,lwd=1,B=TRUE,scale=FALSE)
circos(R=275, cir=db, mapping=Recomb_all_LGs[,c(1:2,4)], type="l",col="purple", W=40 ,lwd=1,B=TRUE,scale=TRUE)
circos(R=240, cir=db, mapping=Recomb_all_LGs[,c(1:2,9)], type="l",col="brown", W=40 ,lwd=1,B=TRUE,scale=TRUE)
circos(R=205, cir=db, mapping=Recomb_all_LGs[,c(1:2,7)], type="l",col="red", W=40 ,lwd=1,B=TRUE,scale=TRUE)
circos(R=170, cir=db, mapping=Recomb_all_LGs[,c(1:2,8)], type="l",col="dark green", W=40 ,lwd=1,B=TRUE,scale=TRUE)
circos(R=165, cir=db, mapping=intra_cons.links, type="link2",B=TRUE, col=col.intra.links)
circos(R=165, cir=db, mapping=intra_cons.links2, type="link2", col="grey35")
circos(R=100, cir=db, mapping=inter_cons.links, type="link",B=TRUE, col=col.inter.links)
circos(R=100, cir=db, mapping=inter_cons.links2, type="link", col="dark blue")
legend("topright",c("Number of scaffolds (scale 0-239)","Genes/scaffold (scale 1-3)","Physical size (scale 0-5.56 Mb)","Genes/Mb (scale 0-376)"),lty=1,bty="n",col=c("purple","brown","red","dark green"),cex=0.7)
dev.off()

write.table(unique(t(as.data.frame(strsplit(c(inter_cons.links$marker1,inter_cons.links$marker2),split=":"))[1,])),"inter_split_scaffolds.txt",quote=F,col.names=F,row.names=F)
write.table(unique(t(as.data.frame(strsplit(c(intra_cons.links$marker1[col.intra.links == "grey35"],intra_cons.links$marker2[col.intra.links == "grey35"]),split=":"))[1,])),"intra_split_scaffolds.txt",quote=F,col.names=F,row.names=F)
