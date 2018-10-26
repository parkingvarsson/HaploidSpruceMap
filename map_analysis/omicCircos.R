######################################
##Make a circular plot of linkage groups and connect markers from the same scaffolds
######################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("OmicCircos")

read.table("../Consensus_maps.txt",head=T)->maps


library (OmicCircos)
options (stringsAsFactors=FALSE) ;
set.seed(1234) ;

##Segment frame
cons.name <- paste("LG",maps$LG,sep=" ")
cons.start <- as.character(round(maps$Consensus),digit=3)
cons.end <- as.character(round(maps$Consensus),digit=3)
cons.f <- data.frame(seg.name=cons.name, seg.start=cons.start,seg.end=cons.end,the.v=runif(length(cons.end)),Note=as.character(maps$Marker))
db <- segAnglePo(cons.f,seg=unique(cons.name),angle.start = 3,angle.end = 357);

##Segment mapping
cons.v<-maps[,c(1,3)]
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
for (i in 1:length(levels(maps$Scaffold))){
  tmp <- maps[maps$Scaffold == levels(maps$Scaffold)[i],]
  if (dim(tmp)[1] > 1){
    for (j in 2:dim(tmp)[1]){
      cons.links$chr1 <- c(cons.links$chr1,as.character(paste("LG",tmp$LG[j-1],sep=" ")))
      cons.links$chr2 <- c(cons.links$chr2,as.character(paste( "LG", tmp$LG[j],sep=" ")))
      cons.links$pos1 <- c(cons.links$pos1,as.character(tmp$Consensus[j-1]))
      cons.links$pos2 <- c(cons.links$pos2,as.character(tmp$Consensus[j]))
      cons.links$marker1 <- c(cons.links$marker1,as.character(tmp$Marker[j-1]))
      cons.links$marker2 <- c(cons.links$marker2,as.character(tmp$Marker[j]))
    }
  }
}

cons.links <- as.data.frame(cons.links)
cons.links <- cons.links[,c(1,3,5,2,4,6)]

intra_cons.links <- cons.links[cons.links$chr1 == cons.links$chr2,]
inter_cons.links <- cons.links[cons.links$chr1 != cons.links$chr2,]

col.intra.links<-NULL
for ( i in 1:dim(intra_cons.links)[1]){
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i])) >5  & abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) <10)){
    col.intra.links[i] <- "grey35"
  }
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) >10)){
    col.intra.links[i] <- "red"
  }
  if(abs(as.numeric(intra_cons.links$pos2[i])-as.numeric(intra_cons.links$pos1[i]) < 5)){
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
intra_cons.links3<-intra_cons.links[which(col.intra.links== "red"),]


pdf("Circosplot_with_split_Scaffolds.pdf")
par(mar=c(2,2,2,2),mfrow=c(1,1))
plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
circos(R=350, cir=db, type="chr",col="grey50", print.chr.lab=FALSE, W=4, scale=FALSE)
circos(R=310, cir=db, mapping=cons.v, type="b3",col="black", W=40 ,lwd=1,B=TRUE,scale=FALSE)
#circos(R=170, cir=db, mapping=Recomb_all_LGs[,c(1:2,5,6,3)], type="ml",col=c(rep("black",2),"red"), W=140 ,lwd=1,B=TRUE,scale=FALSE)
#circos(R=210, cir=db, mapping=Recomb_all_LGs[,c(1:2,13,14,11)], type="ml",col=c(rep("black",2),"purple"), W=50 ,lwd=1,B=TRUE,scale=FALSE)
#circos(R=170, cir=db, mapping=Recomb_all_LGs[,c(1:2,17,18,15)], type="ml",col=c(rep("black",2),"blue"), W=50 ,lwd=1,B=TRUE,scale=FALSE)
circos(R=305, cir=db, mapping=intra_cons.links, type="link2",B=TRUE, col=col.intra.links)
circos(R=305, cir=db, mapping=intra_cons.links2, type="link2", col="grey35")
circos(R=305, cir=db, mapping=intra_cons.links3, type="link2", col="red")
circos(R=185, cir=db, mapping=inter_cons.links, type="link",B=TRUE, col=col.inter.links)
circos(R=185, cir=db, mapping=inter_cons.links2, type="link", col="dark blue")
text(400,730,"A")
text(400,650,"B")
text(400,550,"C")
text(550,750,"LG I")
text(730,600,"LG II")
text(790,430,"LG III")
text(750,240,"LG IV")
text(620,85,"LG V")
text(450,30,"LG VI")
text(280,40,"LG VII")
text(100,150,"LG VIII")
text(20,340,"LG IX")
text(30,500,"LG X")
text(110,650,"LG XI")
text(250,750,"LG XII")
dev.off()

write.table(unique(t(as.data.frame(strsplit(c(inter_cons.links$marker1,inter_cons.links$marker2),split=":"))[1,])),"inter_split_scaffolds.txt",quote=F,col.names=F,row.names=F)
write.table(unique(t(as.data.frame(strsplit(c(intra_cons.links$marker1[col.intra.links != "grey90"],intra_cons.links$marker2[col.intra.links != "grey90"]),split=":"))[1,])),"intra_split_scaffolds.txt",quote=F,col.names=F,row.names=F)

save.image("Circos_plot.RData")