library("data.table")

cluseq<-fread("white-spruce-map-cluseq.txt",head=T)
cluseq$LGn.ws<-as.numeric(substr(cluseq$LG,3,4))


abies<-fread("Consensus-maps.txt")
abies$scaffold<-unlist(lapply(abies$marker,strsplit,split=":"))[seq(1,30010,2)]


blast<-fread("arborea/cluseqs.out",head=F)
names(blast)<-c("query","match","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")

comp<-merge(abies,blast,by.x="scaffold",by.y="match",all=T)
comp<-merge(comp,cluseq,by.x="query",by.y="seq",all=T)
comp<-comp[(!is.na(comp$query))&(!is.na(comp$scaffold)),]

comp<-comp[comp$pident>95,]
pdf("comparative_map_Pavy2017.pdf")
par(mfrow=c(12,12))
par(mar=c(0,0,0,0))
par(omi=c(1,1,0,0))
par(xpd=TRUE)
mismatches<-NULL
for(i in 1:12){
	for(j in 1:12){
		tmp<-subset(comp,LGn.ab==j&LGn.ws==i)
		if(dim(tmp)[1]>150){
			print(cor(tmp$cM.ab,tmp$cM.ws,meth="kend"))
		}else{
			mismatches<-c(mismatches,dim(tmp)[1])
		}
		if(dim(tmp)[1]>0){
			col<-ifelse(dim(tmp)[1]>150,"black","grey65")
			plot(cM.ws~cM.ab,tmp,xlab="",ylab="",col=col,pch=19,xaxt="n",yaxt="n",cex=0.7)
		}else{
			plot(c(0,0),c(0,0),type="n",xlab="",ylab="",col=col,pch=19,xaxt="n",yaxt="n",cex=0.7)
		}
	}
}
mtext(c("LG I","LG II","LG III","LG IV","LG V","LG VI","LG VII","LG VIII","LG IX","LG X","LG XI","LG XII"),side=1,line=1, outer=TRUE, at=seq(0.04,1,0.083))
mtext(c("LG 12","LG 11","LG 10","LG 9","LG 8","LG 7","LG 6","LG 5","LG 4","LG 3","LG 2","LG 1"),side=2,line=1, outer=TRUE, at=seq(0.04,1,0.083),las=1)
mtext(expression(paste(italic("P. abies"),", haploid consensus map")),side=1, line=4,outer=T,at=0.5)
mtext(expression(paste(italic("P. glauca"),", Pavy et al. 2017 map")),side=2, line=6,outer=T,at=0.5)
dev.off()
rm(tmp)