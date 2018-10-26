load("C1_cleaned_maps.RData")

create_recombination_plot<-function(map)
{
    link.phases <- matrix(NA, length(map$seq.num), 2)
    link.phases[1, ] <- rep(1, 2)
    for (j in 1:length(map$seq.phases)) {
      switch(EXPR = map$seq.phases[j], link.phases[j + 1, 
        ] <- link.phases[j, ] * c(1, 1), link.phases[j + 
        1, ] <- link.phases[j, ] * c(1, -1), link.phases[j + 
        1, ] <- link.phases[j, ] * c(-1, 1), link.phases[j + 
        1, ] <- link.phases[j, ] * c(-1, -1), )
    }
    link.phases <- apply(link.phases, 1, function(x) paste(as.character(x), 
                                                           collapse = "."))
    parents <- matrix("", length(map$seq.num), 4)
    for (k in 1:length(map$seq.num)) parents[k, ] <- return.geno(get(map$data.name, 
                                                                        pos = 1)$segr.type[map$seq.num[k]], link.phases[k])
    mother<-t(apply(parents[,1:2],1,function(x)(gsub("a",1,x))))
    mother<-t(apply(mother,1,function(x)(gsub("b",2,x))))
    offspring<-as.data.frame(t(get(map$data.name, pos = 1)$geno[,map$seq.num]))
    cross_over<-as.data.frame(matrix(0,ncol=ncol(offspring),nrow=nrow(offspring)))
    for(i in 1:nrow(mother)){
      cross_over[i,which(offspring[i,] == mother[i,1])]<-1
      cross_over[i,which(offspring[i,] == mother[i,2])]<-2
    }
    
    library(pheatmap)
    pheatmap(as.matrix(cross_over),cluster_cols = FALSE,cluster_rows = FALSE,)
   
}

create_recombination_plot(C1_LG10)
