library(onemap)
load(".RData")
load("rec40_LGs_dropped.RData")
set.map.fun("kosambi")

args <- commandArgs(trailingOnly = TRUE)

LG <- get(paste("LG",args[1],"rec40_dropped",sep="_"))

size_pred <- pick_batch_sizes(LG, size = 50,overlap = 30, around = 10)
full <- map_overlapping_batches(input.seq=LG,
                                size=size_pred,
                                overlap=30,
                                fun.order=NULL,
                                phase.cores=4,
                                max.dist = 50,
                                max.tries = 10,
                                min.tries = 0,
                                verbosity="batch")

out_full <- paste("LG", args[1], "rec40_droppedMap", sep = "_")
assign(out_full,full)

save(list=out_full,file=paste("LG",args[1],"rec40_droppedMap.RData",sep="_"))