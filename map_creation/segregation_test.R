test_segregation_of_a_marker <- function(x, marker) {
  options(warn=-1) #Supress warnings from function chisq.test - when cell have < 5 elements
  ## Segregation pattern for each marker type
  p.a <- rep(1/4, 4); p.b <- c(1/4, 1/2, 1/4); p.c <- c(3/4, 1/4); p.d <- rep(1/2, 2)
  ## Counting each category, removing missing data (coded as 0)
  count <- table(x$geno[,marker], exclude=0)
  ## Do the chisq test, using the appropriate expected segregation
  ## grepl() allows finding the marker type (it has the letter in the argument)
  ## Some markers may not have data for some classes, so fill them with 0
  if (grepl("A",x$segr.type[marker]) & class(x)=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    if (is.element(3,x$geno[,marker])) c3 <- count[names(count)==3] else c3 <- 0
    if (is.element(4,x$geno[,marker])) c4 <- count[names(count)==4] else c4 <- 0
    qui <- chisq.test(as.vector(c(c1,c2,c3,c4)), p=p.a, correct = FALSE)
    H0 <- "1:1:1:1"
  }
  else if (grepl("B",x$segr.type[marker]) & class(x)=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    if (is.element(3,x$geno[,marker])) c3 <- count[names(count)==3] else c3 <- 0
    qui <- chisq.test(as.vector(c(c1,c2,c3)), p=p.b, correct = FALSE)
    H0 <- "1:2:1"
  }
  else if (grepl("C",x$segr.type[marker]) & class(x)=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.c, correct = FALSE)
    H0 <- "3:1"
  }
  else if (grepl("D",x$segr.type[marker]) & class(x)=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.d, correct = FALSE)
    H0 <- "1:1"
  }
  #impossible to test: dominant and co-dominant mixed in the same marker
  #however, it will not work with "qui <- NA"; using NULL instead
  else if (grepl("M.X",x$segr.type[marker])) {
    qui <- NULL
    qui$statistic <- NA
    qui$p.value <- NA
    H0 <- "Check data!"
  }
  
  return(list(Hypothesis=H0, qui.quad=qui$statistic, p.val=qui$p.value,
              perc.genot=100*(sum(table(x$geno[,marker], exclude=0))/x$n.ind)))
}




test_segregation <- function(x) {
  if (is(x,"outcross")) {
    y <- list(Marker=dimnames(x$geno)[[2]],
              Results.of.tests=sapply(1:x$n.mar, function(batchmap.object, marker)
                test_segregation_of_a_marker(batchmap.object, marker),
                batchmap.object=x))
    # sapply iterates from 1 to x$n.mar; x is fixed (onemap object with data)
    class(y) <- c("batchmap_segreg_test")
    invisible(y) #returns y without showing it
  }
  else stop("This is not an BatchMap object with raw data")
}




print.onemap_segreg_test <- function(x,...) {
  Z <- data.frame(Marker=x$Marker,
                  H0=unlist(x$Results.of.tests[1,]),
                  Chi.square=unlist(x$Results.of.tests[2,]),
                  p.value=unlist(x$Results.of.tests[3,]),
                  Perc.genot=round(unlist(x$Results.of.tests[4,]),2))
  colnames(Z) <- c("Marker","H0","Chi-square","p-value","% genot.")
  return(Z)
}




plot.onemap_segreg_test <- function(x, order=TRUE,...) {
  # Create a data frame
  Z <- data.frame(Marker=x$Marker,
                  X.square=unlist(x$Results.of.tests[2,]),
                  p.value=unlist(x$Results.of.tests[3,]))
  Bonf <- -log10(.05/nrow(Z)) #Bonferroni's threshold'
  Z$signif <- factor(ifelse(-log10(Z$p.value)<Bonf,"non sign.","sign."))
  Z$order <- 1:nrow(Z)
  # % of distorted
  perc <- 100*(1-(table(Z$signif)[1]/nrow(Z)))
  # Keeping markers in their original order (not alphanumeric), or by p-values (default)
  if (order!=TRUE) Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
  else Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$p.value, decreasing=TRUE)])
  # Plotting
  par(mar=c(3,2,5,2))
  g <- ggplot(data=Z, aes(x=Marker, y=-log10(p.value)))
  g <- g + ylab(expression(-log[10](p-value)))
  g <- g + geom_point(aes(color=signif), stat="identity",size=1.5)
  g <- g + scale_colour_manual(name=paste("Bonferroni\n","(",round(perc,0),"% distorted)",sep=""),
                               values = c("#46ACC8","#B40F20"))
  g <- g + geom_hline(yintercept = Bonf, colour="#E58601", linetype = "longdash")
  g <- g + coord_flip()
  if (nrow(Z)>30) g <- g + theme(axis.text.y = element_blank())
  g
}


plot_segreg_test <- function(x, order=TRUE,Title) {
  # Create a data frame
  Z <- data.frame(Marker=x$Marker,
                  X.square=unlist(x$Results.of.tests[2,]),
                  p.value=unlist(x$Results.of.tests[3,]))
  Bonf <- -log10(.05/nrow(Z)) #Bonferroni's threshold'
  Z$signif <- factor(ifelse(-log10(Z$p.value)<Bonf,"non sign.","sign."))
  Z$order <- 1:nrow(Z)
  # % of distorted
  perc <- 100*(1-(table(Z$signif)[1]/nrow(Z)))
  # Keeping markers in their original order (not alphanumeric), or by p-values (default)
  if (order!=TRUE) Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
  else Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$p.value, decreasing=TRUE)])
  # Plotting
  if (max(-log10(Z$p.value)) < Bonf){
      lim<- c(0,Bonf)
  }  else{
    lim<-c(0,max(-log10(Z$p.value)))
  }
  par(mar=c(5,3,4,2)+0.1)
  plot.new()
  plot.window(xlim=lim,ylim=c(0,length(Z$Marker)))
  points(-log10(Z$p.value),Z$Marker,col=as.integer(Z$signif),cex=0.5,pch=19)
  axis(1)
  mtext(side = 1,text = expression(-log[10](p-value)),line=3)
  mtext(side = 2,text = "Marker",line=1)
  mtext(side=3,text=Title,line=2)
  abline(v=Bonf, col="orange", lty=2)
  text(x=Bonf-0.5, y=0.5*length(Z$Marker),paste("Bonferroni\n","(",round(perc,0),"% distorted)",sep=""),cex=0.8)
}






Bonferroni_alpha <- function(x, global.alpha=0.05) {
  if (!is(x,"batchmap_segreg_test")) stop("This is not an object of class batchmap_segreg_test")
  alpha.Bonf <- global.alpha/length(x$Marker)
  return(alpha.Bonf)
}




select_segreg <- function(x, distorted=FALSE, numbers=FALSE) {
  if (!is(x,"batchmap_segreg_test")) stop("This is not an object of class batchmap_segreg_test")
  Z <- data.frame(Marker=x$Marker,
                  p.value=unlist(x$Results.of.tests[3,]))
  if (distorted==FALSE) Z <- subset(Z, p.value>=Bonferroni_alpha(x))
  else Z <- subset(Z, p.value<Bonferroni_alpha(x))
  if (numbers==TRUE) return(which(x$Marker %in% as.vector(Z[,1])))
  else return(as.vector(Z[,1]))
}