# ------------------------------------------------------------------------------
#   Estimate degree of dominance
#   See Maurizio et al. (2017) for details
#   P. Maurizio (main scripts) & S. Turner (minor modifications)
#   Ref: Comstock & Robinson (1948); Gardner & Lonnquist (1959); Kacser & Burns (1981)
#   12 September 2016
# ------------------------------------------------------------------------------

# This script uses BayesDiallel results to estimate the degree of dominance
# (Comstock & Robinson) and the dominance index (Kacser & Burns) for crosses in
# a diallel

setwd("~/GitHub/carrot-diallel/data")
source("~/GitHub/carrot-diallel/07_read_AFD_objects.R")

#' @title mcmc.stack.and.burn: Stack MCMC chains, after burning each
#' @description This is a modification of mcmc.stack that also removes specified number of burn-ins from chains.
#' 
#' @param coda.object this is the mcmc object
#' @param burnin this is the common amount to burn off each chain
#' @param ... additional arguments
#' @return returns single chain, stacked and burned MCMC object
#' @examples
#' NULL
#' @export
mcmc.stack.and.burn <- function (coda.object, burnin, ...){
  burnin <- burnin+1
  if (inherits(coda.object, "mcmc")) {
    return(coda.object)
  }
  if (!inherits(coda.object, "mcmc.list")) {
    stop("Non-mcmc object passed to function\n")
  }
  len <- dim(coda.object[[1]])[1]
  chain <- coda.object[[1]][c(burnin:len),]
  for (i in 2:nchain(coda.object)) {
    chain <- rbind(chain, coda.object[[i]][c(burnin:len),])
  }
  as.mcmc(chain)
}

# loop through all phenotypes
# output: mean/median ranks, degree of dominance (CR), dominance index (KB)

for(trait in names(results)){
  # read in data for trait
  AFD <- results[[trait]]
  
  head(AFD$AllDiallelObs[[1]]$centered.chains[,"Sigma:1"])
  
  # Need to remove burn-in samples
  
  ADO <- AFD$AllDiallelObs[[1]]
  MyPostSum <- PosteriorPredSummary(ADO, AFD, burnin = 1, AFD=AFD, keep = TRUE);
  postPred <- mcmc.stack.and.burn(ADO$PostKeeper$FakeCoda, burnin=500)
  
  plot(rank(postPred[1,]), pch=16, col=rgb(0,0,0,alpha=0.2), xaxt='n', ylab="rank",
       xlab="cross")
  axis(side=1, at = c(1:length(colnames(postPred))), labels = colnames(postPred), las=2)
  for(i in c(2:1000)){
    #for(i in c(2:nrow(postPred))){
    points(rank(postPred[i,]), pch=16, col=rgb(0,0,0,alpha=0.005))
  }
  
  # 1. mean/median ranks
  
  ranks <- NULL
  
  for(i in 1:nrow(postPred)){
    ranks <- rbind(ranks, rank(postPred[i,]))
  }
  
  write.csv(colMeans(ranks), paste0("~/Github/carrot-diallel/results/", trait, 
                                    "/cross.mean.ranks_", trait, ".csv", sep=""))
  
  write.csv(colMedians(ranks), paste0("~/GitHub/carrot-diallel/results/", trait,
                                      "/cross.median.ranks_", trait, ".csv", sep=""))
  
  # 2. dominance of strains - Kascer-Burns
  
  dominance.index <- NULL
  j <- 0
  
  for(i in c(1:6)){
    for(k in c(1:6)){
      if(k != i){
        j <- j + 1
        mt.mt <- paste0("i:", i, ";k:", i)
        mt.wt <- paste0("i:", i, ";k:", k)
        wt.mt <- paste0("i:", k, ";k:", i)
        wt.wt <- paste0("i:", k, ";k:", k)
        dominance.index <- cbind(dominance.index, 
                                 (postPred[,wt.wt] - rowMeans(cbind(postPred[,wt.mt], postPred[,mt.wt])))/(postPred[,wt.wt] - postPred[,mt.mt]))
        colnames(dominance.index)[j] <- wt.mt
      }
    }
  }
  
  write.csv(dominance.index, paste0("~/GitHub/carrot-diallel/results/", trait, 
                                    "/cross.dominance.index-KB_", trait, ".csv", sep=""), 
            row.names=FALSE)
  
  # 3. dominance of strains - Comstock & Robinson
  
  degree.of.dominance <- NULL
  j <- 0
  
  for(i in c(1:6)){
    for(k in c(1:6)){
      if(k != i){
        j <- j + 1
        mt.mt <- paste0("i:", i, ";k:", i)
        mt.wt <- paste0("i:", i, ";k:", k)
        wt.mt <- paste0("i:", k, ";k:", i)
        wt.wt <- paste0("i:", k, ";k:", k)
        degree.of.dominance <- cbind(degree.of.dominance, 
                                     1-2*(postPred[,wt.wt] - rowMeans(cbind(postPred[,wt.mt], postPred[,mt.wt])))/(postPred[,wt.wt] - postPred[,mt.mt]))
        colnames(degree.of.dominance)[j] <- wt.mt
      }
    }
  }
  
  write.csv(dominance.index, paste0("~/GitHub/carrot-diallel/results/", trait, 
                                    "/cross.dominance.index-CR_", trait, ".csv", sep=""), 
            row.names=FALSE)
  
  # plot degree of dominance
  
  lower <- -2; upper <- 2
  low <- -2; up <- 2
  lwd <- 3
  ylim <- NULL
  
  j <- 0
  
  pdf(paste0("~/GitHub/carrot-diallel/results/", trait, 
             "/degree-of-dominance-confint-CR_", trait, ".pdf", sep=""), 
      width=12, height=12)
  par(mfrow=c(6,6))
  par(mar=c(3,3,1,1))
  
  for(i in c(1:6)){
    for(k in c(1:6)){
      if(k != i){
        j <- j + 1
        plot(density(degree.of.dominance[,j], from=low, to=up), xlim=c(lower, upper), 
             main="", lwd=lwd, las=2, xlab="", yaxt="n", ylab="", ylim=c(0,2.5))
        abline(v=c(-2, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2), col="grey", lty=2)
        abline(v=mean(degree.of.dominance[,j]), col="red", lty=3)
        abline(v=median(degree.of.dominance[,j]), col="blue")
        conf <- quantile(degree.of.dominance[,j], probs=c(0.025, 0.975))
        thisgrey <- rgb(0,0,0,maxColorValue=1, alpha=0.1)
        rect(xleft=conf[1], ybottom=-1, xright=conf[2], ytop=3.0,
             col=thisgrey, border=thisgrey)
        mtext(colnames(degree.of.dominance)[j], side=2, line=1)
      }else{
        plot(1,1,col="white", xlab="", ylab="", xaxt='n', yaxt='n')
      }
    }
  }
  dev.off()
  
}
