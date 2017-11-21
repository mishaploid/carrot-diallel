# ------------------------------------------------------------------------------
#   Export rankings from BayesDiallel by environment
#   P. Maurizio; edited by S. Turner
#   5 July 2017
# ------------------------------------------------------------------------------

# this script exports means and 95% credibility intervals of rankings for each 
# environment (WI2015, CA2015, and CA2016)

library("BayesDiallel")

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

# ------------------------------------------------------------------------------
# Need to remove burn-in samples for each environment
# ------------------------------------------------------------------------------
source("~/GitHub/carrot-diallel/17_07_18read in AFD objects - separate environments.R")

postPredwi <- vector("list")
postPredca15 <- vector("list")
postPredca16 <- vector("list")

for(i in varNames){
  # WI2015
  AFDwi <- results_wi15[[i]]
  ADOwi <- AFDwi$AllDiallelObs[[1]]
  
  MyPostSum_wi2015 <- PosteriorPredSummary(ADOwi, AFDwi, burnin = 1, AFD=AFDwi, keep = TRUE);
  postPredwi[[i]] <- mcmc.stack.and.burn(ADOwi$PostKeeper$FakeCoda, burnin=500)
  
  # CA2015
  AFDca15 <- results_ca15[[i]]
  ADOca15 <- AFDca15$AllDiallelObs[[1]]
  
  MyPostSum_ca2015 <- PosteriorPredSummary(ADOca15, AFDca15, burnin = 1, AFD=AFDca15, keep = TRUE);
  postPredca15[[i]] <- mcmc.stack.and.burn(ADOca15$PostKeeper$FakeCoda, burnin=500)
  
  # CA2016
  AFDca16 <- results_ca16[[i]]
  ADOca16 <- AFDca16$AllDiallelObs[[1]]
  
  MyPostSum_ca2016 <- PosteriorPredSummary(ADOca16, AFDca16, burnin = 1, AFD=AFDca16, keep = TRUE);
  postPredca16[[i]] <- mcmc.stack.and.burn(ADOca16$PostKeeper$FakeCoda, burnin=500)
}

# ------------------------------------------------------------------------------
# Export ranks - this step takes a long time! 
# ~6 min per trait.. 54 min for all traits
# consider running on remote server
# ------------------------------------------------------------------------------

rankings <- function(postpred, ...){
  ranks <- NULL
  for(i in 1:nrow(postpred)){
    ranks <- rbind(ranks, rank(postpred[i,]))
  }
  mean_rank <- cbind(colMeans(ranks))
  # median_rank <- cbind(colMedians(ranks))
  rank_hpd <- HPDinterval(mcmc(ranks), 0.95)
  cbind(mean = mean_rank, # median = median_rank, 
        rank_hpd)
}

ranks_wi15 <- lapply(postPredwi, FUN = rankings)
ranks_ca15 <- lapply(postPredca15, FUN = rankings)
ranks_ca16 <- lapply(postPredca16, FUN = rankings)

# lapply(ranks_wi15, function(x) write.table( data.frame(x), 
#                                             "~/GitHub/carrot-diallel/results/wi2015/", 
#                                             append= T, sep="," ))

for (i in varNames) {
  filename = paste("~/GitHub/carrot-diallel/results/wi2015/", i, "_ranks_wi2015.csv")
  write.csv(ranks_wi15[[i]], filename)
}

for (i in varNames) {
  filename = paste("~/GitHub/carrot-diallel/results/ca2015/", i, "_ranks_ca2015.csv")
  write.csv(ranks_ca15[[i]], filename)
}

for (i in varNames) {
  filename = paste("~/GitHub/carrot-diallel/results/ca2016/", i, "_ranks_ca2016.csv")
  write.csv(ranks_ca16[[i]], filename)
}

