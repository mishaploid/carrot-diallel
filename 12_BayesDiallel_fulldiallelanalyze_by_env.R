# ------------------------------------------------------------------------------
#   Use 'BayesDiallel' package to create Full Diallel Analyze objects
#   Separated by environment
#   S. Turner 
#   20 June 2017
# ------------------------------------------------------------------------------

# This script runs the DiallelAnalyzer function from BayesDiallel on all traits
# of interest and stores the AFD objects for later use

# Data is split into separate environments and run separately
# Environment is not included as a random effect

try(require(lattice, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE)
try(require(coda, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE)
try(require(R.methodsS3, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE)
try(require(R.oo, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE)
try(require(corpcor, quietly=TRUE, warn.conflicts=FALSE), silent=TRUE)
require(BayesDiallel, quietly = TRUE, warn.conflicts=FALSE)
require(R.oo, warn.conflicts=FALSE, quietly=TRUE)

# ------------------------------------------------------------------------------
# read in data, set factors, and apply transformations
# ------------------------------------------------------------------------------
setwd("~/GitHub/carrot-diallel/data")
dialleltmp <- read.csv("diallelRawData.csv", header = TRUE)

dialleltmp$year <- as.factor(dialleltmp$year)
dialleltmp$rep <- as.factor(dialleltmp$rep)
dialleltmp$geno <- as.factor(dialleltmp$geno)
dialleltmp$male <- as.factor(dialleltmp$male)
dialleltmp$female <- as.factor(dialleltmp$female)
dialleltmp$stand[is.na(dialleltmp$stand)] <- 0 
dialleltmp$flw <- log(dialleltmp$flw)
dialleltmp$dlw <- log(dialleltmp$dlw)
dialleltmp$frw <- log(dialleltmp$frw)
dialleltmp$drw <- log(dialleltmp$drw)
dialleltmp$ratio <- dialleltmp$dlw/dialleltmp$drw
dialleltmp$is.female <- as.logical(0)

# subset by environment
library(plyr)
dialleltmp$env <- revalue(dialleltmp$year, c("1" = "CA2015", "2" = "WI2015", "3" = "CA2016"))
wi2015 <- subset(dialleltmp, env == "WI2015")
ca2015 <- subset(dialleltmp, env == "CA2015")
ca2016 <- subset(dialleltmp, env == "CA2016")


# ------------------------------------------------------------------------------
# load packages, model information, and priors
# ------------------------------------------------------------------------------

# attach Piximus data to extract models
data(PiximusData)

# load tau prior information 
Diallel.tau.Prior.Info <- read.csv("Diallel_tau_prior.csv", header = TRUE, 
                                   row.names = 1)
print(Diallel.tau.Prior.Info)

# load models - see BayesDiallel documentation for details
str(Example.Piximus.Models) 

# vector with names of traits
varNames <- names(dialleltmp[,9:17]) 

# ------------------------------------------------------------------------------
# generate fullDiallelAnalyze objects and save to unique directories
# 'fullu' model w/ stand as fixed and random effect
# run environments separately
# ------------------------------------------------------------------------------
results_wi2015 <- list()

for (i in varNames) {
  dir.create(paste0("~/GitHub/carrot-diallel/results/wi2015", sep = ""))
  setwd(paste0("~/GitHub/carrot-diallel/results/wi2015", sep = ""))
  results_wi2015[[i]] <- DiallelAnalyzer(data = wi2015, father.strain = "male",
                                  mother.strain = "female", phenotype = i,
                                  is.female = "is.female", sep="", 
                                  FixedEffects = "stand",
                                  RandomEffects = "stand",
                                  Models = "fullu",
                                  sigmasq.start = 1,
                                  numChains = 5,
                                  lengthChains = 10000, 
                                  burnin = 1000, 
                                  DIC.Only = FALSE,
                                  tauPriorFile = Diallel.tau.Prior.Info,
                                  SaveAFDFile = paste("AFDbackup_wi2015_", i, 
                                                      ".RData", sep = ""))
}

results_ca2015 <- list()

for (i in varNames) {
  dir.create(paste0("~/GitHub/carrot-diallel/results/ca2015", sep = ""))
  setwd(paste0("~/GitHub/carrot-diallel/results/ca2015", sep = ""))
  results_ca2015[[i]] <- DiallelAnalyzer(data = ca2015, father.strain = "male",
                                         mother.strain = "female", phenotype = i,
                                         is.female = "is.female", sep="", 
                                         FixedEffects = "stand",
                                         RandomEffects = "stand",
                                         Models = "fullu",
                                         sigmasq.start = 1,
                                         numChains = 5,
                                         lengthChains = 10000, 
                                         burnin = 1000, 
                                         DIC.Only = FALSE,
                                         tauPriorFile = Diallel.tau.Prior.Info,
                                         SaveAFDFile = paste("AFDbackup_ca2015_", i, 
                                                             ".RData", sep = ""))
}

results_ca2016 <- list()

for (i in varNames) {
  dir.create(paste0("~/GitHub/carrot-diallel/results/ca2016", sep = ""))
  setwd(paste0("~/GitHub/carrot-diallel/results/ca2016", sep = ""))
  results_ca2016[[i]] <- DiallelAnalyzer(data = ca2016, father.strain = "male",
                                         mother.strain = "female", phenotype = i,
                                         is.female = "is.female", sep="", 
                                         FixedEffects = "stand",
                                         RandomEffects = "stand",
                                         Models = "fullu",
                                         sigmasq.start = 1,
                                         numChains = 5,
                                         lengthChains = 10000, 
                                         burnin = 1000, 
                                         DIC.Only = FALSE,
                                         tauPriorFile = Diallel.tau.Prior.Info,
                                         SaveAFDFile = paste("AFDbackup_ca2016_", i, 
                                                             ".RData", sep = ""))
}

