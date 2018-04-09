# ------------------------------------------------------------------------------
#   Use 'BayesDiallel' package to create Full Diallel Analyze objects
#   S. Turner 
#   12 September 2016
# ------------------------------------------------------------------------------

# This script runs the DiallelAnalyzer function from BayesDiallel on all traits
# of interest and stores the AFD objects for later use

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
setwd("~/GitHub/carrot-diallel")
dialleltmp <- read.csv("data/diallelRawData.csv", header = TRUE)

dialleltmp$year <- as.numeric(dialleltmp$year)
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

# ------------------------------------------------------------------------------
# load model information and priors
# ------------------------------------------------------------------------------

# attach Piximus data to extract models
data(PiximusData)

# load tau prior information 
Diallel.tau.Prior.Info <- read.csv("data/Diallel_tau_prior.csv", header = TRUE, 
                                   row.names = 1)
print(Diallel.tau.Prior.Info)

# load models - see BayesDiallel documentation for details
str(Example.Piximus.Models) 

# vector with names of traits
varNames <- names(dialleltmp[,9:17]) 

# ------------------------------------------------------------------------------
# generate fullDiallelAnalyze objects and save to unique directories
# 'fullu' model
# stand (i.e. planting density) included as fixed effect
# year (i.e. location) and stand included as random effects
# ------------------------------------------------------------------------------
results <- list()

for (i in varNames) {
  dir.create(paste0("~/GitHub/carrot-diallel/results/", i, sep = ""))
  setwd(paste0("~/GitHub/carrot-diallel/results/", i, sep = ""))
  results[[i]] <- DiallelAnalyzer(data = dialleltmp, father.strain = "male",
                                  mother.strain = "female", phenotype = i,
                                  is.female = "is.female", sep="", 
                                  FixedEffects = "stand",
                                  RandomEffects = c("year", "stand"),
                                  Models = "fullu",
                                  sigmasq.start = 1,
                                  numChains = 5,
                                  lengthChains = 10000, 
                                  burnin = 1000, 
                                  DIC.Only = FALSE,
                                  tauPriorFile = Diallel.tau.Prior.Info,
                                  SaveAFDFile = paste("SaveAFDBackUp", i, 
                                                      ".RData", sep = ""))
}
