# ------------------------------------------------------------------------------
#   Plot rankings by year
#   S. Turner 
#   18 July 2016
# ------------------------------------------------------------------------------

# This script reads in AFD objects from the DiallelAnalyzer function and stores
# them in a list ('results')

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
setwd("~/Documents/carrot-diallel-all-files/carrot-diallel-scripts")
dialleltmp <- read.csv("diallelRawData.csv", header = TRUE)

dialleltmp$year <- as.numeric(dialleltmp$year)
dialleltmp$rep <- as.factor(dialleltmp$rep)
dialleltmp$geno <- as.factor(dialleltmp$geno)
dialleltmp$male <- as.factor(dialleltmp$male)
dialleltmp$female <- as.factor(dialleltmp$female)
dialleltmp$stand <- as.factor(dialleltmp$stand)
dialleltmp$ratio <- dialleltmp$dlw/dialleltmp$drw
dialleltmp$flw <- log(dialleltmp$flw)
dialleltmp$dlw <- log(dialleltmp$dlw)
dialleltmp$frw <- log(dialleltmp$frw)
dialleltmp$drw <- log(dialleltmp$drw)
dialleltmp$is.female <- as.logical(0)

# ------------------------------------------------------------------------------
# loop to read in AFD files
# ------------------------------------------------------------------------------

# select traits of interest and store in a vector
varNames <- names(dialleltmp[,9:17])

results_wi15 <- vector("list")

for (i in varNames) {
  setwd("~/Documents/carrot-diallel-all-files/BayesDiallel/wi2015")
  x <- load(paste0("AFDBackUp_wi2015_", i, ".RData", sep = ""))
  results_wi15[[i]] <- get(x)
  setwd("~/Documents/carrot-diallel-all-files/")
}

results_ca15 <- vector("list")

for (i in varNames) {
  setwd("~/Documents/carrot-diallel-all-files/BayesDiallel/ca2015")
  x <- load(paste0("AFDBackUp_ca2015_", i, ".RData", sep = ""))
  results_ca15[[i]] <- get(x)
  setwd("~/Documents/carrot-diallel-all-files/")
}

results_ca16 <- vector("list")

for (i in varNames) {
  setwd("~/Documents/carrot-diallel-all-files/BayesDiallel/ca2016")
  x <- load(paste0("AFDBackUp_ca2016_", i, ".RData", sep = ""))
  results_ca16[[i]] <- get(x)
  setwd("~/Documents/carrot-diallel-all-files/")
}
