# ------------------------------------------------------------------------------
#   Create p2 data frame
#   S. Turner 
#   12 September 2016
# ------------------------------------------------------------------------------

# This script compiles posterior PSq values (VarP in Crowley et al. 2014) into 
# datasets for plotting

# ------------------------------------------------------------------------------
# load AFD objects
# ------------------------------------------------------------------------------
source("~/Documents/carrot-diallel/09 read in AFD objects.R")

# ------------------------------------------------------------------------------
# Diallel Variance Projection (PSq)
# ------------------------------------------------------------------------------
# Proportion of phenotypic variance explained by genetic effects
# Effect of this inheritance component, as expressed in future un-forseen 
#   diallels using new strains sampled independently from a similar population 
# Depends on fitted tau^2^ values   
# See Appendix A in _Lenarcic et al. 2012_
# Output is proportion of heritability attributed to each component
# ot all effects are equal - only inbreds are subject to inbreeding, so 
#   inbreeding component contributes 1/J to heritability in a diallel

# ------------------------------------------------------------------------------
# loop to generate posterior PSq tables for all traits
# store in new directory, 'PSqTables'
# ------------------------------------------------------------------------------
dir.create("~/Documents/carrot-diallel/PSqTables/")

for(i in varNames) {
  setwd("~/Documents/carrot-diallel/PSqTables/")
  psq <- cbind(results[[i]]$AllDiallelObs[[1]]$PosteriorPSqTable[[1]], 
               results[[i]]$AllDiallelObs[[1]]$PosteriorPSqTable[[2]])
  write.csv(psq, file = paste0(i, "psq.csv", sep = ""))
}

# ------------------------------------------------------------------------------
# function to combine posterior PSq tables
# ------------------------------------------------------------------------------
createPsqDF <- function(x) {
  dat <- read.csv(x, header = TRUE)
  id <- tools::file_path_sans_ext(basename(x))
  dat2 <- cbind(dat, id)
}

PSqTmp <- do.call("rbind", lapply(list.files(pattern = "*psq.csv",
                                             full = TRUE), createPsqDF))

# ------------------------------------------------------------------------------
# store output in .csv
# ------------------------------------------------------------------------------
write.csv(PSqTmp, "PSq_all.csv", row.names = FALSE)


