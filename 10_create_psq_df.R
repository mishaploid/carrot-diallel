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
source("~/GitHub/carrot-diallel/07_read_AFD_objects.R")

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

for(i in varNames) {
  setwd(paste0("~/GitHub/carrot-diallel/results/", i, sep = ""))
  psq <- cbind(results[[i]]$AllDiallelObs[[1]]$PosteriorPSqTable[[1]], 
               results[[i]]$AllDiallelObs[[1]]$PosteriorPSqTable[[2]])
  print(psq)
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

psq_files <- dir("~/GitHub/carrot-diallel/results", recursive=TRUE, 
                 full.names=TRUE, pattern="\\psq.csv$")

PSqTmp <- do.call("rbind", lapply(psq_files, createPsqDF))

# ------------------------------------------------------------------------------
# store output in .csv
# ------------------------------------------------------------------------------
write.csv(PSqTmp, "~/GitHub/carrot-diallel/results/PSq_all.csv", row.names = FALSE)


