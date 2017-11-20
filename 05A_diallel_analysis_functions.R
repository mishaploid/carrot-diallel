# ------------------------------------------------------------------------------
#   Edits to plantbreeding::diallele1 function to include multiple environments
#   Source: 'plantbreeding' package - Based on Griffing 1956
#   Note: corrected SCA matrix courtesy of G. Ramstein
#   S. Turner 
#   13 December 2016
# ------------------------------------------------------------------------------

# install.packages("plantbreeding", repos="http://R-Forge.R-project.org")
library(plantbreeding)

# ------------------------------------------------------------------------------
# Full diallel (Method I) - modification only for model I (fixed)
# based on function diallele1 in package 'plantbreeding' but corrects a mistake
# on calculation of SCA's (in diallele1, sca are calculated in such a way that
# the resulting SCA matrix is not symmetric) type ?diallele1 for information
# about syntax

diallelI <- function (dataframe, yvar = "yvar", progeny = "progeny", male = "male", 
                      female = "female", replication = "replication", year = "year") {
  dataframe <- data.frame(yvar = dataframe[, yvar],
                          progeny = dataframe[, progeny],
                          male = dataframe[, male],
                          female = dataframe[, female],
                          replication = dataframe[, replication],
                          year = dataframe[, year])
  dataframe$male = as.factor(dataframe$male)
  dataframe$female = as.factor(dataframe$female)
  dataframe$progeny = as.factor(dataframe$progeny)
  dataframe$replication = as.factor(dataframe$replication)
  dataframe$year <- as.factor(dataframe$year)
  mean.y = mean(as.numeric(dataframe$yvar))
  
  formula = paste("yvar", " ~ ", "progeny", "*", "year", " + ", "replication", 
                  "%in%", "year")
  md1 <- lm(formula, data = dataframe)
  yvar = yvar
  cat("Diallel analysis for trait: ", yvar, "\n\n")
  cat("...........................\n\n")
  anvout <- anova(md1)
  print(anvout)
  if (anvout[1, 5] > 0.05) {
    cat("Note: Progeny effect is not significant at 0.05 p threshold")
  }
  if (anvout[1, 5] > 0.01) {
    cat("Note: Progeny effect is not significant at 0.01 p threshold")
  }
  data1 <- data.frame(aggregate(yvar ~ male:female, data = dataframe, 
                                mean))
  mydf <- aggregate(yvar ~ female, data1, "c")
  myMatrix <- as.matrix(mydf[, -1])
  n = (length(myMatrix))^0.5
  r = nlevels(dataframe$replication)
  y = nlevels(dataframe$year)
# ------------------------------------------------------------------------------
  ## EDIT 1 - modify main effect sum of squares
  # multiply ssgca, sssca, and ssrecp by number of environments and replications
# ------------------------------------------------------------------------------
  acon <- sum((1/(2 * n)) * ((rowSums(myMatrix) + colSums(myMatrix))^2))
  ssgca = r * y * (acon - (2/(n^2)) * (sum(myMatrix)^2))
  sssca = r * y * (sum((1/2) * (myMatrix * (myMatrix + t(myMatrix)))) - 
                     acon + (1/(n^2)) * (sum(myMatrix)^2))
  ssrecp = r * y * (((1/4) * sum((myMatrix - t(myMatrix))^2)))
# ------------------------------------------------------------------------------ 
  ## EDIT 2 - partition sum of squares by environment
# ------------------------------------------------------------------------------
  yearlist <- split(dataframe, dataframe$year, drop = TRUE)
  # function to partition sum of squares for gca, sca, and reciprocal effects by environment
  ssEnv <- function(x) {
    tmp <- aggregate(yvar ~ male:female, data = x, mean)
    mydf <- aggregate(yvar ~ female, tmp, "c")
    myMatrixE <- as.matrix(mydf[, -1])
    n = (length(myMatrixE))^0.5
    acon <- (1/(2 * n)) * sum((rowSums(myMatrixE) + colSums(myMatrixE))^2)
    ssgcaE <- acon - (2/(n^2)) * (sum(myMatrixE)^2)
    ssscaE <- sum((1/2) * (myMatrixE * (myMatrixE + t(myMatrixE)))) -
      acon + (1/(n^2)) * (sum(myMatrixE)^2)
    ssrecpE <- ((1/4) * sum((myMatrixE - t(myMatrixE))^2))
    list(ssgcaE = ssgcaE, ssscaE = ssscaE, ssrecpE = ssrecpE)
  }
  ssqE <- lapply(yearlist, ssEnv)
  ssqEdf <- data.frame(do.call(rbind.data.frame, ssqE))
  ssgcaE <- r * sum(ssqEdf[,1]) - ssgca
  ssscaE <- r * sum(ssqEdf[,2]) - sssca
  ssrecpE <- r * sum(ssqEdf[,3]) - ssrecp
# ------------------------------------------------------------------------------  
  ## EDIT 3 - construct partitioned ANOVA table
  # original code: MSEAD = anvout$`Mean Sq`[3]/r 
# ------------------------------------------------------------------------------
  MSEAD = anvout$`Mean Sq`[5]/r*y
  p = n
  Df <- c((p - 1), (p * (p - 1)/2), (p * (p - 1)/2), (p - 1) * (y - 1),
          (p * (p - 1)/2) * (y - 1), (p * (p - 1)/2) * (y - 1), anvout$Df[5])
  SSS <- c(ssgca, sssca, ssrecp, ssgcaE, ssscaE, ssrecpE)
  # original code:
  # ssq <- c(SSS, anvout$`Sum Sq`[3]/r)
  ssq <- c(SSS, anvout$`Sum Sq`[5])
  MSSS = SSS/Df[1:6]
  MSSS1 = c(MSSS, MSEAD)
  FVAL = c(MSSS1[1:6]/MSEAD, NA)
  pval = c(1 - pf(FVAL[1:6], Df[1:6], Df[7]), NA)
  anovadf.mod1 <- data.frame(Df, `Sum Sq` = ssq, `Mean Sq` = MSSS1, 
                             `F value` = FVAL, `Pr(>F)` = pval, check.names = FALSE)
  rownames(anovadf.mod1) <- c("GCA", "SCA", "Reciprocal", "GCAxE", "SCAxE", 
                              "ReciprocalxE", "error")
  class(anovadf.mod1) <- c("anova", "data.frame")
  
  cat("Anova for combining ability - ModelI", yvar, "\n\n")
  print(anovadf.mod1)
  
  GCAcomp = (MSSS[1] - MSEAD)/(2 * n)
  SCAcomp = (MSSS[2] - MSEAD)
  RCAcomp = (MSSS[3] - MSEAD)/2
  GCASCAratio = GCAcomp/SCAcomp
  components.model1 <- list(GCAcomp = GCAcomp, SCAcomp = SCAcomp, 
                            RCAcomp = RCAcomp, GCASCAratio = GCASCAratio)
  cat("Componets: Model 1", "\n\n")
  cat("GCA :", GCAcomp, "\n\n")
  cat("SCA :", SCAcomp, "\n\n")
  cat("Reciprocal:", RCAcomp, "\n\n")
  cat("GCA to SCA ratio: ", GCASCAratio, "\n\n")
  cat("Note this ratio is biased if parents are included! Consider running Model III")
  
  gcaeff <- ((1/(2 * n)) * (rowSums(myMatrix) + colSums(myMatrix))) - 
    ((1/(n^2)) * (sum(myMatrix)))
  
  rs <- rowSums(myMatrix)
  cs <- colSums(myMatrix)
  M <- matrix(numeric(), nrow=nrow(myMatrix), ncol=ncol(myMatrix))
  
  for (i in 1:nrow(myMatrix)) {
    for (j in 1:ncol(myMatrix)) {
      M[i,j] <- rs[i] + cs[i] + rs[j] + cs[j]
    }
  }
# ------------------------------------------------------------------------------
  # G. Ramstein edit for SCA matrix
  # replaces original code: scaeff <- ((1/2) * (myMatrix + t(myMatrix))) - 
  # ((1/(2 * n)) * (rowSums(myMatrix) + colSums(myMatrix) + colSums(t(myMatrix))
  # + rowSums(t(myMatrix)))) + ((1/(n^2)) * (sum(myMatrix)))
# ------------------------------------------------------------------------------
  scaeff <- ((1/2) * (myMatrix + t(myMatrix))) -
    ((1/(2 * n)) * M) +
    ((1/(n^2)) * (sum(myMatrix)))
  
  recieff <- 0.5 * (myMatrix - t(myMatrix)) 
  effmat <- gcaeff
  
  var.gi = ((n - 1)/(2 * n^2)) * MSEAD
  var.sii = (((n - 1)^2)/n^2) * MSEAD
  var.sij = (1/(2 * n^2)) * ((n^2) - 2 * n + 2) * MSEAD
  var.rij = 0.5 * MSEAD
  var.gi_gj = (1/n) * MSEAD
  var.sij_sji = ((2 * (n - 2))/n) * MSEAD
  var.sii_sij = ((3 * n - 2)/(2 * n)) * MSEAD
  var.sii_sjk = ((3 * (n - 2))/(2 * n)) * MSEAD
  var.sij_sik = ((n - 1)/n) * MSEAD
  var.sij_skl = ((n - 2)/n) * MSEAD
  var.rij_rkl = MSEAD
  varcompare <- list(var.gi = var.gi, var.sii = var.sii, var.sij = var.sij, 
                     var.rij = var.rij, var.gi_gj = var.gi_gj, var.sij_sji = var.sij_sji, 
                     var.sii_sij = var.sii_sij, var.sii_sjk = var.sii_sjk, 
                     var.sij_sik = var.sij_sik, var.sij_skl = var.sij_skl, 
                     var.rij_rkl = var.rij_rkl)
  return(list(anvout = anvout, anova.mod1 = anovadf.mod1, components.model1 = components.model1, 
              gca.effmat = gcaeff, sca.effmat = scaeff, reciprocal.effmat = recieff, 
              varcompare = varcompare))
}

