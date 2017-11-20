# ------------------------------------------------------------------------------
#   Griffing's analysis - method III (w/o parents) to estimate Baker's ratio
#   S. Turner 
#   30 October 2017       
# ------------------------------------------------------------------------------

# calculates and pools Griffing's ANOVA, general combining ability (GCA),
# specific combining ability (SCA), and reciprocal effects for multiply 
# imputed data

library(mice)

setwd("~/GitHub/carrot-diallel/data")

poolingDf <- read.csv("70_iter_full_data.csv", header = TRUE)

# set factor levels
poolingDf$.imp <- as.factor(poolingDf$.imp)
poolingDf$.id <- as.factor(poolingDf$.id)
poolingDf$geno <- as.factor(poolingDf$geno)
poolingDf$rep <- as.factor(poolingDf$rep)
poolingDf$male <- as.factor(poolingDf$male)
poolingDf$female <- as.factor(poolingDf$female)
poolingDf$ratio <- poolingDf$dlw/poolingDf$drw

# ------------------------------------------------------------------------------
# Griffing's ANOVA (Method III, Model I)
# ------------------------------------------------------------------------------

# remove original data w/NAs first (.imp = 0)
poolingDf2 <- poolingDf[!poolingDf$.imp==0,]
poolingDf2 <- droplevels(poolingDf2)

p <- nlevels(poolingDf2$male)
r <- nlevels(poolingDf2$rep)
y <- nlevels(poolingDf2$year)
traits <- colnames(poolingDf2[,9:17])

Baker.ratio <- list()

for(trait in traits){
  for(imp in levels(poolingDf2$.imp)){
    df <- poolingDf2[which(poolingDf2$.imp==imp),]
    form1 <- as.formula(paste(trait, " ~ male:female", sep=""))
    data1 <- data.frame(aggregate(form1, data = df, mean))
    form2 <- as.formula(paste(trait, " ~ female", sep=""))
    mydf <- aggregate(form2, data1, "c")
    myMatrix <- as.matrix(mydf[, -1])
    diag(myMatrix) <- 0
    Xi. <- rowSums(myMatrix)
    X.i <- colSums(myMatrix)
    X.. <- sum(myMatrix)
    acon <- sum((Xi.+X.i)^2)/(2*(p-2))
    SSgca <- (acon - (2/(p*(p-2))) * (X..^2))
    SSsca <- (sum(((myMatrix + t(myMatrix))^2)/2)/2) - acon + 
               (X..)^2/((p -1) * (p - 2))
    MSgca <- SSgca/(p-1)
    MSsca <- SSsca/(p*(p-3)/2)
    GCAtoSCA <- (2*MSgca)/(2*MSgca + MSsca)
    Baker.ratio[[trait]][[imp]] <- GCAtoSCA
  }
}

lapply(Baker.ratio, FUN=function(x) round(mean(x),2))
