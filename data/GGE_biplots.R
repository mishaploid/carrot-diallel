# ------------------------------------------------------------------------------
#   Identification of best testers using GGE biplot method
#   S. Turner 
#   Ref: Yan and Hunt (2002); Frutos et al (2014)
#   19 February 2017        
# ------------------------------------------------------------------------------

library(reshape2)
library(GGEBiplotGUI)

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
poolingDf <- poolingDf[!poolingDf$.imp==0,]
poolingDf <- droplevels(poolingDf)

# Example: height
aggregate(height ~ female + male, data = poolingDf, FUN = mean)

# create a matrix of means with females as rows and males as columns
matrix1 <- acast(poolingDf, female ~ male, value.var = "height", 
                 fun.aggregate = mean, margins = TRUE, na.rm = TRUE)
matrix1
colnames(matrix1)[1:6] <- #c("a", "b", "c", "d", "e", "f")
  c("L6038", "L7550", "P0159", "Nbh2189", "P6139", "B7262")
rownames(matrix1)[1:6] <- #c("A", "B", "C", "D", "E", "F")
  c("L6038", "L7550", "P0159", "Nbh2189", "P6139", "B7262")

matrix2 <- matrix1[-7,-7] # remove row/col means
matrix2

GGEBiplot(Data=matrix2) # select symmetrical, tester centered, & no scaling
