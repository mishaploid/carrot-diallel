# ------------------------------------------------------------------------------
#   Griffing's analysis & combining abilities for imputed data
#   S. Turner 
#   19 February 2017        
# ------------------------------------------------------------------------------

# calculates and pools Griffing's ANOVA, general combining ability (GCA),
# specific combining ability (SCA), and reciprocal effects for multiply 
# imputed data

library(mice)
library(data.table)
library(Hmisc)
library(reshape2)
library(ggplot2)

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
# Griffing's ANOVA
# ------------------------------------------------------------------------------

# remove original data w/NAs first (.imp = 0)
poolingDf2 <- poolingDf[!poolingDf$.imp==0,]
poolingDf2 <- droplevels(poolingDf2)

# source modified dialleleI function (diallelI) from plantbreeding package
source("~/GitHub/carrot-diallel/05A_diallel_analysis_functions.R")

varNames <- c("midHeight", "height", "midWidth", "width", "flw", "dlw",
              "frw", "drw", "ratio")
GriffingResults <- list()
anvFull <- list()
pooled_ANOVA <- list()

# loop to run Griffing's analysis on m imputed data sets
# full output (see plantbreeding documentation) stored in GriffingResults list
# ANOVA output stored in anvFull list

# replace yvar (e.g. "height") with trait of interest
for(trait in varNames){
  for(i in levels(poolingDf2$.imp)) {
    df <- poolingDf2[which(poolingDf2$.imp==i),]
    GriffingResults[[trait]][[i]] <- diallelI(df, yvar = trait, progeny = "geno", 
                                              male = "male", female = "female", 
                                              replication = "rep", year = "year")
    anvFull[[trait]][[i]] <- rbind(GriffingResults[[trait]][[i]]$anvout[1,c(1,3,5)], 
                          GriffingResults[[trait]][[i]]$anova.mod1[1:3,c(1,3,5)],
                          GriffingResults[[trait]][[i]]$anvout[2:3,c(1,3,5)],
                          GriffingResults[[trait]][[i]]$anova.mod1[4:6,c(1,3,5)],
                          GriffingResults[[trait]][[i]]$anvout[4:5,c(1,3,5)])

# ------------------------------------------------------------------------------
# pool F-tests using rules outlined by Raghunathan & Dong (2011, unpublished)
# see http://www-personal.umich.edu/~teraghu/Raghunathan-Dong.pdf
# ------------------------------------------------------------------------------

# use anvFull output above to isolate mean squares for effect classes
MSgeno <- sapply(anvFull[[trait]], "[", 1, 2) #sN
MSGCA <- sapply(anvFull[[trait]], "[", 2, 2)
MSSCA <- sapply(anvFull[[trait]], "[", 3, 2)
MSrecip <- sapply(anvFull[[trait]], "[", 4, 2)
MSyear <- sapply(anvFull[[trait]], "[", 5, 2)
MSgxe <- sapply(anvFull[[trait]], "[", 6, 2)
MSgcaxe <- sapply(anvFull[[trait]], "[", 7, 2)
MSscaxe <- sapply(anvFull[[trait]], "[", 8, 2)
MSrecipxe <- sapply(anvFull[[trait]], "[", 9, 2)
MSrep <- sapply(anvFull[[trait]], "[", 10, 2)
MSE <- sapply(anvFull[[trait]], "[", 11, 2) #sD

# degrees of freedom 
#   p = # of parents
#   r = # of reps
#   y = # of years
p <-  6
r <- 2
y <- 3
Dfgeno <- p^2-1
DfGCA <- p-1
DfSCA <- p*(p-1)/2
Dfrecip <- p*(p-1)/2
Dfyear <- y-1
Dfgxe <- Dfgeno*Dfyear #70
Dfgcaxe <- DfGCA*Dfyear #10
Dfscaxe <- DfSCA*Dfyear #30
Dfrecipxe <- Dfrecip*Dfyear #30
Dfrep <- r
DfE <- anvFull[[trait]][[1]][11,1]

# function to pool F-statistic
# sN = numerator mean squares 
# vN = numerator degrees of freedom
# sD = denominator mean squares
# vD = denominator degrees of freedom
# M = # of imputations

ANVpool <- function(sN, vN, sD, vD, M, trait) {
  sNinv <- 1/sN
  sDinv <- 1/sD
  AN <- sum(sNinv)/M # multiple imputation posterior mean of sigma_N^-2, conditional on observed data
  AD <- sum(sDinv)/M
  F_MI <- AD/AN # compute F-statistic (ratio of harmonic means)
  
  BN <- sum(1/(vN*sN^2))/M
  BD <- sum(1/(vD*sD^2))/M
  
  CN <- sum((sNinv-AN)^2)/(M-1) # variance of reciprocals of mean squares
  CD <- sum((sDinv-AD)^2)/(M-1)
  
  rN <- 2*AN^2/(2*BN+(M+1)*CN/M)
  rD <- 2*AD^2/(2*BD+(M+1)*CD/M)
  pval <- pf(F_MI, rN, rD, lower.tail = FALSE)
  return(data.frame(MS = mean(sN), FVAL = F_MI, rN = rN, rD = rD, pvalue = pval))
}

# pool F-tests for each effect class
geno <- ANVpool(MSgeno, Dfgeno, MSE, DfE, 50) # genotype
gca <- ANVpool(MSGCA, DfGCA, MSE, DfE, 50) # GCA
sca <- ANVpool(MSSCA, DfSCA, MSE, DfE, 50) # SCA
recip <- ANVpool(MSrecip, Dfrecip, MSE, DfE, 50) # Reciprocal
year <- ANVpool(MSyear, Dfyear, MSE, DfE, 50) # Environment
gxe <- ANVpool(MSgxe, Dfgxe, MSE, DfE, 50) # GxE
gcaxe <- ANVpool(MSgcaxe, Dfgcaxe, MSE, DfE, 50) # GCAxE
scaxe <- ANVpool(MSscaxe, Dfscaxe, MSE, DfE, 50) # SCAxE
recipxe <- ANVpool(MSrecipxe, Dfrecipxe, MSE, DfE, 50) # ReciprocalxE
rep <- ANVpool(MSrep, Dfrep, MSE, DfE, 50) # rep%in%year
error <- c(mean(MSE), rep(NA, 4)) # error

# combine effects into data frame
pooled_ANOVA[[trait]] <- data.frame(trait = paste(trait),
                           source = c("genotype", "gca", "sca", "recip", "year", "gxe", 
                                      "gcaxe", "scaxe", "recipxe", "rep:year", "residual"),
                           rbind(geno, gca, sca, recip, year, gxe, gcaxe, scaxe, 
                                 recipxe, rep, error))

print(pooled_ANOVA[[trait]], row.names = FALSE)
  }
}


# ------------------------------------------------------------------------------
# pool general combining ability (GCA) estimates
# ------------------------------------------------------------------------------

parent1 <- poolingDf[,c(1:5, 7:17)] 
parent2 <- poolingDf[,c(1:6, 8:17)]
colnames(parent1)[6] <- "parent"
colnames(parent2)[6] <- "parent"

# convert new data frames to mids objects (input for pooling functions in mice)
p1 <- as.mids(parent1, .imp = 1, .id = 2)
p2 <- as.mids(parent2, .imp = 1, .id = 2)

gcaData <- rbind.mids(p1, p2) # combine male and female poolingDf into single factor

# use with and pool in mice to summarize GCA estimates 
# GCA estimates combined over all environments
# change y to reflect trait of interest (e.g. "height")
gca_i <- with(gcaData, lm(height - mean(height) ~ parent-1)) 
summary(gca_i)
gca_pooled <- summary(pool(gca_i))
gca_pooled

# use with and pool in mice to summarize GCA estimates 
# GCA estimates combined over all environments
# Example shown for canopy height 
# change y to reflect trait of interest (e.g. "height")
gca_i <- with(gcaData, lm(height - mean(height) ~ parent-1)) 

summary(gca_i)
gca_pooled <- summary(pool(gca_i))
gca_pooled

# GCA estimates for each environment (if significant GxE & GCAxE present)
gca_long <- complete(gcaData, action = "long", include = TRUE)

gca_wi2015 <- gca_long[gca_long$year == "WI2015",] # subset by location
gca_ca2015 <- gca_long[gca_long$year == "CA2015" & gca_long$.imp == 0,] # no missing data
gca_ca2016 <- gca_long[gca_long$year == "CA2016",]

gca_wi2015 <- as.mids(gca_wi2015, .imp = 1, .id = 2)
# gca_ca2015 <- as.mids(gca_ca2015, .imp = 1, .id = 2) # no missing data
gca_ca2016 <- as.mids(gca_ca2016, .imp = 1, .id = 2)

gcaxe_wi2015 <- with(gca_wi2015, lm(height - mean(height) ~ parent-1))
gcaxe_ca2015 <- with(gca_ca2015, lm(height - mean(height) ~ parent-1))
gcaxe_ca2016 <- with(gca_ca2016, lm(height - mean(height) ~ parent-1))

gca_wi2015_pooled <- summary(pool(gcaxe_wi2015))
gca_ca2015_pooled <- summary(gcaxe_ca2015)
gca_ca2016_pooled <- summary(pool(gcaxe_ca2016))

# simple output
bound <- cbind(gca_wi2015_pooled[,1], gca_ca2015_pooled$coefficients[,1], gca_ca2016_pooled[,1])

colnames(bound) <- c("WI2015", "CA2015", "CA2016")

rcorr(bound, type = "spearman") # rank correlation among environments

# output with SE
gca_wi2015_pooled <- as.data.frame(gca_wi2015_pooled[,c(1:3,4)])
gca_wi2015_pooled$parent <- row.names(gca_wi2015_pooled)
gca_ca2015_pooled <- as.data.frame(gca_ca2015_pooled[4])
colnames(gca_ca2015_pooled) <- colnames(gca_wi2015_pooled)
gca_ca2015_pooled$parent <- row.names(gca_ca2015_pooled)
gca_ca2016_pooled <- as.data.frame(gca_ca2016_pooled[,c(1:3,4)])
gca_ca2016_pooled$parent <- row.names(gca_ca2016_pooled)

gca_df <- combine(gca_wi2015_pooled, gca_ca2015_pooled, gca_ca2016_pooled)

# export data
write.csv(gca_df, "~/GitHub/carrot-diallel/results/height/gcabyyear_height.csv")

# plot GCAxE interactions
gca_df <- read.csv("~/GitHub/carrot-diallel/results/height/gcabyyear_height.csv", header = TRUE)

interaction.plot(gca_df$parent, gca_df$source, gca_df$est, las = 2)
abline(h = 0, col = "red")

# ------------------------------------------------------------------------------
# pretty plot of GCAxE interactions (Fig 5)
# ------------------------------------------------------------------------------
# recode and specify levels for variables
gca_df$parent <- revalue(gca_df$parent, c("parent1" = "L6038", "parent2" = "L7550", 
                                               "parent3" = "P0159", "parent4" = "Nbh2189", 
                                               "parent5" = "P6139", "parent6" = "B7262"))
gca_df$source <- revalue(gca_df$source, c("gca_wi2015_pooled" = "WI2015", 
                                          "gca_ca2015_pooled" = "CA2015", 
                                          "gca_ca2016_pooled" = "CA2016"))
gca_df$source <- factor(gca_df$source, levels = c("WI2015", "CA2015", "CA2016"))

gca_df$parent <- factor(gca_df$parent, levels = c("L6038", "L7550", "P0159", "Nbh2189", "P6139", "B7262"))

# plot GCA by year interactions with SE bars
p1 <- ggplot(gca_df, aes(x = source, y = est, colour = parent, group = parent)) + 
  geom_point(size = 2, aes(shape = parent)) +
  scale_shape_discrete(name = "Parental Line") + 
  geom_line(size = 0.75) +
  geom_errorbar(aes(ymin = est-se, ymax = est+se), width = 0.05, size = 0.5) + 
  ylab("GCA for Canopy Height (130 DAP)") +
  xlab("") +   
  ylim(-11, 10) + 
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                 "#0072B2", "#D55E00", "#CC79A7"),
                      name = "Parental Line") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank())
p1

ggsave("~/GitHub/carrot-diallel/results/gcaxe_height.eps", height = 80, width = 129, units = "mm")

# ------------------------------------------------------------------------------
# pool specific combining ability (SCA) estimates (e.g. midseason height)
# ------------------------------------------------------------------------------

# create factor for crosses (not ordered; jk = kj)
poolingDf$cross <- paste(pmin(as.numeric(poolingDf$female),
                              as.numeric(poolingDf$male)), 
                         pmax(as.numeric(poolingDf$female),
                              as.numeric(poolingDf$male)))
poolingDf$cross <- as.factor(poolingDf$cross)

# call calcSCA function 
source("~/GitHub/carrot-diallel/05B_calcSCA.R")

# loop to calculate SCA for m imputed data sets
# replace "height" with trait of interest
scaout <- list()
for(i in levels(poolingDf$.imp)){
  df <- poolingDf[which(poolingDf$.imp==i),]
  scaout[[i]] <- calcSCA("midHeight", ".id", "year", "female", "male", "cross", df)
}

# combine outputs for each dataset into single dataframe
scaDF <- rbindlist(scaout) 

# convert data to mids object for pooling in mice
scaMids <- as.mids(scaDF, .imp = 2, .id = 1) 

# SCA estimates combined over all environments
# change y to reflect trait of interest (e.g. "height")
lmsca_ij <- with(scaMids, lm(s_ij ~ cross-1))

sca_pooled <- summary(pool(lmsca_ij))

# SCA estimates for each environment (if significant GxE & SCAxE present)
lmsca_ij_E <- with(scaMids, lm(s_ij ~ cross:year-1))

scaxE_pooled <- summary(pool(lmsca_ij_E))

# export data
write.csv(scaxE_pooled, "~/GitHub/carrot-diallel/results/midHeight/scabyyear_midHeight.csv")

sca_df <- read.csv("~/GitHub/carrot-diallel/results/midHeight/scabyyear_midHeight.csv", header = TRUE)

sca_df <- data.frame(colsplit(sca_df$X, ":", names = c("cross", "year")), 
                     sca = sca_df$est)

# plot SCAxE interactions
interaction.plot(sca_df$cross, sca_df$year, sca_df$sca, las = 2)
abline(h = 0, col = "red")

# ------------------------------------------------------------------------------
# pool reciprocal effect estimates (e.g. height)
# ------------------------------------------------------------------------------

# factor with ordered cross (jk != kj)
poolingDf$orderedCross <- poolingDf$female:poolingDf$male
# factor with reciprocal cross order
poolingDf$reciprocalCross <- poolingDf$male:poolingDf$female 

# call calcRecip function
source("~/GitHub/carrot-diallel/05C_calcRecip.R")

# loop to calculate reciprocal cross differences
# creates new variable with reciprocal cross values
# replace "height" with trait of interest
recipout <- list()
for(i in levels(poolingDf$.imp)){
  df <- poolingDf[which(poolingDf$.imp==i),]
  recipout[[i]] <- calcRecip("height", ".id", "year", "orderedCross", 
                             "reciprocalCross", df)
}

recipDF <- rbindlist(recipout)
recipDF <- as.data.frame(recipDF)

# calculate reciprocal effects
# replace "height" with trait of interest
recipDF$reff <- 1/2*(recipDF$height - recipDF$r_ji) #change variable name

recipMids <- as.mids(recipDF, .imp = 2, .id = 1)

# Reciprocal effect estimates combined over all environments
# change y to reflect trait of interest (e.g. "height")
lmr_ij <- with(recipMids, lm(reff ~ orderedCross-1))

recipEff_pooled <- summary(pool(lmr_ij)) # note: parental selfs should be ~ 0
recipEff_pooled

# Reciprocal effect estimates for each environment 
# (if significant GxE & ReciprocalxE present)
lmr_ij_E <- with(recipMids, lm(reff ~ orderedCross:year-1))

recipbyE_pooled <- summary(pool(lmr_ij_E))

# export data
write.csv(recipbyE_pooled, "~/GitHub/carrot-diallel/results/height/recipbyyear_height.csv")

recip_df <- read.csv("~/GitHub/carrot-diallel/results/height/recipbyyear_height.csv", header = TRUE)

recip_df <- data.frame(colsplit(recip_df$X, ":year", names = c("cross", "year")), 
                       recipEff = recip_df$est)

# plot recipxE interactions
interaction.plot(recip_df$cross, recip_df$year, recip_df$recipEff, las = 2)
abline(h = 0, col = "red")

# ------------------------------------------------------------------------------
# pool ratio of GCA to SCA mean squares using Griffing's Method 3
# ------------------------------------------------------------------------------
source("~/GitHub/carrot-diallel/05D_pool_Griffing_Method3.R", print.eval=TRUE)
