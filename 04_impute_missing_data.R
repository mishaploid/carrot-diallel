# ------------------------------------------------------------------------------
#   Missing data imputation for carrot diallel 
#   S. Turner 
#   29 August 2016        
# ------------------------------------------------------------------------------
library(mice)
library(lattice)
library(VIM)
library(plyr)
library(gridExtra)

setwd("~/GitHub/carrot-diallel/data")
diallelRaw <- read.csv("diallelRawData.csv", header = TRUE)

str(diallelRaw)

# subset with traits and predictors of interest
diallelRaw <- diallelRaw[, c(3:16)]

# apply transformations
diallelRaw$flw <- log(diallelRaw$flw)
diallelRaw$dlw <- log(diallelRaw$dlw)
diallelRaw$frw <- log(diallelRaw$frw)
diallelRaw$drw <- log(diallelRaw$drw)

# specify factors
diallelRaw$geno <- as.factor(diallelRaw$geno)
diallelRaw$year <- as.factor(diallelRaw$year)
diallelRaw$rep <- as.factor(diallelRaw$rep)
diallelRaw$male <- as.factor(diallelRaw$male)
diallelRaw$female <- as.factor(diallelRaw$female)
diallelRaw$stand <- as.factor(diallelRaw$stand)

# ------------------------------------------------------------------------------
# visualize missing data (VIM package) 
# ------------------------------------------------------------------------------
# % missing
pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(diallelRaw, 2, pMiss) # apply pMiss function over columns

# look for patterns in missing data
md.pattern(diallelRaw) 

# visualize missing data - Fig S2
aggr_plot <- aggr(diallelRaw, col = c('navyblue', 'red'), numbers = TRUE, 
                  labels = names(diallelRaw),
                  cex.axis = .7, gap = 3, 
                  ylab = c("Histogram of missing data", "Pattern"))

# check MCAR (missing completely at random) assumption
# distribution of x when y is missing shown in red
marginplot(diallelRaw[c(8,12)]) # dist of midHeight when flw is missing
marginplot(diallelRaw[c(8,14)]) # dist of midHeight when frw is missing

# ------------------------------------------------------------------------------
# impute missing data for diallel data
# Multivariate Imputation by Chained Equations (MICE)
#   m = # of imputations
#   meth = "pmm" = predictive mean matching 
#   maxit = number of iterations (default is 5)
#   default predictor is all other columns in data
# ------------------------------------------------------------------------------

# 50 chains w/ 70 iterations - full data
tempdatafull <- mice(data = diallelRaw, m = 50, method = "pmm", 
                     maxit = 70, seed = 1234)

# ------------------------------------------------------------------------------
# save imputed dataset
# ------------------------------------------------------------------------------
# recode year to represent the locations
tempdatafull$data$year <- revalue(tempdatafull$data$year, 
                                  c("1" = "CA2015", "2" = "WI2015", "3" = "CA2016"))

imputed_data <- mice::complete(tempdatafull, "long", include = TRUE)

write.csv(imputed_data, 
          file = "~/GitHub/carrot-diallel/data/70_iter_full_data.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# plot convergence of MCMC chains - Fig S4
# ------------------------------------------------------------------------------
dir.create("~/GitHub/carrot-diallel/results/imputation/")
png("~/GitHub/carrot-diallel/results/imputation/convergence-50imp-70iter-full.png", 
    height = 9.5, width = 7, units = "in", res = 600)
plot(tempdatafull, y = c("midHeight", "midWidth", "height", "width", 
                         "flw", "dlw", "frw", "drw"), 
     header = "Convergence (70 Iterations)", layout = c(2, 8), lwd = 0.8)
dev.off()

# ------------------------------------------------------------------------------
# export observed vs imputed (combined over years) - Fig S5
# ------------------------------------------------------------------------------
png("~/GitHub/carrot-diallel/results/imputation/50imp-70iter-impvsobs.png", 
    height = 8, width = 6, units = "in", res = 600)
densityplot(tempdatafull, main = "Imputed (pink) and Observed (blue) distributions",
            lwd = 1, layout = c(2,4))
dev.off()

# ------------------------------------------------------------------------------
# plot distribution of imputed (pink) and observed (blue) data by year
# ------------------------------------------------------------------------------
traits <- colnames(tempdatafull$data[,7:14])
pltList <- list()
for (i in traits) {
  pltName <- paste("p", i, sep = "")
  form <- formula(paste("~", i, "|", "year"))
  pltList[[pltName]] <- densityplot(tempdatafull, form, groups = .imp, 
                                    col = c("#006CC2CC", rep("#B61A51CC", 50)), 
                                    lwd = c(2, rep(1, 50)))
  print(pltList[[pltName]])
}

# ------------------------------------------------------------------------------
# export obs vs imputed (by year) - Fig S6
# selected traits where observed/imputed distributions above did not match
# shoot and root biomass (fresh and dry)
# ------------------------------------------------------------------------------
png("~/GitHub/carrot-diallel/results/imputation/impvsobs-year.png", 
    height = 8, width = 8, units = "in", res = 600)
grid.arrange(pltList[[5]], pltList[[6]], pltList[[7]], pltList[[8]])
dev.off()

# ------------------------------------------------------------------------------
# OTHER OPTIONS FOR IMPUTATION
# Run for individual years based on degree of missing data
# Run 500 iterations
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# impute year 3 from year 1, 50 imputations & 50 iterations

# diallelY1Y3 <- subset(diallelRaw, year == "1" | year == "3")
# 
# tempdata2 <- mice(data = diallelY1Y3, meth = "pmm", m = 50, maxit = 50, seed = 1234)
# 
# pdf("~/GitHub/carrot-diallel/results/imputation/Y1Y3_50_iter.pdf", height = 8, width = 12)
# plot(tempdata2, main = "Convergence for Y1 and Y3 - 50 Iterations", lwd = 2)
# densityplot(tempdata2, main = "Imputed (pink) and Observed (blue) distributions",
#             lwd = 3)
# stripplot(tempdata2, pch = 20, cex = 1.2, 
#           main = "Imputed (pink) and Observed (blue) distributions")
# dev.off()
# 
# write.csv(complete(tempdata2, "long", include = TRUE), 
#           file = "~/GitHub/carrot-dialle/results/imputation/Y1Y3_50_iter.csv",
#           row.names = FALSE)
# 
# # ------------------------------------------------------------------------------
# # impute year 2 from Y1Y3 imputed data - 50 imputations & 50 iterations
# 
# Y2 <- subset(diallelRaw, year == "2")
# 
# diallelY2 <- rbind(diallelY1Y3, Y2)
# 
# tempdata3 <- mice(data = diallelY2, meth = "pmm", m = 50, maxit = 50, seed = 1234)
# 
# pdf("~/GitHub/carrot-diallel/results/imputation/Y1Y3-Y2_50_iter.pdf", height = 8, width = 12)
# plot(tempdata3, 
#      main = "Convergence for Y1Y3 (complete) and Y2 (imputed) - 50 iterations",
#      lwd = 2)
# densityplot(tempdata3, main = "Imputed (pink) and Observed (blue) distributions",
#             lwd = 3)
# stripplot(tempdata3, pch = 20, cex = 1.2, 
#           main = "Imputed (pink) and Observed (blue) distributions")
# dev.off()
# 
# write.csv(complete(tempdata3, "long", include = TRUE), 
#           file = "~/GitHub/carrot-diallel/results/imputation/Y1Y3-Y2_50_iter.csv",
#           row.names = FALSE)
# 
# densityplot(completedata, main = "Imputed (pink) and Observed (blue) distributions",
#             lwd = 3)

# # ------------------------------------------------------------------------------
# # increase to 500 iterations
# tempdata <- mice(data = diallelRaw, m = 15, method = "pmm", maxit = 500, seed = 1234)
# 
# summary(tempdata)
# 
# pdf("~/GitHub/carrot-diallel/results/imputation/convergence_500_iter.pdf", height = 12, width = 20)
# plot(tempdata, header = "Convergence (500 Iterations)", lwd = 2)
# densityplot(tempdata, main = "Imputed (pink) and Observed (blue) distributions",
#             lwd = 3)
# stripplot(tempdata, pch = 20, cex = 1.2,
#           main = "Imputed (pink) and Observed (blue) distributions")
# dev.off()
