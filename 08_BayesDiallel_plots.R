# ------------------------------------------------------------------------------
#   BayesDiallel plots - see BayesDiallel vignette for specifics
#   S. Turner 
#   18 July 2016
# ------------------------------------------------------------------------------

library(psych)
library(ggplot2)

# Code to plot observed vs. expected values (TwoDiallelPlot), highest
# posterior density (HPD) intervals (PlotHPD), and strawplots for AFD objects

# ------------------------------------------------------------------------------
# load AFD objects
# ------------------------------------------------------------------------------
setwd("~/GitHub/carrot-diallel")
source("~/GitHub/carrot-diallel/07_read_AFD_objects.R")

# ------------------------------------------------------------------------------
# Observed vs. expected plots - expected shows predictions based on observed data
# ------------------------------------------------------------------------------

# e.g. height
results[["height"]]$AllDiallelObs[[1]]$TwoDiallelPlot(PlotObservedVersusExpectedAll = TRUE,
                                                      HeadTitle = "",
                                                      show.strain.names = TRUE,
                                                      LabelLeft = TRUE)

# loop to generate observed vs expected TwoDiallelPlots for all traits
for (i in varNames) {
  setEPS()
  postscript(paste0("~/GitHub/carrot-diallel/results/", i, "/OvsE", i, ".eps", sep = ""), width = 6, 
      height = 4)
  results[[i]]$AllDiallelObs[[1]]$TwoDiallelPlot(PlotObservedVersusExpectedAll = TRUE, 
                                                 HeadTitle = "", 
                                                 show.strain.names = TRUE,
                                                 LabelLeft = TRUE)
  dev.off()
}

# ------------------------------------------------------------------------------
# highest posterior density (HPD) plots for inheritance classes
# ------------------------------------------------------------------------------

# e.g. height
results[["height"]]$AllDiallelObs[[1]]$PlotHPD(UseDefaultWanted1 = TRUE,
                                               EvenXlims = TRUE, DoMu = FALSE,
                                               xlim = c(-10, 10), main = "height")

# loop to generate HPD plots for all traits
for (i in varNames) {
  setEPS()
  postscript(paste0("~/GitHub/carrot-diallel/results/", i, "HPD_main", i, ".eps", sep = ""), 
             width = 10, height = 6)
  results[[i]]$AllDiallelObs[[1]]$PlotHPD(UseDefaultWanted1 = TRUE, 
                                          EvenXlims = TRUE, DoMu = FALSE,
                                          xlim = c(-10, 10), main = paste(i))
  dev.off()
}

# ------------------------------------------------------------------------------
# HPD plots for random effects
# ------------------------------------------------------------------------------

# loop to generate HPD plots of random effects for all traits
for (i in varNames) {
  postscript(paste0("~/GitHub/carrot-diallel/results/", i, 
                    "/HPD_random", i, ".eps", sep = ""), width = 5, height = 5)
  results[[i]]$AllDiallelObs[[1]]$PlotHPD(wanted = c(5,4,6,3,7:12), 
                                               EvenXlims = TRUE, DoMu = FALSE,
                                               xlim = c(-10, 10), main = paste(i))
  dev.off()
}

# e.g. height & midHeight (for Fig 6)
setEPS()
postscript("~/GitHub/carrot-diallel/results/midheight_random_eff.eps", width = 6, height = 6)
#layout(matrix(c(1,2), 2, 2, byrow=TRUE), widths=c(1.8,1))
results[["midHeight"]]$AllDiallelObs[[1]]$PlotHPD(wanted = c(5,4,6,3,7:12), 
                                                  EvenXlims = TRUE, DoMu = FALSE,
                                                  xlim = c(-15, 15), main = "midHeight")
dev.off()
setEPS()
postscript("~/GitHub/carrot-diallel/results/height_random_eff.eps", width = 3.465, height = 6)
results[["height"]]$AllDiallelObs[[1]]$PlotHPD(wanted = c(5,4,6,3,7:12), 
                                               EvenXlims = TRUE, DoMu = FALSE,
                                               xlim = c(-15, 15), main = "height")
dev.off()

# ------------------------------------------------------------------------------
# straw plots (interactions among inheritance classes)
# see Lenarcic et al. (2012) Genetics and BayesDiallel documentation for details
# ------------------------------------------------------------------------------

# e.g. height
PlotStrawPlot(results[["height"]]$AllDiallelObs[[1]], DoMedian = FALSE, DoMu = FALSE,
              GiveHeadings = TRUE, yline = 2, DoCC = FALSE, lwd = 4,
              col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", 
                      "#D55E00"), ylab = "height")
legend("topleft", c("A", "B", "C", "D", "E", "F"), lty = 1, lwd = 4, 
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"),
       title = "parent")

# loop to generate straw plots for all traits
for (i in varNames) {
  PlotStrawPlot(results[[i]]$AllDiallelObs[[1]], DoMedian = FALSE, DoMu = FALSE,
                GiveHeadings = TRUE, yline = 2, DoCC = FALSE, lwd = 4,
                col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", 
                        "#D55E00"), ylab = paste(i))
  legend("topleft", c("A", "B", "C", "D", "E", "F"), lty = 1, lwd = 4, 
         col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"),
         title = "parent")
}

# ------------------------------------------------------------------------------
# plot correlations for effects
# Fig 2 B and C

summary(results[[1]]$AllDiallelObs[[1]]$cent.chains)

post_means <- list()

for(i in varNames){
  AO <- results[[i]]$AllDiallelObs[[1]]
  means <- summary(AO$cent.chains)
  post_means[[i]] <- means[[1]][13:30,1]
}

post_df <- as.data.frame(post_means)
post_df <- post_df[,c(1:9)]

# rename variables
post_df <- rename(post_df, c("midHeight" = "height (80DAP)", 
                                   "midWidth" = "width (80DAP)",
                                   "height" = "height (130DAP)", 
                                   "width" = "width (130DAP)",
                                   "flw" = "shoot biomass (fresh)",
                                   "dlw" = "shoot biomass (dry)",
                                   "frw" = "root biomass (fresh)",
                                   "drw" = "root biomass (dry)", 
                                   "ratio" = "shoot:root ratio"))

a_cor <- cor(post_df[1:6,]) # additive effect correlations
i_cor <- cor(post_df[13:18,]) # inbred deviation correlations

corrplot(i_cor)

# ------------------------------------------------------------------------------
# create correlation matrix for additive and inbred parameters
# ------------------------------------------------------------------------------
correlation_a <- cor(post_df[1:6,], method = "pearson")
correlation_i <- cor(post_df[13:18,], method = "pearson")

# ------------------------------------------------------------------------------
# function to calculate significance of correlations
# source: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# ------------------------------------------------------------------------------
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

# ------------------------------------------------------------------------------
# export matrix of p-values
# ------------------------------------------------------------------------------
p.mat_a <- cor.mtest(post_df[1:6,])
p.mat_i <- cor.mtest(post_df[13:18,])

# ------------------------------------------------------------------------------
# create matrix for significance notations 
# * = P < 0.05, ** = P < 0.01, *** = P < 0.001
# ------------------------------------------------------------------------------
sigCodes_a <- ifelse(p.mat_a[[1]] > 0.05, "NS", 
                   ifelse(p.mat_a[[1]] > 0.01, "*", 
                          ifelse(p.mat_a[[1]] > 0.001, "**", "***")))

sigCodes_i <- ifelse(p.mat_i[[1]] > 0.05, "NS", 
                     ifelse(p.mat_i[[1]] > 0.01, "*", 
                            ifelse(p.mat_i[[1]] > 0.001, "**", "***")))

# set diagonal to NA
sigCodes_a[upper.tri(sigCodes_a)] <- NA
sigCodes_i[upper.tri(sigCodes_i)] <- NA

# ------------------------------------------------------------------------------
# combine correlation coefficients and significance codes into single matrix
# ------------------------------------------------------------------------------
combinedMat_a <- lowerUpper(upper = sigCodes_a, lower = round(correlation_a, 2), 
                          diff = FALSE)
diag(combinedMat_a) <- ""

combinedMat_i <- lowerUpper(upper = sigCodes_i, lower = round(correlation_i, 2), 
                            diff = FALSE)
diag(combinedMat_i) <- ""

# ------------------------------------------------------------------------------
# construct the plot! 
# ------------------------------------------------------------------------------
# create labels as input for geom_text
labels_a <- melt(combinedMat_a, na.rm = TRUE)
labels_i <- melt(combinedMat_i, na.rm = TRUE)

# melt matrix into vector format for plotting
meltedCormat_a <- melt(correlation_a)
meltedCormat_i <- melt(correlation_i)

ggplot(data = meltedCormat_a, aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#ca0020", high = "#2166AC", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Pearson\nCorrelation") + 
  scale_y_discrete(name = "", limits = rev(levels(meltedCormat_a$Var1))) +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) + 
  geom_text(aes(Var2, Var1, label = labels_a$value), colour = "black", size = 3.15) + 
  coord_fixed()

ggsave("~/GitHub/results/additive_corr.eps", height = 150, width = 150, units = "mm")

ggplot(data = meltedCormat_i, aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#ca0020", high = "#2166AC", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Pearson\nCorrelation") + 
  scale_y_discrete(name = "", limits = rev(levels(meltedCormat_i$Var1))) +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) + 
  geom_text(aes(Var2, Var1, label = labels_i$value), colour = "black", size = 3.15) + 
  coord_fixed()

ggsave("~/GitHub/results/inbred_corr.eps", height = 150, width = 150, units = "mm")

