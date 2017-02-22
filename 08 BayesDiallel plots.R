# ------------------------------------------------------------------------------
#   BayesDiallel plots - see BayesDiallel vignette for specifics
#   S. Turner 
#   18 July 2016
# ------------------------------------------------------------------------------

# This script plots observed vs. expected values (TwoDiallelPlot), highest
# posterior density (HPD) intervals (PlotHPD), and strawplots for AFD objects

# ------------------------------------------------------------------------------
# load AFD objects
# ------------------------------------------------------------------------------
setwd("~/Documents/carrot-diallel/")
source("~/Documents/carrot-diallel/09 read in AFD objects.R")

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
  png(paste0("~/Documents/carrot-diallel/OvsE", i, ".png", sep = ""), width = 6, 
      height = 4, units = "in", res = 300)
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
  results[[i]]$AllDiallelObs[[1]]$PlotHPD(UseDefaultWanted1 = TRUE, 
                                          EvenXlims = TRUE, DoMu = FALSE,
                                          xlim = c(-10, 10), main = paste(i))
}

# ------------------------------------------------------------------------------
# HPD plots for random effects
# ------------------------------------------------------------------------------

# e.g. height
results[["height"]]$AllDiallelObs[[1]]$PlotHPD(wanted = c(4:6, 3, 7:12), 
                                        EvenXlims = TRUE, DoMu = FALSE,
                                        xlim = c(-10, 10), main = "height")

# loop to generate HPD plots of random effects for all traits
for (i in varNames) {
  results[[i]]$AllDiallelObs[[1]]$PlotHPD(wanted = c(4:6, 3, 7:12), 
                                               EvenXlims = TRUE, DoMu = FALSE,
                                               xlim = c(-10, 10), main = paste(i))
}

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