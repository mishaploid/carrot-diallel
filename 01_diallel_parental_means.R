# ------------------------------------------------------------------------------
#   Plot parental phenotypes (mean + 95% confidence interval)
#   S. Turner 
#   06 January 2017
# ------------------------------------------------------------------------------

# This script plots parental means and 95% confidence intervals 

library(plyr)
library(tidyr)
library(ggplot2)

setwd("~/GitHub/carrot-diallel/data/")

# ------------------------------------------------------------------------------
# read in parental data, subset, and calculate shoot:root ratio
# ------------------------------------------------------------------------------
pMeans <- read.csv("parentalData.csv", header = TRUE) 

# subset to exclude extra info (plot, root source, etc.)
pMeans <- pMeans[,c(2, 9:14)]

# calculate raito of shoot to root biomass 
# note: no ln(x) transformations in this data
pMeans$ratio <- pMeans$dlw/pMeans$drw 

# optional: 
# subset new data frames to plot A (sterile) and B (fertile) lines separately
bLines <- subset(pMeans, pedigree %in% c("L6038B", "L7550B", "P0159B", 
                                         "Nbh2189B", "P6139B", "B7262B"))
aLines <- subset(pMeans, !(pedigree %in% c("L6038B", "L7550B", "P0159B", 
                                           "Nbh2189B", "P6139B", "B7262B")))
# ------------------------------------------------------------------------------
# recode inbred line names to combine A and B lines
# ------------------------------------------------------------------------------
pMeans$pedigree <- revalue(pMeans$pedigree, c("B7262A" = "B7262", 
                                              "B7262B" = "B7262", 
                                              "L6038A" = "L6038",
                                              "L6038B" = "L6038", 
                                              "L7550A" = "L7550", 
                                              "L7550B" = "L7550",
                                              "Nbh2189B" = "Nbh2189", 
                                              "P.S.C. x P0159B" = "P0159",
                                              "P0159B" = "P0159", 
                                              "P6139A" = "P6139", 
                                              "P6139B" = "P6139",
                                              "S.C. x Nbh2189B" = "Nbh2189"))

# ------------------------------------------------------------------------------
# restructure data
# ------------------------------------------------------------------------------
# set factor levels for parental lines
pMeans$pedigree <- factor(pMeans$pedigree, levels = c("L6038", "L7550", "P0159", 
                                                      "Nbh2189", "P6139", "B7262"))

# extract trait names
traits <- names(pMeans[2:8]) # extract trait names 

# stack data by trait (long format)
pMeans2 <- gather(pMeans, -pedigree, key = "var", value = "value") 

# ------------------------------------------------------------------------------
# calculate mean and standard error of each trait for each parental line
# ------------------------------------------------------------------------------
pMeans3 <- aggregate(.~pedigree + var, data = pMeans2, 
                     FUN = function(x) c(mu = mean(x), 
                                         se = sd(x)/sqrt(length(x))))

# ------------------------------------------------------------------------------
# set formatting for plot
# ------------------------------------------------------------------------------
# format results as a dataframe 
pMeans4 <- do.call(data.frame, pMeans3) 

# set factor levels for plotting
pMeans4$var <- factor(pMeans4$var, levels = c("height", "width", "flw", "dlw",
                                              "frw", "drw", "ratio"))
# recode variables with pretty names
pMeans4$var <- revalue(pMeans4$var, c("height" = "height (130 DAP)", 
                                      "width" = "width (130 DAP)", 
                                      "flw" = "shoot biomass (fresh)", 
                                      "dlw" = "shoot biomass (dry)",
                                      "frw" = "root biomass (fresh)",
                                      "drw" = "root biomass (dry)"))

# ------------------------------------------------------------------------------
# plot parental means + 95% confidence interval
# ------------------------------------------------------------------------------
parents <- ggplot(pMeans4, aes(x = pedigree, y = value.mu, colour = pedigree)) +
  geom_errorbar(aes(ymin = value.mu-2*value.se, ymax = value.mu + 2*value.se), 
                width = 0.3, size = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                 "#0072B2", "#D55E00", "#CC79A7")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  facet_wrap(~ var, scales = "free", ncol = 3)

parents

dir.create("~/GitHub/carrot-diallel/results")
ggsave("~/GitHub/carrot-diallel/results/parental_means.eps", height = 8, width = 6, units = "in")

# ------------------------------------------------------------------------------
# plot height only (e.g. Fig 1 in Turner et al.)
# ------------------------------------------------------------------------------
ht <- subset(pMeans4, var == "height (130 DAP)")

ggplot(ht, aes(x = pedigree, y = value.mu, colour = pedigree)) +
  geom_errorbar(aes(ymin = value.mu-2*value.se, ymax = value.mu+2*value.se), 
                width = 0.3, size = 1) +
  ylab("Canopy height (130 DAP)") + ylim(c(25,60)) +
  geom_point(size = 2) +
  #scale_y_continuous(limits = c(25, 70)) + 
  scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                 "#0072B2", "#D55E00", "#CC79A7")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

# export plot as .eps file
ggsave("~/GitHub/carrot-diallel/results/parental_means_ht.eps", 
       height = 84, width = 84, units = "mm")

