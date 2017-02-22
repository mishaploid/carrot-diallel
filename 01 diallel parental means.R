# ------------------------------------------------------------------------------
#   Diallel - plot parental means for height (130 DAP) - Fig1
#   S. Turner 
#   06 January 2017
# ------------------------------------------------------------------------------

library(plyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/carrot-diallel-all-files/")

# ------------------------------------------------------------------------------
# read in parental data, subset, and calculate shoot:root ratio
# ------------------------------------------------------------------------------
pMeans <- read.csv("parentalData.csv", header = TRUE) 

# subset to exclude extra info (plot, root source, etc.)
pMeans <- pMeans[,c(2, 9:14)]

# calculate raito of shoot to root biomass 
# note: no ln(x) transformations in this data
pMeans$ratio <- pMeans$dlw/pMeans$drw 

# subset new data frames to plot A and B lines separately
bLines <- subset(pMeans, pedigree %in% c("L6038B", "L7550B", "P0159B", "Nbh2189B", "P6139B", "B7262B"))
aLines <- subset(pMeans, !(pedigree %in% c("L6038B", "L7550B", "P0159B", "Nbh2189B", "P6139B", "B7262B")))

# recode inbred line names to combine A and B lines
pMeans$pedigree <- revalue(pMeans$pedigree, c("B7262A" = "B7262", "B7262B" = "B7262", "L6038A" = "L6038",
                                              "L6038B" = "L6038", "L7550A" = "L7550", "L7550B" = "L7550",
                                              "Nbh2189B" = "Nbh2189", "P.S.C. x P0159B" = "P0159",
                                              "P0159B" = "P0159", "P6139A" = "P6139", "P6139B" = "P6139",
                                              "S.C. x Nbh2189B" = "Nbh2189"))

# set factor levels for parental lines
pMeans$pedigree <- factor(pMeans$pedigree, levels = c("L6038", "L7550", "P0159", "Nbh2189", "P6139", "B7262"))

# extract trait names
traits <- names(pMeans[2:8]) # extract trait names 

# stack data by trait (long format)
pMeans2 <- gather(pMeans, -pedigree, key = "var", value = "value") 

# calculate mean and standard error of each trait for each parental line
pMeans3 <- aggregate(.~pedigree + var, data = pMeans2, 
                     FUN = function(x) c(mu = mean(x), 
                                         se = sd(x)/sqrt(length(x))))

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

# plot parental means + 95% confidence interval
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

ggsave("FigS1.eps", height = 8, width = 6, units = "in")

# plot height only (Fig 1)
ht <- subset(pMeans4, var == "height (130 DAP)")

ggplot(ht, aes(x = pedigree, y = value.mu, colour = pedigree)) +
  geom_errorbar(aes(ymin = value.mu-2*value.se, ymax = value.mu+2*value.se), 
                width = 0.3, size = 1) +
  ylab("Canopy Height (cm)") + 
  geom_point(size = 2) +
  scale_y_continuous(limits = c(25, 60)) + 
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
ggsave("Fig1.eps", height = 84, width = 84, units = "mm")

