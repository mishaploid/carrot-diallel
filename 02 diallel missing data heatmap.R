# ------------------------------------------------------------------------------
#   Diallel heat maps (missing data visualization)
#   S. Turner 
#   06 July 2016
# ------------------------------------------------------------------------------

# Plots missing data by location, rep, and hybrid combination

library(ggplot2)
library(plyr)

setwd("~/Documents/carrot-diallel")

# read in raw data
diallelRaw <- read.csv("diallelRawData.csv", header = TRUE)

str(diallelRaw)

# ------------------------------------------------------------------------------
# set, reorder, and rename factors 
# ------------------------------------------------------------------------------
diallelRaw$female <- ordered(diallelRaw$female, levels = c("1", "2", "3", "4", "5", "6"))
diallelRaw$male <- ordered(diallelRaw$male, levels = c("6", "5", "4", "3", "2", "1"))
diallelRaw$year <- as.factor(diallelRaw$year)
diallelRaw$year <- revalue(diallelRaw$year, c("1" = "CA2015", "2" = "WI2015", "3" = "CA2016"))
diallelRaw$rep <- as.factor(diallelRaw$rep)

# ------------------------------------------------------------------------------
# take natural log for biomass
# ------------------------------------------------------------------------------
diallelRaw$lnFLW <- log(diallelRaw$flw)
diallelRaw$lnDLW <- log(diallelRaw$dlw)
diallelRaw$lnFRW <- log(diallelRaw$frw)
diallelRaw$lnDRW <- log(diallelRaw$drw)

# shoot:root ratio
diallelRaw$ratio <- diallelRaw$lnDLW/diallelRaw$lnDRW

# ------------------------------------------------------------------------------
# function to draw heatmaps
# ------------------------------------------------------------------------------

plotNA <- function(df, xname) {
  p <- ggplot(df, aes(female, male)) + 
    geom_tile(aes_string(fill = xname), colour = "black") + 
    scale_fill_gradient(low = "light grey", high = "black", na.value = "white") + 
    geom_text(data = df[is.na(df[xname]),], aes(label = "NA", fontface = "bold")) + 
    facet_grid(year ~ rep, space = "free", scales = "free", labeller = label_both)
  print(p)
}

# ------------------------------------------------------------------------------
# loop for heatmap output
# ------------------------------------------------------------------------------

# extract trait names
traits <- c("midHeight", "midWidth", "height", "width", "flw", 
            "dlw", "frw", "drw", "ratio")

for (i in traits) {
  plotNA(diallelRaw, i)
}


