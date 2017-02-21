# ------------------------------------------------------------------------------
#   Plot phenotypic correlations - Fig2
#   S. Turner 
#   06 January 2017       
# ------------------------------------------------------------------------------
library(plyr)
library(corrplot)
library(psych)
library(reshape2)
library(ggplot2)

setwd("~/Documents/carrot-diallel")
# read in raw data
diallelRaw <- read.csv("diallelRawData.csv", header = TRUE)

# apply transformations for biomass
diallelRaw$dlw <- log(diallelRaw$dlw)
diallelRaw$drw <- log(diallelRaw$drw)
diallelRaw$ratio <- diallelRaw$dlw/diallelRaw$drw

# rename variables
diallelRaw <- rename(diallelRaw, c("midHeight" = "height (80DAP)", 
                                       "midWidth" = "width (80DAP)",
                                       "height" = "height (130DAP)", 
                                       "width" = "width (130DAP)", 
                                       "dlw" = "shoot biomass", 
                                       "drw" = "root biomass", 
                                       "ratio" = "shoot:root ratio"))

# create correlation matrix for selected traits
correlation <- cor(diallelRaw[,c(9:12, 14, 16:17)], method = "pearson", 
                   use = "complete.obs")

# function to calculate significance of correlations
# source: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
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

# export matrix of p-values
p.mat <- cor.mtest(diallelRaw[,c(9:12, 14, 16:17)])

# create matrix for significance notations 
# * = P < 0.05, ** = P < 0.01, *** = P < 0.001
sigCodes <- ifelse(p.mat[[1]] > 0.05, "NS", 
                   ifelse(p.mat[[1]] > 0.01, "*", 
                          ifelse(p.mat[[1]] > 0.001, "**", "***")))

sigCodes[upper.tri(sigCodes)] <- NA

# combine correlation coefficients and significance codes into single matrix
combinedMat <- lowerUpper(upper = sigCodes, lower = round(correlation, 2), 
                          diff = FALSE)
diag(combinedMat) <- ""

# create labels as input for geom_text
labels <- melt(combinedMat, na.rm = TRUE)
labels$value <- revalue(labels$value, c("0.7" = "0.70", "0.1" = "0.10"))

# melt matrix into vector format for plotting
meltedCormat <- melt(correlation)

ggplot(data = meltedCormat, aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#ca0020", high = "#2166AC", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Pearson\nCorrelation") + 
  scale_y_discrete(name = "", limits = rev(levels(meltedCormat$Var1))) +
  theme(axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) + 
  geom_text(aes(Var2, Var1, label = labels$value), colour = "black", size = 3.15) + 
  coord_fixed()

ggsave("Fig2.eps", height = 129, width = 129, units = "mm")


