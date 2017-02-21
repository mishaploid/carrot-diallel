# ------------------------------------------------------------------------------
#   Plot diallel variance projection (VarP)
#   See Crowley et al. (2014) Genetics for details
#   S. Turner 
#   12 September 2016
# ------------------------------------------------------------------------------

# This script uses posterior P^2 values to plot the diallel variance projection
# (VarP) as described by Crowley et al. (2014) 
# e.g. Fig 6 in the carrot diallel manuscript

library(plyr)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------------------
# load compiled PSq table (from create psq data frame.R)
# ------------------------------------------------------------------------------
setwd("~/Documents/carrot-diallel/PSqTables/")

PSqTmp <- read.csv("PSq_all.csv", header = TRUE)

# ------------------------------------------------------------------------------
# formatting details
# ------------------------------------------------------------------------------
# specify color palette (colorblind friendly)
cbPalette <- c("#009E73", "#E69F00", "#56B4E9", "#0072B2", "#999999", "#F0E442",
               "#D55E00", "#CC79A7")

# select inheritance classes (i.e. exclude random effects) and traits
PSq <- subset(PSqTmp, X %in% c("aj", "ASymCrossjkDkj", "BetaInbred:Av",
                               "dominancej", "motherj", "SymCrossjk", 
                               "ASymCrossjkDkj", "Noise") &
                id %in% c("midHeightpsq", "midWidthpsq", "heightpsq", 
                          "widthpsq", "dlwpsq", "drwpsq", "ratiopsq"))

# set factor levels for inheritance classes
PSq$X <- factor(PSq$X, levels = c("aj", "motherj", "BetaInbred:Av", 
                                  "dominancej", "SymCrossjk", 
                                  "ASymCrossjkDkj", "Noise"))

# rename inheritance classes 
PSq$X <- revalue(PSq$X, c("aj" = "Additive (a)", "motherj" = "Maternal (m)",
                          "BetaInbred:Av" = "Inbred Penalty (B)",
                          "dominancej" = "Inbreeding (b)",
                          "SymCrossjk" = "Symmetric (v)", 
                          "ASymCrossjkDkj" = "Asymmetric (w)"))
# rename traits
PSq$id <- revalue(PSq$id, c("widthpsq" = "width\n(130DAP)", 
                            "midWidthpsq" = "width\n(80DAP)",
                            "midHeightpsq" = "height\n(80DAP)", 
                            "heightpsq" = "height\n(130DAP)",
                            "drwpsq" = "root biomass", 
                            "dlwpsq" = "shoot biomass",
                            "ratiopsq" = "shoot:root ratio"))
# set levels for traits
PSq$id <- factor(PSq$id, levels = c("shoot:root ratio", 
                                    "root biomass",
                                    "shoot biomass", 
                                    "width\n(130DAP)", 
                                    "width\n(80DAP)",
                                    "height\n(130DAP)",
                                    "height\n(80DAP)"))

# reverse order of traits and inheritance classes for plotting
PSq$revX <- factor(PSq$X, levels = rev(levels(PSq$X)))
PSq$revid <- factor(PSq$id, levels = rev(levels(PSq$id)))

# ------------------------------------------------------------------------------
# Plot PSq mean and 95% credibility intervals
# ------------------------------------------------------------------------------
p1 <- ggplot(data = PSq[order(PSq$revX),], aes(x = revX, y = Mean, ymin = X2.5., ymax = X97.5.,
                                               group = revX, color = revX, order = revX)) +
  geom_hline(yintercept = 0, col = "black", lwd = 1, linetype = 2) +
  geom_pointrange(stat = "identity", position = position_dodge(1), size = 0.5) + 
  coord_flip() +
  scale_color_manual(values = cbPalette[c(5, 8, 7, 3, 4, 2, 1)], name = "") +
  facet_grid(revid ~ ., switch = "y", labeller = label_wrap_gen(width = 9.5)) +
  scale_y_continuous(breaks = seq(-0.2, 1, by = .2), limits = c(-0.2, 1)) + 
  ylab("Diallel Variance Projection") +
  theme(legend.position = "none", panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
p1

# ------------------------------------------------------------------------------
# Stacked barplot - proportion of variation attributed to each inheritance class
# ------------------------------------------------------------------------------
# create a dummy variable for plotting
PSq$dummy <- 1

# ggplot doesn't like negative values for geom_bar
if(PSq$Mean < 0) PSq$Mean <- 0

p2 <- ggplot(PSq[PSq$Mean > 0,], aes(x = dummy, y = Mean, fill = revX, order = revX)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  ylab("Diallel Variance Projection") +
  xlab("") +
  scale_fill_manual(values = cbPalette[c(5, 8, 7, 3, 4, 2, 1)], name = "") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  facet_grid(revid ~., labeller = label_wrap_gen(width = 9.5)) + 
  theme(legend.position = "none", panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        text = element_text(size = 12), 
        axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

p2


# ------------------------------------------------------------------------------
# create and extract legend
# ------------------------------------------------------------------------------
legend <- ggplot(PSq[order(PSq$X),], aes(x = id, y = Mean, fill = X, order = X)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cbPalette[c(1, 2, 4, 3, 7, 8, 5)], name = "") +
  theme(legend.position = "top", text = element_text(size = 12))

# function to extract legend
# source: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(legend)

# ------------------------------------------------------------------------------
# Combine plots!
# ------------------------------------------------------------------------------
grid.arrange(mylegend, p1, p2, ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,1), c(2,3)),
             widths = c(3, 2.7), heights = c(0.2, 2))

