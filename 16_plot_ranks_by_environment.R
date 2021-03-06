# ------------------------------------------------------------------------------
#   Plot rankings from BayesDiallel by environment
#   P. Maurizio; edited by S. Turner
#   6 July 2017
# ------------------------------------------------------------------------------

# this script plots mean rankings and 95% credibility intervals for each 
# environment. Hybrids with 95% intervals that always fall in the bottom 15 or
# top 15 ranks are color coded green and red, respectively.

library(data.table)
library(ggplot2)

varNames <- c("dlw", "drw", "flw", "frw", "height", "midHeight", "midWidth",
                       "ratio", "width")

cross_names <- c("L6038", "L7550 x L6038", "P0159 x L6038", "Nbh2189 x L6038",
                 "P6139 x L6038", "B7262 x L6038", "L6038 x L7550", "L7550",
                 "P0159 x L7550", "Nbh2189 x L7550", "P6139 x L7550",
                 "B7262 x L7550", "L6038 x P0159", "L7550 x P0159", "P0159",
                 "Nbh2189 x P0159", "P6139 x P0159", "B7262 x P0159", 
                 "L6038 x Nbh2189", "L7550 x Nbh2189", "P0159 x Nbh2189",
                 "Nbh2189", "P6139 x Nbh2189", "B7262 x Nbh2189", 
                 "L6038 x P6139", "L7550 x P6139", "P0159 x P6139",
                 "Nbh2189 x P6139", "P6139", "B7262 x P6139", 
                 "L6038 x B7262", "L7550 x B7262", "P0159 x B7262", 
                 "Nbh2189 x B7262", "P6139 x B7262", "B7262")

# read in rankings
setwd("~/GitHub/carrot-diallel/results/wi2015/")
temp_wi15 <- list.files(pattern="*_ranks_wi2015.csv")
ranks_wi15 <- lapply(temp_wi15, read.csv)
names(ranks_wi15) <- c("dlw", "drw", "flw", "frw", "height", "midHeight", "midWidth",
                       "ratio", "width")
ranks_wi15 <- mapply(cbind, ranks_wi15, "loc"="wi15", SIMPLIFY=F)
ranks_wi15 <- lapply(ranks_wi15, cbind, cross = cross_names)

setwd("~/GitHub/carrot-diallel/results/ca2015/")
temp_ca15 <- list.files(pattern="*_ranks_ca2015.csv")
ranks_ca15 <- lapply(temp_ca15, read.csv)
names(ranks_ca15) <- c("dlw", "drw", "flw", "frw", "height", "midHeight", "midWidth",
                       "ratio", "width")
ranks_ca15 <- mapply(cbind, ranks_ca15, "loc"="ca15", SIMPLIFY=F)
ranks_ca15 <- lapply(ranks_ca15, cbind, cross = cross_names)

setwd("~/GitHub/carrot-diallel/results/ca2016/")
temp_ca16 <- list.files(pattern="*_ranks_ca2016.csv")
ranks_ca16 <- lapply(temp_ca16, read.csv)
names(ranks_ca16) <- c("dlw", "drw", "flw", "frw", "height", "midHeight", "midWidth",
                       "ratio", "width")
ranks_ca16 <- mapply(cbind, ranks_ca16, "loc"="ca16", SIMPLIFY=F)
ranks_ca16 <- lapply(ranks_ca16, cbind, cross = cross_names)

setwd("~/GitHub/carrot-diallel/results/")

# merge data
ranks_merged <- mapply(rbind, ranks_wi15, ranks_ca15, ranks_ca16, SIMPLIFY=FALSE)

# reorder crosses by mean rank
ordered_ranks <- lapply(ranks_merged, transform, ordered_ranks = reorder(cross, -X.1))

for(i in varNames){
  ordered_ranks[[i]]$col <- 
    ifelse(ordered_ranks[[i]]$lower > 20, "top 15", 
           ifelse(ordered_ranks[[i]]$upper < 16, "bottom 15", "middle"))
}

pd <- position_dodge(0.5) # move them .05 to the left and right

for(i in varNames){
  ggplot(ordered_ranks[[i]], aes(x=ordered_ranks, y=X.1, shape = loc, colour=col)) + 
    ylim(0, 36) + ylab("Rank") +
    coord_flip() +
    scale_color_manual(values=c("#1b9e77", "gray", "#d95f02")) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=pd) +
    geom_point(position=pd) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          panel.background = element_blank(),
          axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_blank())
  
  ggsave(paste("~/GitHub/carrot-diallel/results/ranks_by_loc_", i, ".eps", sep = ""), 
         width = 6, height = 9, units = "in")
}
