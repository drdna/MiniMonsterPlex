## Plots trees and colors taxa according to lineage affinities

library(ape)
#library(phangorn)
#library(phytools)
#library(RColorBrewer)
#library(tidyverse)
library(ggtree)
library(ggrepel)
#library(viridisLite)
#library(data.table)
library(ggtext)

# read arguments
args = commandArgs(trailingOnly=TRUE)

#read in tree data
Tree1 <- read.tree(args[1])

Tree1 <- drop.tip(Tree1, tip = "GN0001_Z")

# assign the groups
groupInfo <- split(Tree1$tip.label, gsub(".*_", "", Tree1$tip.label))

# code to write out IDs of strains used
#labels <- Tree1$tip.label
#write.csv(list(labels), "~/MLTreeStrains.csv")

groupIDs <- gsub(".*_", "", Tree1$tip.label)
groupLabels <- unique(lapply(sort(groupIDs), function(x) paste0("*", x, "*")))
groupInfo <- lapply(groupInfo, function(x) gsub("_.*", "", x))

# group isolates by host genus
Tree1$tip.label <- gsub("_.*", "", Tree1$tip.label)
Tree1 <- groupOTU(Tree1, groupInfo)

# adjust branch lengths for plotting
Tree1$edge.length <- Tree1$edge.length * 500

# set a color palette
c25 <- c(
  "yellow3",
  "#E31A1C",
  "green4",
  "#6A3D9A",
  "#FF7F00",
  "black",
  "gold1",
  "skyblue2",
  "#FB9A99",
  "palegreen2",
  "#FDBF6F",
  "#CAB2D6",
  "gray70",
  "khaki2",
  "maroon",
  "orchid1",
  "deeppink1",
  "dodgerblue2",
  "steelblue4",
  "darkturquoise",
  "green1",
  "yellow4",
  "blue1",
  "darkorange4",
  "brown"
)

# plot the tree
t <- ggtree(Tree1,
       layout = "rectangular", 
       branch.length="branch.length",
       aes(label = gsub(".*_", "", label))) +
  
       # the following block holds key tricks for creating decent legends
       geom_tiplab(aes(color = group), size = 2, linesize = 0.25, align = TRUE, offset = 2, show.legend=F) +
       scale_color_manual(values = c25) +
       geom_polygon(aes(fill = group, x = 0, y = 0)) +
       scale_fill_manual(name = "Host Genus",
                         labels = groupLabels,
                         values = c25,
                         guide=guide_legend(nrow=25, keywidth=1,
                                            keyheight=1,
                                            order=1,
                                            override.aes=list(size=5,alpha=1))) +
       # end of block
  
       geom_treescale(x= 5, y = -4) + xlim_tree(12) +
       theme(legend.position = c(0.1, 0.88), legend.text = element_markdown())

# save the tree
ggsave(paste0(args[2], ".pdf"), device="pdf", width =8, height = 24)

