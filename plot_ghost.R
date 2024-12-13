#!/usr/bin/env Rscript
# by Theo Tricou

library(ape)
library(ggtree)

sim_dir <- Sys.glob("sim_*")

# plot number of transfer by branch lenth
sapply(sim_dir, function(x){
  d <- read.table(paste(x, "/SAMPLE_1/", "stat_ale_ghost_", x, sep = ""), h = T)
  # create variable for coloring tips
  d$color = "black"
  # color tips with ghost red
  d[which(d$L_ghost_branch > 0), "color"] = "#ff0000"
  pdf(file = paste("plot_transfer_",x , ".pdf", sep = ""), width = 8, height = 8)
  png(filename = paste("plot_transfer_", x, ".png", sep = ""), width = 800, height = 800, res = 100)
  plot(d$br_length, d$N_transfers_donor,col = d$color, pch=19 ,
       cex = .8, xlab = "Branch length", ylab = "Sum of inferred transfers")
  lm_big = lm(d$N_transfers_donor ~ d$br_length)
  abline(lm_big, col = "blue")
  dev.off()
})

# plot phylegenetic tree
sapply(sim_dir, function(x){
  #read complete tree
  tree = read.tree(paste(x,'/T/CompleteTree.nwk', sep = ""))
  #read sample tree
  stree = read.tree(paste(x,'/SAMPLE_1/aletree', sep = ""))
  ghost_tips = tree$tip.label[!(tree$tip.label %in% stree$tip.label)]
  sample_tip = tree$tip.label[(tree$tip.label %in% stree$tip.label)]
  # create list of ghost and no ghsot tips
  grp = list(ghost = ghost_tips, samp = sample_tip)
  # plot the tree
  if (length(grp$ghost) > 0){
    ttree = groupOTU(tree, grp, 'Species')
    ggtree(tree, size = 1, ladderize=FALSE) %>% groupOTU( grp, 'Species' ) + aes(color=Species) +
      scale_color_manual(values = c("#B5C0D0", "#432E54")) + theme(legend.position="none") # nolint
  }else{
    ggtree(tree, size = 1, ladderize=FALSE)
  }
  ggsave(filename = paste("plot_plot_",x ,".pdf", sep = ""))
  ggsave(filename = paste("plot_plot_",x ,".png", sep = ""))
})

# GNU Ghost
