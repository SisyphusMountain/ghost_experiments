# April 2025
# CODE USED FOR THE ANALYSIS OF THE CYANOBACTERIA TREE USABILITY

require(ape)
require(phangorn)
require(ggtree)
require(ggplot2)

tr<-read.tree("cyano/analysis_results/cyano_renamed.nwk")
tab<-read.csv("cyano/analysis_results/combined_spr_classifications.csv")
tipnodes<-as.numeric(c(tr$tip.label, tr$node.label))
namecorresp<-read.csv("cyano/analysis_results/cyano_correspondence.csv")
tree<-read.tree("cyano/analysis_results/cyano.nwk")


###################
#### ANALYSIS #####
###################

SCORE<-NULL
for (n in unique(tab$removed_node)) {
  # We look for transfers for the ghost clade (n) ablation experiment only
  smalltab <- tab[which(!is.na(match(tab$removed_node, as.numeric(n)))),]
  # we extract the list of ghost nodes
  allghost <- tipnodes[Descendants(tr, match(n, tipnodes),type="all")]
  # we get the induced node
  induced_n<-smalltab$induced_donor[1]
  # we get only the transfers from the ghost clade (its induced node in fact, but it's the same)
  smalltab_donorghost <- smalltab[smalltab$donor==induced_n,]
  # Small script to get the id of the topology that is identical to the no-transfer one
  pruned_tr <- drop.tip(tr, Descendants(tr, match(n, tipnodes))[[1]])
  tipnodes2 <- as.numeric(c(pruned_tr$tip.label, pruned_tr$node.label))
  sib_n <- tipnodes2[Siblings(pruned_tr, match(induced_n, tipnodes2))]
  if (length(sib_n)==0) { #we are treating one of the two descendants of the root (node 2 and 3)
    if (n==2) topo_identique <- smalltab$topology_id[smalltab$donor==6&smalltab$receiver==7]
    if (n==3) topo_identique <- smalltab$topology_id[smalltab$donor==4&smalltab$receiver==5]
    }
  else {
    topo_identique <- smalltab_donorghost$topology_id[smalltab_donorghost$receiver==sib_n]
  }
  # we remove transfers that do not change the topology
  smalltab_donorghost_rmidentical<-smalltab_donorghost[smalltab_donorghost$topology_id!=topo_identique,]
  # the smalltab_donorghost_rmidentical matrix contains all possible transfers from the ghost clade (to all possible recipients outside the ghost clade) that are not identical to the no transfer topology
  # Now we get ALL possible transfers that do not have the ghost clade as donor and that change the topology
  smalltab_donor_no_ghost <- smalltab[smalltab$donor!=induced_n,]

  #we look for each topology obtained by transfers from the ghost whether it exists in the list of transfers from other clades
  computescore <- match(smalltab_donorghost_rmidentical$topology_id, smalltab_donor_no_ghost$topology_id)
  # we compute the proportion.
  SCORE <- rbind(SCORE, c(n, sum(!is.na(computescore))/length(computescore)))
}

colnames(SCORE)<-c("node", "score")
SCORE<-as.data.frame(SCORE)

###################
#### PLOT    ######
###################

pdf("NodesTree.pdf", height=9, width=7)
p <- ggtree(tr) + geom_tiplab(size = 3, align = TRUE, offset=0.05, linetype="blank")
p$data$score1<-SCORE$score[match(p$data$label, SCORE$node)]
#rename tip and node labels
p$data$label<-namecorresp$name[as.numeric(p$data$label)]
p + geom_point(aes(fill = score1), shape=21, size = 7.5, color="black", stroke=0.5) + geom_text(aes(label = round(score1, 2)), size = 2.5, color = "white", fontface="bold") + scale_fill_gradient(low = "blue", high = "red", limits = c(0, NA)) + labs(fill = "Proportion of\nnon-specific\ntopologies\n") + xlim(0,1.1)

dev.off()

