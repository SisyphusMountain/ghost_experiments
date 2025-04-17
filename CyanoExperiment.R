# April 2025
# CODE USED FOR THE ANALYSIS OF THE CYANOBACTERIA EXPERIMENT

require(ape)
require(phangorn)
require(scales)
require(ggplot2)
require(ggtree)

### READ DATASETS
## table with all scores, all amputations, and name  mapping before and after amputation 
tab <- read.csv("cyano/analysis_results/uTs_df_with_original_names.csv")
## Table with original transfers scores (before pruning experiment)
tabinit <- read.csv("cyano/analysis_results/gene_uts.csv")

## GET GENE LIST AND TOTAL SUM OF TRANSFERS FROM EACH BRANCH (BEFORE PRUNING)
tab2 <- read.csv("cyano/analysis_results/gene_uml.csv")
tab2 <- tab2[tab2$Gene_Index != "",]
nbtransfpergene2 <- tab2$Total_Transfers
genes2 <- tab2$Gene_Index
names(nbtransfpergene2) <- genes2
genes2remove <- c('HBG635057', 'HBG634854', 'HBG636434')
genes2 <- setdiff(genes2, genes2remove)
nbtransfpergene2 <- nbtransfpergene2[genes2]
#rm genes with many duplications
nbtransfpergene2 <- nbtransfpergene2[!is.na(match(genes2, tab2$Gene_Index[tab2$Total_Duplications < 1]))]


## READ TREE AND GET INFORMATION
sptr <- read.tree("cyano/gene_tree/aletree")
tipnodes <- c(sptr$tip.label, sptr$node.label) #tip and nodes in the right order

##
MatchNodesBeforeAfter<-function(genetree, speciestree, ghosts) {
  #ghosts is the list of all ghost nodes. 
  ghoststips <- ghosts[!is.na(match(ghosts, speciestree$tip.label))]
  speciestree2 <- drop.tip(speciestree,ghoststips)
  tipnodesgn <- c(genetree$tip.label, genetree$node.label)
  tipnodessp <- c(speciestree2$tip.label, speciestree2$node.label)
  descgn <- unlist(
    lapply(Descendants(genetree), function(x,sp) paste(sort(sp[x]), collapse="-"), sp=genetree$tip.label)
  )
  descsp <- unlist(
    lapply(Descendants(speciestree2), function(x,sp) paste(sort(sp[x]), collapse="-"), sp=speciestree2$tip.label)
  )
  df <- data.frame(ablatnode=tipnodesgn, orignode=tipnodessp[match(descgn, descsp)])
  df
}

ComputeScores <- function(genetree, matching,matgn, matsp, indnode, dadnod, ghosts) {

	focalnode <- matching$ablatnode[matching$orignode==indnode]
	tipnodesgn <- c(genetree$tip.label, genetree$node.label)
	matgn2 <- data.frame(fromto=paste0(matgn$from_ale, "-",matgn$to_ale), freq=matgn$freq)
	matgn3 <- aggregate(freq ~ fromto, data = matgn2, sum)
	matgn4 <- cbind(t(apply(matgn3, 1, function(x) strsplit(x[1],"-")[[1]])),matgn3$freq)

	#Rename the dadnode donor and receiver as induced node donor to not lose them.
	matsp[which(matsp$from==dadnod),1] <- indnode
	matsp[which(matsp$to==dadnod),2] <- indnode

	#add new names to matsp
	matsp$from_ale <- matching$ablatnode[match(matsp$from, matching$orignode)]
	matsp$to_ale <- matching$ablatnode[match(matsp$to, matching$orignode)]
	matsp <- matsp[!is.na(matsp$from_ale)&!is.na(matsp$to_ale),] #remove NA
	matsp2 <- data.frame(fromto=paste0(matsp$from_ale, "-",matsp$to_ale), freq=matsp$freq)
	matsp3 <- aggregate(freq ~ fromto, data = matsp2, sum)
	matsp4 <- cbind(t(apply(matsp3, 1, function(x) strsplit(x[1],"-")[[1]])),matsp3$freq)
	original.source <- sapply(tipnodesgn, function(x,mat) sum(as.numeric(mat[mat[,1]==x,3])),mat=matsp4)
	ablated.source <- sapply(tipnodesgn, function(x,mat) sum(as.numeric(mat[mat[,1]==x,3])),mat=matgn4)
	original.target <- sapply(tipnodesgn, function(x,mat) sum(as.numeric(mat[mat[,2]==x,3])),mat=matsp4)
	ablated.target <- sapply(tipnodesgn, function(x,mat) sum(as.numeric(mat[mat[,2]==x,3])),mat=matgn4)
	resu <- data.frame(original.source, ablated.source, original.target, ablated.target)
	return(resu)
}

DOIT<-function(ablation_node, gene_filter=4, transfer_filter=0.5, plot=TRUE, printDF=FALSE) {
	#filter gene list based on gene_filter score
	genelist<-names(nbtransfpergene2[nbtransfpergene2<=gene_filter])
	#filter df of ablation experiment (`tab`) to select (i) good genes based on genelist (see above), 
	#(ii) good transfers (transfer scores >= transfer_filter) and (iii) scores 
	#for the current pruning experiment (ablation_node)
	tab.small<-tab[!is.na(match(tab$gene_name, genelist)),] #keep only good genes
	tab.small2<-tab.small[tab.small$freq>=transfer_filter,] #keep only good transfers
	tab.small3<-tab.small2[tab.small2$removed_clade==ablation_node, ]
	if (nrow(tab.small3) == 0) {
	    return(c(NA,NA,NA))
	}
	#filter df of initial experiment (`tabinit`, no ablation) to select (i) good genes based on genelist (see above) and 
	#(ii) good transfers (transfer scores >= transfer_filter) 
	tabinit.small<-tabinit[!is.na(match(tabinit$identifier, genelist)),]
	tabinit.small2<-tabinit.small[tabinit.small$freq>=transfer_filter,]

	#remove from tabinit.small2 the transfers TO nodes belonging to the ghost clade. 
	#this is required for comparing what is comparable before and after ablation.
	ghostclade <- if (is.na(match(ablation_node, sptr$node.label))) ablation_node else c(ablation_node, tipnodes[Descendants(sptr, match(ablation_node, sptr$node.label) + Ntip(sptr), "all")])

	#remove these nodes from tabinit.small2 to create `tabinit.small3`
	tabinit.small3<-tabinit.small2[is.na(match(tabinit.small2$to, ghostclade)),]

	#GET THE TRANSFERS SCORES BEFORE AND AFTER, MATCHING CORRECTLY THE NODES.
	#NB: AFTER means AFTER ABLATION, BEFORE means BEFORE ABLATION

	ablatedtree<-read.tree(paste("cyano/pruned/",ablation_node,"/species_tree.spTree", sep=""))
	tipnodesgn<-c(ablatedtree$tip.label, ablatedtree$node.label)
	corrnodes<-MatchNodesBeforeAfter(ablatedtree, sptr, ghostclade) #correspondence between nodes

	sibnode<-tipnodes[Siblings(sptr,match(ablation_node,tipnodes))]
	dadnode<-tipnodes[Ancestors(sptr,match(ablation_node,tipnodes), "parent")]
	##COMPUTE ALL SCORES
	DF<-ComputeScores(ablatedtree, corrnodes, tab.small3, tabinit.small3, sibnode, dadnode, ghostclade)
	##GET THE ghost lengthmatsp$freq[match(ghostclade,matsp$from)
	ghostlength<-sum(sptr$edge.length[match(ghostclade,tipnodes[sptr$edge[,2]])])
	ghosttransfers<-sum(tabinit.small3$freq[!is.na(match(tabinit.small3$from, ghostclade))], na.rm=TRUE)

	#get index of the focal branch the one initially carrying the ghost
	focal<-corrnodes$ablatnode[corrnodes$orignode==sibnode]
	focal.granddad<-tipnodesgn[Ancestors(ablatedtree,match(focal,tipnodesgn), "parent")]
	focal.granddad<-focal.granddad[!is.na(focal.granddad)] #remove NA
	focal.nephewnieces<-tipnodesgn[Descendants(ablatedtree,match(focal,tipnodesgn), "children")]
	focal.nephewnieces<-focal.nephewnieces[!is.na(focal.nephewnieces)] #remove NA
	non.focal<-setdiff(rownames(DF), c(focal, focal.granddad, focal.nephewnieces))
	#temp: add non.focal.all for all rows except the focal one. For keeping choice to make a simpler plot.
	non.focal.all<-setdiff(rownames(DF), focal)

	diff.focal.source<-DF$ablated.source[match(focal, rownames(DF))]-DF$original.source[match(focal, rownames(DF))]
	diff.focal.target<-DF$ablated.target[match(focal, rownames(DF))]-DF$original.target[match(focal, rownames(DF))]
	diff.focal.granddad.source<-DF$ablated.source[match(focal.granddad, rownames(DF))]-DF$original.source[match(focal.granddad, rownames(DF))]
	diff.focal.granddad.target<-DF$ablated.target[match(focal.granddad, rownames(DF))]-DF$original.target[match(focal.granddad, rownames(DF))]
	diff.focal.nephewnieces.source<-mean(DF$ablated.source[match(focal.nephewnieces, rownames(DF))]-DF$original.source[match(focal.nephewnieces, rownames(DF))])
	diff.focal.nephewnieces.target<-mean(DF$ablated.target[match(focal.nephewnieces, rownames(DF))]-DF$original.target[match(focal.nephewnieces, rownames(DF))])
	diff.non.focal.source<-mean(DF$ablated.source[match(non.focal, rownames(DF))]-DF$original.source[match(non.focal, rownames(DF))])
	diff.non.focal.target<-mean(DF$ablated.target[match(non.focal, rownames(DF))]-DF$original.target[match(non.focal, rownames(DF))])
	diff.non.focal.all.source<-mean(DF$ablated.source[match(non.focal.all, rownames(DF))]-DF$original.source[match(non.focal.all, rownames(DF))])
	diff.non.focal.all.target<-mean(DF$ablated.target[match(non.focal.all, rownames(DF))]-DF$original.target[match(non.focal.all, rownames(DF))])

  RES1 <- data.frame(
        ghostlength, ghosttransfers, diff.focal.source,diff.focal.target,diff.focal.granddad.source,diff.focal.granddad.target,diff.focal.nephewnieces.source,diff.focal.nephewnieces.target,diff.non.focal.source,diff.non.focal.target, diff.non.focal.all.source, diff.non.focal.all.target
        )
  RES2 <- DF$ablated.source[match(non.focal.all, rownames(DF))]-
      DF$original.source[match(non.focal.all, rownames(DF))]
  names(RES2) <- non.focal.all
	return(list(RES1=RES1, RES2=RES2))
}


###################
#### ANALYSIS #####
###################

nodes2visit<-setdiff(tipnodes, c(tipnodes[Ntip(sptr)+1],tipnodes[Children(sptr, Ntip(sptr)+1)]))
res<-NULL
res2<-list()
cpt<-0
for (s in nodes2visit) {
  cpt<-cpt+1
  print(s)
  temp<-DOIT(s, gene_filter=6, transfer_filter=0.5, plot=FALSE)
  res<-rbind(res, temp$RES1)
  res2[[cpt]]<-temp$RES2 
}

rownames(res)<-nodes2visit
names(res2)<-nodes2visit

# NODES REMOVED FROM THE ANALYSIS (see Supp Material)
res.rmn<-res[is.na(match(rownames(res), c("56","68","63","54"))),]
res2.rmn<-res2[is.na(match(names(res2), c("56","68","63","54")))]


###################
#### PLOTS    #####
###################

##### 1. THE TREE
require(ggtree)

pdf("thetree.pdf", height=7, width=5)
p <- ggtree(sptr)
p$data$label_color <- ifelse(p$data$label %in% c("70", "65", "69","54","63","56","68"), "grey60", "black")
p + geom_tiplab(offset=0.02, size=3) +
  geom_label2(aes(subset = !isTip, label=label, color=label_color), size=3) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(5, 30, 5, 5, "pt")) + 
  scale_color_identity()
dev.off()

##### 2. THE LARGE BOXPLOT
# Prepare data for "non-induced branches"
boxplot_data <- data.frame(
  Branch = rep(names(res2.rmn), sapply(res2.rmn, length)),
  Difference = unlist(res2.rmn)
)
# Prepare data for induced branches
induced_data <- data.frame(
  Branch = rownames(res.rmn),
  Difference = res.rmn$diff.focal.source
)
# order branches by increasing induced branch difference
boxplot_data$Branch <- factor(boxplot_data$Branch, levels = rownames(res.rmn)[order(res.rmn$diff.focal.source)])
# plot
pdf("resPerBranch.pdf", width=7, height=12)
ggplot(boxplot_data, aes(x = Branch, y = Difference)) + 
  geom_boxplot(size=0.3, outlier.shape=21, outlier.fill=NA, outlier.size=0.8) + 
  coord_flip() + 
  geom_point(induced_data, mapping=aes(x=Branch, y=Difference), color="red", size=1.5) + 
  labs(x = "Pruned clade", y = "Difference in number of transfers per branch (after pruning - before pruning)") +
   theme_minimal()
dev.off()

##### 3. THE NORMALIZED RANKS
# Get and format normalized ranks
NormalizedRanks <- sapply(1:length(res.rmn$diff.focal.source), function(x,ind, nonind) (rank(c(ind[x],nonind[[x]]))[1]-1)/(length(nonind[[x]])), ind=res.rmn$diff.focal.source, nonind=res2.rmn)
# Prepare data for ranks
rank_data <- data.frame(
  Branch = rownames(res.rmn),
  NormalizedRank = NormalizedRanks
)
# Create boxplot + jitter
pdf("NormalizedRanksBoxplot.pdf", width=5, height=5)
ggplot(rank_data, aes(x = "", y = NormalizedRank)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill="grey90") +
  geom_jitter(width = 0.2, size = 1.5, color="red") +
  labs(x = "", y = "Normalized Rank") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
  #+ stat_summary(fun = mean, geom = "crossbar", width = 0.1,
  #             color = "grey10", linetype = "12", size = 0.3)
dev.off()

##### 4. THE WILCOXON TESTS

# Perform Wilcoxon test 1000 times for 1000 samplings
PVALUES_RESULTS<-replicate(1000, wilcox.test(unlist(lapply(res2.rmn, function(x) sample(x,1))), res.rmn$diff.focal.source, alternative="l")$p.value)
DF_PVALUES<-data.frame(pval=PVALUES_RESULTS)
pdf("wilcox_test_pvalues.pdf")
ggplot(DF_PVALUES, aes(x = pval)) +
  geom_histogram(bins=30,fill = "steelblue", alpha = 0.6, color="black") + 
  geom_vline(xintercept = 0.05, color = "darkgrey", linetype = "dashed", linewidth = 1) + 
  xlab("P-values from Wilcoxon test") + ylab("Counts")
dev.off()

