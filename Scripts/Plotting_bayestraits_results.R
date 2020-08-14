#importing bayestraits results
library(BTRTools)
library(phytools)
library(BTprocessR)
library(tidyverse)
library(pbapply)
library(ggplot2)
library(ggtree)
library(coda)
library(plotrix)

############

##This code is used to import the results of BayesTraits RJ-MCMC analysaes
##first run BayesTraits sing the supplied trees and PC scores. See details in file "BT_Variable_Rates_Readme.txt"
##This can can be used to generate figures 1, and S1-S36



color3<-colorRampPalette(c("#0c2c84","#225ea8","#31a354","#ffff00","#fe9929","#fc4e2a","red","darkred"))

#load in dated phylogenetic trees for each module
#Files marked with MR are dated with FBD model in Mr. Bayes
#files marked with T are traditional dinosauria topology dated with MBL method
#files marked with B are the Baron et al Ornithoscelida hypothesis dated with MBL method
#Confirm that time trees imported here match the BT RJMCMC results imported below (ie, MR, T, or B)

tree_rostrum<-read.nexus("Data/Trees/tree_MR_rostrum.nex")
tree_palate<-read.nexus("Data/Trees/tree_MR_palate.nex")
tree_vault<-read.nexus("Data/Trees/tree_MR_vault.nex")
tree_pterygoid<-read.nexus("Data/Trees/tree_MR_pterygoid.nex")
tree_occipital<-read.nexus("Data/Trees/tree_MR_occipital.nex")
tree_quadrate<-read.nexus("Data/Trees/tree_MR_quadrate.nex")
tree_sphenoid<-read.nexus("Data/Trees/tree_MR_sphenoid.nex")

##########
#Utility function for plotting trees with branches colored by rate
#modified from [phytools] function "plotBranchbyTrait"
mytreebybranch<-function (tree, x, mode = c("edges", "tips", "nodes"), palette = "rainbow", nbin=15,
                          legend = TRUE, xlims = NULL, ...) 
{
  mode <- mode[1]
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (mode == "tips") {
    x <- c(x[tree$tip.label], fastAnc(tree, x))
    names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  else if (mode == "nodes") {
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  if (hasArg(tol)) 
    tol <- list(...)$tol
  else tol <- 1e-06
  if (hasArg(prompt)) 
    prompt <- list(...)$prompt
  else prompt <- FALSE
  if (hasArg(type)) 
    type <- list(...)$type
  else type <- "phylogram"
  if (hasArg(show.tip.label)) 
    show.tip.label <- list(...)$show.tip.label
  else show.tip.label <- TRUE
  if (hasArg(show.node.label)) 
    show.node.label <- list(...)$show.node.label
  else show.node.label <- FALSE
  if (hasArg(edge.width)) 
    edge.width <- list(...)$edge.width
  else edge.width <- 4
  if (hasArg(edge.lty)) 
    edge.lty <- list(...)$edge.lty
  else edge.lty <- 1
  if (hasArg(font)) 
    font <- list(...)$font
  else font <- 3
  if (hasArg(cex)) 
    cex <- list(...)$cex
  else cex <- par("cex")
  if (hasArg(adj)) 
    adj <- list(...)$adj
  else adj <- NULL
  if (hasArg(srt)) 
    srt <- list(...)$srt
  else srt <- 0
  if (hasArg(no.margin)) 
    no.margin <- list(...)$no.margin
  else no.margin <- TRUE
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (hasArg(label.offset)) 
    label.offset <- list(...)$label.offset
  else label.offset <- 0.01 * max(nodeHeights(tree))
  if (hasArg(underscore)) 
    underscore <- list(...)$underscore
  else underscore <- FALSE
  if (hasArg(x.lim)) 
    x.lim <- list(...)$x.lim
  else x.lim <- NULL
  if (hasArg(y.lim)) 
    y.lim <- list(...)$y.lim
  else y.lim <- if (legend && !prompt && type %in% c("phylogram", 
                                                     "cladogram")) 
    c(1 - 0.06 * length(tree$tip.label), length(tree$tip.label))
  else NULL
  if (hasArg(direction)) 
    direction <- list(...)$direction
  else direction <- "rightwards"
  if (hasArg(lab4ut)) 
    lab4ut <- list(...)$lab4ut
  else lab4ut <- NULL
  if (hasArg(tip.color)) 
    tip.color <- list(...)$tip.color
  else tip.color <- "black"
  if (hasArg(plot)) 
    plot <- list(...)$plot
  else plot <- TRUE
  if (hasArg(rotate.tree)) 
    rotate.tree <- list(...)$rotate.tree
  else rotate.tree <- 0
  if (hasArg(open.angle)) 
    open.angle <- list(...)$open.angle
  else open.angle <- 0
  if (is.function(palette)) 
    cols <- palette(n = nbin)
  else {
    if (palette == "heat.colors") 
      cols <- heat.colors(n = nbin)
    if (palette == "gray") 
      cols <- gray(nbin:1/nbin)
    if (palette == "rainbow") 
      cols <- rainbow(nbin, start = 0.7, end = 0)
  }
  if (is.null(xlims)) 
    xlims <- range(x) + c(-tol, tol)
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
        1
    cols[i]
  }
  colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  par(lend = 2)
  xx <- plot.phylo(tree, type = type, show.tip.label = show.tip.label, 
                   show.node.label = show.node.label, edge.color = colors, 
                   edge.width = edge.width, edge.lty = edge.lty, font = font, 
                   cex = cex, adj = adj, srt = srt, no.margin = no.margin, 
                   root.edge = root.edge, label.offset = label.offset, underscore = underscore, 
                   x.lim = x.lim, y.lim = y.lim, direction = direction, 
                   lab4ut = lab4ut, tip.color = tip.color, plot = plot, 
                   rotate.tree = rotate.tree, open.angle = open.angle, lend = 2, 
                   new = FALSE)
  if (legend == TRUE && is.logical(legend)) 
    legend <- round(0.3 * max(nodeHeights(tree)), 2)
  if (legend) {
    if (hasArg(title)) 
      title <- list(...)$title
    else title <- "trait value"
    if (hasArg(digits)) 
      digits <- list(...)$digits
    else digits <- 1
    if (prompt) 
      add.color.bar(legend, cols, title, xlims, digits, 
                    prompt = TRUE)
    else add.color.bar(legend, cols, title, xlims, digits, 
                       prompt = FALSE, x = par()$usr[1] + 0.05 * (par()$usr[2] - 
                                                                    par()$usr[1]), y = par()$usr[3] + 0.05 * (par()$usr[4] - 
                                                                                                                par()$usr[3]))
  }
  invisible(xx)
}


#Import Bayestraits data and plot
#########

#define the folder that contains the output of Bayestraits runs
resultsfolder <- ""
#Rostrum region
#########




{
  BTraits_rostrum<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_rostrum_MR.txt.VarRates.txt", sep=""),
                                  rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_rostrum_MR.txt.Output.trees", sep=""),
                                  tree = tree_rostrum,
                                  burnin = 1250)
  
  cophylotrees <- cophylo(tree_rostrum,BTraits_rostrum$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
  tree_rostrum.a <- cophylotrees$trees[[1]]
}

label %in%  c(cladedefs$TipLabel[c(9993:10029)],"Gallus_gallus")))



mytreebybranch(tree = tree_rostrum.a, x=log(BTraits_rostrum$data$meanRate)[-1], 
               mode="edges",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("rostrum rates")


#########
#occipital region
#########

  BTraits_occipital<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_occipital_MR.txt.VarRates.txt", sep=""),
                                    rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_occipital_MR.txt.Output.trees", sep=""),
                                    tree = tree_occipital,
                                    burnin = 0)
  
  cophylotrees <- cophylo(tree_occipital,BTraits_occipital$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
  tree_occipital.a <- cophylotrees$trees[[1]]

mytreebybranch(tree = tree_occipital.1, x=log(BTraits_rostrum$data$meanRate)[-1], 
               mode="edges",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("rostrum rates")


#########
#Palate region
#########

BTraits_palate<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_palate_MR.txt.VarRates.txt", sep=""),
                               rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_palate_MR.txt.Output.trees", sep=""),
                               tree = tree_palate,
                               burnin = 0)


cophylotrees <- cophylo(tree_palate,BTraits_palate$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_palate.a <- cophylotrees$trees[[1]]

mytreebybranch(tree = tree_palate.a, x=log(BTraits_palate$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("palate rates")


#########
#Vault region
#########

BTraits_vault<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_vault_MR.txt.VarRates.txt", sep=""),
                              rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_vault_MR.txt.Output.trees", sep=""),
                              tree = tree_vault,
                              burnin = 0)

cophylotrees <- cophylo(tree_vault,BTraits_vault$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_vault.a <- cophylotrees$trees[[1]]

mytreebybranch(tree = tree_vault.a, x=log(BTraits_vault$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("vault rates")




#########
#Sphenoid region
#########

BTraits_sphenoid<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_sphenoid_MR.txt.VarRates.txt", sep=""),
                                 rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_sphenoid_MR.txt.Output.trees", sep=""),
                                 tree = tree_sphenoid,
                                 burnin = 0)

cophylotrees <- cophylo(tree_sphenoid,BTraits_sphenoid$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_sphenoid.a <- cophylotrees$trees[[1]]

mytreebybranch(tree = tree_sphenoid.a, x=log(BTraits_sphenoid$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("sphenoid rates")




#########
#Pterygoid region
#########

#######
BTraits_pterygoid<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_pterygoid_MR.txt.VarRates.txt", sep=""),
                                  rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_pterygoid_MR.txt.Output.trees", sep=""),
                                  tree = tree_pterygoid,
                                  burnin = 0)

cophylotrees <- cophylo(tree_pterygoid,BTraits_pterygoid$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_pterygoid.a <- cophylotrees$trees[[1]]

mytreebybranch(tree = tree_pterygoid.a, x=log(BTraits_sphenoid$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("pterygoid rates")



#########
#Quadrate region
#########

#######
BTraits_quadrate<-BTRTools::rjpp(rjlog = paste(resultsfolder,"Phylo_PC_SCORES_quadrate_MR.txt.VarRates.txt", sep=""),
                                 rjtrees = paste(resultsfolder, "Phylo_PC_SCORES_quadrate_MR.txt.Output.trees", sep=""),
                                 tree = tree_quadrate,
                                 burnin = 0)

cophylotrees <- cophylo(tree_quadrate,BTraits_quadrate$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_quadrate.a <- cophylotrees$trees[[1]]


mytreebybranch(tree = tree_quadrate.a, x=log(BTraits_sphenoid$data$meanRate)[-1], 
               mode="edges",type="fan",
               palette = color3,cex=.3,
               edge.width=3,
               show.tip.label= TRUE)
title("quadrate rates")

