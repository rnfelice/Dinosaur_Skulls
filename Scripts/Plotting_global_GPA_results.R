#Summarizing the global GPA results:
library(tidyverse)
library(ape)
library(ggpubr)
library(BTprocessR)
library(phytools)



#load time trees:
tree_rostrum<-read.nexus("./Data/tree_rostrum.nex")
tree_vault<-read.nexus("./Data/tree_vault.nex")
tree_occipital<-read.nexus("./Data/tree_occiput.nex")
tree_quadrate<-read.nexus("./Data/tree_quadrate.nex")

  
#summarize the posterior of the bayestraits reversable jump mcmc 
#NOTE: these files are archived in .zip format and need to be extracted before use
#NOTE: this is very slow!
#First we need to summarize one of each output tree:
rostrum_rjpp <- rjpp(rjlog = "/Data/Phylo_PC_SCORES_rostrum_full.txt.VarRates.txt",
                     rjtrees =  "/Data/Phylo_PC_SCORES_rostrum_full.txt.Output.trees",
                     tree=tree_rostrum, burnin = 0,verbose = FALSE)

vault_rjpp <- rjpp(rjlog = "/Data/Phylo_PC_SCORES_vault_full.txt.VarRates.txt",
                   rjtrees =  "/Data/Phylo_PC_SCORES_vault_full.txt.Output.trees",
                   tree=tree_vault, burnin = 0,verbose = FALSE)

occipital_rjpp <- rjpp(rjlog = "/Data/Phylo_PC_SCORES_occiput_full.txt.VarRates.txt",
                       rjtrees =  "/Data/Phylo_PC_SCORES_occiput_full.txt.Output.trees",
                       tree=tree_occiput, burnin = 0,verbose = FALSE)

quadrate_rjpp <- rjpp(rjlog = "/Data/Phylo_PC_SCORES_quadrate_full.txt.VarRates.txt",
                      rjtrees =  "/Data/Phylo_PC_SCORES_quadrate_full.txt.Output.trees",
                      tree=tree_quadrate, burnin = 0,verbose = FALSE)

#load workspace:
load("/Data/global gpa rjpp outs4.RData")

#load the taxonomic data
lookup <- read_csv("/Data/taxonomy.lookup.table.csv")

#load the rate-scaled trees from the posterior
post_trees_rostrum <- read.nexus(paste(resultsfolder, "Phylo_PC_SCORES_rostrum_Tnew.txt.Output.trees", sep=""))
post_trees_vault<- read.nexus(paste(resultsfolder, "Phylo_PC_SCORES_vault_Tnew.txt.Output.trees", sep=""))
post_trees_occipital <- read.nexus(paste(resultsfolder, "Phylo_PC_SCORES_occipital_Tnew.txt.Output.trees", sep=""))
post_trees_quadrate <- read.nexus(paste(resultsfolder, "Phylo_PC_SCORES_quadrate_Tnew.txt.Output.trees", sep=""))

###
#Define a function to calculate mean rate scalars

rate_summary <- function(rjpp_out, time_tree, lookup.table){
  #select the species from  the lookup table that are part of the current dataset
  current.clades <- filter(lookup, newnames %in% rjpp_out$meantree$tip.label)
  #select the sub tree for each group
  bird_time_tree <- keep.tip(time_tree, current.clades$newnames[which(current.clades$group.id2=="bird")])
  thero_time_tree <- keep.tip(time_tree, current.clades$newnames[which(current.clades$group.id2=="theropod")])
  nonthero_time_tree <- keep.tip(time_tree, current.clades$newnames[which(current.clades$group.id2=="other")])
  #get the time represented by each sub tree
  bird_tree_time <- sum(bird_time_tree$edge.length)
  thero_tree_time <- sum(thero_time_tree$edge.length)
  nonthero_tree_time <- sum(nonthero_time_tree$edge.length)
  #make blank vectors for the rates and scaled rates in each group
  Mean_Rate_Bird<-Mean_Scaled_Rate_Bird<-Mean_Rate_Thero<-Mean_Scaled_Rate_Thero<-Mean_Rate_NonThero<-Mean_Scaled_Rate_NonThero<-c()
  #select just the posterior distribution of rate scalars from the BayesTraits output
  #the first column is the rate on the root and our time tree has no root edge so this is removed with the "[-1,]"
  rate_table <- rjpp_out$scalars$rates[-1,]
  #make copy the mean tree x times where x is the number of trees in the posterior distribution
  treelist<-rep(list(rjpp_out$meantree),dim(rate_table)[2])
  #now convert the edge lengths of each of those trees to be the rate scalar for that edge
  #this allows us to keep rates associated with the appropriate edges as we go to the next step
  for (k in 1:length(treelist)){
    treelist[[k]]$edge.length<-rate_table[,k]
  }
  #now subset those trees based on group identity
  for (i in 1:length(treelist)){
    bird_tree <- keep.tip(treelist[[i]], current.clades$newnames[which(current.clades$group.id2=="bird")])
    thero_tree <- keep.tip(treelist[[i]], current.clades$newnames[which(current.clades$group.id2=="theropod")])
    nonthero_tree <- keep.tip(treelist[[i]], current.clades$newnames[which(current.clades$group.id2=="other")])
    #get the mean of each group
    mean_rate_bird <- mean(bird_tree$edge.length)
    #and scale that mean to "tree_time", AKA the time represented that sub tree
    mean_rate_bird_scaled <- mean_rate_bird/bird_tree_time
    mean_rate_thero <- mean(thero_tree$edge.length)
    mean_rate_thero_scaled <- mean_rate_thero/thero_tree_time
    mean_rate_nonthero <- mean(nonthero_tree$edge.length)
    mean_rate_nonthero_scaled <- mean_rate_nonthero/nonthero_tree_time
    Mean_Rate_Bird <- c(Mean_Rate_Bird,  mean_rate_bird)
    Mean_Scaled_Rate_Bird <- c(Mean_Scaled_Rate_Bird, mean_rate_bird_scaled)
    Mean_Rate_Thero <- c(Mean_Rate_Thero, mean_rate_thero)
    Mean_Scaled_Rate_Thero <- c(Mean_Scaled_Rate_Thero, mean_rate_thero_scaled)
    Mean_Rate_NonThero <- c(Mean_Rate_NonThero, mean_rate_nonthero)
    Mean_Scaled_Rate_NonThero <- c(Mean_Scaled_Rate_NonThero, mean_rate_nonthero_scaled)
  }
  #comnbine into two tables, one for unscaled one for scaled 
  results <- tibble(Mean_Rate_Bird, Mean_Rate_Thero,  Mean_Rate_NonThero) %>% pivot_longer(cols = starts_with("Mean_Rate_"),
                                                                                           names_to = "Group",
                                                                                           names_prefix = "Mean_Rate_",
                                                                                           values_to = "MeanRate")
  results_scaled <- tibble(Mean_Scaled_Rate_Bird, Mean_Scaled_Rate_Thero, Mean_Scaled_Rate_NonThero)%>% pivot_longer(cols = starts_with("Mean_Scaled_Rate_"),
                                                                                                                     names_to = "Group",
                                                                                                                     names_prefix = "Mean_Scaled_Rate_",
                                                                                                                     values_to = "MeanRate")
  
  return(list(results, results_scaled))
}



#execute this function for each module
rostrum_res <- rate_summary(rjpp_out = rostrum_rjpp, time_tree =tree_rostrum, lookup.table= lookup) 
# and rename the group labels for prettier plotting
rostrum_res2 <- rostrum_res[[2]] %>% mutate(., Group = recode(Group, Bird = "Birds", Thero = "Non-Avian\nTheropods", NonThero = "Non-Theropod\nDinosaurs"))
occipital_res <- rate_summary(rjpp_out = occipital_rjpp, time_tree =tree_occipital, lookup.table= lookup)
occipital_res2 <- occipital_res[[2]] %>% mutate(., Group = recode(Group, Bird = "Birds", Thero = "Non-Avian\nTheropods", NonThero = "Non-Theropod\nDinosaurs"))
quadrate_res <- rate_summary(rjpp_out = quadrate_rjpp, time_tree =tree_quadrate, lookup.table= lookup)
quadrate_res2 <- quadrate_res[[2]] %>% mutate(., Group = recode(Group, Bird = "Birds", Thero = "Non-Avian\nTheropods", NonThero = "Non-Theropod\nDinosaurs"))
vault_res <- rate_summary(rjpp_out = vault_rjpp, time_tree =tree_vault, lookup.table= lookup)
vault_res2 <- vault_res[[2]] %>% mutate(., Group = recode(Group, Bird = "Birds", Thero = "Non-Avian\nTheropods", NonThero = "Non-Theropod\nDinosaurs"))


#plot
my_comparisons <- list( c("Birds", "Non-Avian\nTheropods"), c("Birds", "Non-Theropod\nDinosaurs"), c("Non-Avian\nTheropods", "Non-Theropod\nDinosaurs") )
g.rostrum <- ggboxplot(rostrum_res2, x = "Group", y = "MeanRate", fill = "Group", palette = "jco",
                       ylab = "Mean Rate Scalar") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 25,hjust = 1),
        text=element_text(family = "Arial"))+
  ggtitle("Rostrum")
g.occiput <- ggboxplot(occipital_res2, x = "Group", y = "MeanRate", fill = "Group", palette = "jco",
                       ylab = "Mean Rate Scalar")+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 25,hjust = 1),
        text=element_text(family = "Arial"))+
  ggtitle("Occipital")
g.quadrate <- ggboxplot(quadrate_res2, x = "Group", y = "MeanRate", fill = "Group", palette = "jco",
                        ylab = "Mean Rate Scalar")+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 25,hjust = 1),
        text=element_text(family = "Arial"))+
  ggtitle("Quadrate")
g.vault <- ggboxplot(vault_res2, x = "Group", y = "MeanRate", fill = "Group", palette = "jco",
                     ylab = "Mean Rate Scalar")+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 25,hjust = 1),
        text=element_text(family = "Arial"))+
  ggtitle("Vault")

library(patchwork)
g.combined_global<-(g.rostrum | g.quadrate | g.vault | g.occiput)+ plot_layout(guides = "collect") & theme(legend.position = 'bottom',
                                                                              axis.title.x = element_blank(),
                                                                              axis.ticks.x = element_blank(),
                                                                              axis.text.x = element_text(size=9),
                                                                              plot.title = element_text(vjust=0.5))
g.combined_global
