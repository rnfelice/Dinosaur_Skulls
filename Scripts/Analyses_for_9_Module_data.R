library(geomorph)
library(hot.dots)

####Code for analysing 9 module dataset and calculating rates and disparity 
####Loading data####
#Load slid and procrustes aligned coordinate data for 9 module data set
nine.module.coords <- read.csv("/Dinosaur_Skulls/Data/Coordinate_data/nine_modules_GPA.csv",header= TRUE, row.names = 1)
nine.module.coords3d <- arrayspecs(nine.module.coords, p = 1189, k = 3) 

#load module definitions
module.id.nine <- read.csv("/Data/Coordinate_data/nine_modules_moduledefs_rightside.csv")
module.id.nine.1<-module.id.nine$region_number
#load phylogeny
traditional_tree<-read.tree("/Data/Trees/full_traditional_tree.tre")
#trim phylogeny to match data:

nine.module.tree <- name.check(traditional_tree, two.d.array(nine.module.coords3d))

ninemodule.tree <- drop.tip(fulltree,topandsides.nc1$tree_not_data)

####Rates and Disparity####
#calculate per module rates
evoratesgeomorph<-compare.multi.evol.rates(A = nine.module.coords3d , gp = module.id.nine.1, phy = ninemodules.tree, Subset = TRUE, iter = 9999) 
evoratesgeomorph
modnamesrates<-c("Maxilla","Premaxilla","Nasal","Frontal","Parietal","Postorbital", "Jugal+Quadratojugal", "Squamosal", "Lacrimal")
pvalmat<-as.matrix(evoratesgeomorph$pairwise.pvalue)
rownames(pvalmat)<-colnames(pvalmat)<-modnamesrates

###########per region variance
premax.points<-nine.module.coords3d[which(module.id.nine.1==4),,] %>% two.d.array(.)
max.points<-nine.module.coords3d[which(module.id.nine.1==3),,] %>% two.d.array(.)
nasal.points<-nine.module.coords3d[which(module.id.nine.1==5),,] %>% two.d.array(.)
frontal.points<-nine.module.coords3d[which(module.id.nine.1==6),,] %>% two.d.array(.)
parietal.points<-nine.module.coords3d[which(module.id.nine.1==7),,] %>% two.d.array(.)
postorbital.points<-nine.module.coords3d[which(module.id.nine.1==15),,] %>% two.d.array(.)
squamosal.points<-nine.module.coords3d[which(module.id.nine.1==17),,] %>% two.d.array(.)
jugal.points<-nine.module.coords3d[which(module.id.nine.1==16),,] %>% two.d.array(.)
lacrimal.points<-nine.module.coords3d[which(module.id.nine.1==19),,] %>% two.d.array(.)


premax.disp<-morphol.disparity(premax.points~1)#n=143
max.disp<-morphol.disparity(max.points~1) #n=140
nasal.disp<-morphol.disparity(nasal.points~1) #n=129
frontal.disp<-morphol.disparity(frontal.points~1) #n=137
parietal.disp<-morphol.disparity(parietal.points~1) #n=136
postorbital.disp<-morphol.disparity(postorbital.points~1) #n=95
squamosal.disp<-morphol.disparity(squamosal.points~1) #n=83
jugal.disp<-morphol.disparity(jugal.points~1) #n=186
lacrimal.disp<-morphol.disparity(lacrimal.points~1) #n=140


#### Fig 3 ####

landmark_rates<-hot.dots::per_lm_rates(shape.data = nine.module.coords3d, phy = ninemodules.tree)
landmark_variance <- hot.dots::per_lm_variance(shape.data = nine.module.coords3d)
  