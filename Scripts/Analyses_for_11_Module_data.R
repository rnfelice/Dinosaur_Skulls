library(geomorph)
library(hot.dots)


####Code for analysing 9 module dataset and calculating rates and disparity 
####Loading data####
#Load slid and procrustes aligned coordinate data for 9 module data set
eleven.module.coords <- read.csv("/Dinosaur_Skulls/Data/Coordinate_data/eleven_modules_GPA.csv",header= TRUE, row.names = 1)
eleven.module.coords3d <- arrayspecs(eleven.module.coords, p = 1515, k = 3) 

#load module definitions
module.id.eleven <- read.csv("/Data/Coordinate_data/eleven_modules_moduledefs_rightside.csv")
module.id.eleven.1<-module.id.eleven$region_number
#load phylogeny
traditional_tree<-read.tree("/Data/Trees/full_traditional_tree.tre")
#trim phylogeny to match data:

eleven.module.tree <- name.check(traditional_tree, two.d.array(eleven.module.coords3d))

elevenmodule.tree <- drop.tip(fulltree,topandsides.nc1$tree_not_data)

####Rates and Disparity####
#calculate per module rates
evoratesgeomorph<-compare.multi.evol.rates(A = eleven.module.coords3d , gp = module.id.eleven.1, phy = elevenmodules.tree, Subset = TRUE, iter = 9999) 
evoratesgeomorph
modnamesrates<-c("Maxilla","Premaxilla","Nasal","Frontal","Parietal","Postorbital", "Jugal+Quadratojugal", "Squamosal", "Lacrimal")
pvalmat<-as.matrix(evoratesgeomorph$pairwise.pvalue)
rownames(pvalmat)<-colnames(pvalmat)<-modnamesrates

###########per region variance
premax.points<-eleven.module.coords3d[which(module.id.eleven.1==4),,] %>% two.d.array(.)
max.points<-eleven.module.coords3d[which(module.id.eleven.1==3),,] %>% two.d.array(.)
nasal.points<-eleven.module.coords3d[which(module.id.eleven.1==5),,] %>% two.d.array(.)
frontal.points<-eleven.module.coords3d[which(module.id.eleven.1==6),,] %>% two.d.array(.)
parietal.points<-eleven.module.coords3d[which(module.id.eleven.1==7),,] %>% two.d.array(.)
postorbital.points<-eleven.module.coords3d[which(module.id.eleven.1==15),,] %>% two.d.array(.)
squamosal.points<-eleven.module.coords3d[which(module.id.eleven.1==17),,] %>% two.d.array(.)
jugal.points<-eleven.module.coords3d[which(module.id.eleven.1==16),,] %>% two.d.array(.)
lacrimal.points<-eleven.module.coords3d[which(module.id.eleven.1==19),,] %>% two.d.array(.)


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
#calculate per-landmark rates and variance (Fig 3)

landmark_rates<-hot.dots::per_lm_rates(shape.data = eleven.module.coords3d, phy = elevenmodules.tree)
landmark_variance <- hot.dots::per_lm_variance(shape.data = eleven.module.coords3d)




###Fig 4###

#load in the mean shape from the procustes alignment:
eleven_module_meanshape<-read.csv("/Data/Coordinate_data/eleven_module_meanshape.csv")

library(devtools)
install_github("TGuillerme/landvR")
library(landvR)

n.points=1515
differences_from_mean <- coordinates.difference(coordinates = eleven.module.coords3d,
                                                reference = eleven_module_meanshape,
                                                type = "vector",
                                                absolute.distance= FALSE)

get.color.spectrum<-function(data, nbin=100,limits=NULL, colorlist = c("#6e016b","#0c2c84","#225ea8","#005a32","#ffff00","#fe9929","#fc4e2a","red")){
  x=data
  xlims<-NULL
  tol <- 1e-06
  cols1<-colorRampPalette(colorlist)
  cols<-cols1(100)
  if (is.null(limits)==TRUE){
    xlims <- range(x) + c(-tol, tol)}
  else{
    xlims <- limits}
  nbin=nbin
  breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
        1
    cols[i]
  }
  spec1colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
}

#open an rgl window, plot the mean specimen and align it in the position you want
open3d()
spheres3d(eleven_modules_aligned_bilateral$consensus[c(1:n.points),],radius = .002)


#save rgl window parameters

zoom<-par3d()$zoom
userMatrix<-par3d()$userMatrix
windowRect<-par3d()$windowRect

folderA<-"~/outputfolder/"
for (i in 1:length(differences_from_mean)){
  datcol1<-get.color.spectrum(differences_from_mean[[i]][,1])
  name<-names(differences_from_mean)[i]
  open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
  spheres3d(dinos.only.f.v.occ.jj.bilateral.GPA.2$consensus[c(1:n.points),],radius = .002,col=datcol1)
  rgl.snapshot(filename = paste(folderA, names(differences_from_mean)[i],".png",sep=""))
  rgl.close()
}

  