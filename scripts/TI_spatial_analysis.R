## TEXTURE IMAGE SPATIAL ANALYSIS ##

# 1. Import data: point pattern, covariates: texture image principal component maps, DEM
# 2. Autocorrelation of covariates
# 3. Isotropy of covariates
# 4. Parametric heterogenous Poisson model
# 5. Maximum entropy model

# Libraries
library(spatstat)
library(terra)
library(tidyverse)
library(pgirmess)
library(geoR)
library(dismo)
library(rgeos)
library(rJava)
library(tidyverse)
library(ecospat)

source("./raster_processing_functions.R") # functions: 'as.im.SpatRaster()', 'as.mat.SpatRaster()', 'remove_outliers()'
source("./TI_ordination.R") #data: pc_rasts

#-------------------------------------------------------------------------------
# (1) Import data
#-------------------------------------------------------------------------------

# A. frond point pattern
df <- read.csv("./data/pc_fronds23_all.csv")[, -1]

# B. principal component maps
pc_rasts <- lapply(pc_rasts, terra::rast) # output from ordination

#smoothing 
pc_rasts2 <- lapply(pc_rasts, spatialEco::raster.gaussian.smooth, type = "mean", na.rm = T, n = 7) # varied n

# C. DEM
pc_dem <- terra::rast("./data/pc_all_dems.tif")
pc_dem <- terra::aggregate(pc_dem, 20, na.rm = T) # equivalent to processing of principal component maps
pc_dem <- remove_outliers(pc_dem, 0.999) 


# D. Height map
pc_h <- terra::rast("./data/lidar_c3.tif")

# mask with using texture image outline
r <- classify(pc_rasts[[1]], cbind(-Inf, Inf, 256))
ext(pc_h) <- ext(r)
pc_h <- change_res(max(res(r)), pc_h)
pc_h2 <- mask(pc_h, r)

#smooth
pc_h3 <- lapply(pc_h2, spatialEco::raster.gaussian.smooth, type = "mean", na.rm = T, n = 7) 

# E. convert rasters to spatstat image

pc_ims <- lapply(pc_rasts2, as.im.SpatRaster) # convert to spatstat format
names(pc_ims) <- c("pc1", "pc2", "pc3", "pc5", "pc6")

dem_im <- list(as.im.SpatRaster(pc_dem))
names(dem_im) <- "dem"

h_im <- list(as.im.SpatRaster(pc_h3))
names(h_im) <- "height"

# F. data.frames of sampled rasters

df_dem <- as.df.Spatraster(pc_dem, 30)
df_pc <- lapply(pc_rasts, as.df.Spatraster, 30)
df_pc2 <- lapply(pc_rasts2, as.df.Spatraster, 30)

#-------------------------------------------------------------------------------
# (2) Spatial autocorrelation
#-------------------------------------------------------------------------------

# Maximum distance to consider
coords <- crds(aggregate(pc_dem, 10))
distmat <- as.matrix(dist(coords))
maxdist <- 2/3*max(distmat)

# Morans's correlogram

clg <- pgirmess::correlog(df_dem[,1:2], df_dem[,3], nbclass = 30)
clg_1 <- pgirmess::correlog(df_pc2[[1]][,1:2], df_pc2[[1]][,3], nbclass = 30)
clg_2 <- pgirmess::correlog(df_pc2[[2]][,1:2], df_pc2[[2]][,3], nbclass = 30)
clg_3 <- pgirmess::correlog(df_pc2[[3]][,1:2], df_pc2[[3]][,3], nbclass = 30)

# Plot correlog: function from pgirmess, edited to account for bonferroni correction
bpv <- 0.05/30
plot.correlog<-function (x,type,xlab,ylab,main,...) {
  if (!inherits(x, "correlog")) stop("Object must be of class 'correlog'")
  if (missing(main)) main<-paste(attributes(x)$Method," = f(distance classes)",sep="")
  if (missing(type)) type<-"b"
  if (missing(ylab)) ylab<-attributes(x)$Method
  if (missing(xlab)) xlab<-"distance classes"
  plot(x[,1:2,drop=FALSE],type=type,xlab=xlab,ylab=ylab,main=main,xaxt="n",...)
  inc<-(x[2,1]-x[1,1])/2
  breaks <- pretty(c(x[1,1]-inc,x[length(x[,1]),1]+inc), n = length(x[,1]), min.n = 2)
  axis(1,at=breaks,...)
  points(x[x[,3]<bpv,1:2,drop=FALSE],pch=19,col="red",cex=2)
}

plot.correlog(clg, xlim = c(0, maxdist), main = "DEM") %>% abline(h = 0)
plot.correlog(clg_1, xlim = c(0, maxdist), main = "PC1") %>% abline(h = 0)
plot.correlog(clg_2, xlim = c(0, maxdist), main = "PC2") %>% abline(h = 0)
plot.correlog(clg_3, xlim = c(0, maxdist), main = "PC3") %>% abline(h = 0)

#-------------------------------------------------------------------------------
# (3) Isotropy - partial correlogram by angle categories
#-------------------------------------------------------------------------------

# Convert data to geoR format
geo.dem <- geoR::as.geodata(df_dem)
geo.pc <- lapply(df_pc2, geoR::as.geodata)

# partial correlogram split by angle

clgi <- geoR::variog4(geo.dem, maxdist)
clgi_pc <- lapply(geo.pc, geoR::variog4, maxdist)

# plot correlograms
geoR::plot(clgi, lwd = 2, max.dist = maxdist)
abline(v = maxdist)

lapply(clgi_pc, plot, lwd = 2, max.dist = maxdist)
abline(v = maxdist)

#-------------------------------------------------------------------------------
# (3) Parametric heterogeneous Poisson model
#-------------------------------------------------------------------------------

# point pattern
r <- classify(pc_rasts2[[1]], cbind(-Inf, Inf, 256))
win_pc <- as.owin(as.im.SpatRaster(r))
pp <- ppp(df$x, df$y, window = win_pc)
den <- density(pp, kernel = 'epanechnikov')

# heterogeneous poisson models
hpp1 <- ppm(pp ~ pc3, covariates = pc_ims)
hpp2 <- ppm(pp ~ pc1 + pc2 + pc3, covariates = pc_ims)
hpp3 <- ppm(pp ~ den + pc1 + pc2 + pc3, covariates = pc_ims)
hpp4 <- ppm(pp ~ den)
hpp5 <- ppm(pp ~ den + pc3, covariates = pc_ims)
hpp6 <- ppm(pp ~ den* pc3, covariates = pc_ims)
hpp7 <- ppm(pp ~ pc1 + pc3, covariates = pc_ims)
hpp8 <- ppm(pp ~ pc1 + pc2, covariates = pc_ims)
hpp9 <- ppm(pp ~ pc1 , covariates = pc_ims)
hpp10 <- ppm(pp ~ pc2, covariates = pc_ims)
hpp11 <- ppm(pp ~ pc3 + pc2, covariates = pc_ims)
hpp12 <- ppm(pp ~ pc1 + den, covariates = pc_ims)
hpp13 <- ppm(pp ~ pc2 + den, covariates = pc_ims)
hpp14 <- ppm(pp ~ pc1*pc2*pc3*den, covariates = pc_ims)
hpp15 <- ppm(pp ~ pc1*den, covariates = pc_ims)
hpp16 <- ppm(pp ~ pc2*den, covariates = pc_ims)
hpp17 <- ppm(pp ~ pc1*pc2*den, covariates = pc_ims)
hpp18 <- ppm(pp ~ pc1*pc3*den, covariates = pc_ims)
hpp19 <- ppm(pp ~ pc2*pc3*den, covariates = pc_ims)
hpp20 <- ppm(pp ~ pc1*pc2, covariates = pc_ims)
hpp21 <- ppm(pp ~ pc1*pc3, covariates = pc_ims)
hpp22 <- ppm(pp ~ pc2*pc3, covariates = pc_ims)

# model selection
models <- list(hpp1, hpp2, hpp3, hpp4, hpp5, hpp6, hpp7, hpp8, hpp9, hpp10, hpp11, hpp12, hpp13, hpp14, hpp15, hpp16, hpp17, hpp18, hpp19, hpp20, hpp21, hpp22)
aic.values <- lapply(models, AIC)

write.csv(unlist(aic.values), "./results/aic.csv")

# model evaluation
plot(predict.ppm(hpp5))
plot(pp, add = T)
diagnose.ppm(hpp6)

#-------------------------------------------------------------------------------
# (4) Maximum entropy model
#-------------------------------------------------------------------------------

pc_rasts2 <- lapply(pc_rasts, spatialEco::raster.gaussian.smooth, type = "mean", na.rm = T, n = 15)

# Create a raster stack of first 3 principle components from texture_ordination.R
pc_rasts2 <- lapply(pc_rasts2, raster::raster)
tex15 <- raster::stack(pc_rasts2[[1]], pc_rasts2[[2]], pc_rasts2[[3]])
plot(tex15)

# Specimen locations
fronds <- df

# Split into testing and training datasets
frondocc <- cbind.data.frame(fronds$x, fronds$y)
frondocc$fold <- kfold(frondocc, k=5)

frondtest <- subset(frondocc, fold == 1)[1:2]
frondtrain <- subset(frondocc, fold != 1)[1:2]


# Fit MaxEnt model
me <- maxent(tex, frondtrain)

# Habitat suitability map
par(mfrow = c(1,1))
hs <- predict(me, tex)
plot(hs)
points(frondocc, pch = 20, cex = 1.5, col = "chartreuse3")

# Response curves
plot(me)
response(maxent.beta3)
str(me)

# Mean AUC
auc <- list()
for(i in 1:100){
  set.seed(sample(1:100000, 1))
  bg <- randomPoints(tex, 68) # background pseudoabsences
  e1 <- dismo::evaluate(me, p=frondtest, a=bg, x=tex)
  auc[[i]]<- e1@auc
}

summary(unlist(auc))

# Boyce index
obs <- cbind(x = fronds$x, y = fronds$y)
raster::extract(pc_rasts[[1]], obs)
obs <- obs[-c(31,34),]

alt.ecospat.boyce(hs, obs)$cor

#get over a number of test and training splits
split.fronds <- function(df){
  fronds <- df
  
  frondocc <- cbind.data.frame(fronds$x, fronds$y)
  frondocc$fold <- kfold(frondocc, k=5)
  
  frondtest <- subset(frondocc, fold == 1)[1:2]
  frondtrain <- subset(frondocc, fold != 1)[1:2]
  return(list(frondtest, frondtrain))
}

eval <- list()
for(i in 1:100){
  set.seed(sample(100000,1))
 
  tt <- split.fronds(df)
  me <- maxent(tex, frondtrain)
  hs <- predict(me, tex)
  
  me5 <- maxent(tex5, frondtrain)
  hs5 <- predict(me5, tex5)
  
  me15 <- maxent(tex15, frondtrain)
  hs15 <- predict(me15, tex15)
  
  auc <- list()
  for(j in 1:100){
    set.seed(sample(100000,1))
    bg <- randomPoints(tex, 68) # background pseudoabsences
    e1 <- dismo::evaluate(me, p=frondtest, a=bg, x=tex)
    e5 <- dismo::evaluate(me5, p=frondtest, a=bg, x=tex5)
    e15 <- dismo::evaluate(me15, p=frondtest, a=bg, x=tex15)
    auc[[j]]<- c(e1@auc,e5@auc, e15@auc)
  }
  auc_df <- plyr::ldply(auc)
  mauc <- mapply(mean, auc_df)
 
  
  obs <- cbind(x = fronds$x, y = fronds$y)
  obs <- obs[-c(31,34),]
  
  boyce1 <- alt.ecospat.boyce(hs, obs)$cor
  boyce5 <- alt.ecospat.boyce(hs5, obs)$cor
  boyce15 <- alt.ecospat.boyce(hs15, obs)$cor
  
  boyce <-c(boyce1, boyce5, boyce15)
  ev <- rbind(mauc, boyce)
  
  eval[[i]] <- ev

}

eval_df <- ldply(eval)
eval_auc <- eval_df[ c(T,F), ]
eval_boyce <- eval_df[ c(F,T), ]

mapply(summary, eval_auc)
summary(eval_df$V2)
hist(eval_df$V1)
sd(eval_df$V1)

#-------------------------------------------------------------------------------
# (5) Impact of dip towards the coast
#-------------------------------------------------------------------------------
hpp1 <- ppm(pp ~ height, covariates = h_im)
dclf.test(hpp1)
cdf.test(hpp1, h_im[[1]])
diagnose.ppm(hpp1)
plot(h_im[[1]])
plot(pp, pch = 20, cols = "skyblue", add = T)



