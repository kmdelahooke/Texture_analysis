### SPECIES DISTRIBUTION MODELLING ###

# Libraries
library(raster)
library(dismo)
library(rgeos)
library(rJava)
library(tidyverse)
library(ecospat)

#------------------------------------------------------------------------------------
# (1) Read in and prepare data
#------------------------------------------------------------------------------------

# Create a Raster stack of first 3 principle components from TI_post_processing.R
pc_rasts <- lapply(pc_rasts, raster)
tex <- raster::stack(pc_rasts[1:3])

# Read in specimen locations
fronds <- read.csv("E:/pc_fronds23.csv")[-1]

# Split into testing and training datasets
frondocc <- cbind.data.frame(fronds$x, fronds$y)
frondocc$fold <- kfold(frondocc, k=5)
frondtest <- frondocc[fold ==1, 1:2]
frondtrain <- frondocc[fold != 1, 1:2]

#-------------------------------------------------------------------------------------
# (2) Fit MaxEnt model
#-------------------------------------------------------------------------------------

me <- maxent(tex, frondtrain)

# Habitat suitibility map
r <- predict(me, tex)
plot(r)
points(frondocc)

#-------------------------------------------------------------------------------------
# (3) Testing the model
#-------------------------------------------------------------------------------------

# Response curves
plot(me)
response(me)


# Mean AUC
auc <- list()
for(i in 1:100){
  set.seed(sample(1:100000, 1))
  bg <- randomPoints(tex, 57) # background pseudoabsences
  e1 <- evaluate(me, p=frondtest, a=bg, x=tex)
  auc[[i]]<- e1@auc
  
}

summary(unlist(auc))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.6029  0.7010  0.7177  0.7181  0.7388  0.8070

# Boyce model 
obs = cbind(x = fronds$x, y = fronds$y)
obs <- obs[-32,]

ecospat.boyce(r, obs) ## good
## write as loop like above#

#------------------------------------------------------------------------------------

## bioclim ## - NPS, not properly evaluated

bioc1 <- bioclim(tex, frondtrain)
pairs(bioc1, pa = "p")
hsm1 <- predict(tex, bioc1)
plot(hsm1)

par(mfrow=c(1, 2))
plot(hsm1)
plot(r)

diffs <- hsm1 - r
plot(abs(diffs))

ecospat.boyce(hsm1, obs) # bioclim is shit, and probably overly sensitive
