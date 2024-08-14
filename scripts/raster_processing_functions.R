## RASTER PROCESSING FUNCTIONS ##

#includes:
#1. get_masks()
#2. apply_masks()
#3. set_zero()
#4. change_res()
#5. remove_outliers()
#6. interpolate_raster()
#7. as.im.Spatraster()
#8. as.df.Spatraster()

#Libraries
library(terra)
library(fields)
library(spatstat)

#--------------------------------------------------------------------------
#(1) Import a mask from a .png file
#--------------------------------------------------------------------------
# filepaths: filepaths to .png of mask (shaded regions to be removed)

get_masks <- function(filepaths){
  m <- rast(filepaths)
  mc <- classify(m[[1]], cbind(-Inf, 254, NA)) #make all shaded regions NA
  return(mc)
}

#--------------------------------------------------------------------------
#(2) Apply mask to raster
#--------------------------------------------------------------------------
# dems: terra raster to be masked
# masks: mask obtained by 'get_masks()'

apply_masks <- function(dems, masks){
  ext(masks) <- ext(dems)
  m <- masks
  dems_m <- terra::mask(x = dems, mask = m[[1]])
  return(dems_m)
}

#--------------------------------------------------------------------------
#(3) Transform raster coordinates to begin at origin
#--------------------------------------------------------------------------
# rast: terra raster

set_zero <- function(rast){
  xext <- ext(rast)$xmax - ext(rast)$xmin
  yext <- ext(rast)$ymax - ext(rast)$ymin
  
  terra::ext(rast) <- c(0, xext, 0, yext)
  
  return(rast)
}

#--------------------------------------------------------------------------
#(4) Change resolution of a raster
#--------------------------------------------------------------------------
# resol: resolution to change the raster to
# dem: terra raster

change_res <- function(resol, dem){
  dem2 <- dem 
  res(dem2) <- resol
  dem2 <- resample(dem, dem2)
  return(dem2)
}


#--------------------------------------------------------------------------
#(5) Remove outliers
#--------------------------------------------------------------------------
# Replace raster values that are beyond a specified threshold with NA
# dem: terra raster
# quantile: threshold for value replacement eg. 0.95

remove_outliers <- function(dem, quantile){
  q <- global(dem,  \(i) quantile(i, quantile, na.rm = T))
  q2 <- global(dem, \(i) quantile(i, (1-quantile), na.rm = T))
  map <- classify(dem, cbind(q, Inf, NA))[[1]]
  map <- classify(map, cbind(-Inf, q2, NA))[[1]]
  return(map)
}

#--------------------------------------------------------------------------
#(6) Interpolate raster
#--------------------------------------------------------------------------
# interpolate NAs in raster (e.g. masked cracks) using a fitted thin plate spline model (fields package)
# dem: terra raster
# aggregate_by: factor to aggregate raster by for fitting thin plate spline surface

interpolate_raster <- function(dem, aggregate_by){
  ra <- terra::aggregate(dem, aggregate_by, na.rm = T)
  xy <- data.frame(xyFromCell(ra, 1:ncell(ra)))
  v <- values(ra) 
  i <- !is.na(v)
  xy <- xy[i,]
  v <- v[i]
  #tps <- fields::fastTps(xy, v, aRange = 20)
  tps <- fields::Tps(xy, v)
  p <- terra::rast(dem)
  p <- interpolate(p, tps)
  return(p)
}

#--------------------------------------------------------------------------
#(7) Convert terra raster to spatstat image
#--------------------------------------------------------------------------
# Thanks to stack exchange (user)
# X: terra raster to be converted to spatstat image

as.im.SpatRaster <- function(X) {
  X <- X[[1]]
  rs <- terra::res(X)
  e <- as.vector(terra::ext(X))
  out <- list(
    v = as.matrix(X, wide=TRUE)[nrow(X):1, ],
    dim = dim(X)[1:2],
    xrange = e[1:2],
    yrange = e[3:4],
    xstep = rs[1],
    ystep = rs[2],
    xcol = e[1] + (1:ncol(X)) * rs[1] + 0.5 * rs[1],
    yrow = e[4] - (nrow(X):1) * rs[2] + 0.5 * rs[2],
    type = "real",
    units  = list(singular=units(X), plural=units(X), multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}

#--------------------------------------------------------------------------
#(8) Convert raster to data frame of downsampled coordinates and values
#--------------------------------------------------------------------------
# raster: terra raster
# every: factor to downsample by, e.g. every = 5, every 5th value kept in data frame

as.df.Spatraster <- function(raster, every){
  
  xy <- data.frame(xyFromCell(raster, 1:ncell(raster)))
  v <- values(raster) 
  print(length(v))
  i <- !is.na(v)
  xy <- xy[i,]
  newrows <- seq(1, length(xy$x), every)
  xy <- xy[newrows,]
  v <- v[i]
  v <- v[seq(1, length(v), every)]
  print(paste("new length:", length(v)))
  xy$v <- v
  return(xy)
}

#-------------------------------------------------------------------------------
# (9) fixed ecospat boyce
#-------------------------------------------------------------------------------
alt.ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                              PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman') {
  
  #### internal function calculating predicted-to-expected ratio for each class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }
  
  if (inherits(fit,"RasterLayer")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- raster::extract(fit, obs)
    }
    fit <- getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  #if (inherits(fit, "SpatRaster")) {
  #  if (is.data.frame(obs) || is.matrix(obs)) {
  #    obs <- terra::extract(fit, as.data.frame(obs),ID=FALSE)
  #  }
  #  fit <- terra::values(fit,na.rm=T)
  # }
  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  if(length(nclass)==1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}
