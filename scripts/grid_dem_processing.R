## PRE-PROCESSING OF DEMS PRIOR TO TEXTURE ANALYSIS (WHOLE SURFACE) ##

# 1. Mask cracks
# 2. Crop DEMs based on filename
# 3. Mask fronds
# 4. Match resolution across grid square
# 5. Export processed grid squares
# 6. Merge grid squares and export


# Libraries
library(raster)
library(terra)
library(plyr)

source("./raster_processing_functions.R") # functions: 'get_masks()', 'apply_masks()', 'change_res()'

#-------------------------------------------------------------------------------
# (1) Mask cracks
#-------------------------------------------------------------------------------

#Import .tif files of DEMs as rasters
filelist <- list.files(path = "./data/grid_dem", pattern = ".*.tif")
filepaths <- unlist(lapply("./data/grid_dem/", paste0, filelist))

dems <- lapply(filepaths, terra::rast)

#Import .png files of masked cracks as rasters
maskfiles <- list.files(path = "./data/grid_dem_masks", pattern = ".*.png")
maskfps <- unlist(lapply("./data/grid_dem_masks/", paste0, maskfiles))

masks <- lapply(maskfps, FUN = get_masks)

#Mask DEMs
dems_m <- mapply(FUN = apply_masks, dems = dems,  masks = masks)

#Rescale no. 21 to be the same as the others
ext(dems_m[[21]]) <- ext(dems_m[[21]])[1:4]*1000
values(dems_m[[21]]) <- values(dems_m[[21]])*1000

#Check masks: Find skewing pts
#r <- dems_m[[36]]
#plot(classify(r, cbind(4, 10, 10))[[1]]) #very positive (>4) regions will appear as dark green
#plot(classify(r, cbind(-10, -4, 10))[[1]]) #very negative(<-4) regions will appear as dark green

#-------------------------------------------------------------------------------
# (2) Crop DEMs to grid based on filename
#-------------------------------------------------------------------------------

#Get filenames
filenames <- gsub(".tif", "", filelist)

#Get x and y values of the bottom left corner
get_x <- function(filename){
  x <- unlist(strsplit(filename, "[^[:digit:]]"))[2]
  return(x)
}

get_y <- function(filename){
  y <- unlist(strsplit(filename, "[^[:digit:]]"))[3]
  return(y)
}

x <- unlist(lapply(filenames, get_x))
y <- unlist(lapply(filenames, get_y)) 

add_decimal <- function(X){
  m <- unlist(strsplit(X, ""))
  ifelse(length(m)== 2, {X <- as.numeric(paste0(m[1], ".", m[2]))}, {X <- as.numeric(X)})
  return(X)
}

x <- unlist(lapply(x, add_decimal))*1000
y <- unlist(lapply(y, add_decimal))*1000

#Crop DEM based on filename
crop_to_grid <- function(dem, x, y){
  dem_c <- crop(dem, ext(x, x + 500, y, y + 500))
}

dems_c <- mapply(crop_to_grid, dems_m, x, y )

#-------------------------------------------------------------------------------
# (3) Mask fronds
#-------------------------------------------------------------------------------

#Import frond map as raster
fronds <- terra::rast("./data/pc_fronds.png")
fronds <- classify(fronds, cbind(-Inf, 254, NA))[[1]]

#change extent to 0, 4500, 0, 3500 by adding extra space
xmx <- xmax(fronds)
ymx <- ymax(fronds)

c1 <- rast(ext(0, xmx, ymx, 3500), resolution = res(fronds))
values(c1) <- 255

c2 <- rast(ext(xmx, 4500, 0, 3500), resolution = res(fronds))
values(c2) <- 255

fronds2 <- merge(fronds, c1, c2)

#Replicate frond raster for number of grid squares
fronds_r <- list()
for(i in 1:length(x)){
  fronds_r[[i]]<- fronds2
}


#Crop frond raster to grid cells
fronds_c <- mapply(crop_to_grid, fronds_r, x, y )

#Apply frond masks to crack masked gridded DEMs
apply_masks2 <- function(dems, masks){
  resmasks <- resample(masks, dems)
  ext(resmasks) <- ext(dems)
  m <- resmasks
  dems_m <- terra::mask(x = dems, mask = m)
  return(dems_m)
}

fronds_m <- mapply(FUN = apply_masks2, dems = dems_c,  masks = fronds_c)

#------------------------------------------------------------------------------
# (4) Make resolution uniform across grid-squares
#------------------------------------------------------------------------------

max_res <- max(unlist(lapply(dems_c, res))) #find the coarsest resolution

dems_fc <- lapply(fronds_m, change_res, resol = max_res)

#-------------------------------------------------------------------------------
# (5) Export processed DEM grid-squares
#-------------------------------------------------------------------------------

outputfiles <- paste0("./data/processed_grid/",filenames, ".tif")
mapply(writeRaster, dems_fc, outputfiles, overwrite = T)

#-------------------------------------------------------------------------------
# (6) Merge all rasters to get map of whole surface
#-------------------------------------------------------------------------------

dems_merged <- do.call(merge, dems_fc)
plot(dems_merged)

#Export
writeRaster(dems_merged, "./data/pc_all_dems.tif", overwrite = T)
plot(rast("./data/pc_all_dems.tif"))
