## POST-PROCESSING OF TEXTURE IMAGES ##

# 1. Import texture images
# 2. Remove outliers
# 3. OPTIONAL: Reduce resolution
# 4. OPTIONAL: Interpolate cracks
# 5. Export processed texture images as .tif

# Libraries 
library(terra)
library(geodiv)
library(raster)
library(spatialEco)

source("./raster_processing_functions.R") # functions: 'interpolate_raster()'

#-------------------------------------------------------------------------------
# (1) Read in texture images
#-------------------------------------------------------------------------------

filelist <- list.files("./data/texture_images_w2/", ".*.tif")
filenames <- gsub(".tif", "", filelist)
filepaths <- unlist(lapply("./data/texture_images_w2/", paste0, filelist))

maps <- lapply(filepaths, terra::rast)

#-------------------------------------------------------------------------------
# (2) Remove outliers 
#-------------------------------------------------------------------------------
maps95 <- lapply(maps, remove_outliers, 0.95)

#-------------------------------------------------------------------------------
# (3) Reduce resolution 
#-------------------------------------------------------------------------------

maps95 <- lapply(maps95, terra::aggregate, 20, na.rm = T)

#-------------------------------------------------------------------------------
# (3) OPTIONAL - interpolate cracks
#-------------------------------------------------------------------------------

#maps95 <- lapply(maps95, interpolate_raster, 5) # note that this fails for edge squares

#-------------------------------------------------------------------------------
# (4) OPTIONAL - smoothing
#-------------------------------------------------------------------------------

#maps95 <- lapply(maps95, spatialEco::raster.gaussian.smooth, type = "mean", na.rm = T, n = 5)

#-------------------------------------------------------------------------------
# (5) Export processed texture images as .tif
#-------------------------------------------------------------------------------

outputfiles <- paste0("./data/postprocessed_ti_w2/", filelist)
mapply(writeRaster, maps95, outputfiles, overwrite = T)

