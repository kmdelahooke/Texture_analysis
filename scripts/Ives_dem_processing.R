##  PRE-PROCESSING OF IVESHEADIOMORPH DEMS PRIOR TO TEXTURE ANALYSIS ##

# 1. Mask cracks
# 2. Set scale to zero
# 3. Make resolution consistent
# 4. Remove outliers
# 5. Remove plane
# 6. Export processed dems

# Libraries
library(imager)
library(magick)

source("./raster_processing_functions.R") # functions: 'get_masks()', 
#'apply_masks()', 'change_res()', 'remove_outliers()'

# Import .tif as rasters
filelist <- list.files(path = "./Pigeon Cove/Ivesheadiomorphs/ives_dem", pattern = ".*.tif")
filepaths <- unlist(lapply("./Pigeon Cove/Ivesheadiomorphs/ives_dem/", paste0, filelist))
filenames <- gsub(".tif", "", filelist)

ives_dems <- lapply(filepaths, terra::rast)
ext(ives_dems[23][[1]]) <- ext(ives_dems[23][[1]])*10  #need to check geomagic file to see how this has occurred
#ives_dems <- lapply(ives_dems, terra::aggregate, 10, na.rm = T)

#-------------------------------------------------------------------------------
# (1) Mask cracks
#-------------------------------------------------------------------------------

# Import masks
maskfiles <- list.files(path = "./data/ives_masks", pattern = ".*.png")
maskfps <- unlist(lapply("./data/ives_masks/", paste0, maskfiles))

masks <- lapply(maskfps, FUN = get_masks)

# Apply masks
ives2 <- ives_dems[2:28] # not all have masks

ives_masked <- mapply(FUN = apply_masks, ives = ives2,  masks = masks)
ives_masked <- c(list(ives_dems[[1]]), ives_masked, list(ives_dems[[29]]))

#-------------------------------------------------------------------------------
# (2) Set scale to zero
#-------------------------------------------------------------------------------

ives_masked <- lapply(ives_masked, set_zero)

#-------------------------------------------------------------------------------
# (3) Make DEM resolutions consistent
#-------------------------------------------------------------------------------

max_res <- max(unlist(lapply(ives_masked, res)))
print(paste('maximum resolution:', max_res))

ives_masked <- lapply(ives_masked, change_res, resol = max_res)

# Further resultion changes

#ives5 <- lapply(ives_masked, terra::aggregate, 5, na.rm = T)
#ives10 <- lapply(ives_masked, terra::aggregate, 10, na.rm = T)
#ives20 <- lapply(ives_masked, terra::aggregate, 20, na.rm = T)

#-------------------------------------------------------------------------------
# (4) Remove outliers
#-------------------------------------------------------------------------------

ives_masked <- lapply(ives_masked, remove_outliers, quantile = 0.99)

#-------------------------------------------------------------------------------
# (5) Remove plane
#-------------------------------------------------------------------------------

ives_masked <- lapply(ives_masked, geodiv::remove_plane) # n.b. if old geodiv, convert to raster first

#-------------------------------------------------------------------------------
# (6) Export processed dems as .tif and .png
#-------------------------------------------------------------------------------

# Export .tif for GSMs
#---------------------

outputfiles <- paste0("./processed_ives/", filenames, ".tif")
mapply(writeRaster, ives_masked, outputfiles, overwrite = T)


#  Export greyscale .png for PH
#------------------------------

# Create greyscale image
ims <- lapply(ives_masked, raster)
ims <- lapply(ims, imager::as.cimg)

# Make png
ims <- lapply(ims, magick::image_read)
ims <- lapply(ims, magick::image_convert, "png")
str(ims)

# Export pngs
for(i in seq_along(ims)){
  
  outputfile <- paste0("./data/processed_ives/ives0/ives_m_",i,".png")
  image_write(ims[[i]], path = outputfile, format = "png")
}


