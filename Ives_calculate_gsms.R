#### IVESHEADIOMORPH SURFACE METRICS ####

#Libraries
library(terra)
library(plyr)
library(imager)
library(magick)
library(geodiv)

#--------------------------------------------------------------------------------------------------------
## (1) Import DEMS and mask
#--------------------------------------------------------------------------------------------------------

# Import .tif as rasters

filelist <- list.files(path = "D:/ives_dem", pattern = ".*.tif")
filepaths <- unlist(lapply("D:/ives_dem/", paste0, filelist))

ives_dems <- lapply(filepaths, terra::rast)

# Plot all rasters
#par(mfrow = c(5,6))
#lapply(ives_dems, plot)


# Import masks

maskfiles <- list.files(path = "D:/ives_masks", pattern = ".*.png")
maskfps <- unlist(lapply("D:/ives_masks/", paste0, maskfiles))

get_masks <- function(filepaths){
  m <- terra::rast(filepaths)
  mc <- classify(m, cbind(-Inf, 256, NA))
  mc <- mc[[1]]
  return(mc)
}

masks <- lapply(maskfps, FUN = get_masks)

# Apply masks

apply_masks <- function(ives, masks){
  ext(masks) <- ext(ives)
  m <- masks
  ives_m <- terra::mask(x = ives, mask = m[[1]])
  return(ives_m)
}

ives2 <- ives_dems[2:28] # not all have masks

ives_masked <- mapply(FUN = apply_masks, ives = ives2,  masks = masks)
ives_masked <- c(list(ives_dems[[1]]), ives_masked, list(ives_dems[[29]]))

# Plot masked rasters
#par(mfrow = c(5,6))
#lapply(ives_masked, plot)

#--------------------------------------------------------------------------------------------------------
## (2) Preprocessing of rasters
#--------------------------------------------------------------------------------------------------------

# Set scale from zero

set_zero <- function(rast){
  xext <- ext(rast)$xmax - ext(rast)$xmin
  yext <- ext(rast)$ymax - ext(rast)$ymin
  
  terra::ext(rast) <- c(0, xext, 0, yext)
  
  return(rast)
}

ives_masked <- lapply(ives_masked, set_zero)


# Convert to raster::raster and remove best fit polynomial surface

ives_masked <- lapply(ives_masked, raster)
ives_masked <- lapply(ives_masked, geodiv::remove_plane)


#---------------------------------------------------------------------------------------------------------
## (3) Export PNG images for PH analysis
#---------------------------------------------------------------------------------------------------------

# Create greyscale image
ims <- lapply(ives_masked, imager::as.cimg)

# Make png
ims <- lapply(ims, magick::image_read)
ims <- lapply(ims, magick::image_convert, "png")
str(ims)

# Export pngs
for(i in seq_along(ims)){
  
  outputfile <- paste0("D:/ives_png/ives_m_",i,".png")
  image_write(ims[[i]], path = outputfile, format = "png")
}

# Resize pngs

#pngs <- list.files(path = "E:/Pigeon Cove/Ivesheadiomorphs/ives_png", pattern = ".*.png")
#pngpath <- unlist(lapply("E:/Pigeon Cove/Ivesheadiomorphs/ives_png/", paste0, pngs))
#ims <- lapply(pngpath, magick::image_read)
#ims <-lapply(ims, magick::image_resize, "300x300")

#for(i in seq_along(ims)){

#  outputfile <- paste0("E:/Pigeon Cove/Ivesheadiomorphs/ives_png2/ives_m_",i,".png")
#  magick::image_write(ims[[i]], path = outputfile, format = "png")
#}


#------------------------------------------------------------------------------------------------------
## (4) Calculate surface metrics
#------------------------------------------------------------------------------------------------------

cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(geodiv)})

# List of Surface Metrics to calculate
m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
               'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
               'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')

m_list <- m_list[-c(16, 24)] # package errors at these positions (scl, stxr)

# Run over all rasters calculating surface metrics globally

metrics <- list()

for(i in 1:length(m_list)){
  metrics[[i]] <- parLapply(cl = cl, X = ives_masked, fun = m_list[[i]])
  cat(m_list[[i]], sep = "\n")
}

stopCluster(cl)

# Create dataframe and export as .csv

df <- lapply(metrics, ldply)
df <- as.data.frame(matrix(unlist(df), ncol = 25))
names(df) <- c('Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr', 'Sbi',
               'Sci', 'Ssk', 'Sku', 'Sds', 'Sfd', 'Srw', 'Srwi', 'Shw',
               'Std', 'Stdi', 'Svi', 'Ssc', 'Sv', 
               'Sp', 'Sk', 'Smean', 'Spk', 'Svk')

write.csv(df, "D:/ives_gsms.csv")
