## CALCULATE SURFACE METRICS VIA A MOVING WINDOW ##

#libraries
library(terra)
library(geodiv)
library(raster)

#Read in .tiff files as rasters
filelist <- list.files(path = "D:/processed_grid/", pattern = ".*.tif")
filepaths <- unlist(lapply("D:/processed_grid/", paste0, filelist))
dems <- lapply(filepaths, terra::rast)

#-------------------------------------------------------------------------------
# (Option 1) calculate every metric for each square
#-------------------------------------------------------------------------------

create_texture_image <- function(dem, filename){
  
  dem <- raster::raster(dem)
  
  #remove trends
  dem <- remove_plane(dem)
  
  #List metrics to calculate
  m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
                 'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
                 'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')
  
  m_list <- m_list[-c(16, 24)] # package errors at these positions (scl, stxr)
  
  #calculate metrics
  outrasts <- list()
  
  system.time(for (i in 20:length(m_list)) {
    outrasts[[i]] <- texture_image(dem, window_type = 'square', 
                                   size = 2, in_meters = FALSE, 
                                   metric = m_list[[i]], parallel = TRUE,
                                   nclumps = 100)})

  outrasts <- stack(unlist(outrasts))
  
  data <- data.frame(x = coordinates(outrasts)[, 1], 
                         y = coordinates(outrasts)[, 2])
  
  for (i in 1:25) {
    data[, i + 2] <- outrasts[[i]][]
  }
  names(data) <- c('x', 'y', 'Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr', 'Sbi',
                       'Sci', 'Ssk', 'Sku', 'Sds', 'Sfd', 'Srw', 'Srwi', 'Shw',
                       'Std', 'Stdi', 'Svi', 'Ssc', 'Sv', 
                       'Sp', 'Sk', 'Smean', 'Spk', 'Svk')
  
  write.csv(data, paste0("D:/grid_tex_images/", filename, ".csv"))
  print("Finished one square!")
  

  return(data)
  
}


filenames <- gsub(".tif", "", filelist)

tex_images<- mapply(create_texture_image, dems, filenames)

#-------------------------------------------------------------------------------
# (Option 2) Run over all squares for an individual metric
#-------------------------------------------------------------------------------

create_texture_image2 <- function(dem, metric){
  
  dem <- raster::raster(dem)
  
  #remove trends
  dem <- remove_plane(dem)
  
  #calculate metric
  outrast <- texture_image(dem, window_type = 'square', 
                                   size = 2, in_meters = FALSE, 
                                   metric = metric, parallel = TRUE,
                                   nclumps = 100)
  outrast <- terra::rast(outrast)
  return(outrast)
  
}

#List metrics to calculate
m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
               'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
               'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')

m_list <- m_list[-c(16, 24)] # package errors at these positions (scl, stxr)


for(i in 2: length(m_list)){
  
  outrasts <- lapply(dems, create_texture_image2, m_list[[i]])
  
  #merge rasters
  rasts_merged <- do.call(merge, outrasts)
  plot(rasts_merged, main = m_list[[i]])
  
  #export
  writeRaster(rasts_merged, paste0("D:/pc_",m_list[[i]],".tif"), overwrite = T)
}




