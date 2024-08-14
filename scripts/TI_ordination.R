## ORDINATION OF TEXTURE IMAGES ##

# Follows the Smith et al., 2021 vignette to create maps of principal components

# 1. Input: processed texture images
# 2. Clean data
# 3. Calculate principal components
# 4. Create maps of principal components
# 5. Export maps  as .tif
# 6. Evaluate PCA

# Libraries:
library(geodiv)
library(terra)
library(raster)
library(rasterVis)
library(tidyverse)

#-------------------------------------------------------------------------------
# (1) Read in processed texture images and convert to data.frame
#-------------------------------------------------------------------------------

# A. Read in Texture Images
filelist <- list.files("./data/postprocessed_ti_w2/", ".*.tif")
filenames <- gsub(".tif", "", filelist)
filepaths <- unlist(lapply("./data/postprocessed_ti_w2/", paste0, filelist))

maps95 <- lapply(filepaths, terra::rast)
maps95 <- lapply(maps95, terra::aggregate, 20, na.rm = T)

# B. Make list into a data.frame
maps95 <- lapply(maps95, raster::raster) # make into raster
maps95 <- raster::stack(unlist(maps95))

data_pc <- data.frame(x = raster::coordinates(maps95)[, 1], 
                      y = raster::coordinates(maps95)[, 2])

for (i in 1:length(filenames)) {
  data_pc[, i + 2] <- maps95[[i]][]
}
names(data_pc) <- c('x', 'y', 'S10z', 'Sa', 'Sbi', 'Sci', 'Sdq', 'Sdq6', 'Sdr',
                    'Sds', 'Sfd', 'Sku', 'Smean', 'Sp', 'Sq', 'Srw', 'Srwi',
                    'Shw', 'Ssc', 'Ssk', 'Std', 'Stdi', 'Sv', 'Svi')

#-------------------------------------------------------------------------------
# (2) Clean data
#-------------------------------------------------------------------------------

# Remove infinite values
data_pc <- do.call(data.frame,lapply(data_pc, function(x){replace(x, is.infinite(x), NA)}))

# Remove nas
clean_data <- function(df) {
  NAs <- sapply(df, function(x) sum(is.na(x)))
  rm_cols <- which(NAs >= 200000) # depends on resolution, check maps
  df <- df[, -rm_cols]
  df <- na.omit(df)
  return(df)
}

for(i in 1:24){
  test <- data_pc[,i]
  print(length(test[is.na(test)]))
}

pc_noNA <- clean_data(data_pc)


#-------------------------------------------------------------------------------
# (3) Principal Component Analysis
#-------------------------------------------------------------------------------

pc_prc <- prcomp(pc_noNA[,3:length(pc_noNA)], center = TRUE, scale = TRUE)
summary(pc_prc)


#-------------------------------------------------------------------------------
# (4) Create maps of principal components
#-------------------------------------------------------------------------------

pc_dem <- rast("./pc_all_dems.tif")
pc_dem <- terra::aggregate(pc_dem, 20, na.rm = T) # make sure same resolution as Texture images
pc_dem <- raster::raster(pc_dem)

eviCols <- colorRampPalette(c('lightyellow1', 'darkgreen'))(100)
eviTheme <- rasterVis::rasterTheme(region = eviCols)

map_comps <- function(pc_dat, noNA_df, full_df, r, theme) {
  # Add pc values to no-NA dataframe.
  for (i in 1:5) {
    colname <- paste0('prc', i)
    noNA_df[, colname] <- pc_dat$x[, i]
  }
  
  # Add PCA results back to full raster dataframe.
  full_dat <- full_df %>% left_join(noNA_df)
  # Cut to only the prc columns.
  full_dat <- full_dat[, grep('prc', names(full_dat))]
  
  # Create rasters and maps with principle component values.
  out_maps <- list()
  out_rasts <- list()
  for (i in 1:5) {
    new_rast <- setValues(r, full_dat[, i])
    pc_map <- rasterVis::levelplot(new_rast, margin = F,
                                   par.settings = theme,
                                   ylab = NULL, xlab = NULL,
                                   main = paste0('PC', i))
    pc_map$par.settings$layout.heights[c( 'bottom.padding',
                                          'top.padding',
                                          'key.sub.padding',
                                          'axis.xlab.padding',
                                          'key.axis.padding',
                                          'main.key.padding') ] <- 1
    pc_map$aspect.fill <- TRUE
    out_maps[[i]] <- pc_map
    out_rasts[[i]]<- new_rast
  }
  
  # Plot in a grid.
  gridExtra::grid.arrange(grobs = out_maps, nrow = 2, ncol = 3)
  return(out_rasts)
}

pc_rasts <- map_comps(pc_prc, pc_noNA, data_pc, pc_dem, eviTheme)

#-------------------------------------------------------------------------------
# (4) Export principal component maps
#-------------------------------------------------------------------------------

pc_rasts <- lapply(pc_rasts, terra::rast)
pcfiles <- c("pc1","pc2","pc3","pc4","pc5")
outputfiles <- paste0("D:/pc_maps_20/", pcfiles, ".tif")

mapply(writeRaster, pc_rasts, outputfiles, overwrite = T)

#---------------------------------------------------------------------------------
# (5) Assess PCA
#---------------------------------------------------------------------------------

# A. Scree plot
plot_scree <- function(pc_dat) {
  screeplot(pc_dat, type = "l", npcs = 15,
            main = "Screeplot of the first 10 PCs")
  abline(h = 1, col = "red", lty = 5)
  legend("topright", legend = c("Eigenvalue = 1"),
         col = c("red"), lty = 5, cex = 0.6)
}

plot_scree(pc_prc)

# B. Cumulative Variance
plot_cvar <- function(pc_dat) {
  # Get cumulative variance explained.
  cumpro <- summary(pc_dat)$importance[3, ][1:16]
  
  # Create plot of cumulative variance, marking the 5th component as the cutoff.
  plot(cumpro, xlab = "PC #", ylab = "Amount of explained variance",
       main = "Cumulative variance plot")
  abline(v = 5, col = "blue", lty = 5)
  abline(h = cumpro[5], col = "blue", lty = 5)
  legend("topleft", legend = c("Cut-off @ PC5"),
         col = c("blue"), lty = 5, cex = 0.6)
}

plot_cvar(pc_prc)

# E. Plot loadings.
plt_names <- data.frame(old = names(data_pc)[3:ncol(data_pc)],
                        new = c('S10z', 'Sa', 'Sbi', 'Sci', 'Sdq', 'Sdq6', 'Sdr',
                                'Sds', 'Sfd', 'Sku', 'Smean', 'Sp', 'Sq', 'Srw', 'Srwi',
                                'Shw', 'Ssc', 'Ssk', 'Std', 'Stdi', 'Sv', 'Svi'))

plot_loadings <- function(pc_dat) {
  # Get rotation for top 5 components.
  loadings <- pc_dat$rotation[, 1:5]
  
  # Figure out the relative loadings.
  aload <- abs(loadings)
  rel <- sweep(aload, 2, colSums(aload), "/")
  
  # Convert relative loadings to dataframe.
  rel <- as.data.frame(rel)
  # Get good variable names (from dataframe created earlier).
  rel$var <- plt_names$new[match(rownames(rel), plt_names$old)]
  
  # Create importance plots.
  imp_plts <- list()
  for (i in 1:5) {
    temp <- rel
    # Determine whether component loading is postive or negative.
    temp$sign <- factor(sapply(loadings[, i], FUN = function(x) x / abs(x)),
                        levels = c(-1, 1))
    
    # Order loadings by value.
    temp <- temp[order(temp[, i]),]
    
    temp$var <- factor(temp$var, levels = temp$var)
    
    temp_plt <- ggplot(temp, aes(x = temp[, i], y = var)) +
      geom_point(size = 3, aes(pch = sign)) +
      scale_shape_manual(name = element_blank(),
                         breaks = c(1, -1),
                         values = c(19, 8),
                         labels = c("Positive", "Negative")) +
      xlab(paste0('PC', i)) +
      ylab('Metric') +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            legend.justification = c(1, 0),
            legend.position = c(1, 0),
            legend.background = element_blank(),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 12))
    
    imp_plts[[i]] <- temp_plt
  }
  
  # Return grid of first three components.
  gridExtra::grid.arrange(grobs = imp_plts[1:3], ncol = 3)
}

plot_loadings(pc_prc)
