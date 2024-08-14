## SURFACE METRIC SENSITIVITY ANALYSIS ##

# 1. Effects of .tif resolution
# 2. Preprocessing vs post processing DEMs

source("./raster_processing_functions.R") # functions: 'set_zero()', 'change_res()'

#-------------------------------------------------------------------------------
# (1) Effects of .tif resolution
#-------------------------------------------------------------------------------

# Import data
filelist <- list.files(path = "./data/hi_res_grid", pattern = ".*.tif")
filepaths <- unlist(lapply("./data/hi_res_grid/", paste0, filelist))
sens_dems <- lapply(filepaths, terra::rast)

# Process
sens_dems<- lapply(sens_dems, set_zero)

# Convert to raster and remove best fit polynomial surface
sens_dems <- lapply(sens_dems, raster::raster)
sens_dems <- lapply(sens_dems, geodiv::remove_plane)

# Calculate texture metrics
m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
               'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
               'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')

m_list <- m_list[-c(16, 24)]

metrics <-list()

for(i in 1:length(m_list)){
  metrics[[i]] <- lapply(sens_dems, FUN = m_list[[i]])
  cat(m_list[[i]], sep = "\n")
}

str(metrics)

df <- lapply(metrics, ldply)
df <- as.data.frame(matrix(unlist(df), ncol = 25))
names(df) <- c('Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr', 'Sbi',
               'Sci', 'Ssk', 'Sku', 'Sds', 'Sfd', 'Srw', 'Srwi', 'Shw',
               'Std', 'Stdi', 'Svi', 'Ssc', 'Sv', 
               'Sp', 'Sk', 'Smean', 'Spk', 'Svk')
df <- df[,-12]


# Compare 
x <- c(0.0005, 0.0005, 0.001)
par(mfrow = c(5,5))
for(i in 1:length(df)){
  plot(df[,i] ~ log(x))
}

#-------------------------------------------------------------------------------
# (2) Preprocessing vs post processing DEMs
#-------------------------------------------------------------------------------

#remove outliers
test95 <- remove_outliers(dems[[8]], 0.95)
test99 <- remove_outliers(dems[[8]], 0.99)

#whole surface 
dem <- do.call(merge, dems)
dem95 <- remove_outliers(dem, 0.95)
dem99 <- remove_outliers(dem, 0.99)
plot(dem99)

#aggregate
test2 <- terra::aggregate(dems[[8]], 2, na.rm = T)
test5 <- terra::aggregate(dems[[8]], 5, na.rm = T) 

# remove outliers, aggregate
test99_2 <- terra::aggregate(test99, 2, na.rm = T) 
test99_5 <- terra::aggregate(test99, 5, na.rm = T) 

# aggregate, remove outliers
test2_99 <- remove_outliers(test2, 0.99) 
test5_99 <- remove_outliers(test5, 0.99) 
test5_95 <- remove_outliers(test5, 0.95)

#interpolate
test5_99_i <- interpolate_dem(test5_99, 2) 
test5_95_i <- interpolate_dem(test5_95, 2) 
test5_i <- interpolate_dem(test5, 2) 
test5_5i <- interpolate_dem(test5, 5)


##TEXTURE IMAGES
sq <- texture_image(x = dems[[8]], window = 'square',size = 2, metric = 'sq', parallel = TRUE)
test95_sq <- texture_image(x = test95, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test99_sq <- texture_image(x = test99, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()

test2_sq <- texture_image(x = test2, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test5_sq <- texture_image(x = test5, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()

test99_2_sq <- texture_image(x = test99_2, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test99_5_sq <- texture_image(x = test99_5, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()

test2_99_sq <- texture_image(x = test2_99, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test5_99_sq <- texture_image(x = test5_99, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()

test5_99_i_sq <- texture_image(x = test5_99_i, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test5_95_i_sq <- texture_image(x = test5_95_i, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test5_i_sq <- texture_image(x = test5_i, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()
test5_5i_sq <- texture_image(x = test5_5i, window = 'square',size = 2, metric = 'sq', parallel = TRUE) %>% plot()


##window size

testw5 <- texture_image(x = test5_99_i, window = 'square',size = 5, metric = 'sq', parallel = TRUE) %>% plot()
testw10 <- texture_image(x = test5_99_i, window = 'square',size = 10, metric = 'sq', parallel = TRUE) %>% plot()
testw15 <- texture_image(x = test5_99_i, window = 'square',size = 15, metric = 'sq', parallel = TRUE) %>% plot()
testw25 <- texture_image(x = test5_99_i, window = 'square',size = 25, metric = 'sq', parallel = TRUE) %>% plot()
testw50 <- texture_image(x = test5_99_i, window = 'square',size = 50, metric = 'sq', parallel = TRUE) %>% plot()
testw100 <- texture_image(x = test5_99_i, window = 'square',size = 100, metric = 'stxr', parallel = TRUE) %>% plot()


##COMPARE POST AND PRE
psq <- remove_outliers(sq, 0.99) %>% plot
psq2 <- aggregate(sq, 5, na.rm = T) 
psq3 <- aggregate(psq, 5, na.rm = T)
psq4 <- interpolate_dem(psq3, 2)
spatialEco::raster.gaussian.smooth(psq4, type = "mean", na.rm = T, n = 9)%>% plot()

##post and pre-processed

ppsq <- remove_outliers(test5_99_i_sq, 0.99) %>% plot

##value distributions
hist(values(ppsq))
boxplot(values(ppsq))
