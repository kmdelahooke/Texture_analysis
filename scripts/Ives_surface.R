## DISTINGUISHING IVESHEADIOMORPHS FROM BACKGROUND SURFACE TEXTURE ##

# 1. Import processed ivesheadiomorph DEMs
# 2. Sample non-ivesheadiomorph textures
# 3. Process non-ives texture samples
# 4. Surface metrology/PCA
# 5. Persistent homology/tNSE

# Libraries
library(Rtsne)
library(parallel)

source("./raster_processing_functions.R") # functions: 'set_zero()', 'change_res()'
source("./Ives_gsm_clustering.R") # functions: 'clean_data()', 'scale_mets()'

#-------------------------------------------------------------------------------
# (1) Import processed ivesheadiomorph DEMs
#-------------------------------------------------------------------------------
filelist <- list.files(path = "./data/processed_ives/", pattern = ".*.tif")
filepaths <- unlist(lapply("./data/processed_ives/", paste0, filelist))
filenames <- gsub(".tif", "", filelist)

ives_masked <- lapply(filepaths, terra::rast)

par(mfrow = c(5,6))
lapply(ives_masked, plot)
#-------------------------------------------------------------------------------
# (2) Sample non-ivesheadiomorph textures
#-------------------------------------------------------------------------------

# Import grid
filelist <- list.files(path = "F:/Pigeon Cove/processed_grid/", pattern = ".*.tif")
filepaths <- unlist(lapply("F:/Pigeon Cove/processed_grid/", paste0, filelist))
dems <- lapply(filepaths, terra::rast)

#crop out non-ives bits ext ( 50 x 50 )
s1 <- crop(dems[[2]], ext(800, 900, 1000, 1100))
s2 <- crop(dems[[8]], ext(1500, 1600, 1200, 1300))
s3 <- crop(dems[[8]], ext(1500, 1600, 1400, 1500))
s4 <- crop(dems[[8]], ext(1800, 1900, 1300, 1400))
s5 <- crop(dems[[9]], ext(1700, 1800, 1550, 1650))
s6 <- crop(dems[[9]], ext(1900, 2000, 1600, 1700))
s7 <- crop(dems[[9]], ext(1700, 1800, 1700, 1800))
s8 <- crop(dems[[17]], ext(2700, 2800, 1300, 1400))
s9 <- crop(dems[[18]], ext(2600, 2700, 1850, 1950))
s10 <- crop(dems[[18]], ext(2700, 2800, 1500, 1600))

surface <- list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)

#-------------------------------------------------------------------------------
# (3) Process non-ives texture samples
#-------------------------------------------------------------------------------

# Set zero
surface <- lapply(surface, set_zero)

# Make resolution consistent
max_res <- max(unlist(lapply(surface, res))) #should all be the same
print(paste('maximum resolution:', max_res))

ives_masked <- lapply(ives_masked, change_res, resol = max_res)

print(paste('Ivesheadiomorph resolution:', res(ives_masked[[1]])[1]))
print(paste('Surface resolution:', res(surface[[1]])[1]))

#Export .png for PH analysis

ims <- lapply(ives_masked, raster)
ims <- lapply(ims, imager::as.cimg)
ims <- lapply(ims, magick::image_read)
ims <- lapply(ims, magick::image_convert, "png")

for(i in seq_along(ims)){
  
  outputfile <- paste0("./data/processed_ives/ives_comp/ives_",i,".png")
  image_write(ims[[i]], path = outputfile, format = "png")
}


#-------------------------------------------------------------------------------
# (4) Surface metrology/PCA
#-------------------------------------------------------------------------------

# Remove plane
surface <- lapply(surface, geodiv::remove_plane)

ives_masked <- lapply(ives_masked, raster)
ives_masked <- lapply(ives_masked, geodiv::remove_plane)

# Calculate global surface metrics
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(geodiv)})

m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
               'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
               'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')

m_list <- m_list[-c(16, 24)] # package errors at these positions (scl, stxr)

metrics <- list()
for(i in 1:length(m_list)){
  metrics[[i]] <- parLapply(cl = cl, X = surface, fun = m_list[[i]])
  cat(m_list[[i]], sep = "\n")
}

stopCluster(cl)

# Create data.frame and export as .csv

df <- lapply(metrics, ldply)
df <- as.data.frame(matrix(unlist(df), ncol = 25))
names(df) <- c('Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr', 'Sbi',
               'Sci', 'Ssk', 'Sku', 'Sds', 'Sfd', 'Srw', 'Srwi', 'Shw',
               'Std', 'Stdi', 'Svi', 'Ssc', 'Sv', 
               'Sp', 'Sk', 'Smean', 'Spk', 'Svk')

write.csv(df, "./results/surface_gsms.csv")

# Read back in
surface <- read.csv("./results/surface_gsms.csv", header = T)[,-1]
ives <- read.csv("./results/ives_masked_gsms_comp.csv", header = T)[,-1]

# Combine
df <- data.frame(rbind(ives, surface))

#clean & scale metrics
df <- clean_data(df) 
df <- scale_mets(df)

# add 'surface' or 'ives' labels
type <- c(rep('ives', length(ives[,1])), rep('surface', length(surface[,1])))

# Calculate principal components
ives_prc <- prcomp(df, center = TRUE, scale = TRUE)

# Evaluate PCs
par(mfrow = c(1,1))
plot_scree(ives_prc)
plot_cvar(ives_prc)

plt_names <- data.frame(old = names(df)[1:ncol(df)],
                        new = c('Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr',
                                'Sbi', 'Sci', 'Ssk', 'Sku', 'Sds',
                                'Srw', 'Srwi', 'Shw', 'Stdi',
                                'Svi', 'Ssc', 'Sv', 
                                'Sp', 'Sk', 'Smean', 'Spk', 'Svk'))
plot_loadings(ives_prc)

# Plot PCA
df2 <- as.data.frame(-ives_prc$x[,1:2])

new_plot <- data.frame(x = df2$PC1, y = df2$PC2)%>% as_tibble()
colnames(new_plot) <- c("PC1", "PC2")
new_plot <- new_plot  %>% mutate(groups = type %>% as.factor())

ggpubr::ggscatter(new_plot, x = "PC1", y = "PC2",
                  color = "groups",
                  palette = 'pal5',
                  size =3,
                  repel = TRUE)

#-------------------------------------------------------------------------------
# (5) Persistent homology/tNSE
#-------------------------------------------------------------------------------

# Input data
wssc<- read.csv("./results/wasserstein_surface.csv", header = FALSE)[-1,-1]

# tSNE visualisation
tsne <- Rtsne(wssc, perplexity = 5, is_distance = T)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2])%>% as_tibble()
colnames(tsne_plot) <- c("Dim.1", "Dim.2")
tsne_plot <- tsne_plot  %>% mutate(groups = as.factor(type))
str(tsne_plot)
ggscatter(tsne_plot, x = "Dim.1", y = "Dim.2",
          color = "groups",
          palette = 'pal5',
          size = 2,
          repel = TRUE)
