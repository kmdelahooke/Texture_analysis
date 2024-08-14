## CALCULATE GLOBAL SURFACE METRICS FOR IVESHEADIOMORPHS

# 1. Import processed ivesheadiomorphs
# 2. Calculate global suface metrics for each ivesheadiomorph
# 3. Export surface metrics as .csv

# Libraries
library(plyr)
library(parallel)

# N.B. geodiv now runs with terra rasters, I will update this code

#-------------------------------------------------------------------------------
# (1) Import processed ivesheadiomorph DEMs
#-------------------------------------------------------------------------------

filelist <- list.files(path = "./data/processed_ives/", pattern = ".*.tif")
filepaths <- unlist(lapply("./data/processed_ives/", paste0, filelist))
filenames <- gsub(".tif", "", filelist)

ives_masked <- lapply(filepaths, terra::rast)

#-------------------------------------------------------------------------------
# (2)  Run over all rasters calculating surface metrics globally
#-------------------------------------------------------------------------------

cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(geodiv)})

# List of Surface Metrics to calculate
m_list <- list('sa', 'sq', 's10z', 'sdq', 'sdq6', 'sdr', 'sbi', 'sci', 'ssk', 
               'sku', 'sds', 'sfd', 'srw', 'std', 'svi', 'stxr', 'ssc', 'sv', 
               'sph', 'sk', 'smean', 'spk', 'svk', 'scl', 'sdc')

m_list <- m_list[-c(16, 24)] # package errors at these positions (scl, stxr) - check update

metrics <- list()

for(i in 1:length(m_list)){
  metrics[[i]] <- parLapply(cl = cl, X = ives_masked, fun = m_list[[i]])
  cat(m_list[[i]], sep = "\n")
}

stopCluster(cl)

#-------------------------------------------------------------------------------
# (3) Create data frame and export as .csv
#-------------------------------------------------------------------------------

df <- lapply(metrics, ldply)
df <- as.data.frame(matrix(unlist(df), ncol = 25))
names(df) <- c('Sa', 'Sq', 'S10z', 'Sdq', 'Sdq6', 'Sdr', 'Sbi',
               'Sci', 'Ssk', 'Sku', 'Sds', 'Sfd', 'Srw', 'Srwi', 'Shw',
               'Std', 'Stdi', 'Svi', 'Ssc', 'Sv', 
               'Sp', 'Sk', 'Smean', 'Spk', 'Svk')

write.csv(df, "./results/ives_masked_gsms_new.csv")
