## CLUSTERING OF IVESHEADIOMORPHs ##

# 1. Import data - gsm, ph
# 2. Determine optimal number of clusters
# 3. Clustering
# 4. Visualisation
# 5. Cluster evaluation

# Libraries
library(fpc)
library(Rtsne)
library(missForest)
library(patchwork)

#-------------------------------------------------------------------------------
# (1) Import data
#-------------------------------------------------------------------------------

# GSMs
gsms <- read.table("./results/ives_masked_gsms_new.csv", sep = ",", header = T)[,-1]

# PH
wss <- read.csv("./results/wasserstein_lowerstar.csv", header = FALSE)[-1,-1]
wssc <- read.csv("./results/wasserstein_cubical.csv", header = FALSE)[-1,-1]
wssb <- read.csv("./results/wasserstein_binarised.csv", header = FALSE)[-1,-1]


#-------------------------------------------------------------------------------
# (2) Transform for analysis, PCA follows geodiv vignette
#-------------------------------------------------------------------------------
##post processing
clean_data <- function(df) {
  # Remove columns with very large numbers of NAs.
  NAs <- sapply(df, function(x) sum(is.na(x)))
  NANs <- sapply(df, function(x) sum(x == 0))
  rm_cols <- which(NAs >= 10 | NANs >=10)
  df <- df[, -rm_cols]
  # impute NAs from remaining columns.
  df <-missForest(df)$ximp
  print <- missForest(df)$OOBerror
  
  return(df)
}

gsms2 <- clean_data(gsms)

# Get all metrics on same scale (0-1).
scale_mets <- function(df) {
  for (i in 1:ncol(df)) {
    df[,i] <- (df[, i] - min(df[, i], na.rm = TRUE)) /
      (max(df[, i], na.rm = TRUE) - min(df[, i], na.rm = TRUE))
  }
  return(df)
}
gsms2 <- scale_mets(gsms2)s

# PCA GSMs
ives_prc <- prcomp(gsms2, center = TRUE, scale = TRUE)
df <- as.data.frame(-ives_prc$x[,1:2])
summary(ives_prc)

#df0 <- df
#df10 <- df
#df20 <- df

# PH: specify as distance matrix
#df <- as.dist(m = wss)
#df2 <- as.dist(m = wssc)


#-------------------------------------------------------------------------------
# (2) Determine optimal number of clusters
#-------------------------------------------------------------------------------

plot_optimal_clusters <- function(dist, clust){
  p1 <- fviz_nbclust(dist, clust, method = 'wss') 
  p2 <- fviz_nbclust(dist, clust, method = 'silhouette')
  p3 <- fviz_nbclust(dist, clust, method = 'gap_stat')
  
  (p1|p2|p3) # requires patchwork
}


plot_optimal_clusters(df, kmeans) # replace accordingly

#-------------------------------------------------------------------------------
# (3) Clustering
#-------------------------------------------------------------------------------

hc.res <- eclust(df, "hclust", k = 2, graph = FALSE)
km.res <- eclust(df, "kmeans", k = 2, nstart = 25, graph = FALSE)

# K-medians
pam.res <- pam(df, 2)
clust <- pam.res$clustering %>% as.factor()

#-------------------------------------------------------------------------------
# (4) Visualisation
#-------------------------------------------------------------------------------

#GSMs
fviz_dend(hc.res, show_labels = TRUE,
          palette = "jco", as.ggplot = TRUE)
fviz_cluster(hc.res, data = df)
fviz_cluster(pam.res, data = df)

new_plot <- data.frame(x = df$PC1, y = df$PC2)%>% as_tibble()
colnames(new_plot) <- c("PC1", "PC2")
new_plot <- new_plot  %>% mutate(groups = pam.res$cluster %>% as.factor())

ggscatter(new_plot, x = "PC1", y = "PC2",
          color = "groups",
          label = rownames(gsms2),
          palette = 'pal5',
          size = 2)

# PH
tsne <- Rtsne(wssc, perplexity = 5, is_distance = T)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2])%>% as_tibble()
colnames(tsne_plot) <- c("Dim.1", "Dim.2")
tsne_plot <- tsne_plot  %>% mutate(groups = pam.res$cluster %>% as.factor())

ggscatter(tsne_plot, x = "Dim.1", y = "Dim.2",
          color = "groups",
          palette = 'pal5',
          label = as.character(seq(1,29,1)),
          size = 2)

#-------------------------------------------------------------------------------
# (5) Cluster evaluation
#-------------------------------------------------------------------------------

# Dunn Index
km_stats <- cluster.stats(dist(new_plot),  km.res$cluster)
km_stats$dunn

pam.res1 <- eclust(df, "pam", k = 2)
pam.res2 <- eclust(df2, "pam", k = 2)

fviz_cluster(pam.res2)

#pam.res <- eclust(df, "pam", k = 2)
clust_stats <- cluster.stats(d = df, km.res$cluster)
clust_stats$dunn

#agreement between clustering: Rand index
table(pam.res1$cluster, pam.res2$cluster)
clust_stats <- cluster.stats(d = df, clustering = pam.res1$cluster, alt.clustering = pam.res2$cluster)
clust_stats$corrected.rand

# Jitter
cbj <- clusterboot(wssc, bootmethod = "jitter", k= 2, clustermethod = pamkCBI)
cbj$jittermean

# Bootstrap
cbb <- clusterboot(df2, bootmethod = "boot", k=2, clustermethod = pamkCBI, method = "average")
#cbb <- clusterboot(wssc, bootmethod = "boot", k=2, clustermethod = kmeansCBI)
cbb$bootmean


