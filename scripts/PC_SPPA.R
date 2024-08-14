 ## SPATIAL POINT PATTERN ANALYSIS OF PIGEON COVE ##

# 1. Import point pattern data
# 2. Create binary mask window
# 3. Create point pattern
# 4. Orientations
# 5. Summary statistics - pcf, l-function, log nnk
# 6. Non-parametric heterogeneous poisson model
# 7. Directional changes in intensity

# Libraries
library(spatstat)
library(circular)
library(CircMLE)
library(ggplot2)
library(tidyverse)

source("./raster_processing_functions.R") # functions: 'as.im.SpatRaster()'

#-------------------------------------------------------------------------------
# (1) Import point pattern data
#-------------------------------------------------------------------------------
#Import data
df <- read.csv("./data/pc_fronds23_all.csv")[, -1]

hist(log(df$length), breaks = 10)
summary(df$length)

shapiro.test(log(df$length))
#-------------------------------------------------------------------------------
# (2) Create binary mask window
#-------------------------------------------------------------------------------

r <- rast("./data/texture_images_w2/pc_sa.tif")
r <- aggregate(r, 20, na.rm = TRUE)
r <- classify(r, cbind(-Inf, Inf, 256))

win_pc <- as.owin(as.im.SpatRaster(r))

#-------------------------------------------------------------------------------
# (3) Create point pattern
#-------------------------------------------------------------------------------

pp <- ppp(df$x, df$y, window = win_pc)
plot(pp, cols = "indianred", pch = 20, cex = 2.5)

#-------------------------------------------------------------------------------
# (4) Orientations
#-------------------------------------------------------------------------------

# Rose diagram
angles <- circular(df$angle, units='degrees', rotation='clock', zero=pi/2)
plot.circular(angles, col='lightblue', stack=TRUE,sep =0.05, shrink = 1.5, cex = 1.5)
arrows.circular(mean(angles, na.rm = TRUE)) 

# tests for multimodality
circular::rayleigh.test(angles) 
CircMLE::HR_test(angles)

#angles plotted on map
pc <- ppp(df$x, df$y, marks = as.character(df$angle), window = win_pc)

# convert so direction is angle in degree acw from the x-axis
convert_angles <- function(x){
  ifelse(x < 90, {x = -x + 90}, {x = -x + 450})}

d <- unlist(lapply(as.numeric(marks(pc)), convert_angles))

plot(pc, shape = 'arrows', direction = d, size=100, cols="#B40F20")

#-------------------------------------------------------------------------------
# (5) Summary statistics - pcf, l-function, log nnk
#-------------------------------------------------------------------------------

# A. PCF
E_pcf <- envelope(pp, fun = 'pcf', nrank = 50, global = F, nsim = 999, bw = 35) 
plot(E_pcf, ylim = c(0, 2.5))

bw = 0.2/(intensity(pp)^(0.5))

# B. L-Function
E_l <- envelope(pp, 'Lest',  nsim = 999, nrank = 50)
plot(E_l)

# transformed: L(r) -r
p1 <- ggplot(E_l, aes(r))+
  theme_classic()+
  xlab('Distance (mm)') + ylab('L(r) - r') +
  geom_ribbon(aes(ymin=lo -r, ymax=hi -r), fill ='grey70')+
  geom_line(aes(y=obs - r), linewidth=1.1) +
  geom_line(aes(y=0), lty=2, col='red')+
  expand_limits(x = 0, y = 0) 
p1

# C: log nnk ~ log k
mean_nndist <- function(pp, k){
  mean(nndist(pp, k = k))
}

k <- seq(1, 20, 1)
nnk <- unlist(lapply(k, mean_nndist, pp = pp))

lambda <- intensity(pp)
CSR <- (k/(lambda*pi))^0.5

plot(log(nnk)~ log(k), type = "l")
lines(log(CSR) ~log(k), lty = 2, col = "red")

#-------------------------------------------------------------------------------
# (6) Non-parametric heterogeneous poisson model
#-------------------------------------------------------------------------------
# kernel smoothed intensity
den <- density.ppp(pp, kernel = "epanechnikov", sigma = 450)  #here default sigma
plot(den)

# fit heterogenous poisson  model to smoothed intensity function
hpp <- ppm(pp ~ den) 
plot(predict.ppm(hpp))

# goodness of fit
dclf.test(hpp)

# envelope of fitted model
E <- envelope(hpp, fun = 'pcf', nsim = 999, nrank = 50, bw = 35)
plot(E, ylim = c(0, 2.5))
abline(v = 450) #sigma

# spatial structure of residuals
diagnose.ppm(hpp)

#-------------------------------------------------------------------------------
# (7) Directional changes in intensity
#-------------------------------------------------------------------------------

# Point of max/min erosion
dr <- raster(den)
max <- xyFromCell(dr, which.max(dr)) 
min <- xyFromCell(dr, which.min(dr)) 

a <- max[1]
b <- max[2]
c <- min[1]
d <- min[2]

# Add to intensity map
plot(den, main = "Fossil density")
plot(pp, add = TRUE)
maxp <- ppp(x = max[1], y = max[2], window = pp$window)
plot(maxp, add = TRUE, cex = 2, pch = 8)
minp <- ppp(x = min[1], y = min[2], window = pp$window)
plot(minp, add = TRUE, col = 'red', cex = 2, pch = 8)
legend("top", box.col = FALSE, legend = c("maximum density", "minimum density"), pch = 8, col = c("black", "red"))

# Fitting heterogenous Poisson models
csr <- ppm(pp)
hx <- ppm(pp ~ x)
hy <- ppm(pp ~ y)
hxy <- ppm(pp ~x + y)
hmax <- ppm(pp ~ sqrt((x - a)^2 +(y - b)^2))
hmin <- ppm(pp ~ sqrt((x - c)^2 +(y - d)^2))

predict(hxy) %>% plot() #for example

# Assess model fit
dclf.test(csr)
dclf.test(hx)
dclf.test(hy)
dclf.test(hxy)
dclf.test(hmax)
dclf.test(hmin)

#Assess model fit using spatial Kolmogorov-Smirnov tests 
# - if p>0.05, the difference between two samples is not significant enough to say that they have different distribution
fxy <- function(x,y) {x + y}
fmax <- function(x,y) {sqrt((x-a)^2 +(y-b)^2)}
fmin <- function(x,y) {sqrt((x-c)^2 +(y-d)^2)}

csr.test <- cdf.test(csr, 'x', test='ks') 
x.test <- cdf.test(hx, 'x', test='ks')
y.test <- cdf.test(hy, 'y', test='ks')
xy.test <- cdf.test(hxy, fxy, test='ks')
max.test <- cdf.test(hmax, fmax, test='ks')
min.test <- cdf.test(hmin, fmin, test='ks')

#-------------------------------------------------------------------------------
# (8) Anisotropy of point pattern
# ------------------------------------------------------------------------------

# pair-orientation distribution
fv <- pairorient(pp, r1=100, r2=800, sigma = 10)
f <- convert_angles(f)

rose(fv)
rose(f, start = "N", clockwise = F)



