
## Demo example illustrating the methodology of Wadoux and Heuvelink, Methods in Ecology
## and Evolution 14, 1320-1332, DOI: 10.1111/2041-210X.14106

############
############  load libraries
############

library(ggplot2)
library(gridExtra)
library(ranger)
library(sp)
library(raster)
library(matrixStats)
library(spatstat)
library(gstat)

############
############  Initialize script
############

root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

############
############ Simulate fields 
############

# Reference field 
## grid of size 100*100 cells
## we simulate a realization of a map composed of a linear spatial trend superimposed on a Gaussian random field
## trend has an intercept of 5 and slope parameter of 0.1 for the x-axis and 0.05 for the y-axis
## the Gaussian random field has a mean of zero and a covariance given by C(h) = 5 exp(-h/25), where h is the lag distance
## We take a sample of 200 points for the analysis. 

# Define discretization grid of 100 times 100
grid <- expand.grid(x1 = seq(0, 100, length.out = 100),
                    x2 = seq(0, 100, length.out = 100))

# Compute spatial trend; x1 and x2 ares used as covariates, the linear model has an intercept
grid$mu <- 5 +  0.1*grid$x1 + 0.05*grid$x2

# Define covariance function for simulation of residuals
covfun <- function(sill, range, Dist) {
  sill * exp(-Dist / range)
}

# Compute matrix with distances between simulation nodes
distx<-outer(grid$x1,grid$x1,FUN="-")
disty<-outer(grid$x2,grid$x2,FUN="-")
Dist<-sqrt(distx^2+disty^2)

# Compute matrix with mean covariances
sill <- 5
range <- 25
C <- covfun(sill, range, Dist = Dist)

# Simulate values for residuals by Cholesky decomposition
set.seed(31415)
Upper <- chol(C)
G <- rnorm(n=nrow(grid),0,1) #simulate random numbers from standard normal distribution
grid$residuals <- crossprod(Upper,G)

# Add the trend and the Cholesky decomposed values
grid$y <- grid$mu+grid$residuals

# Add residuals to trend
grid <- grid[,c('x1', 'x2', 'y')]

# Plot reference map
ggplot(grid) + geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + theme_bw() + coord_fixed()

# Mean and variance of the reference map
mean(grid$y)
var(grid$y)

# Sample 200 points for fitting the model
set.seed(1234)
pt_all <- grid[sample(1:nrow(grid), 200),]

# Plot the two covariates and the sample of 200 observations
plot1 <- ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = x1)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

plot2 <- ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = x2)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

plot3 <- ggplot(pt_all) + 
  geom_tile(aes(x1, x2, fill = y)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

# Arrange the plots next to each other
grid.arrange(plot1, plot2, plot3, ncol = 2)

############
############ Model cross-validation 
############

# Prediction at points with CV. Note that this could already be obtained from the ranger output but better make our own cross-validation. 
# Number of folds
nfolds <- 10

# Ceate a data.frame to store the prediction of each fold
pt_all$predRF <- NA
pt_all$predVarRF <- NA
pt_all$interval90 <- NA

# Samples that will be used for validation
sp <- sample(1:nfolds, nrow(pt_all), replace = T)
sp <- sample(sp)
pt_all$sp <- sp

# make the for loop over the validation folds 
for(i in 1:nfolds){
  # Extracting the training and validation subsets
  dat.train <- as.data.frame(pt_all[pt_all$sp != i,][c('y', c('x1', 'x2'))])
  dat.val <- as.data.frame(pt_all[pt_all$sp == i,][c('x1', 'x2')])
  
  # Train QRF on calibration set
  m.rf <- ranger(y ~ x1 + x2, data = dat.train, num.trees = 250, quantreg = TRUE) 
  
  # Make 500 simulations from the model on the validation set
  Nsim <- 500
  simQRF <- predict(m.rf, data = dat.val, type = "quantiles", what = function(x) sample(x, Nsim, replace = TRUE))
  
  # Calculate mean and variance of the realizations
  pt_all$predRF[pt_all$sp == i] <- rowMeans(simQRF$predictions)
  pt_all$predVarRF[pt_all$sp == i] <- rowVars(simQRF$predictions) * Nsim/(Nsim - 1)
  
  # check the 90% coverage, return TRUE or FALSE if covered or not covered by the interval
  cov90QRF <- predict(m.rf, data = dat.val, type = "quantiles", quantiles = c(0.05, 0.95))$predictions
  pt_all$interval90[pt_all$sp == i] <- cov90QRF[,1] <= pt_all[pt_all$sp == i,][,'y'] & pt_all[pt_all$sp == i,][,'y'] <= cov90QRF[,2]
  
  # Print the fold number
  print(paste0("Fold ", i, " done"))
}  

# Scatter plot of the observaed vs predicted obtained by CV
plot(pt_all$y, pt_all$predRF) ; abline(0, 1)

# Calculate mean error (bias)
mean(pt_all$predRF - pt_all$y)

# Calculate RMSE
sqrt(mean((pt_all$predRF - pt_all$y)^2))

# Calculate modelling efficiency (i.e. R2)
ss_res <- sum((pt_all$y - pt_all$predRF)^2)
ss_tot <- sum((pt_all$y - mean(pt_all$y))^2)
1 - (ss_res / ss_tot)

# check the 90% coverage probability (PICP90)
sum(pt_all$interval90 == TRUE)/nrow(pt_all)*100 # should be close to 90

############
############ Final model fitting 
############

# Fit a ranger quantile regression forest model with 250 trees 
model.rf <- ranger(y ~x1+x2, data = pt_all, num.trees = 250, quantreg = TRUE)
model.rf

# Simulate realizations of the map prediction with QRF
Nsim <- 500
simQRF <- predict(model.rf, data = grid, type = "quantiles", what = function(x) sample(x, Nsim, replace = TRUE))
str(simQRF)

# Take the mean and variance of the simulations
grid$predVarRF <- rowVars(simQRF$predictions)
grid$predRF <- rowMeans(simQRF$predictions)

# Plot the prediction mean and variance maps
plot11 <- ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = predRF)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

plot22 <- ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = predVarRF)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

# Arrange the plots next to each other
grid.arrange(plot11, plot22, ncol = 2)

############
############ Uncertainty of spatial averages
############
# follow the steps described in Section 4.1 in https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14106

############ Step 1: prediction uncertainty at point support for all locations
# this is grid$predVarRF 
ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = predVarRF)) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

############ Step 2: compute the point support prediction error standard deviation
# the standard deviation is obtained from the sqrt of the prediction variance
ggplot(grid) + 
  geom_tile(aes(x1, x2, fill = sqrt(predVarRF))) +
  scale_fill_distiller(palette="Spectral") + 
  theme_bw() + 
  coord_fixed()

############ Step 3: at every observation location, divide the error (observed minus predicted) by the standard deviation (square root of variance) to get a standardized error.
pt_all$stdErrRF <- (pt_all$y - pt_all$predRF)/sqrt(pt_all$predVarRF)

############ Step 4: Compute correlation function from the standardized errors. Note that here we assume a stationary correlation function, that is correlation only depends on Euclidean distance.
# Sample variogram
lzn.vgm.RF = variogram(stdErrRF ~ 1, pt_all, locations = ~ x1 + x2)
lzn.vgm.RF
plot(lzn.vgm.RF, pl = T)

# Fit the variogram model
# lzn.fit.RF = fit.variogram(lzn.vgm.RF, fit.method = 6, model = vgm(0.7, "Ste", 25, 0.1, kappa = 0.5), fit.kappa = F)
# lzn.fit.RF
# plot(lzn.vgm.RF, lzn.fit.RF)

# We attempt to improve the fit with a nested variogram to better account for the range. Note that this is more difficult to fit
lzn.fit.RF.nested = fit.variogram(lzn.vgm.RF, vgm(psill = 0.46, model = "Exp", range = 2, nugget = 0.05, add.to = 
                                                    vgm(psill = 0.08, model = "Sph", range = 15, nugget = 0, add.to = 
                                                          vgm(psill = 0.06, model = "Exp", range = 60, nugget = 0))),
                                  fit.sills = FALSE, 
                                  fit.ranges = FALSE, 
                                  fit.method = 2)
lzn.fit.RF.nested
plot(lzn.vgm.RF, lzn.fit.RF.nested)

# The correlation function is derived from the variogram
corRF <- variogramLine(lzn.fit.RF.nested, dist_vector = seq(1, 100, 1), covariance = TRUE)
corRF$gamma <- (sum(lzn.fit.RF.nested$psill) - variogramLine(lzn.fit.RF.nested, dist_vector = seq(1, 100, 1), covariance = FALSE)$gamma)/sum(lzn.fit.RF.nested$psill)

# plot the correlogram
ggplot(data = corRF, aes(x = dist, y=gamma))+
  geom_line(linewidth=2, alpha = 0.4)

############ Step 5: Obtain the variance of the prediction error of the spatial average 

# Function that return Eq. 7 from the paper for a pair of any two points in the area
f <- function(pt1x, pt1y, pt2x, pt2y){
  # Any two spatial points in the area
  sp.grid <- SpatialPoints(coords = rbind(c(pt1x, pt1y),c(pt2x, pt2y)))
  
  # Compute euclidean distance between points
  dist.mat <- pointDistance(sp.grid, lonlat = FALSE)[1,][2]
  
  # Calculate rho(|s-u|) from Eq. 7
  gamma.pt <- (sum(lzn.fit.RF.nested$psill) - variogramLine(lzn.fit.RF.nested, dist_vector = dist.mat, covariance = FALSE)$gamma)/sum(lzn.fit.RF.nested$psill)
  
  # Extract standard deviations at points 
  sd.pt <- sqrt(raster::extract(rasterFromXYZ(grid), sp.grid)[,'predVarRF'])
  
  # Return Eq. 7
  return(sd.pt[1]*sd.pt[2]*gamma.pt)
  
}

# Vectorize the function so it takes a input vector of points
ff <- Vectorize(f)

############ Step 6: Evaluation of the double integral with MC integration

# Evaluation of the double integral 
mcint2.mod <- function (f, shape, m = m) {
  # Select points from an uniform distribution 
  pairs1 <- as.data.frame(runifpoint(m, win = shape))
  pairs2 <- as.data.frame(runifpoint(m, win = shape))
  x1 <- pairs1[,1]
  y1 <- pairs1[,2]
  x2 <- pairs2[,1]
  y2 <- pairs2[,2]
  
  # Run the function for calculating the variance of the spatial average
  z.hat <- ff(pt1x = x1, pt1y = y1, pt2x = x2, pt2y = y2)
  return(mean(z.hat))
}

# Run the evaluation of the integral with MC integration and m samples
var_average <- mcint2.mod(ff, shape = owin(xrange=c(0,100), yrange=c(0,100)), m = 1000)

# Report the results for mean and standard deviation of the field
mean(grid$predRF)
sqrt(var_average)

# true value inside 90 percent prediction interval?
mean(grid$y)
mean(grid$predRF) - 1.64*sqrt(var_average); mean(grid$predRF) + 1.64*sqrt(var_average)

# Optional: test for convergence with various m values. It takes some time to run. 
# test convergence for various MC sample size
# val <- c()
# for (i in seq(10, 10000, 500)){
#   val <- append(val, mcint2.mod(ff, shape = owin(xrange=c(0,100), yrange=c(0,100)), m = i))
# }
# plot(seq(10, 10000, 500), val, type='l')

############ Step 7: To obtain block total you multiply the results of the average by the square of the block area
# block area, here simple because the area of a rectangle = length × width 
B = 100*100
var_total <- var_average*(B^2)

# report the results for mean and standard deviation of the field total
sum(grid$predRF)
sqrt(var_total)