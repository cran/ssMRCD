## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE, 
  fig.dim = c(7, 4.5),
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ssMRCD)
library(ggplot2)
library(dplyr)

## ----metadata, eval = FALSE---------------------------------------------------
# # get meta data for the data set
# ? weatherAUT2021

## ----load data, eval = TRUE---------------------------------------------------
# load the data
data("weatherAUT2021")

# inspect the data
head(weatherAUT2021)

# select variables, station names and number of observations
data = weatherAUT2021 %>% select(p:rel)
stations = weatherAUT2021$name
n = dim(data)[1]

## ----build groups-------------------------------------------------------------
# build 5 groups of observations based on spatial proximity and geography
cut_lon = c(min(weatherAUT2021$lon)-0.2, 12, 16, max(weatherAUT2021$lon) + 0.2)
cut_lat = c(min(weatherAUT2021$lat)-0.2, 48, max(weatherAUT2021$lat) + 0.2)
groups = ssMRCD::groups_gridbased(weatherAUT2021$lon, 
                                  weatherAUT2021$lat, 
                                  cut_lon, 
                                  cut_lat)
N = length(unique(groups))
table(groups)

## ----run model----------------------------------------------------------------
# calculate MG-GMM
model = cellMGGMM(X = data, groups = groups,
                  nsteps = 100, alpha = 0.5,
                  maxcond = 100)

## ----mixture probabilities----------------------------------------------------
# mixture probabilities
cat("Pi (in %):\n")
round(model$pi_groups*100, 2)

## ----percentage of outlier----------------------------------------------------
# percentage of outliers
cat("% Outliers per group and variable:\n")
round(sapply(1:N, function(x) colMeans(1-model$W[groups == x, ]))*100, 2)

## ----residuals----------------------------------------------------------------
# calculate residuals
res = residuals_mggmm(X = data, 
                groups = groups,
                Sigma = model$Sigma,
                mu = model$mu, 
                probs = model$probs,
                W = model$W)

