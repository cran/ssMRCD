## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  warning = FALSE, 
  fig.dim = c(7, 4.5),
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ssMRCD)
library(ggplot2)

## ---- eval = FALSE------------------------------------------------------------
#  ? weatherAUT2021

## ----setup--------------------------------------------------------------------
data("weatherAUT2021")
head(weatherAUT2021)

## -----------------------------------------------------------------------------
cut_lon = c(9:16, 18)
cut_lat = c(46, 47, 47.5, 48, 49)
N = ssMRCD::N_structure_gridbased(weatherAUT2021$lon, weatherAUT2021$lat, cut_lon, cut_lat)
table(N)
N[N == 2] = 1
N[N == 3] = 4
N[N == 5] = 4
N[N == 6] = 7
N[N == 11] = 15
N = as.numeric(as.factor(N))

## -----------------------------------------------------------------------------
g1 = ggplot() + 
  geom_text(data = weatherAUT2021, aes(x = lon, y = lat, col = as.factor(N), label = N)) + 
  geom_hline(aes(yintercept = cut_lat), lty = "dashed", col = "gray") +
  geom_vline(aes(xintercept = cut_lon), lty = "dashed", col = "gray") +
  labs(x = "Longitude", y = "Latitude", title = "Austrian Weather Stations: Neighborhood Structure") +
  coord_fixed(1.4) + 
  theme_classic() +
  theme(legend.position = "none")
g1

## -----------------------------------------------------------------------------
GW = geo_weights(coordinates = weatherAUT2021[, c("lon", "lat")], 
                 N_assignments = N)
GW$W[4, ]

g1 + 
  labs(title = "Austrian Weather Stations: Weighting Matrix W") +
  geom_segment(aes(x = GW$centersN[4, 1], 
                   y = GW$centersN[4, 2], 
                   xend = GW$centersN[-4, 1]-0.1, 
                   yend = GW$centersN[-4, 2]-0.05,
                   alpha = GW$W[4, -4]), 
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 2,
               col = "blue")+
  geom_text(aes(x = GW$centersN[, 1], y = GW$centersN[, 2], label = 1:length(GW$centersN[, 2])))

## ---- eval = FALSE------------------------------------------------------------
#  ? parameter_tuning

## ---- eval = TRUE, message = FALSE--------------------------------------------
set.seed(123)
parameter_tuning(data = weatherAUT2021[, 1:6], 
                 coords = weatherAUT2021[, c("lon", "lat")], 
                 N_assignments = N, 
                 repetitions = 3, 
                 k = c(5, 10, 15), 
                 lambda = c(0, 0.25, 0.5, 0.75))$plot

## ---- eval = FALSE------------------------------------------------------------
#  ? local_outliers_ssMRCD

## ---- message = FALSE---------------------------------------------------------
res = local_outliers_ssMRCD(data = weatherAUT2021[, 1:6],
                            coords = weatherAUT2021[, c("lon", "lat")],
                            N_assignments = N,
                            lambda = 0.5,
                            k = 10)
summary(res)

## ---- message = FALSE---------------------------------------------------------
cat(weatherAUT2021$name[res$outliers], sep = ",\n")

## -----------------------------------------------------------------------------
covariance_res = res$ssMRCD
plot(covariance_res, 
     centersN = res$centersN, 
     manual_rescale = 0.5, 
     type = c("convergence", "ellipses"), 
     colour_scheme = "regularity", 
     legend = TRUE,
     xlim = c(9, 19))

## -----------------------------------------------------------------------------
biplot(stats:: princomp(weatherAUT2021[, 1:6],
                        cor  = TRUE, 
                        covmat = robustbase::covMcd(weatherAUT2021[, 1:6])), 
       col = c("grey", "black"))

## ---- eval = FALSE------------------------------------------------------------
#  ? plot.locOuts

## ---- fig.dim = c(5,3.5)------------------------------------------------------
plot(res, type = "hist", pos = 4)

## -----------------------------------------------------------------------------
plot(res, type = "spatial", colour = "outScore", xlim = c(9.5, 19))

## ---- fig.dim = c(7,3)--------------------------------------------------------
# SCHOECKL
plot(res, type = "lines", focus = res$outliers[3])

## ---- fig.dim = c(7,4)--------------------------------------------------------
plot(res, type = "3D", colour = "outScore", theta = 0, phi = 10)

## ---- eval = FALSE, include = T-----------------------------------------------
#  library(animation)
#  
#  # fineness of movement
#  n = 90
#  param = rbind(cbind(rep(0, 2*n + n/2), c(seq(90, 0, length.out = n), seq(0, 90, length.out = n), seq(90, 30, length.out = n/2))),
#                cbind(c(seq(0, 90, length.out = n), seq(90, 0, length.out = n)), rep(30, 2*n)),
#                cbind(rep(0, n/2), seq(30, 90, length.out = n/2)))
#  
#  # use function saveGIF
#  saveGIF(interval = 0.2,
#          movie.name = "local_outliers_3D.gif",
#          expr = {
#            for (i in 1:dim(param)[1]) {
#              plot(outs, type = "3D", colour = "outScore", theta = param[i, 1], phi = param[i, 2])
#            }}
#          )

