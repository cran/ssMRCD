## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  warning = FALSE, 
  fig.dim = c(7, 4.5),
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ssMRCD)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

## ----eval = FALSE-------------------------------------------------------------
# ? weatherAUT2021

## ----setup--------------------------------------------------------------------
data("weatherAUT2021")
head(weatherAUT2021)

rownames(weatherAUT2021) = weatherAUT2021$name

## -----------------------------------------------------------------------------
# Load Austria as sf object
austria <- ne_countries(scale = "medium", country = "Austria", returnclass = "sf")

g_boundary = ggplot() + 
  geom_sf(data = austria, fill = "transparent", color = "black") +
  theme_classic()

## -----------------------------------------------------------------------------
# group by spatial grid
cut_lon = c(9:16, 18)
cut_lat = c(46, 47, 47.5, 48, 49)
groups = groups_gridbased(x = weatherAUT2021$lon, 
                     y = weatherAUT2021$lat, 
                     cutx = cut_lon, 
                     cuty = cut_lat)
table(groups)

# join particularly small groups together
groups[groups == 2] = 1
groups[groups == 3] = 4
groups[groups == 5] = 4
groups[groups == 6] = 7
groups[groups == 11] = 15
groups = as.numeric(as.factor(groups))
table(groups)

## -----------------------------------------------------------------------------
g_groups = g_boundary + 
  geom_text(data = weatherAUT2021, aes(x = lon, y = lat, col = as.factor(groups), label = groups)) + 
  geom_hline(aes(yintercept = cut_lat), lty = "dashed", col = "gray") +
  geom_vline(aes(xintercept = cut_lon), lty = "dashed", col = "gray") +
  labs(x = "Longitude", y = "Latitude", title = "Austrian Weather Stations: Neighborhood Structure") +
  theme_classic() +
  theme(legend.position = "none")

g_groups

## -----------------------------------------------------------------------------
GW = geo_weights(coordinates = weatherAUT2021[, c("lon", "lat")], 
                 groups = groups)
g_weights = g_groups + 
  labs(title = "Austrian Weather Stations: Weighting Matrix W") +
  geom_segment(aes(x = GW$centersN[4, 1], 
                   y = GW$centersN[4, 2], 
                   xend = GW$centersN[-4, 1]-0.1, 
                   yend = GW$centersN[-4, 2]-0.05,
                   alpha = GW$W[4, -4]), 
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 2,
               col = "blue")+
  geom_text(aes(x = GW$centersN[, 1],
                y = GW$centersN[, 2], 
                label = 1:length(GW$centersN[, 2])))

g_weights

## ----eval = FALSE-------------------------------------------------------------
# ? ssMRCD

## ----eval = TRUE, message = FALSE---------------------------------------------
set.seed(123)
out = ssMRCD(X = weatherAUT2021[, 1:6], 
             groups = groups, 
             weights = GW$W, 
             lambda = 0.5,
             tuning = NULL)
class(out)

## ----eval = TRUE, message = FALSE---------------------------------------------

# plot the tolerance ellipses in the geographical space
plots = plot(x = out, 
             type = c("convergence","ellipses", "ellipses_geo"),
             geo_centers = GW$centersN, 
             variables = c("s", "t"),
             manual_rescale = 0.001)

plots$plot_geoellipses +
  geom_sf(data = austria, fill = "transparent", color = "black")

plots$plot_ellipses
plots$plot_convergence

## ----eval = TRUE, message = FALSE---------------------------------------------
set.seed(123)
out = ssMRCD(X = weatherAUT2021[, 1:6], 
             groups = groups, 
             weights = GW$W, 
             tuning = list(method = "local contamination", 
                           plot = TRUE,
                           k = 10, 
                           coords = weatherAUT2021[, c("lon", "lat")],
                           cont = 0.05, 
                           repetitions = 3), 
             lambda = c(0.25, 0.5, 0.75))

## ----eval = TRUE, message = FALSE---------------------------------------------
set.seed(123)
out = ssMRCD(X = weatherAUT2021[, 1:6], 
             groups = groups, 
             weights = GW$W, 
             tuning = list(method = "residuals", 
                           plot = TRUE), 
             lambda = seq(0, 1, 0.1))

## ----eval = FALSE-------------------------------------------------------------
# ? locOuts

## ----eval = FALSE, message = FALSE--------------------------------------------
# set.seed(123)
# res = locOuts(data = weatherAUT2021[, 1:6],
#                             coords = weatherAUT2021[, c("lon", "lat")],
#                             groups = groups,
#                             lambda = 0.5,
#                             k = 10)
# summary(res)

## ----eval = FALSE,  message = FALSE-------------------------------------------
# cat(weatherAUT2021$name[res$outliers], sep = ",\n")

## ----eval = FALSE-------------------------------------------------------------
# ? plot.locOuts

## ----eval = FALSE,fig.dim = c(5,3.5)------------------------------------------
# plot(res, type = "hist")$p_hist

## ----eval = FALSE-------------------------------------------------------------
# plot(res, type = "spatial")$p_spatial +
#   geom_sf(data = austria, fill = "transparent", color = "black")

## ----eval = FALSE, fig.dim = c(7,3)-------------------------------------------
# # SCHOECKL
# plot(res, type = "pcp", observation = "SCHOECKL", scale = "zscore")$p_pcp

