## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  warning = FALSE, 
  fig.dim = c(7, 4.5),
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ssMRCD)
library(dplyr)
library(ggplot2)
library(robustbase)
library(tidyr)
library(ggridges)

## ----eval = FALSE-------------------------------------------------------------
# ? weatherHoheWarte

## ----load data----------------------------------------------------------------
data("weatherHoheWarte")
head(weatherHoheWarte)

## ----data preparation---------------------------------------------------------
varnames_short= c("cl", "rad", "vp", "wmax", "ap", "hum", "prec", "sight", "sun", 
               "tmax", "tmin", "t", "w")

# filter data
weather =  weatherHoheWarte %>% 
  dplyr::filter(year > 2000)

# select variables
X_data = weather %>% 
  select(cloud_cover:wind_velocity)  %>%
  as.matrix
colnames(X_data) = varnames_short

# get number of variables
p = dim(X_data)[2]

## ----scale data---------------------------------------------------------------
X_data = scale(X_data, center = robustbase::colMedians(X_data), scale = apply(X = X_data, MARGIN = 2, FUN = mad))

## ----groups and weights-------------------------------------------------------
N_groups = as.numeric(as.factor(weather$year))
N = max(N_groups)
W = time_weights(max(N_groups), c(5:1))

## ----hyperparameter tuning, eval = FALSE--------------------------------------
# set.seed(1234)
# 
# lambda = seq(0, 1, by = 0.05)
# opt_lambda = ssMRCD(X = X_data,
#                     groups = N_groups,
#                     weights = W,
#                     lambda = lambda,
#                     tuning = list(method = "residuals", plot = TRUE),
#                     alpha = 0.75)

## -----------------------------------------------------------------------------
lambda_grid = seq(0, 1, by = 0.05)
residual_values = c(3.199500,3.188801,3.182297,3.175784,3.171462,3.167695,3.165151,3.162758,
                    3.161896,3.161929,3.162562, 3.163831,3.165651,3.168228,3.171867,3.175790,
                    3.180868,3.186437,3.192634,3.200292,3.208845)
lambda_opt = 0.4

ggplot() + 
  geom_line(aes(x = lambda_grid,
                y = residual_values)) +
  geom_point(aes(x = lambda_grid,
                y = residual_values)) +
  geom_point(aes(x = lambda_opt,
                y = residual_values[lambda_grid == lambda_opt]),
             col = "red") +
  labs(x = expression(lambda),
       y = "Residual Value",
       title = "Optimal Smoothing Based on Residuals") +
  theme_minimal()

## ----ssMRCD-------------------------------------------------------------------
# get optimal smoothing parameter and estimator
set.seed(123)
ssMRCD_weather = ssMRCD(X = X_data,
                        groups = weather$year,
                        weights = W,
                        lambda = 0.4,
                        alpha = 0.75)

# collect covariance matrices
COVS = ssMRCD_weather$MRCDcov

plot(ssMRCD_weather, type = c("ellipses", "convergence"), variables = c("hum", "sun"))

## ----dynamic means------------------------------------------------------------
# collect means
mu_timeseries = do.call(cbind, ssMRCD_weather$MRCDmu)
colnames(mu_timeseries) = 2001:2023
rownames(mu_timeseries) = colnames(X_data)

# plot means dependent on year
ggplot(mu_timeseries %>%
         data.frame %>%
         mutate(Var1 = rownames(.)) %>%
         tidyr::pivot_longer(cols = X2001:X2023) %>% 
         group_by(Var1) %>% 
         mutate(value = (value-min(value))/(max(value) - min(value)),
                name = as.numeric(gsub("X", "", name)))) +
  geom_line(aes(x = name, y = value, group = Var1)) + 
  geom_smooth(aes(x = name, y = value, group = Var1), se = F) + 
  facet_wrap(vars(Var1)) + 
  theme_minimal() +
  labs(x = "", y = "")

## ----pca tuning, eval = TRUE--------------------------------------------------
set.seed(12345)
pca = msPCA(eta = seq(0, 3, 0.25), gamma = 0.5, COVS = COVS, k = p)
pca$plot_tpo

## -----------------------------------------------------------------------------
summary(pca, k = 3)

## ----scree plot---------------------------------------------------------------
# plot scree plot
screeplot(x = pca, 
                 text = TRUE,
                 size = 3, 
                 gnames = 2001:2023,
                 cutoff = 0.8)

## ----optimal k----------------------------------------------------------------
# optimal number of components
optimal_k = 3

## ----align loadings-----------------------------------------------------------
?align
pca$PC[, 1:optimal_k, ] = align(PC = pca$PC[, 1:optimal_k, ],type = "mean")

## ----plot loadings------------------------------------------------------------
# plot aligned loadings
plot(x = pca, 
     type = "loadings",
     gnames = 2001:2023, 
     k = 1) 

plot(x = pca, 
     type = "loadings",
     gnames = 2001:2023, 
     k = 2) 

plot(x = pca, 
     type = "loadings", 
     gnames = 2001:2023, 
     k = 3) 

## ----biplot-------------------------------------------------------------------
biplot(pca)

## ----scores-------------------------------------------------------------------
sc = scores(X = X_data, groups = weather$year, PC = pca$PC, mu = ssMRCD_weather$MRCDmu, Sigma = ssMRCD_weather$MRCDcov)
sc2 = scores(PC = pca$PC, ssMRCD = ssMRCD_weather)

## -----------------------------------------------------------------------------
plot(x = pca, ssMRCD = ssMRCD_weather, type = "scores")

## ----distance distribution----------------------------------------------------
# distances
score_sd = sc$SD
score_od = sc$OD

dat = data.frame(time = weather$time,
                 year = weather$year,
                 OD = score_od,
                 SD = score_sd,
                 groups = weather$year)
dat_melt = dat %>% tidyr::pivot_longer(cols = c("OD", "SD"))

ggplot(dat_melt %>% 
         dplyr::filter(name %in% c( "SD", "OD"))) +
  ggridges::geom_density_ridges_gradient(aes(x = value, 
                                   y = groups, 
                                   fill = as.factor(groups), 
                                   group = groups),
                               scale = 10,
                               rel_min_height = 0.01, 
                               gradient_lwd = 0.5) +
  scale_x_continuous(expand = c(0.03, 0)) +
  scale_y_discrete(expand = c(0.0, 0, 0.1, 0)) +
  scale_fill_discrete(name = "") +
  ggridges::theme_ridges(font_size = 13, grid = TRUE) + 
  theme(axis.title.y = element_blank(),
        legend.position = "right",
        panel.grid.major.y = element_line()
        ) + 
  labs(x = "") +
  facet_grid(cols = vars(name), scales = "free") 

## ----distance distance plot---------------------------------------------------
plot(pca, type = "score_distances", ssMRCD = ssMRCD_weather, k = 3)

## ----gif, eval = FALSE--------------------------------------------------------
# # use function saveGIF
# library(animation)
# param = 2001:2023
# 
# saveGIF(interval = 0.1,
#         movie.name = "PCA_scoreDist_bivariate.gif",
#         expr = {
#           for (i in param ){
#            g = ggplot(dat %>% dplyr::filter (year %in% seq(i-5, i+5, 1))) +
#                  geom_density_2d_filled(aes(x = SD,
#                      y = OD),
#                      alpha = 1) +
#                   scale_x_sqrt(limits = c(0,40)) +
#                   ylim(c(0, 6.2)) +
#                   labs(title = i)
# 
#             print(g)
#           }}
#         )

