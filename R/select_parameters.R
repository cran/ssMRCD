# TUNING FUNCTIONS:
# Two different methods:
#     - for local outlier detection based on swapping
#     - for PCA based on residuals

##########################################################################################
#' Optimal Smoothing Parameter for ssMRCD based on Local Outliers
#'
#' This function provides insight into the effects of different parameter settings.
#'
#' @param data matrix with observations.
#' @param coords matrix of coordinates of these observations.
#' @param groups numeric vector, the neighborhood structure that should be used for \code{\link[ssMRCD]{ssMRCD}}.
#' @param lambda scalar, the smoothing parameter.
#' @param weights weighting matrix used in \code{\link[ssMRCD]{ssMRCD}}.
#' @param k vector of possible k-values to evaluate.
#' @param dist vector of possible dist-values to evaluate.
#' @param cont level of contamination, between 0 and 1.
#' @param repetitions number of repetitions wanted to have a good picture of the best parameter combination.
#'
#' @return Returns a matrix of average false-negative rate (FNR) values and the total number of outliers found by the method as aproxy for the false-positive rate.
#' Be aware that the FNR does not take into account that there are also natural outliers included in the data set that might or might not be found.
#' Also a plot is returned representing these average.
#' The best parameter selection depends on the goal of the analysis.
#' @export
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \donttest{
#' # get data set
#' data("weatherAUT2021")
#'
#' # make neighborhood assignments
#' cut_lon = c(9:16, 18)
#' cut_lat = c(46, 47, 47.5, 48, 49)
#' N = ssMRCD::groups_gridbased(weatherAUT2021$lon, weatherAUT2021$lat, cut_lon, cut_lat)
#' table(N)
#' N[N == 2] = 1
#' N[N == 3] = 4
#' N[N == 5] = 4
#' N[N == 6] = 7
#' N[N == 11] = 15
#' N = as.numeric(as.factor(N))
#'
#' # tune parameters
#' set.seed(123)
#' parameter_tuning(data = weatherAUT2021[, 1:6 ],
#'                  coords = weatherAUT2021[, c("lon", "lat")],
#'                  groups = N,
#'                  lambda = c(0.5, 0.75),
#'                  k = c(10),
#'                  repetitions = 1)
#'}

parameter_tuning = function(data,
                            coords,
                            groups,
                            lambda = c(0, 0.25, 0.5, 0.75, 0.9),
                            weights = NULL,
                            k = NULL,
                            dist = NULL,
                            cont = 0.05,
                            repetitions = 5){

  data = as.matrix(data)
  coords = as.matrix(coords)
  check_input(groups, "vector")
  check_input(lambda, "vector")
  if(!is.null(k)) check_input(k, "vector")
  if(!is.null(dist)) check_input(dist, "vector")


  # distance method
  distance_method = "k"
  parameter_combinations = expand.grid(lambda, k, 1:repetitions)
  if(is.null(k) & is.null(dist)){
    stop("Neither k nor dist are given. Please specify a value range for one of them.")
  }
  if(is.null(k) & !is.null(dist)) {
    distance_method = "dist"
    parameter_combinations = expand.grid(lambda, dist, 1:repetitions)
  }
  if(!is.null(k) & !is.null(dist)) {
    stop("Both k and dist are given. Please set the one you do not like to use to NULL.")
  }
  colnames(parameter_combinations) = c("lambda", "k_dist", "reps")


  # parameter settings
  n_par = dim(parameter_combinations)[1]

  n_outliers_found = rep(NA, n_par)
  FNR = rep(NA, n_par)
  reps = -1

  for(i in 1:n_par){

    # data construction
    if(reps != parameter_combinations[i, 3]){
      reps = parameter_combinations[i, 3]

      # contamination
      contaminated = contamination_random(cont = cont, data = data)
      data_contam = contaminated$data
      outliers = contaminated$out
    }

    # performance
    lam = parameter_combinations[i, 1]
    k_dist = parameter_combinations[i, 2]
    if(distance_method == "k"){
      tmp = local_outliers_ssMRCD(data_contam,
                                  coords,
                                  groups,
                                  lambda = lam,
                                  k = k_dist,
                                  dist = NULL,
                                  weights = weights)
    }
    if(distance_method == "dist"){
      tmp = local_outliers_ssMRCD(data_contam,
                                  coords,
                                  groups,
                                  lambda = lam,
                                  dist = k_dist,
                                  k = NULL,
                                  weights = weights)
    }
    outs_method = tmp$outliers
    n_outliers_found[i] = length(outs_method)
    FNR[i] = 1 - sum(outliers %in% outs_method)/length(outliers)
  }

  param_results = cbind(parameter_combinations, FNR, n_outliers_found) %>%
    dplyr::group_by(lambda, k_dist) %>%
    dplyr::summarise(mean_FNR  = mean(FNR ),
                     mean_outs = mean(n_outliers_found))

  plot = ggplot() +
    geom_point(aes(x = param_results$mean_FNR,
                   y = param_results$mean_outs,
                   col = as.factor(param_results$lambda),
                   shape = as.factor(param_results$k_dist)),
               size = 2,
               alpha = 0.6) +
    scale_colour_discrete(name = "\U03BB") +
    scale_shape_discrete(name = "k/dist") +
    theme_bw() +
    labs(x = "FNR",
         y = "# outliers detected",
         title = paste("Parameter Tuning for", repetitions, "Repetitions"))

  return(list(plot = plot, values = param_results))

}


#' Contamination Through Swapping
#'
#' This function swaps observations completely random in order to introduce contamination
#' in the data. Used in \code{\link[ssMRCD]{parameter_tuning}}.
#'
#' @param cont numeric, amount of contamination in data.
#' @param data data whose observations should be switched.
#'
#' @return A matrix with switched observations.
#' @export
#'
#' @examples
#' # set seed
#' set.seed(1)
#'
#' # get data
#' data(weatherAUT2021)
#'
#' # switch 5% of observations
#' contamination_random(cont = 0.05, data = weatherAUT2021[,1:6])

contamination_random = function(cont, data){

  check_input(cont, "scalar")
  data = as.matrix(data)

  n = dim(data)[1]
  n_cont = 2*round(n*cont/2)

  # select outliers
  outliers = sample(1:n, size = n_cont , replace = FALSE)

  # switch first outlier observation with last outlier and so forth
  for (i in 1:(n_cont/2)){
    tmp = data[outliers[i], ]
    data[outliers[i],] =  data[outliers[n_cont - (i-1)],]
    data[outliers[n_cont - (i-1)],] = tmp
  }
  data_switched = data

  return(list(data = data_switched, out = outliers))

}



##########################################################################################
#' Optimal Smoothing Parameter for ssMRCD based on Residuals
#'
#' The optimal smoothing value for the ssMRCD estimator is based on the residuals and the
#' trimmed mean of the norm.
#'
#' @param X data matrix containing observations.
#' @param groups grouping vector corresponding to \code{X}.
#' @param weights weight matrix for groups, see \code{\link[ssMRCD]{rescale_weights}}, and \code{\link[ssMRCD]{geo_weights}}.
#' @param lambda vector of parameter values for smoothing, between 0 and 1.
#' @param TM target matrix, if not given MCD (or MRCD if non regular) is used with default values and \code{alpha}.
#' @param alpha percentage of outliers to be expected.
#' @param seed seed for ssMRCD calculations.
#' @param return_all logical, if FALSE the function returns only the optimal lambda.
#' @param cores integer, number of cores used for parallel computing.
#'
#' @return \tabular{ll}{
#'    \code{lambda_opt} \tab optimal lambda for smoothing. \cr
#'    \tab \cr
#'    \code{COVS} \tab \code{ssMRCD} object with optimal parameter setting. \cr
#'    \tab \cr
#'    \code{plot} \tab plot for optimal parameter setting. \cr
#'    \tab \cr
#'    \code{residuals} \tab mean of norm of residuals for varying lambda. \cr
#' }
#' @export
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom grDevices recordPlot
#' @importFrom graphics points
#' @importFrom rrcov CovMcd
#' @importFrom rrcov CovMrcd
#'
#' @examples
#' # create data set
#' x1 = matrix(runif(200), ncol = 2)
#' x2 = matrix(rnorm(200), ncol = 2)
#'
#' # create weighting matrix
#' W = matrix(c(0, 1, 1, 0), ncol = 2)
#'
#'select_smoothing (X = rbind(x1, x2),
#'                  groups = rep(c(1,2), each = 100),
#'                  weights = W,
#'                  lambda = seq(0, 1, 0.1),
#'                  return_all = TRUE,
#'                  cores = 1)

select_smoothing = function(X,
                            groups,
                            weights,
                            lambda = seq(0, 1, 0.1),
                            TM = NULL,
                            alpha = 0.75,
                            seed = 123436,
                            return_all = TRUE,
                            cores = 1){

  X = as.matrix(X)
  groups = as.numeric(groups)
  weights = as.matrix(weights)

  x_data = restructure_as_list(data = X, groups)

  if(is.null(TM)) {
    TM = tryCatch({
      round(rrcov::CovMcd(X, alpha = alpha)$cov, 7)
    }, error = function(cond){
      coloured_print("Use MRCD as target matrix.")
      round(rrcov::CovMrcd(X, alpha = alpha)$cov, 7)
    })
  }
  maxcond = eigen(TM)$values[1]/ eigen(TM)$values[dim(TM)[1]] + 100

  i = 0
  COV_list = list()
  if(cores !=  1) {
    cl = parallel::makeCluster(cores)
    COV_list = foreach::foreach(i = 1:length(lambda),
                                .combine = "c",
                                .packages = c("ssMRCD")) %dopar% {
                                  set.seed(seed)

                                  COVS = ssMRCD(x = x_data,
                                                weights = weights,
                                                maxcond = maxcond,
                                                TM = TM,
                                                lambda = lambda[i],
                                                alpha = alpha)
                                  list(COVS)
                                }
    parallel::stopCluster(cl)
  }

  if(cores == 1){
    for(i in 1:length(lambda)){
      COVS = ssMRCD(x = x_data,
                    weights = weights,
                    maxcond = maxcond,
                    TM = TM,
                    lambda = lambda[i],
                    alpha = alpha)
      COV_list = c(COV_list, list(COVS))
    }
  }

  names(COV_list) = paste0("lambda", lambda)
  resid = sapply(X = COV_list,
                 function(x) residuals.ssMRCD(object = x,
                                              remove_outliers = FALSE,
                                              mean = TRUE))
  names(resid) = paste0("lambda", lambda)

  lambda_opt = lambda[which.min(resid)]

  plot(lambda,
       resid,
       col = "black",
       type = "l",
       main = "Optimal smoothing for ssMRCD",
       ylab = "trimmed mean of residuals")
  graphics::points(lambda_opt,
                   min(resid),
                   col = "red",
                   pch = 43)
  plot_optimallambda = grDevices::recordPlot()

  if(return_all) return(list(lambda_opt = lambda_opt,
                             COVS = COV_list[[which(lambda == lambda_opt)]],
                             plot = plot_optimallambda,
                             residuals = resid))
  return(lambda_opt)
}


