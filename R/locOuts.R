# OUTLIER DETECTION

#' Local Outlier Detection Technique based on ssMRCD
#'
#' This function applies the local outlier detection method based on the spatially
#' smoothed MRCD estimator developed in Puchhammer and Filzmoser (2023).
#'
#' @param data data matrix with measured values.
#' @param coords matrix of coordinates of observations.
#' @param groups vector of neighborhood assignments.
#' @param lambda scalar used for spatial smoothing (see also \code{\link[ssMRCD]{ssMRCD}}).
#' @param weights weight matrix used in \code{\link[ssMRCD]{ssMRCD}}.
#' @param k integer, if given the \code{k} nearest neighbors per observations are used to calculate next distances. Default value is \code{k = NULL}.
#' @param dist scalar, if given the neighbors closer than given distance are used for next distances. If \code{dist} is given, \code{dist} is used, otherwise \code{k} is used.
#'
#' @return Returns an object of class \code{"locOuts"} with following components:\tabular{ll}{
#'    \code{outliers} \tab indices of found outliers. \cr
#'    \tab \cr
#'    \code{next_distance} \tab vector of next distances for all observations. \cr
#'    \tab \cr
#'    \code{cutoff} \tab upper fence of adjusted boxplot (see \code{\link[robustbase]{adjbox}}) used as cutoff value for next distances. \cr
#'    \tab \cr
#'    \code{coords} \tab matrix of observation coordinates.\cr
#'    \tab \cr
#'    \code{data} \tab matrix of observation values. \cr
#'    \tab \cr
#'    \code{groups} \tab vector of neighborhood assignments. \cr
#'    \tab \cr
#'    \code{k, dist} \tab specifications regarding neighbor comparisons. \cr
#'    \tab \cr
#'    \code{centersN} \tab coordinates of centers of neighborhoods. \cr
#'    \tab \cr
#'    \code{matneighbor} \tab matrix storing information which observations where used to calculate next distance for each observation (per row). 1 indicates it is used. \cr
#'    \tab \cr
#'    \code{ssMRCD} \tab object of class \code{"ssMRCD"} and output of \code{\link[ssMRCD]{ssMRCD}} covariance estimation. \cr
#' }
#'
#' @seealso See also functions \code{\link[ssMRCD]{ssMRCD}, \link[ssMRCD]{plot.locOuts}, \link[ssMRCD]{summary.locOuts}}.
#'
#' @references Puchhammer P. and Filzmoser P. (2023): Spatially smoothed robust covariance estimation for local outlier detection. \doi{10.48550/arXiv.2305.05371}
#' @export
#' @importFrom dbscan kNN
#'
#' @examples
#' # data construction
#' data = matrix(rnorm(2000), ncol = 4)
#' coords = matrix(rnorm(1000), ncol = 2)
#' groups = sample(1:10, 500, replace = TRUE)
#' lambda = 0.3
#'
#' # apply function
#' outs = local_outliers_ssMRCD(data = data,
#'                              coords = coords,
#'                              groups = groups,
#'                              lambda = lambda,
#'                              k = 10)
#' outs
local_outliers_ssMRCD = function(data, coords, groups, lambda, weights = NULL, k = NULL, dist = NULL){

  # check inputs
  data = as.matrix(data)
  coords = as.matrix(coords)
  groups = as.vector(groups)
  check_input(lambda, "lambda")
  if(!is.null(k)) check_input(k, "scalar")
  if(!is.null(dist)) check_input(dist, "scalar")

  n <- nrow(data)
  p <- ncol(data)
  N <- length(unique(groups))

  # weights
  weights_geo = geo_weights(coordinates = coords, groups = groups)
  centersN = weights_geo$centersN
  if(is.null(weights)){
    message("Weights for ssMRCD are not given. Default setting of geographic inverse-distance weights is used.")
    weights = weights_geo$W
  } else{
    check_input(weights, "W")
  }

  # neighborhood comparison matrix (matneighbors)
  distance_method = "k"
  if(is.null(k) & is.null(dist)){
    k = 10
    message("No k or dist given. Use default value of k = 10.")
  }
  if(is.null(k) & !is.null(dist)) distance_method = "dist"
  if(!is.null(k) & !is.null(dist)) {
    message("Both k and dist are given. Distance dist will be used.")
    distance_method = "dist"
    k = NULL
  }

  matneighbor = matrix(nrow = n, ncol = n)
  if(distance_method == "k"){
    knn_sorted = dbscan::kNN(coords, k = k)
    for(i in 1:n){
      matneighbor[i,] = as.numeric((1:n) %in% knn_sorted$id[i,])
    }
  }
  if(distance_method == "dist"){
    knn_sorted = dbscan::kNN(coords, k = n-1)
    for(i in 1:n){
      ni = sum(knn_sorted$dist[i,] < dist)
      matneighbor[i,] = as.numeric((1:n) %in% knn_sorted$id[i,1:ni])
    }
  }

  # ssMRCD covariances
  data_list = restructure_as_list(data, groups)
  cov_object = ssMRCD(x = data_list, weights = weights, lambda = lambda)

  center <- matrix(NA, N, p)
  cov <- matrix(NA, N, p^2)
  icov <- matrix(NA, N, p^2)
  for(i in 1:N) {
    center[i,] <- cov_object$MRCDmu[[i]]
    cov[i,] <- cov_object$MRCDcov[[i]]
    icov[i,] <- cov_object$MRCDicov[[i]]
  }

  # next distances
  next_distance <- rep(NA, n)
  for(i in 1:n){
    Ni = groups[i]
    cinv <- matrix(icov[Ni,], p, p)
    distn <- rep(NA, sum(matneighbor[i, ] != 0))
    kk = 1
    for(j in (1:n)[matneighbor[i, ] != 0]) {
      distn[kk] <- t(data[i, ] - data[j, ]) %*% cinv %*% (data[i, ] - data[j, ])
      kk = kk + 1
    }
    next_distance[i] <- sort(distn)[1]
  }

  # outlier classification
  box <- robustbase::adjbox(next_distance, horizontal=TRUE, plot=FALSE, range = 1.5)
  outliers = which(next_distance > box$fence[2])

  out = list(outliers = outliers,
             next_distance = next_distance,
             cutoff = box$fence[2],
             coords = coords,
             data = data,
             groups = groups,
             k = k,
             dist = dist,
             centersN = centersN,
             matneighbor = matneighbor,
             ssMRCD = cov_object)
  class(out) = c("locOuts", "list")
  return(out)
}




