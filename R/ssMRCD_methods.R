# METHODS FOR cov_ssMRCD OBJECT

##########################################################################################
#' Plot Method for ssMRCD Object
#'
#' Produces diagnostic plots for an object of class \code{"ssMRCD"} including convergence behavior and visualizations of the covariance matrices.
#'
#' @param x An object of class \code{"ssMRCD"}.
#' @param type Character string or vector specifying the type of plot(s) to produce. Available options are \code{"convergence"}, \code{"ellipses"}, and \code{"ellipses_geo"}. See Details.
#' @param variables A character vector of length 2 specifying the variable names (columns of the data) used to compute and plot the covariance ellipses.
#' @param geo_centers A matrix specifying the spatial/geographical coordinates of the group centers. Required when \code{type = "ellipses_geo"}.
#' @param manual_rescale Numeric scaling factor to adjust the size of ellipses in \code{"ellipses_geo"} plots.
#' @param tolerance Numeric value (between 0 and 1) specifying the quantile used to define tolerance ellipses. Default is \code{0.95}.
#' @param ... Further arguments passed to plotting functions.
#'
#' @details
#' \strong{type = "convergence"}:
#' Displays the convergence behavior of the objective function across C-step iterations for each initialization. Red lines indicate non-monotonic convergence.
#'
#' \strong{type = "ellipses"}:
#' Shows Mahalanobis tolerance ellipses for each group based on their estimated covariance matrix. Only the two variables specified in \code{variables} are visualized. The global MCD ellipse may be shown for comparison (if included elsewhere).
#'
#' \strong{type = "ellipses_geo"}:
#' Projects group-wise covariance ellipses onto a geographical coordinate system (e.g., spatial map), using the positions given in \code{geo_centers}. The ellipses are scaled using \code{manual_rescale} and drawn using the same two variables as for \code{type = "ellipses"}.
#'
#' @return A named list of ggplot2 plot objects:
#' \item{plot_convergence}{Plot showing convergence diagnostics (if \code{"convergence"} selected).}
#' \item{plot_ellipses}{Plot of covariance ellipses in variable space (if \code{"ellipses"} selected).}
#' \item{plot_geoellipses}{Plot of covariance ellipses in geographical space (if \code{"ellipses_geo"} selected).}
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(2000), ncol = 4)
#' colnames(data) <- paste0("V", 1:4)
#' coords <- matrix(rnorm(1000), ncol = 2)
#' groups <- sample(1:10, 500, replace = TRUE)
#' lambda <- 0.3
#'
#' outs <- locOuts(data = data,
#'                               coords = coords,
#'                               groups = groups,
#'                               lambda = lambda,
#'                               k = 10)
#'
#' # Plot convergence
#' plot(x = outs$ssMRCD, type = "convergence")
#'
#' # Plot ellipses in variable space
#' plot(x = outs$ssMRCD, type = "ellipses", variables = c("V1", "V2"))
#'
#' # Plot ellipses in geographical space
#' centers <- matrix(rnorm(20), ncol = 2)  # example centers for 10 groups
#' plot(x = outs$ssMRCD, type = "ellipses_geo", geo_centers = centers, variables = c("V1", "V2"))
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}},
#'   \code{\link[ssMRCD]{locOuts}},
#'   \code{\link[ssMRCD]{plot.locOuts}}
#'
#' @exportS3Method plot ssMRCD
#' @importFrom ellipse ellipse
#' @import ggplot2
plot.ssMRCD = function(x,
                       type = c("convergence", "ellipses", "ellipses_geo"),
                       variables = NULL,
                       geo_centers = NULL,
                       manual_rescale = 1,
                       tolerance = 0.95,
                       ...){



  N = length(x$MRCDcov)
  p = ncol(x$mT)

  g_conv = NULL
  if("convergence" %in% type){

    # prepare data
    n_starts = nrow(x$obj_fun_values)
    max_steps = ncol(x$obj_fun_values)

    decr = which(apply(X = x$obj_fun_values,
                       MARGIN = 1,
                       FUN = monotonic_decreasing))

    tmp = cbind(values = c(x$obj_fun_values),
                start = rep(1:n_starts, times = max_steps),
                step = rep(1:max_steps, each = n_starts))
    tmp = na.omit(tmp)

    # plot
    g_conv = ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = tmp[,"step"],
                    y = tmp[, "values"],
                    group = tmp[, "start"],
                    col = tmp[, "start"] %in% decr),
                alpha = 0.3) +
      ggplot2::scale_color_manual(name = "",
                         values = c("FALSE" = "red", "TRUE"= "grey"),
                         labels = c("FALSE" = "not descreasing", "TRUE" = "decreasing")) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Iteration Step",
                    y = "Objective Value",
                    title = "Convergence (tolerance 1e-08)")
  }

  g_ell = NULL
  if("ellipses" %in% type){

    if(is.null(variables)) variables = colnames(x$MRCDcov[[1]])[1:2]

    ell = matrix(NA, nrow = 0, ncol = 3)
    for(i in 1:N){
      elltmp = data.frame(ellipse::ellipse(x$MRCDcov[[i]][variables, variables],
                             centre = x$MRCDmu[[i]][variables],
                             level = tolerance))
      ell = rbind(ell, cbind(elltmp, x$gnames[i]))
    }

    g_ell = ggplot2::ggplot() +
      ggplot2::geom_polygon(aes(x = ell[, 1],
                    y = ell[, 2],
                    group = ell[, 3],
                    fill = as.factor(ell[,3]),
                    col = as.factor(ell[,3])),
                   alpha = 0.05) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_discrete(name = "Groups") +
      ggplot2::scale_fill_discrete(name = "Groups") +
      ggplot2::labs(x = variables[1],
           y = variables[2],
           title = "Tolerance Ellipses (95%)")

  }

  g_ellgeo = NULL
  if("ellipses_geo" %in% type){

    if(is.null(variables)) variables = colnames(x$MRCDcov[[1]])[1:2]

    if(is.null(geo_centers)) {
      stop("You need to specify the centers of the neighborhoods geo_centers to plot covariance ellipses in the geographical space.")
    }

    ellgeo = matrix(NA, nrow = 0, ncol = 3)
    for(i in 1:N){
      elltmp = data.frame(ellipse::ellipse(x$MRCDcov[[i]][variables, variables],
                                centre = geo_centers[i, ],
                                level = 0.5*manual_rescale))
      ellgeo = rbind(ellgeo, cbind(elltmp, x$gnames[i]))
    }

    g_ellgeo = ggplot2::ggplot() +
      ggplot2::geom_path(aes(x = ellgeo[,1],
                    y = ellgeo[,2],
                    group = ellgeo[, 3],
                    col = as.factor(ellgeo[, 3]))) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_discrete(name = "Groups") +
      ggplot2::labs(x = paste0("x | ",  variables[1]),
           y = paste0("y | ",  variables[2]),
           title = paste0("Tolerance Ellipses"))
  }

  return(list("plot_convergence" = g_conv,
              "plot_ellipses" = g_ell,
              "plot_geoellipses" = g_ellgeo))
}


##########################################################################################
#' Locally Center and/or Scale or Data Using an ssMRCD Object
#'
#' Applies local standardization (scaling and/or centering) of either the original data
#' from an \code{ssMRCD} object or new data provided via the \code{X} argument,
#' using group-wise robust means and variances from the ssMRCD estimation.
#'
#' @param x An object of class \code{"ssMRCD"}. See \code{\link[ssMRCD]{ssMRCD}}.
#' @param ... List of additional arguments including:
#' \describe{
#'   \item{\code{X}}{A numeric matrix or data frame containing new observations to be scaled. If not provided, the data stored in the \code{ssMRCD} object is used.}
#'   \item{\code{groups}}{An integer vector from 1 to number of groups of group assignments corresponding
#'    to the rows of \code{X}. If \code{X} is not provided, defaults to the group
#'    assignments used in the original ssMRCD estimation.}
#'   \item{\code{center_only}}{Logical. If \code{TRUE}, only centering is
#'   applied; if \code{FALSE}, both centering and scaling are applied.
#'   Default is \code{FALSE}.}
#' }
#'
#' @return A numeric matrix of the same dimension as \code{X}, where each observation has
#' been standardized (or centered) using the corresponding group-wise robust mean and
#' (if applicable) variance from the ssMRCD model.
#' If \code{X = NULL}, the original data from the ssMRCD object is returned in scaled form,
#' sorted according to group labels.
#'
#' @details
#' For each group, the function applies scaling (or just centering) using the robust
#' location and scale (square root of the diagonal of the covariance) estimates obtained
#' during ssMRCD estimation.
#'
#' @examples
#' # Simulated example
#' x1 <- matrix(runif(200), ncol = 2)
#' x2 <- matrix(rnorm(200), ncol = 2)
#' x <- list(x1, x2)
#'
#' W <- matrix(c(0, 1, 1, 0), ncol = 2)
#' localCovs <- ssMRCD(x, weights = W, lambda = 0.5)
#'
#' # Scale original data
#' sc = scale(localCovs)
#'
#' # Scale new observations
#' sc = scale(localCovs,
#'            list(X = matrix(rnorm(20), ncol = 2, nrow = 10),
#'            groups = rep(2, 10)))
#'
#' # Center only
#' sc = scale(localCovs,
#'            list(X = matrix(rnorm(20), ncol = 2, nrow = 10),
#'            groups = rep(2, 10),
#'            center_only = TRUE))
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}}
#'
#' @exportS3Method scale ssMRCD
scale.ssMRCD = function(x, ...){

  args = list(...)
  N = x$N

  if(is.null(args$X) | is.null(args$groups)) {
    args$X = do.call(rbind, x$mX)
    args$groups = rep(1:N, times = sapply(X = x$mX, nrow))
  }
  if(is.null(args$center_only)) args$center_only = FALSE

  X = as.matrix(args$X)

  for(i in 1:N){
    ind = which(args$groups == i)
    if(length(ind)!= 0){
      if(!args$center_only){
        X[ind,] = scale(x = X[ind,], center = x$MRCDmu[[i]], scale = sqrt(diag(x$MRCDcov[[i]])))
      } else {
        X[ind,] = scale(x = X[ind,], center = x$MRCDmu[[i]])
      }
    }
  }

  return(X)
}


##########################################################################################
# #' Get Indices of Non-Local Outliers
# #'
# #' This function extracts indices of non-local outliers from the `ssMRCD` object. The indices are returned as a vector, indicating the positions of the non-local outliers in the sorted groups.
# #'
# #' @param ssMRCD An object containing the results of the local outlier detection, which includes:
# #'   \itemize{
# #'     \item `hset`: A list of indices for each neighborhood, indicating the positions of local outliers.
# #'     \item `mX`: A list of matrices, where each matrix represents the data for a neighborhood.
# #'     \item `N`: The number of neighborhoods.
# #'   }
# #'
# #' @return A numeric vector containing the indices of non-local outliers. The indices are sorted according to the original order of the groups.
# #'
# #' @details
# #' The function aggregates the indices of non-local outliers across all neighborhoods. It constructs the indices based on the structure of the `ssMRCD` object, taking into account the sizes of the neighborhoods.
#
# #' @seealso \code{\link[ssMRCD]{ssMRCD}}
# #' @keywords internal
#
# outliers = function(ssMRCD){
#   hsets = ssMRCD$hset
#   size_ns = sapply(ssMRCD$mX, function(x) dim(x)[1])
#   N = ssMRCD$N
#
#   ind = c()
#   size_sum = 0
#   for(i in 1:N){
#     ind = c(ind, hsets[[i]] + size_sum)
#     size_sum = size_sum + size_ns[i]
#   }
#
#   return(ind)
# }

##########################################################################################
#' Residual Method from an ssMRCD Object
#'
#' Computes group-wise Mahalanobis residuals (standardized distances) using the robust local covariance and location estimates from an \code{ssMRCD} object. Residuals can be computed for the fitted data or for new data, and optionally summarized as a trimmed mean.
#'
#' @param object An object of class \code{"ssMRCD"}, typically the result of \code{\link[ssMRCD]{ssMRCD}}.
#' @param ... Additional arguments, see Details.
#'
#' @details
#' The function supports several modes of use, controlled by the \code{type} argument in \code{...}:
#' \describe{
#'   \item{\code{type}}{\code{"residuals"} (default), \code{"trimmed_mean"}, or \code{"additional_data"}.}
#'   \item{\code{X}}{A numeric matrix of new observations to compute residuals for. Required if \code{type = "additional_data"}.}
#'   \item{\code{groups}}{A vector of group assignments for the new data in \code{X}. Required if \code{type = "additional_data"}.}
#'   \item{\code{alpha}}{A numeric value (default taken from the \code{ssMRCD} object if missing) indicating the quantile for trimmed mean calculation. Only used if \code{type = "trimmed_mean"}.}
#' }
#'
#' Notes:
#' \itemize{
#'   \item If \code{type = "residuals"}, residuals are computed for the original data stored
#'    in the \code{ssMRCD} object.
#'   \item If \code{type = "additional_data"}, both \code{X} and \code{groups} must be
#'   provided. All residuals of \code{X} are returned (i.e., \code{alpha = 1} is used internally).
#'   \item If \code{type = "trimmed_mean"}, the mean of the \code{alpha} proportion of
#'   smallest residual norms is returned. This is also used for parameter tuning.
#' }
#'
#' @return
#' Depending on the \code{type}:
#' \describe{
#'   \item{\code{"residuals"} or \code{"additional_data"}}{A numeric matrix of residuals.}
#'   \item{\code{"trimmed_mean"}}{A single numeric value: the trimmed mean of residual norms.}
#' }
#'
#' @examples
#' # Create data
#' x1 <- matrix(runif(200), ncol = 2)
#' x2 <- matrix(rnorm(200), ncol = 2)
#' x <- list(x1, x2)
#'
#' # Define neighborhood weights
#' W <- matrix(c(0, 1, 1, 0), ncol = 2)
#'
#' # Compute ssMRCD
#' localCovs <- ssMRCD(x, weights = W, lambda = 0.5)
#'
#' # Residuals for original data (all)
#' residuals(localCovs, type = "residuals")
#'
#' # Trimmed mean of residual norms
#' residuals(localCovs, type = "trimmed_mean", alpha = 0.8)
#'
#' # Residuals for new data
#' newX <- matrix(rnorm(20), ncol = 2, nrow = 10)
#' newGroups <- rep(2, 10)
#' residuals(localCovs, type = "additional_data", X = newX, groups = newGroups)
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}}
#'
#' @exportS3Method residuals ssMRCD
#' @importFrom expm sqrtm
residuals.ssMRCD = function(object, ...){

  args = list(...)

  # type: "trimmed_mean", "residuals", "additional_data"
  if(is.null(args$type)) args$type = "residuals"
  if(args$type == "additional_data" & (is.null(args$X) | is.null(args$groups))) {
    stop("Type is 'additional data' but no additional data X or additional groups assignemnts in groups are provided.")
  }

  if(args$type == "trimmed_mean" & is.null(args$alpha)) args$alpha = object$alpha
  if(args$type == "residuals" | args$type == "additional_data")  args$alpha = 1

  N = length(object$MRCDcov)
  p = ncol(object$MRCDcov[[1]])
  n = sum(sapply(object$mX, nrow))

  ind = 1:n

  # calculate residuals
  if(args$type != "additional_data"){
    residuals_model = matrix(NA, nrow = n, ncol = p)
    groups = rep(1:N, times = sapply(object$mX, nrow))
    for(i in 1:N){
      ind = which(groups == i)
      if(length(ind)!= 0){
        centered = sweep(object$mX[[i]], 2, object$MRCDmu[[i]], "-")
        residuals_model[ind,] = t(apply(centered, 1, function(x) expm::sqrtm(object$MRCDicov[[i]]) %*% x))
      }
    }
  }

  if(args$type == "additional_data"){
    residuals_model = matrix(NA, nrow = nrow(args$X), ncol = p)
    for(i in 1:N){
      ind = which(args$groups == i)
      if(length(ind)!= 0){
        centered = sweep(args$X[ind, ,drop = FALSE], 2, object$MRCDmu[[i]], "-")
        residuals_model[ind, ] = t(apply(centered, 1, function(x) expm::sqrtm(object$MRCDicov[[i]]) %*% x))
      }
    }
  }

  if(args$type %in% c("residuals", "additional_data")) return(residuals_model)

  # calculate trimmed mean of norm
  if(args$type == "trimmed_mean"){
    res_norm = sqrt(diag(residuals_model %*% t(residuals_model)))
    ind = sort.int(res_norm, index.return = T)$ix[1:round(n*args$alpha)]
    res_trimmed = res_norm[ind]

    return(mean(res_trimmed))
  }

}


