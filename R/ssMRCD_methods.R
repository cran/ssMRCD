# METHODS FOR cov_ssMRCD OBJECT


##########################################################################################
#' Summary Method for ssMRCD Object
#'
#' Summarises most important information of output \code{\link[ssMRCD]{ssMRCD}}.
#'
#' @param object object of class \code{"ssMRCD"}, output of \code{\link[ssMRCD]{ssMRCD}}.
#' @param ... further parameters.
#'
#' @return Prints a summary of the \code{ssMRCD} object.
#'
#' @seealso See also \code{\link[ssMRCD]{ssMRCD}}, \code{\link[ssMRCD]{plot.ssMRCD}}.
#'
#' @exportS3Method summary ssMRCD
summary.ssMRCD = function(object, ...){

  # Additional info
  cat("ss-MRCD with", object$numiter,"many iterations. \nParameter Setting:", "\n")
  cat("\t * N =", object$N, "\n")
  cat("\t * n = (")
  cat(unlist(lapply(X = object$mX,  FUN = function(x) dim(x)[1])), sep = ", ")
  cat(")\n")
  cat("\t * \U03BB =", object$lambda, "\n")
  cat("\t * \U03C1 = (")
  cat(object$rho, sep = ", ")
  cat(")\n\n")

  # Covariance
  cat("Covariance Estimator:\n")
  names(object$MRCDcov) = paste0("N", 1:object$N)
  print(object$MRCDcov)

  # Mean
  cat("Mean Estimator:\n")
  names(object$MRCDmu) = paste0("N", 1:object$N)
  print(object$MRCDmu)
}



##########################################################################################
#' Plot Method for ssMRCD Object
#'
#' Plots diagnostics for function output of \code{\link[ssMRCD]{ssMRCD}} regarding convergence behavior
#' and the resulting covariances matrices.
#'
#' @param x object of class \code{"ssMRCD"}.
#' @param type type of plot, possible values are \code{"convergence"} and \code{"ellipses"}. See details.
#' @param centersN for plot type \code{"ellipses"} a matrix specifying the positions of
#'                 the centers of the covariance estimation centers, see also \code{\link[ssMRCD]{geo_weights}}.
#' @param colour_scheme coloring scheme used for plot type \code{"ellipses"}, either \code{"trace"} or \code{"regularity"} or \code{"none"}.
#' @param xlim_upper numeric giving the upper x limit for plot type \code{"convergence"}.
#' @param manual_rescale for plot type \code{"ellipses"} numeric used to re-scale ellipse sizes.
#' @param legend logical, if color legend should be included.
#' @param xlim vector of xlim (see \code{\link{par}}).
#' @param ylim vector of ylim (see \code{\link{par}}).
#' @param ... further plotting parameters.
#'
#' @details For \code{type = "convergence"} a plot is produced displaying the convergence behaviour.
#' Each line represents a different initial value used for the c-step iteration. On the x-axis the
#' iteration step is plotted with the corresponding value of the objective function. Not monotonically
#' lines are plotted in red. \cr
#'
#' For \code{type = "ellipses"} and more than a 2-dimensional data setting plotting the exact tolerance ellipse is
#' not possible anymore. Instead the two eigenvectors with highest eigenvalue from the
#' MCD used on the full data set without neighborhood assignments are taken and used as axis for
#' the tolerance ellipses of the ssMRCD covariance estimators. The tolerance ellipse for the global MCD
#' covariance is plotted in grey in the upper left corner. It is possible to set the colour scheme
#' to \code{"trace"} to see the overall amount of variabilty and compare the plotted covariance and
#' the real trace to see how much variance is not plotted. For \code{"regularity"} the regularization of each
#' covariance is shown.
#'
#' @return Returns plots of the ssMRCD methodology and results.
#'
#' @examples
#' # set seed
#' set.seed(1)
#'
#' # create data set
#' data = matrix(rnorm(2000), ncol = 4)
#' coords = matrix(rnorm(1000), ncol = 2)
#' groups = sample(1:10, 500, replace = TRUE)
#' lambda = 0.3
#'
#' # calculate ssMRCD by using the local outlier detection method
#' outs = local_outliers_ssMRCD(data = data,
#'                              coords = coords,
#'                              groups = groups,
#'                              lambda = lambda,
#'                              k = 10)
#'
#' # plot ssMRCD object included in outs
#' plot(x = outs$ssMRCD,
#'      centersN = outs$centersN,
#'      colour_scheme = "trace",
#'      legend = FALSE)
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}, \link[ssMRCD]{summary.ssMRCD},
#' \link[ssMRCD]{local_outliers_ssMRCD}, \link[ssMRCD]{plot.locOuts}}
#'
#' @exportS3Method plot ssMRCD
#' @importFrom robustbase covMcd
#' @importFrom  car ellipse
#' @importFrom scales alpha
plot.ssMRCD = function(x,
                       type = c("convergence", "ellipses"),
                       centersN = NULL,
                       colour_scheme = "none",
                       xlim_upper = 9,
                       manual_rescale = 1,
                       legend = TRUE,
                       xlim = NULL,
                       ylim = NULL,
                       ...){


  if("convergence" %in% type){
    color_incr = "darkred"
    ymax = NA

    p = dim(x$mT)[1]
    lambda = x$lambda

    # prepare data
    tmp = x$obj_fun_values
    max_steps = length(tmp[1,])

    # get colours
    decr = which(apply(X = tmp, MARGIN = 1, FUN = monotonic_decreasing))
    col = rep("grey", length(tmp[,1]))
    col[!decr] = color_incr

    alpha_val = rep(0.3, length(tmp[,1]))
    alpha_val[!decr] = 0.6

    graphics::matplot(t(tmp),
            type ="l",
            xlim = c(1, xlim_upper),
            col = scales::alpha(col, alpha_val), lwd  = 3, lty = 1,
            xlab = "Iterations",
            ylab = "Objective function value",
            main = "Convergence",
            ...)
    graphics::legend("topright",
           col = c("grey", color_incr),
           c("decreasing", "non-decreasing"),
           lty = c(1,1),
           lwd = c(3,3) )
  }

  if("ellipses" %in% type){

    if(is.null(centersN)) {
      stop("You need to specify the centers of the neighborhoods (centersN) to plot covariance ellipses.")}
    if(is.null(dim(centersN))){
      centersN = cbind(centersN, 1) # one dimensional space
    }

    N = length(x$MRCDcov)
    p = dim(x$mT)[1]
    lambda = x$lambda

    # base plot
    diffx = max(centersN[,1]) - min(centersN[,1])
    diffy = max(centersN[,2]) - min(centersN[,2])
    scale_lim = 0.3

    if(is.null(xlim)){
      xlim = c(min(centersN[,1]) - diffx * scale_lim, max(centersN[,1]) + diffx * scale_lim)
    }
    if(is.null(ylim)){
      ylim = c(min(centersN[,2]) - diffy * scale_lim, max(centersN[,2]) + diffy * scale_lim)
    }

    plot(min(centersN[,1]), min(centersN[,2]),
         type = "l",
         #asp = 1,
         xlim = xlim,
         ylim = ylim,
         main = bquote(Covariances~(lambda == .(lambda) )),
         xlab = "Coordinate 1",
         ylab = "Coordinate 2",
         ...)

    graphics::text(centersN[,1], centersN[,2], 1:N,
                   cex = 1, col = "black")

    # colorizing
    legend_labels = rep(NA, N)
    color_scale = rep(NA, N)
    if(colour_scheme == "none"){
      min = 0
      max = 0
      color_scale = rep("black", N)
    }
    if(colour_scheme == "regularity"){
      min = min(x$rho)
      max = max(x$rho)
      if(min == max) {
        color_scale = rep("black", N)
      } else {
        color_scale = cut(x$rho,
                          breaks = seq(0, max, length.out = 6),
                          labels = grDevices::heat.colors(8)[5:1])
        legend_labels = cut(x$rho,
                           breaks = seq(0, max, length.out = 6))
      }
    }
    if(colour_scheme == "trace"){
      tr = rep(NA, N)
      for(i in 1:N){
        tr[i] = sum(diag(x$MRCDcov[[i]]))
      }
      min = min(tr)
      max = max(tr)
      if(min == max) {
        color_scale = rep("black", N)
      } else {
        color_scale = cut(tr,
                          breaks = round(seq(min-0.01, max + 0.01, length.out = 6), 2),
                          labels = grDevices::heat.colors(8)[5:1])
        legend_labels = cut(tr,
                            breaks = round(seq(min-0.01, max+0.01, length.out = 6), 2))
      }
    }


    # MCD ellipse
    MCD = robustbase::covMcd(do.call(rbind,x$mX),
                             nsamp = "deterministic")$cov
    S = eigen(MCD)$vectors
    D = diag(eigen(MCD)$values)
    Dsqrti = diag(sqrt(diag(D)^(-1)))
    S = S %*% Dsqrti
    scale_ell = (max(centersN[,1]) - min(centersN[,1]) + 2*diffx*scale_lim)*0.25/(sqrt(N))
    ell = car::ellipse(c(min(centersN[,1]), max(centersN[,2])),
                       shape=diag(1,2),
                       radius=scale_ell*manual_rescale,
                       segments=1e3,
                       draw = FALSE)
    traceMCD = round(sum(diag(MCD)), 2)
    graphics::lines(ell[,"x"], ell[,"y"], col = "lightgray", lty= "dashed")
    graphics::text(min(centersN[,1]), max(centersN[,2]), "glob", col="lightgray", cex = 0.9)

    # neighborhood ellipses
    for(i in 1:N){
      eigtraf = t(S) %*% x$MRCDcov[[i]] %*% S
      ell_tmp = car::ellipse(centersN[i,],
                             shape=eigtraf[1:2, 1:2],
                             radius=scale_ell*manual_rescale,
                             segments=1e3,
                             draw = FALSE)
      graphics::lines(ell_tmp[, "x"], ell_tmp[, "y"],
                      col = as.character(color_scale[i]),
                      lwd = 1.5)
    }

    # legend
    if (legend & colour_scheme == "regularity" & min != max){
      graphics::legend("topright",
                       col = c(grDevices::heat.colors(8)[5:1]),
                       legend = c(levels(legend_labels)),
                       lwd = 1.5,
                       cex = 0.8,
                       title = "Regularity",
                       lty = c(rep("solid", 5)))
    }
    if (legend & colour_scheme == "trace" & min != max){
      graphics::legend("topright",
                       col = c(grDevices::heat.colors(8)[5:1], "grey"),
                       legend = c(levels(legend_labels), paste(traceMCD, "(global)")),
                       lwd = 1.5,
                       cex = 0.8,
                       title = "Trace",
                       lty = c(rep("solid", 5)))
    }

  }
}


##########################################################################################
#' Scale Data Locally
#'
#' @param ssMRCD \code{ssMRCD} object, see \code{\link[ssMRCD]{ssMRCD}}
#' @param X matrix, new data to scale with ssMRCD estimation.
#' @param groups vector, group assignments of new data \code{X}.
#' @param multivariate logical, \code{TRUE} if multivariate structure should be used.
#' Otherwise, univariate variances from the ssMRCD estimator is used.
#' @param center_only logical, if \code{TRUE} observations are only centered.
#'
#' @return Returns matrix of observations. If \code{X = NULL} X from the ssMRCD object is
#' used and sorted according to group numbering.
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}}
#'
#' @export
#' @importFrom expm sqrtm
#'
#' @examples
#'# create data set
#' x1 = matrix(runif(200), ncol = 2)
#' x2 = matrix(rnorm(200), ncol = 2)
#' x = list(x1, x2)
#'
#' # create weighting matrix
#' W = matrix(c(0, 1, 1, 0), ncol = 2)
#'
#' # calculate ssMRCD
#' localCovs = ssMRCD(x, weights = W, lambda = 0.5)
#'
#' # scale used data
#' scale_ssMRCD(localCovs,
#'       multivariate = TRUE)
#'
#' # scale new data
#' scale_ssMRCD(localCovs,
#'       X = matrix(rnorm(20), ncol = 2, nrow = 10),
#'       groups = rep(2, 10),
#'       multivariate =TRUE)

scale_ssMRCD = function(ssMRCD,
                        X = NULL,
                        groups = NULL,
                        multivariate = FALSE,
                        center_only = FALSE){

  if(is.null(X) | is.null(groups)) {
    X = do.call(rbind, ssMRCD$mX)
    groups = rep(1:length(ssMRCD$MRCDcov), times = sapply(X = ssMRCD$mX,
                                                          FUN = function(x) dim(x)[1]))
  }
  X = as.matrix(X)
  N = ssMRCD$N

  if(!multivariate){
    for(i in 1:N){
      ind = which(groups == i)
      if(length(ind)!= 0){
        if(!center_only){
          X[groups == i,] = scale(x = X[groups == i,],
                                  center = ssMRCD$MRCDmu[[i]],
                                  scale = sqrt(diag(ssMRCD$MRCDcov[[i]])))
        } else {
          X[groups == i,] = scale(x = X[groups == i,],
                                  center = ssMRCD$MRCDmu[[i]])
        }
      }
    }
  }

  if(multivariate){
    for(i in 1:N){
      ind = which(groups == i)
      if(length(ind)!= 0){
        centered = t(X[groups == i,]) - matrix(ssMRCD$MRCDmu[[i]],
                                               ncol = sum(groups == i),
                                               nrow = dim(X)[2])
        if(!center_only){
          X[groups == i,] = t(expm::sqrtm(ssMRCD$MRCDicov[[i]]) %*% centered)
        } else{
          X[groups == i,]  = centered
        }
      }
    }
  }

  return(X)
}


##########################################################################################
outliers = function(ssMRCD){
  # get indices of non local outliers in data

  hsets = ssMRCD$hset
  size_ns = sapply(ssMRCD$mX, function(x) dim(x)[1])
  N = ssMRCD$N

  ind = c()
  size_sum = 0
  for(i in 1:N){
    ind = c(ind, hsets[[i]] + size_sum)
    size_sum = size_sum + size_ns[i]
  }

  return(ind)   # returns index vector for sorted groups!
}

##########################################################################################
#' Extracting Residuals from Local Fit
#'
#' @param object \code{ssMRCD} object, see \code{\link[ssMRCD]{ssMRCD}}.
#' @param ... see details
#'
#' @return Returns either all residuals or the mean of the residual norms lower than the \code{alpha}- Quantile.
#'
#' @details Other input variables are: \tabular{ll}{
#'    \code{remove_outliers} \tab logical (default \code{FALSE}). If TRUE, only residuals
#'    from not outlying observations are calculated. If FALSE, trimmed residuals are used (see \code{alpha}). \cr
#'    \tab \cr
#'    \code{X} \tab matrix of new data, if data from the \code{ssMRCD} object is used. \cr
#'    \tab \cr
#'    \code{groups} \tab vector of groups for new data, if \code{NULL} data from the \code{ssMRCD} object is used. \cr
#'    \tab \cr
#'    \code{mean} \tab logical (default \code{FALSE}), specifying if mean of trimmed
#'    observations is returned or all residuals. \cr
#' }
#'
#
#' If \code{X} and \code{groups} are provided, \code{alpha} is set to one and all residuals are used.
#' If \code{remove_outliers} is TRUE, \code{alpha} is set to 1 automatically.
#'
#'
#' @exportS3Method residuals ssMRCD
#'
#' @examples
#'# create data set
#' x1 = matrix(runif(200), ncol = 2)
#' x2 = matrix(rnorm(200), ncol = 2)
#' x = list(x1, x2)
#'
#' # create weighting matrix
#' W = matrix(c(0, 1, 1, 0), ncol = 2)
#'
#' # calculate ssMRCD
#' localCovs = ssMRCD(x, weights = W, lambda = 0.5)
#'
#' # residuals of model
#' residuals(localCovs, remove_outliers = TRUE, mean = FALSE)
#'
#' # residuals of new data
#' residuals(localCovs,
#'       X = matrix(rnorm(20), ncol = 2, nrow = 10),
#'       groups = rep(2, 10),
#'       mean =TRUE)
#'
residuals.ssMRCD = function(object, ...){

  args = list(...)
  if(is.null(args$remove_outliers)) args$remove_outliers = FALSE
  if(is.null(args$mean)) args$mean = FALSE

  N = length(object$MRCDcov)
  if(is.null(args$X) | is.null(args$groups)){
    args$X = do.call(rbind, object$mX)
    args$groups = rep(1:N, times = sapply(X = object$mX, FUN = function(x) dim(x)[1]))
  }
  n = dim(args$X)[1]

  ind = 1:n
  if(args$remove_outliers) ind = outliers(ssMRCD = object)  # only sensible for X, groups not new data
  if(!args$remove_outliers)  alpha = object$alpha
  if(args$remove_outliers)  alpha = 1

  # calculate residuals
  residuals = scale_ssMRCD(object, X = NULL, groups = NULL, multivariate = TRUE)
  if(!args$mean) return(residuals)

  # calculate mean of norm
  res_norm = sqrt(diag(residuals[ind, ] %*% t(residuals[ind, ])))
  ind = sort.int(res_norm, index.return = T)$ix[1:round(n*alpha)]
  res_trimmed = res_norm[ind]

  return(mean(res_trimmed))
}


