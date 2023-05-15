# METHODS FOR cov_ssMRCD OBJECT



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
#' N_assignments = sample(1:10, 500, replace = TRUE)
#' lambda = 0.3
#'
#' # calculate ssMRCD by using the local outlier detection method
#' outs = local_outliers_ssMRCD(data = data,
#'                              coords = coords,
#'                              N_assignments = N_assignments,
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
plot.ssMRCD = function(x, type = c("convergence", "ellipses"),
                           centersN = NULL, colour_scheme = "none",
                           xlim_upper = 9, manual_rescale = 1,
                           legend = TRUE, xlim = NULL, ylim = NULL,  ...){

  # Check object
  check_ssMRCD(x)

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

    check_input(centersN, "matrix")
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

