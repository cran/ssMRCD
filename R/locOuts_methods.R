# PLOTTING FUNCTIONS OUTLIER DETECTION AND METHODS FOR res_local_outlier

#' Summary of Local Outlier Detection
#'
#' Prints a summary of the locOuts object obtained by the function \code{\link[ssMRCD]{local_outliers_ssMRCD}}.
#'
#' @param object a locOuts object.
#' @param ... further parameters passed on.
#'
#' @exportS3Method summary locOuts
#'
#' @seealso \code{\link[ssMRCD]{plot.locOuts}}
#'
#' @return Prints a summary of the \code{locOuts} object.
#'
#' @examples
#' # set seed
#' set.seed(1)
#'
#' # make locOuts object
#' data = matrix(rnorm(2000), ncol = 4)
#' coords = matrix(rnorm(1000), ncol = 2)
#' N_assignments = sample(1:10, 500, replace = TRUE)
#' lambda = 0.3
#'
#'# local outlier detection
#' outs = local_outliers_ssMRCD(data = data,
#'                              coords = coords,
#'                              N_assignments = N_assignments,
#'                              lambda = lambda,
#'                              k = 10)
#'
#'# summary method
#' summary(outs)
summary.locOuts = function(object, ...){

  n_out = length(object$outliers)
  n_all = length(object$next_distance)
  cut_off = object$cutoff
  k = object$k
  dist = object$dist

  cat("Local outliers: ", n_out, " of ", n_all, " observations (", round(n_out/n_all*100, 2),"%).\n", sep = "")
  cat("Used cut-off value: ", cut_off, ".\n", sep = "")
  if(!is.null(k)) cat("Observations are compared to their k=", k, " neighbors.\n", sep = "")
  if(!is.null(dist)) cat("Observations are compared to neighbors closer than dist = ", dist, ".\n", sep = "")

}


#' Diagnostic Plots for Local Outlier Detection
#'
#' This function plots different diagnostic plots for local outlier detection.
#' It can be applied to an object of class \code{"locOuts"} which is the output of the function \code{\link[ssMRCD]{local_outliers_ssMRCD}}.
#'
#'
#' @param x a locOuts object obtained by the function \code{\link[ssMRCD]{local_outliers_ssMRCD}}.
#' @param type vector containing the types of plots that should be plotted, possible values \code{c("hist", "spatial", "lines", "3D")}.
#' @param colour character specifying the color scheme (see details). Possible values \code{"all", "onlyOuts", "outScore"}.
#' @param focus an integer being the index of the observation whose neighborhood should be analysed more closely.
#' @param pos integer specifying the position of the text "cut-off" in the histogram (see \code{\link{par}}).
#' @param alpha scalar specifying the transparancy level of the points plotted for plot type \code{"spatial", "3D"} and \code{"lines"}.
#' @param data optional data frame or matrix used for plot of type \code{"line"}. Will be used to plot lines based scaled \code{data} instead of the data used for local outlier detection.
#' @param add_map TRUE if a map should be plotted along the line plot (\code{type = "lines"}).
#' @param ... further parameters passed on to base-R plotting functions.
#'
#' @details Regarding the parameter \code{type} the value \code{"hist"} corresponds to a plot of the
#' histogram of the next distances together with the used cutoff-value.
#' When using \code{"spatial"} the coordinates of each observation are plotted and colorized according to the color setting.
#' The \code{"lines"} plot is used with the index \code{focus} of one observation whose out/inlyingness to its neighborhood
#' should by plotted. The whole data set is scaled to the range [0,1] and the scaled value of the selected observation and
#' its neighbors are plotted. Outliers are plotted in orange.
#' The \code{"3D"} setting leads to a 3D-plot using the colour setting as height.
#' The view can be adapted using the parameters \code{theta} and \code{phi}. \cr \cr
#' For the \code{colour} setting possible values are \code{"all"} (all next distances are
#' used and colored in an orange palette), \code{"onlyOuts"} (only outliers are
#' plotted in orange, inliers are plotted in grey) and \code{"outScore"} (the next
#' distance divided by the cutoff value is used to colourize the points; inliers are colorized in blue, outliers in orange).
#'
#' @return Returns plots regarding next distances and spatial context.
#'
#' @seealso \code{\link[ssMRCD]{local_outliers_ssMRCD}}
#'
#' @importFrom plot3D scatter3D
#' @importFrom graphics plot
#' @exportS3Method plot locOuts
#'
#' @examples
#' # set seed
#' set.seed(1)
#'
#' # make locOuts object
#' data = matrix(rnorm(2000), ncol = 4)
#' coords = matrix(rnorm(1000), ncol = 2)
#' N_assignments = sample(1:10, 500, replace = TRUE)
#' lambda = 0.3
#'
#'# local outlier detection
#' outs = local_outliers_ssMRCD(data = data,
#'                              coords = coords,
#'                              N_assignments = N_assignments,
#'                              lambda = lambda,
#'                              k = 10)
#'
#' # plot results
#' plot(outs, type = "hist")
#' plot(outs, type = "spatial", colour = "outScore")
#' plot(outs, type = "3D", colour = "outScore", theta = 0)
#' plot(outs, type ="lines", focus = outs$outliers[1])
plot.locOuts = function(x,
                        type = c("hist", "spatial", "lines", "3D"),
                        colour = "all",
                        focus = NULL,
                        pos = NULL,
                        alpha = 0.3,
                        data = NULL,
                        add_map = TRUE,
                        ...){

  # histogram
  if("hist" %in% type){
    if(missing(pos)) pos = 2

    h = graphics::hist(x$next_distance,
                       xlab = "Next distances",
                       main = "Histogram of Next Distances", ...)
    graphics::abline(v = x$cutoff, col = "darkorange", lty = 2, lwd = 2)
    graphics::text(x$cutoff, max(h$counts), "cut-off", col = "darkorange", pos = pos)
  }

  # color scheme
  if(colour == "all"){
    breaks = round(seq(0, max(x$next_distance), length.out = 6), 2)
    z = x$next_distance
    color_scale = cut(x$next_distance,
                      breaks = breaks,
                      labels = grDevices::heat.colors(8)[5:1])
    legend_text = cut(x$next_distance,
                      breaks =breaks)
  }
  if(colour == "onlyOuts"){
    breaks = round(c(0, seq(x$cutoff, max(x$next_distance), length.out = 6)), 2)
    z = x$next_distance
    color_scale = cut(x$next_distance,
                      breaks = breaks,
                      labels = c("grey", grDevices::heat.colors(8)[5:1]))
    legend_text = cut(x$next_distance,
                      breaks = breaks)
  }
  if(colour == "outScore"){
    breaks = round(unique(c(seq(0, 1, length.out = 4), seq(1, max(x$next_distance/x$cutoff), length.out = 4))),2)
    z = x$next_distance/x$cutoff
    color_scale = cut(x$next_distance/x$cutoff,
                      breaks = breaks,
                      labels = c( grDevices::hcl.colors(3, palette = "Teal"), grDevices::heat.colors(6)[c(4, 3, 2)]))
    legend_text = cut(x$next_distance/x$cutoff,
                      breaks = breaks)
  }


  # spatial plot next distances
  if("spatial" %in% type){

    # input check
    stopifnot(dim(x$coords)[2] == 2)

    graphics::plot(x = x$coords[, 1],
                   y = x$coords[, 2] ,
                   col = scales::alpha(as.character(color_scale), alpha),
                   type = "p",
                   pch = 43,
                   xlab = "Coordinate 1",
                   ylab = "Coordinate 2",
                   main = "Spatial Distribution: Next distances",
                   ...)
    graphics::text(x = x$centersN[,1], y = x$centersN[, 2], labels = 1:x$ssMRCD$N, col = scales::alpha("black", 0.5))
    graphics::legend("topright", legend = levels(legend_text), col = as.character(levels(color_scale)), lwd = 2, title = "Next distance", cex = 0.8)
  }


  if("3D" %in% type){
    plot3D::scatter3D(x = x$coords[, 1],
              y = x$coords[, 2],
              z = z,
              col = alpha(as.character(levels(color_scale)), 1),
              breaks = breaks,
              type = "h",
              ...)
  }

  if("lines" %in% type){

    if(is.null(data)){
      data = as.matrix(x$data)
    }
    columns = colnames(data)
    if(is.null( columns)) columns = 1:dim(data)[2]
    p =length(columns)
    n = dim(data)[1]
    data_scaled = matrix(NA, n, p)

    # scale for plotting
    for(i in 1:p){
      min = min(data[, i])
      max = max(data[, i])
      scale_factor = 1/(max - min)
      data_scaled[,i] = (data[,i] - min)*scale_factor
    }


    if(is.null(focus)) {

      # empty plot
      xvar = 1:p
      graphics::plot(xvar, rep(NA, p), ylim = c(-0.1,1.1), xaxt = "n", yaxt = "n", ylab = "", xlab = "", las = 2)
      graphics::axis(1, at=1:p, labels=colnames(data), las = 1)
      graphics::axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), las = 1)
      graphics::title("Multivariate Values")

      # fill with lines
      for( i in 1:n){
        if (i %in% x$outliers){
          color = "darkorange"
        } else{
          color = alpha("grey", 0.3 )
        }
        graphics::lines(xvar, data_scaled[i, ], col = color, lwd = 2)
      }
    }
    if(!is.null(focus)){
      stopifnot(length(focus)== 1)

      # empty plot
      xvar = 1:p
      graphics::plot(xvar, rep(NA, p), ylim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "", las = 2, ...)
      graphics::axis(1, at=1:p, labels=colnames(data), las = 2)
      graphics::axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), las = 1)
      graphics::title(paste("Multivariate Values (i=", focus, ")"))
      subset = which(x$matneighbor[focus, ]== 1)

      # fill with lines
      for( i in c(subset, focus) ){
        if (i %in% x$outliers){
          color = "darkorange"
        }
        if(i == focus){
          color = "darkred"
          graphics::text(p, data_scaled[i, p], focus, col = color, pos = 3)
        }
        if(i != focus & ! i %in% x$outliers) {
          color = alpha("grey", 0.3 )
        }
        graphics::lines(xvar, data_scaled[i, ], col = color, lwd = 2)
      }

      if(add_map){
        # second plot: spatial context
        colour = rep("grey", n)
        colour[subset] = "black"
        colour[(1:n) %in% x$outliers & (1:n) %in% subset] = "orange"
        colour[focus] = "darkred"
        graphics::plot(x = x$coords[, 1],
                       y = x$coords[, 2] ,
                       col = scales::alpha(colour, alpha),
                       type = "p",
                       pch = 43,
                       xlab = "Coordinate 1",
                       ylab = "Coordinate 2",
                       main = "Spatial Context",
                       ...)
      }
    }
  }
}





