#' Local Outlier Detection using Spatially Smoothed MRCD
#'
#' Identifies local multivariate outliers using spatially smoothed robust covariance estimation (ssMRCD) as proposed by Puchhammer and Filzmoser (2023). For each observation, the Mahalanobis distance to its nearest neighbor is computed, and an adjusted boxplot is used to detect outliers.
#'
#' @param data A numeric matrix of observations (rows = observations, columns = variables).
#' @param coords A numeric matrix of spatial coordinates corresponding to the observations.
#' @param groups A vector assigning each observation to a neighborhood/group.
#' @param lambda Smoothing parameter for the \code{\link[ssMRCD]{ssMRCD}} estimator.
#' @param weights Optional weighting matrix for spatial smoothing. If omitted, inverse-distance weights are computed automatically.
#' @param k Integer. Number of nearest neighbors to use if \code{dist} is not provided.
#' @param dist Numeric. Use neighbors within this distance instead of \code{k}-nearest neighbors. If both are provided, \code{dist} is used.
#'
#' @return An object of class \code{"locOuts"} containing:
#' \describe{
#'   \item{\code{outliers}}{Indices of detected outliers.}
#'   \item{\code{next_distance}}{Vector of Mahalanobis next distances (min distance to neighbors).}
#'   \item{\code{cutoff}}{Upper fence of the adjusted boxplot used as outlier threshold.}
#'   \item{\code{coords}}{Matrix of observation coordinates.}
#'   \item{\code{data}}{Original data matrix.}
#'   \item{\code{groups}}{Group assignments.}
#'   \item{\code{k, dist}}{Neighborhood comparison parameters used.}
#'   \item{\code{centersN}}{Centers of neighborhoods.}
#'   \item{\code{matneighbor}}{Binary matrix indicating which neighbors were used for each observation.}
#'   \item{\code{ssMRCD}}{The fitted \code{ssMRCD} object.}
#' }
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}}, \code{\link[ssMRCD]{plot.locOuts}}
#'
#' @references
#' Puchhammer, P. and Filzmoser, P. (2023). Spatially Smoothed Robust Covariance Estimation for Local Outlier Detection. *Journal of Computational and Graphical Statistics*, 33(3), 928–940. \doi{10.1080/10618600.2023.2277875}
#'
#' @export
#' @importFrom dbscan kNN
#' @importFrom robustbase adjbox
#'
#' @examples
#' data <- matrix(rnorm(2000), ncol = 4)
#' coords <- matrix(runif(1000), ncol = 2)
#' groups <- sample(1:10, 500, replace = TRUE)
#' result <- locOuts(data, coords, groups, lambda = 0.3, k = 10)
locOuts <- function(data, coords, groups, lambda, weights = NULL, k = NULL, dist = NULL) {
  # --- Input Checks ---
  data <- as.matrix(data)
  coords <- as.matrix(coords)
  groups <- as.vector(groups)

  n <- nrow(data)
  p <- ncol(data)
  N <- length(unique(groups))

  # --- Spatial Weights ---
  if (is.null(weights)) {
    message("Using default inverse-distance weights.")
    geo <- geo_weights(coordinates = coords, groups = groups)
    weights <- geo$W
    centersN <- geo$centersN
  } else {
    geo <- geo_weights(coordinates = coords, groups = groups)
    centersN <- geo$centersN
  }

  # --- Neighborhood Matrix ---
  distance_method <- "k"
  if (is.null(k) && is.null(dist)) {
    k <- 10
    message("No 'k' or 'dist' given. Defaulting to k = 10.")
  } else if (!is.null(dist)) {
    distance_method <- "dist"
    k <- NULL
    if (!is.null(k)) {
      message("Both 'k' and 'dist' are given. Using 'dist'.")
    }
  }

  matneighbor <- matrix(0, nrow = n, ncol = n)
  knn_all <- dbscan::kNN(coords, k = n - 1)

  for (i in 1:n) {
    if (distance_method == "k") {
      neighbor_ids <- dbscan::kNN(coords, k = k)$id[i, ]
    } else {
      dists <- knn_all$dist[i, ]
      valid_idx <- which(dists < dist)
      neighbor_ids <- knn_all$id[i, valid_idx]
    }
    matneighbor[i, neighbor_ids] <- 1
  }

  # --- Robust Covariance Estimation (ssMRCD) ---
  cov_object <- ssMRCD(X = data, groups = groups, weights = weights, lambda = lambda)

  center_list <- cov_object$MRCDmu
  icov_list <- cov_object$MRCDicov

  # --- Compute Next Distances ---
  next_distance <- numeric(n)
  for (i in 1:n) {
    group_i <- groups[i]
    xi <- data[i, ]
    icov <- icov_list[[group_i]]
    neighbors <- which(matneighbor[i, ] == 1)
    if (length(neighbors) == 0) {
      next_distance[i] <- NA
      next
    }
    distances <- sapply(neighbors, function(j) {
      xj <- data[j, ]
      diff <- xi - xj
      t(diff) %*% icov %*% diff
    })
    next_distance[i] <- min(distances)
  }

  # --- Outlier Detection via Adjusted Boxplot ---
  adjbox_res <- robustbase::adjbox(next_distance, plot = FALSE, horizontal = TRUE, range = 1.5)
  cutoff <- adjbox_res$fence[2]
  outliers <- which(next_distance > cutoff)

  # --- Output ---
  structure(list(
    outliers = outliers,
    next_distance = next_distance,
    cutoff = cutoff,
    coords = coords,
    data = data,
    groups = groups,
    k = k,
    dist = dist,
    centersN = centersN,
    matneighbor = matneighbor,
    ssMRCD = cov_object
  ), class = c("locOuts", "list"))
}




# PLOTTING FUNCTIONS OUTLIER DETECTION AND METHODS FOR locOuts
#' Diagnostic Plots for Local Outlier Detection (`locOuts`)
#'
#' Produces diagnostic plots for local outlier detection results
#' returned by \code{\link[ssMRCD]{locOuts}}.
#' Available visualizations include a histogram of next distances, spatial distribution of next distances,
#' and a parallel coordinate plot (PCP) for a selected observation and their neighborhood.
#'
#' @param x An object of class \code{"locOuts"} obtained from \code{\link[ssMRCD]{locOuts}}.
#' @param type A character vector indicating which plots to generate. Options are:
#'   \describe{
#'     \item{"hist"}{Histogram of next distances with cutoff visualized.}
#'     \item{"spatial"}{Spatial distribution of observations, colored by relative next distance.}
#'     \item{"pcp"}{Parallel coordinate plot for an observation and its neighbors.}
#'   }
#' @param scale Character indicating how variables are scaled in the parallel coordinate plot. One of:
#'   \describe{
#'     \item{"none"}{Use raw values (no scaling).}
#'     \item{"minmax"}{Min-max scaling to [0, 1].}
#'     \item{"zscore"}{Standardization: mean 0, standard deviation 1.}
#'   }
#' @param bins Integer, number of histogram bins (default = 30).
#' @param observation Integer or character; index or name of a specific observation to analyze in the PCP plot.
#'   Used only when \code{type} includes \code{"pcp"}.
#' @param ... Additional parameters passed to low-level plotting functions (currently unused in ggplot versions).
#'
#' @details
#' The function visualizes outlier behavior in different ways:
#' \itemize{
#'   \item \strong{Histogram}: Shows the distribution of next distances across observations. The cutoff is shown as a dashed line.
#'   \item \strong{Spatial Plot}: 2D plot of observation coordinates. Color encodes the ratio of next distance to cutoff.
#'   \item \strong{Parallel Coordinate Plot (PCP)}: Shows scaled values across all variables for a selected observation
#'     (in red) and its neighbors (in blue or grey). The type of scaling can be controlled via the \code{scale} parameter.
#' }
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{p_hist}}{ggplot object of the histogram (or \code{NULL} if not requested).}
#'     \item{\code{p_spatial}}{ggplot object of the spatial plot (or \code{NULL}).}
#'     \item{\code{p_pcp}}{ggplot object of the parallel coordinate plot (or \code{NULL}).}
#'   }
#'
#' @seealso \code{\link[ssMRCD]{locOuts}}
#'
#' @import ggplot2
#' @importFrom scales muted
#' @importFrom graphics hist
#'
#' @exportS3Method plot locOuts
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(2000), ncol = 4)
#' coords <- matrix(rnorm(1000), ncol = 2)
#' groups <- sample(1:10, 500, replace = TRUE)
#' outs <- locOuts(data = data,
#'                 coords = coords,
#'                 groups = groups,
#'                 lambda = 0.3,
#'                 k = 10)
#'
#' # Generate all plots
#' plots <- plot(outs,
#'               type = c("hist", "spatial", "pcp"),
#'               observation = outs$outliers[1],
#'               scale = "minmax")
#' plots$p_hist
#' plots$p_spatial
#' plots$p_pcp

plot.locOuts = function(x,
                        type = c("hist", "spatial", "pcp"),
                        scale = c("none", "minmax", "zscore"),
                        bins = 30,
                        observation = 1,
                        ...){

  p_hist = NULL
  if ("hist" %in% type) {

    # Create histogram using ggplot2
    p_hist <- ggplot2::ggplot() +
      ggplot2::geom_histogram(aes(x =  x$next_distance),
                              color = "black",
                              fill = "lightblue",
                              bins = bins) +
      ggplot2::geom_vline(xintercept = x$cutoff,
                 color = "black",
                 linetype = "dashed",
                 linewidth = 1) +
      ggplot2::annotate("label",
               x = x$cutoff,
               y = max(graphics::hist(x$next_distance, plot = FALSE, breaks = bins)$counts),
               label = "cut-off",
               color = "black") +
      ggplot2::labs(x = "Next distances",
           y = "Count",
           title = "Histogram of Next Distances") +
      ggplot2::theme_minimal()
  }

  # spatial plot next distances
  p_spatial = NULL
  if("spatial" %in% type){

    # Input check
    stopifnot(ncol(x$coords) == 2)

    # Create plotting data
    df_plot <- data.frame(
      x = x$coords[, 1],
      y = x$coords[, 2],
      group = as.factor(x$groups),
      color = x$next_distance/x$cutoff
    )

    # Plot with ggplot2
    p_spatial <- ggplot2::ggplot() +
      ggplot2::geom_point(aes(x = df_plot$x,
                     y = df_plot$y,
                     color = df_plot$color),
                 size = 2) +
      ggplot2::scale_color_gradient2(name = "Next Distance/Cut-Off",
                            midpoint = 1,
                            low = scales::muted("blue"),
                            high = scales::muted("red")) +
      ggplot2::labs(
        x = "Coordinate 1",
        y = "Coordinate 2",
        title = "Spatial Distribution: Next Distances"
      ) +
      ggplot2::theme_minimal()
  }


  p_pcp = NULL
  if("pcp" %in% type){

    if(is.null(observation)) stop("For local PCP please select one observation by name or index.")
    if(is.character(observation) & !is.null(x$ssMRCD$rnames) ){
      if(! (observation %in% x$ssMRCD$rnames)) stop(paste("There is no observation with name", observation,
                                                 "in the rownames given by x$rnames."))
    }
    observation_i = observation # if numeric
    if(is.character(observation)) {
      observation_i = which(observation == x$ssMRCD$rnames)
    }

    cnames = x$ssMRCD$cnames
    if(is.null(cnames)) cnames = 1:ncol(x$ssMRCD$MRCDcov[[1]])

    ind_neighbors = which(x$matneighbor[observation_i, ] == 1)
    p = length(cnames)

    # Scaling
    df = x$data
    scale_data <- switch(scale,
                         "none" = df,
                         "minmax" = as.data.frame(apply(df, 2, function(x) (x - min(x)) / (max(x) - min(x)))),
                         "zscore" = as.data.frame(scale(df))
    )
    scale_data = scale_data[c(observation_i, ind_neighbors), ]

    # Convert to long format manually (no tidyr)
    id <- rep(1:length(c(observation_i, ind_neighbors)), each = p)
    variable <- rep(cnames, times = length(c(observation_i, ind_neighbors)))
    value <- as.vector(t(scale_data))
    df_long <- data.frame(id = id,
                          variable = factor(variable, levels = cnames),
                          value = value,
                          group = rep(c("a", rep("b", length(ind_neighbors))), each = p) )

    p_pcp = ggplot2::ggplot() +
      ggplot2::geom_line(
                aes(x = df_long$variable,
                    y = df_long$value,
                    group = df_long$id,
                    color = as.factor(df_long$group)),
                alpha = 0.5,
                linewidth = 0.75) +
      ggplot2::scale_color_manual(values = c("a" = scales::muted("red"), "b"= scales::alpha("grey", 0.2)), guide = "none") +
      ggplot2::labs(title = paste("Parallel Coordinates Plot - Scaling:", scale),
           subtitle = paste("Observation", observation),
           x = "Variables", y = "Scaled Value") +
      ggplot2::theme_minimal()
    }

  return(list(p_hist = p_hist, p_spatial = p_spatial, p_pcp = p_pcp))
}





