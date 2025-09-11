##########################################################################################
#### PLOTTING METHOD                                                                  ####
##########################################################################################

#' Plot Method for msPCA Objects
#'
#' Generic plotting interface for objects of class \code{"msPCA"}. Depending on the \code{type} argument,
#' this function visualizes loadings, scores, screeplots, biplots, or score distances.
#'
#' @param x An object of class \code{"msPCA"}.
#' @param type Type of plot to produce. One of: \code{"loadings"} (default), \code{"screeplot"}, \code{"scores"},
#' \code{"biplot"}, or \code{"score_distances"}.
#' @param ... Additional arguments passed to internal functions.
#'
#' @details
#' If \code{type = "loadings"} a heatmap of loadings across groups is displayed. Optional input arguments are:
#' \tabular{ll}{
#'    \code{k} \tab Integer. The k-th principal component will be plotted.  \cr
#'    \code{text} \tab Boolean, whether the loading values should also be displayed as text. Default \code{FALSE}. \cr
#'    \code{size} \tab Numeric. Text size when \code{text} is \code{TRUE}. \cr
#'    \code{gnames} \tab Character vector. Names for groups (shown in plots). \cr
#'    \code{cnames} \tab Character vector. Names for variables (shown in plots). \cr
#'    \code{textrotate} \tab Rotation angle of text in the heatmap.\cr
#'    \code{tolerance} \tab Numeric, default 1e-04. Specifies the band of white color values around zero. \cr
#' }
#'
#' If \code{type = "screeplot"} boxplots of the explained variance per component and cumulative variance per group are plotted. Optional input arguments are:
#' \tabular{ll}{
#'    \code{text} \tab Boolean, whether the loading values should also be displayed as text. Default \code{TRUE}. \cr
#'    \code{size} \tab Numeric. Text size when \code{text} is \code{TRUE}. \cr
#'    \code{gnames} \tab Character vector. Names for groups (shown in plots). \cr
#'    \code{cutoff} \tab Scalar with default value 0.8. The cumulative percentage cutoff value for the overall explained variance. \cr
#'    \code{textrotate} \tab Rotation angle of text in the heatmap.\cr
#' }
#'
#' If \code{type = "biplot"} the loadings are visualized in the first and second component over all groups. Optional input arguments are:
#' \tabular{ll}{
#'    \code{color} \tab Character. Either \code{"variable"} (default) when the color should be connected to the variables or \code{"groups"} if the color should correspond to the groups. \cr
#'    \code{size} \tab Numeric. Text size when \code{text} is \code{TRUE}. \cr
#'    \code{alpha} \tab Alpha value for the loading points, default is \code{0.7} .\cr
#' }
#'
#' If \code{type = "scores"} a histogram of the k-th scores per group are shown. Optional input arguments are:
#' \tabular{ll}{
#'    \code{k} \tab Integer. The k-th principal component will be plotted. Default value is one.  \cr
#'    \code{ssMRCD} An object of class \code{ssMRCD} including list elements \code{MRCDcov, MRCDmu, mX}.  Alternatively, \code{X, groups, mu, Sigma} can be provided..
#'    \code{X, groups, mu, Sigma} The original data matrix, a vector of groups as for \code{ssMRCD}, a list of location vectors and a list of covariance matrices. \cr
#'   Alternatively, a \code{ssMRCD} object can be given. \cr
#'    \code{alpha} \tab Alpha value for histogramm, default is \code{0.7} .\cr
#' }
#'
#' If \code{type = "score_distances"} a distance-distance plot of score and orthogonal distances is shown for each group. of the k-th scores per group are shown. Optional input arguments are:
#' \tabular{ll}{
#'    \code{k} \tab Integer. Using the first k PCs for the distances. Default is the number of provided PCs.
#'    \code{ssMRCD} An object of class \code{ssMRCD} including list elements \code{MRCDcov, MRCDmu, mX}.  Alternatively, \code{X, groups, mu, Sigma} can be provided..
#'    \code{X, groups, mu, Sigma} The original data matrix, a vector of groups as for \code{ssMRCD}, a list of location vectors and a list of covariance matrices. \cr
#'   Alternatively, a \code{ssMRCD} object can be given. \cr
#'    \code{shape} \tab Point shape. \cr
#'    \code{size} \tab Numeric. Point size. \cr
#'    \code{alpha} \tab Alpha value for the points, default is \code{0.7} .\cr
#' }
#'
#' @return A ggplot2 object or list of plots, depending on the type.
#'
#' @seealso \code{\link[ssMRCD]{msPCA}}, \code{\link[ssMRCD]{align}},
#' \code{\link[ssMRCD]{screeplot.msPCA}}, , \code{\link[ssMRCD]{biplot.msPCA}},
#' \code{\link[ssMRCD]{scores}},
#'
#' @exportS3Method plot msPCA
#'
#' @examples
# \donttest{
#' # set seed
#' set.seed(236)
#'
#' # create data and setup
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariances
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # calculate sparse PCA
#' pca = msPCA(eta = 1.3, gamma = 0.7, COVS = covs$MRCDcov)
#'
#' # plot screeplot
#' plot(x = pca, type = "screeplot")
#'
#' # align and plot loadings
#' pca$PC = align(PC = pca$PC, type = "mean")
#' plot(x = pca, type = "loadings", k = 1)
#' pca$PC = align(PC = pca$PC, type = "maxvar")
#' plot(x = pca, type = "loadings", k = 1)
#' pca$PC = align(PC = pca$PC, type = "largest")
#' plot(x = pca, type = "loadings", k = 1)
#'
#' # plot different PCA plots
#' plot(x = pca, type = "score_distances", k = 2,
#' groups = groups, X = data, mu = covs$MRCDmu, Sigma = covs$MRCDcov)
#' plot(x = pca, type = "biplot", color = "variable")
#' plot(x = pca, type = "scores", ssMRCD = covs, k = 1)
#' plot(x = pca, type = "loadings", k = 1)
# }

plot.msPCA = function(x,
                       type = c("loadings"),
                       ...){

  type = match.arg(type, choices = c("loadings", "screeplot", "scores", "biplot", "score_distances"), several.ok = TRUE)

  g_all = list()
  if("loadings" %in% type){
    g = plot_loadings(object = x, ...)
    g_all = c(list(g_all, g))
  }
  if("screeplot" %in% type){
    g = screeplot(x = x, ...)
    g_all = c(list(g_all, g))
  }
  if("scores" %in% type){
    g = plot_scores(PC = x$PC, ...)
    g_all = c(list(g_all, g))
  }
  if("score_distances" %in% type){
    g = plot_score_distances(PC = x$PC, ...)
    g_all = c(list(g_all, g))
  }
  if("biplot" %in% type){
    g = biplot(x = x, ...)
    g_all = c(list(g_all, g))
  }

  names(g_all) <- type
  if (length(g_all) == 1) {
    return(g_all[[1]])
  }
  return(g_all)
}


#' Scree Plot and Cumulative Explained Variance for msPCA Objects
#'
#' Generates a scree plot displaying the explained variance of principal components
#' and a heatmap of cumulative explained variance per group for an object of class \code{msPCA}.
#'
#' @param x An object of class \code{msPCA}.
#' @param ... Additional arguments to customize the plot:
#' \describe{
#'   \item{\code{text}}{Logical; whether to display numeric values on the heatmap (default: \code{TRUE}).}
#'   \item{\code{size}}{Numeric; text size for labels in the heatmap (default: 5).}
#'   \item{\code{cutoff}}{Numeric; cutoff threshold for explained variance lines and color midpoint, as a proportion (default: 0.8).}
#'   \item{\code{gnames}}{Character vector; optional custom group names (default: \code{N1}, \code{N2}, ..., \code{Nn}).}
#'   \item{\code{textrotate}}{Numeric; rotation angle of text labels on the heatmap (default: 90 degrees).}
#' }
#'
#' @return A list containing two \code{ggplot2} plot objects:
#' \itemize{
#'   \item A boxplot-based scree plot showing explained variance per principal component across groups.
#'   \item A tile plot showing cumulative explained variance per group and principal component.
#' }
#'
#' @examples
#' set.seed(236)
#' data <- matrix(rnorm(1500), ncol = 5)
#' groups <- sample(1:5, 300, replace = TRUE)
#' W <- time_weights(N = 5, c(3,2,1))
#' covs <- ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#' pca <- msPCA(eta = 0.3, gamma = 0.7, COVS = covs$MRCDcov)
#' screeplot(pca, text = TRUE, cutoff = 0.8, size = 4, textrotate = 0)
#' screeplot(pca, text = FALSE, cutoff = 0.6)
#'
#' @exportS3Method screeplot msPCA
#' @import ggplot2
#' @importFrom stats screeplot
screeplot.msPCA = function(x, ...){

  args = list(...)
  if(is.null(args$size)) args$size = 5
  if(is.null(args$text)) args$text = TRUE
  if(is.null(args$cutoff)) args$cutoff = 0.8
  if(is.null(args$gnames)) args$gnames = paste0("N", 1:x$N)
  if(is.null(args$textrotate)) args$textrotate = 90

  # prepare data: explained variance per neighborhood and component
  explvar = matrix(NA, x$N, x$k)
  for(i in 1:x$k){
    for(j in 1:x$N) {
      ind = jth_col(j = j,p = x$p)
      explvar[j,i] =  (t(x$PC[, i, j]) %*%  x$COVS[[j]] %*% x$PC[, i, j])/sum(diag(x$COVS[[j]]))
    }
  }
  colnames(explvar) = paste0(1:x$k)
  rownames(explvar) = args$gnames

  # cumulative explained variance (CPV)
  meanvar = sapply(1:x$k, function(y) sum(explvar[, 1:y])/x$N)

  # scree plot and global explained variance

  dat_scree = cbind(data.frame(Var1 = rep(args$gnames, times = x$k)),
                    Var2 = rep(1:x$k, each = x$N),
                    value = c(explvar))
  g_screeplot = ggplot2::ggplot() +
    ggplot2::geom_hline(aes(yintercept = args$cutoff),
               col ="grey",
               lty = 2) +
    ggplot2::geom_hline(aes(yintercept = 1),
               col ="grey",
               lty = 1)+
    ggplot2::geom_boxplot(data = dat_scree,
                 aes(y = dat_scree$value,
                     x = as.factor(dat_scree$Var2),
                     group = dat_scree$Var2),
                 fill = "grey90",
                 col = "black") +
    ggplot2::geom_line(aes(x = as.factor(1:x$k),
                  y = meanvar,
                  group = 1),
              col = "black") +
    ggplot2::geom_point(aes(x = as.factor(1:x$k),
                   y = meanvar,
                   group = 1),
               col = "black") +
    ggplot2::scale_y_continuous("Explained Variance",
                       limits = c(0, 1.05),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
    ggplot2::labs(x = "Number of PCs") +
    ggplot2::theme_bw()


  # cumulative explained variance per neighborhood
  explvar_cumu = t(apply(X = explvar,
                         FUN = cumsum,
                         MARGIN = 1))

  dat1 = data.frame(value = c(explvar_cumu),
                    Var1 = factor(rep(args$gnames, times = x$k),
                                  levels = args$gnames,
                                  ordered = TRUE),
                    Var2 = rep(1:x$k, each = length(args$gnames)))
  g_explvar = ggplot2::ggplot() +
    ggplot2::geom_tile(aes(x = dat1$Var1,
                  y = dat1$Var2,
                  fill = dat1$value)) +
    ggplot2::scale_fill_gradient2(name = "Cumulative Explained \nVariance (%)", midpoint = args$cutoff) +
    ggplot2::geom_text(aes(x = dat1$Var1,
                  y = dat1$Var2,
                  label = 100*round(dat1$value, 2)),
              angle = args$textrotate,
              col = ifelse(args$text, "black", "transparent"),
              size = args$size) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 45)) +
    ggplot2::labs(x = "", y = "Number of PCs",
         title = "Cumulative Explained Variance per Group (%)")

  return(list(g_screeplot, g_explvar))
}



#' Biplot Visualization for msPCA Objects
#'
#' Creates a biplot of the first two principal components from an \code{msPCA} object,
#' displaying variables or groups in the component space.
#'
#' @param x An object of class \code{msPCA}.
#' @param ... Additional graphical parameters to customize the plot:
#' \describe{
#'   \item{\code{shape}}{Integer or character; shape of the points (default: 43).}
#'   \item{\code{size}}{Numeric; size of the points/text (default: 3).}
#'   \item{\code{alpha}}{Numeric between 0 and 1; transparency of the points/text (default: 0.7).}
#'   \item{\code{color}}{Character; determines the coloring scheme of points, either \code{"variable"} or \code{"groups"} (default: \code{"variable"}).}
#' }
#'
#' @return A \code{ggplot2} object representing the biplot of the first two principal components.
#'
#' @examples
#' set.seed(236)
#' data <- matrix(rnorm(2000), ncol = 4)
#' groups <- sample(1:10, 500, replace = TRUE)
#' W <- time_weights(N = 10, c(3,2,1))
#' covs <- ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#' pca <- msPCA(eta = 0.3, gamma = 0.5,COVS = covs$MRCDcov, k = 2)
#' pca$PC = align(PC = pca$PC, type = "largest")
#' biplot(pca, alpha = 1, shape = 16, size = 5, color = "variable")
#'
#' @exportS3Method biplot msPCA
#' @import ggplot2
#' @importFrom stats biplot

biplot.msPCA = function(x, ...){

  stopifnot(x$k >= 2)
  args = list(...)

  if(is.null(args$color)) args$color = "variable"
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 3
  if(is.null(args$alpha)) args$alpha = 0.7

  pc = data.frame(apply(x$PC, 2, rbind))
  colnames(pc) = paste0("PC", 1:x$k)
  pc$groups = rep(dimnames(x$PC)[[3]], each = x$p)
  pc$variable = rep(rownames(x$PC), times = x$N)

  g = ggplot2::ggplot() +
    ggplot2::geom_text(aes(x = pc$PC1,
                           y = pc$PC2,
                           col = as.factor(pc[,args$color]),
                           label = pc$groups),
                       size = args$size,
                       alpha = args$alpha) +
    ggplot2::scale_colour_discrete(args$color) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "PC 1", y = "PC 2", title = "Biplot") +
    ggplot2::lims(x = c(-1, 1), y = c(-1, 1))
  g
  return(g)
}



#' @import ggplot2
#' @importFrom scales muted
#' @importFrom scales rescale
plot_loadings = function(object, ...){

  args = list(...)
  if(is.null(args$k)) args$k = 1
  if(is.null(args$size)) args$size = 3
  if(is.null(args$text)) args$text = FALSE

  if(is.null(args$gnames)) args$gnames = paste0("N", 1:object$N)
  if(is.null(args$cnames)) args$cnames = paste0("Var", 1:object$p)
  if(!is.null(names(object$COVS))) args$gnames = names(object$COVS)
  if(!is.null(colnames(object$COVS[[1]])))  args$cnames = colnames(object$COVS[[1]])

  if(is.null(args$textrotate)) args$textrotate = 90
  if(is.null(args$tolerance)) args$tolerance = 0.0001

  PC_flat = cbind(data.frame(value = c(object$PC[,args$k,])),
                  N = factor(rep(args$gnames, each = object$p),
                             levels = args$gnames, ordered = TRUE),
                  variable = factor(rep(args$cnames, times = object$N),
                                    levels = args$cnames, ordered = TRUE))

  g = ggplot2::ggplot() +
    ggplot2::geom_tile(aes(x = PC_flat$N,
                  y = PC_flat$variable,
                  fill = PC_flat$value,
                  col = PC_flat$value)) +
    ggplot2::scale_fill_gradientn(name = "",
                         values = scales::rescale(c(1, args$tolerance,  0,  -args$tolerance, -1)),
                         colours = c(scales::muted("darkred"), "white", "white", "white", scales::muted("darkblue")),
                         limits = c(-1, 1)) +
    ggplot2::scale_colour_gradientn(name = "",
                           values = scales::rescale(c(1, args$tolerance,  0,  -args$tolerance, -1)),
                           colours = c(scales::muted("darkred"), "white", "white", "white", scales::muted("darkblue")),
                           limits = c(-1, 1)) +
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = element_text(angle = 90)) +
    ggplot2::labs(x = "",
         y = "Groups",
         title = paste0("Loadings (PC ", args$k, ")")) +
    ggplot2::geom_text(aes(x = PC_flat$N,
                  y = PC_flat$variable,
                  label = round(PC_flat$value, 2)),
              size = args$size,
              angle = args$textrotate,
              col = ifelse(args$text, "black", "transparent"))

  return(g)
}

#' @import ggplot2
plot_scores = function(PC, ssMRCD = NULL, X = NULL, groups = NULL, mu = NULL, Sigma = NULL, ...){

  args = list(...)
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 5
  if(is.null(args$alpha)) args$alpha = 0.7
  if(is.null(args$k)) args$k = 1

  if(!is.null(ssMRCD) & is.null(groups)) {
    groups = rep(ssMRCD$gnames, times= sapply(ssMRCD$mX, nrow))
  }

  sc = scores( PC[, args$k, , drop = FALSE], ssMRCD = ssMRCD,
               X = X, groups = groups, mu = mu, Sigma = Sigma)$scores
  sc = data.frame(sc, groups = groups)

  g = ggplot2::ggplot(data = sc) +
    ggplot2::geom_histogram(aes(x = sc$Score1,
                       fill = as.factor(sc$groups)),
                   alpha = args$alpha) +
    ggplot2::facet_wrap(vars(sc$groups)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Score value",
         y = "Count",
         title = paste0("Scores for PC ", args$k, " per Group"))
  return(g)
}


#' @import ggplot2
plot_score_distances = function(X = NULL, PC, groups = NULL, ssMRCD = NULL, mu = NULL, Sigma = NULL, k, ...){

  args = list(...)
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 3
  if(is.null(args$alpha)) args$alpha = 0.7
  if(is.null(args$k)) args$k = dim(PC)[2]

  if(!is.null(ssMRCD) & is.null(groups)) {
    groups = rep(ssMRCD$gnames, times= sapply(ssMRCD$mX, nrow))
  }

  sc_all = scores( PC[, args$k, , drop = FALSE], ssMRCD = ssMRCD,
               X = X, groups = groups, mu = mu, Sigma = Sigma)

  sd = sc_all$SD
  od = sc_all$OD

  sc = data.frame(sd = sd,
                  od = od,
                  groups = groups)


  g = ggplot2::ggplot(sc) +
    ggplot2::geom_point(aes(x = sd,
                   y = od,
                   col = as.factor(groups)),
               shape = args$shape,
               size = args$size,
               alpha = args$alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Score Distance",
         y = "Orthogonal Distance",
         title = "Distance-Distance Plot") +
    ggplot2::scale_colour_discrete("Group") +
    ggplot2::facet_wrap(vars(groups))
  g

  return(g)
}


##########################################################################################
#### SUMMARY METHOD                                                                   ####
##########################################################################################


#' Summary Method for msPCA Objects
#'
#' Provides a summary of a sparse multi-source PCA (`msPCA`) result, including group-wise sparsity
#' and explained variance.
#'
#' @param object An object of class \code{"msPCA"} as returned by the main msPCA fitting function.
#' @param ... Currently ignored.
#'
#' @return Prints a summary to the console. No return value.
#' @exportS3Method summary msPCA
summary.msPCA <- function(object, ...) {

  args = list(...)
  if(is.null(args$k)) args$k = object$k

  cat("\n")
  cat("Sparse Multi-Source PCA Summary\n")
  cat("==============================================\n")

  cat("Number of groups (N):              ", object$N, "\n")
  cat("Number of variables (p):           ", object$p, "\n")
  cat("Number of components (k):          ", object$k, "\n")
  cat("Sparsity parameters (gamma, eta):  ", object$gamma, ",", object$eta, "\n")

  group_names <- names(object$COVS)
  if (is.null(group_names)) {
    group_names <- paste0("Group", seq_len(object$N))
  }

  cat("\nComponent-wise Summary:\n")
  for (j in seq_len(args$k)) {
    cat(paste0("\nPC", j, ":\n"))
    cat("Group    | Explained Var (%) | Non-zero Loadings\n")
    cat("-----------------------------------------------\n")
    for (g in seq_len(object$N)) {
      cov_g <- object$COVS[[g]]
      loading <- object$PC[, j, g]
      var_total <- sum(diag(cov_g))
      var_explained <- as.numeric(t(loading) %*% cov_g %*% loading)
      var_percent <- round(100 * var_explained / var_total, 2)
      nonzero <- sum(abs(loading) > 1e-8)
      cat(sprintf("%-9s| %7.2f%%          | %3d / %d\n",
                  group_names[g], var_percent, nonzero, length(loading)))
    }
  }


  cat("\nConvergence:\n")
  cat("  Converged: ", ifelse(object$converged, "Yes", "No"), "\n")
  cat("  Steps:     ", object$n_steps, "\n")

  cat("\nUse `plot()` to visualize components.\n")
}
