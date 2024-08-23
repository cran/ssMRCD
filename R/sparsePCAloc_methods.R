# METHOD FOR PCAloc OBJECTS AND PLOTS




#' Summary method for PCAloc
#'
#' @param object object of class PCAloc
#' @param ... other input variables.
#'
#' @seealso \code{\link[ssMRCD]{sparsePCAloc}}
#'
#' @return Summary for PCAloc
#' @exportS3Method summary PCAloc
#'
#' @examples#'
#' C1 = diag(c(1.1, 0.9, 0.6))
#' C2 = matrix(c(1.1, 0.1, -0.1,
#'               0.1, 1.0, -0.2,
#'              -0.1, -0.2, 0.7), ncol = 3)
#' C3 = (C1 + C2)/2
#'
#' pca = sparsePCAloc(eta = 1, gamma = 0.5, cor = FALSE, COVS = list(C1, C2, C3),
#'              n_max = 100, increase_rho = list(FALSE, 100, 1), trace = FALSE)
#'
#' summary(pca)

summary.PCAloc = function(object, ...){

  rownames(pca$PC) = rep(paste0("Var", 1: (dim(pca$PC)[1]/pca$N)), times = pca$N)
  cat("Loadings: \n")
  print(pca$PC)
  cat("\nParameters: \n")
  cat("Gamma =", pca$gamma, "\nEta =", pca$eta)
}




#' Plotting method PCAloc object
#'
#' @param x object of class PCAloc
#' @param type character indicating the type of plot, see details.
#' @param ... further arguments passed down.
#'
#' @return Returns plots in ggplot2.
#' @exportS3Method plot PCAloc
#'
#' @examples
#' \donttest{
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
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # align loadings
#' pca$PC = align_PC(PC = pca$PC, N = pca$N, p = pca$p, type = "mean")
#'
#' # plot different PCA plots
#' plot(x = pca, type = "score_distances", groups = groups, X = data, ssMRCD = covs, k = 2)
#' plot(x = pca, type = "biplot", color = "variable")
#' plot(x = pca, type = "scores", groups = groups, X = data, ssMRCD = covs, k = 1)
#' plot(x = pca, type = "screeplot")
#' plot(x = pca, type = "loadings", k = 1)
#' }

plot.PCAloc = function(x,
                       type = c("loadings", "screeplot", "scores", "score_distances", "biplot"),
                       ...){
  # additional inputs: X, groups, ssMRCD, k
  # color, shape, size, alpha
  # cutoff, groupnames, textrotate

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

  return(g_all)
}




#' Plots of loadings of PCAloc object
#'
#' @param object object of class PCAloc
#' @param ... other input arguments, see details.
#'
#' @return Returns loading heatmap for component \code{k}.
#' @export
#'
#' @details
#' Additional parameters that can be given to the function are: \tabular{ll}{
#'    \code{text} \tab logical if values should be added as text.  \cr
#'    \tab \cr
#'    \code{size} \tab point size.\cr
#'    \tab \cr
#'    \code{tolerance} \tab tolerance for rounding to zero.\cr
#'    \tab \cr
#'    \code{k} \tab integer, which component scores should be plotted. \cr
#'    \tab \cr
#'    \code{groupnames} \tab names of groups. \cr
#'    \tab \cr
#'    \code{varnames} \tab names of variables. \cr
#'    \tab \cr
#'    \code{textrotate} \tab angle of text rotation, if included.\cr
#'    \tab \cr
#' }
#'
#'
#' @examples
#' # set seed
#' set.seed(236)
#'
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariance matrices
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # sparse PCA
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # plot score distances
#' plot_loadings(object = pca,
#'             k = 1,
#'             size = 2)
#' @import ggplot2
#' @importFrom scales muted
#' @importFrom scales rescale
plot_loadings = function(object, ...){

  args = list(...)
  if(is.null(args$k)) args$k = 1
  if(is.null(args$size)) args$size = 2
  if(is.null(args$text)) args$text = TRUE
  if(is.null(args$groupnames)) args$groupnames = paste0("N", 1:object$N)
  if(is.null(args$varnames)) args$varnames = paste0("Var", 1:object$p)
  if(is.null(args$textrotate)) args$textrotate = 90
  if(is.null(args$tolerance)) args$tolerance = 0.005

  PC_flat = cbind(data.frame(object$PC[,args$k]),
                  N = rep(args$groupnames, each = object$p),
                  variable = rep(args$varnames, times = object$N)) %>%
    reshape2::melt(id.vars= c("N", "variable"))
  PC_flat$variable = factor(PC_flat$variable, levels = sort(args$varnames), ordered = T)
  PC_flat$N = factor(PC_flat$N, levels = sort(args$groupnames), ordered = T)

  g = ggplot() +
    geom_tile(aes(x = PC_flat$N,
                  y = PC_flat$variable,
                  fill = PC_flat$value,
                  col = PC_flat$value)) +
    scale_fill_gradientn(name = "",
                         values = scales::rescale(c(1, args$tolerance,  0,  -args$tolerance, -1)),
                         colours = c(scales::muted("darkred"), "white", "white", "white", scales::muted("darkblue")),
                         limits = c(-1, 1)) +
    scale_colour_gradientn(name = "",
                           values = scales::rescale(c(1, args$tolerance,  0,  -args$tolerance, -1)),
                           colours = c(scales::muted("darkred"), "white", "white", "white", scales::muted("darkblue")),
                           limits = c(-1, 1)) +
    scale_y_discrete(limits=rev) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "",
         y = "",
         title = paste0("Loadings of PC ", args$k)) +
    geom_text(aes(x = PC_flat$N,
                  y = PC_flat$variable,
                  label = round(PC_flat$value, 2)),
              size = args$size,
              angle = args$textrotate,
              col = ifelse(args$text, "black", "transparent"))

  return(g)
}



#' Plots of score distribution
#'
#' @param X data matrix.
#' @param PC loadings from PCA.
#' @param groups vector containing group assignments.
#' @param ssMRCD ssMRCD object.
#' @param ... other input arguments, see details.
#'
#' @return Returns histograms of scores for component \code{k}.
#' @export
#'
#' @details
#' Additional parameters that can be given to the function are: \tabular{ll}{
#'    \code{shape} \tab point shape  \cr
#'    \tab \cr
#'    \code{size} \tab point size \cr
#'    \tab \cr
#'    \code{alpha} \tab transparency  \cr
#'    \tab \cr
#'    \code{k} \tab integer, which component scores should be plotted \cr
#'    \tab \cr
#' }
#'
#'
#' @examples
#' \donttest{
#' # set seed
#' set.seed(236)
#'
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariance matrices
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # sparse PCA
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # plot score distances
#' plot_scores(PC = pca$PC,
#'             groups = groups,
#'             X = data,
#'             ssMRCD = covs,
#'             k = 1,
#'             alpha = 0.4,
#'             shape = 16,
#'             size = 2)
#'}
#' @import ggplot2
plot_scores = function(X, PC, groups, ssMRCD, ...){

  args = list(...)
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 5
  if(is.null(args$alpha)) args$alpha = 0.7
  if(is.null(args$k)) args$k = 1

  sc = scores(X = X,
              PC = PC[, args$k],
              groups = groups,
              ssMRCD = ssMRCD)$scores
  sc = data.frame(sc, groups = groups)

  g = ggplot(sc) +
    geom_histogram(aes(x = sc,
                       fill = as.factor(groups)),
                   alpha = args$alpha) +
    facet_wrap(vars(groups)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Score value",
         y = "Count",
         title = paste0("Scores for PC ", args$k, " per group"))

  return(g)
}


#' Biplot for PCAloc
#'
#' @param x object of class PCAloc.
#' @param ... other input arguments, see details.
#'
#' @return Returns version of biplot for PCAloc object.
#' @exportS3Method biplot PCAloc
#'
#' @details
#' Additional parameters that can be given to the function are: \tabular{ll}{
#'    \code{shape} \tab point shape  \cr
#'    \tab \cr
#'    \code{size} \tab point size \cr
#'    \tab \cr
#'    \code{alpha} \tab transparency  \cr
#'    \tab \cr
#'    \code{color} \tab either \code{"variable"} or \code{"groups"}
#'    indication how points should be coloured.   \cr
#'    \tab \cr
#' }
#'
#'
#' @examples
#' # set seed
#' set.seed(236)
#'
#' # make data
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariance matrices
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # sparse PCA
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # plot biplot
#' biplot(pca, alpha = 0.4, shape = 16, size = 2, color = "variable")
#' @import ggplot2
#' @importFrom stats biplot
biplot.PCAloc = function(x, ...){

  stopifnot(x$k >= 2)
  args = list(...)

  if(is.null(args$color)) args$color = "variable"
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 5
  if(is.null(args$alpha)) args$alpha = 0.7

  pc = data.frame(x$PC)
  colnames(pc) = paste0("PC", 1:x$k)
  pc$groups = rep(1:x$N, each = x$p)
  pc$variable = rep(1:x$p, times = x$N)

  g = ggplot(pc) +
    geom_point(aes(x = pc$PC1,
                   y = pc$PC2,
                   col = as.factor(get(args$color))),
               shape = args$shape,
               size = args$size,
               alpha = args$alpha) +
    scale_colour_discrete(args$color) +
    theme_minimal() +
    lims(x = c(-1, 1), y = c(-1, 1))

  return(g)
}


#' Screeplot for PCAloc
#'
#' @param x object of class PCAloc.
#' @param ... other input arguments, see details.
#'
#' @return Returns version of scree plot and cumulative explained variance per group for PCAloc object.
#' @exportS3Method screeplot PCAloc
#'
#' @details
#' Additional parameters that can be given to the function are: \tabular{ll}{
#'    \code{text} \tab logical if text should be plotted  \cr
#'    \tab \cr
#'    \code{size} \tab text size \cr
#'    \tab \cr
#'    \code{cutoff} \tab cutoff line for scree plot  \cr
#'    \tab \cr
#'    \code{groupnames} \tab name of groups   \cr
#'    \tab \cr
#'     \code{textrotate} \tab angle of text, if text is plotted.   \cr
#'    \tab \cr
#' }
#'
#'
#' @examples
#' # set seed
#' set.seed(236)
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariance matrices
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # sparse PCA
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # plot biplot
#' screeplot(pca, text = TRUE, cutoff = 0.8, size = 2)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats screeplot
screeplot.PCAloc = function(x, ...){

  args = list(...)
  if(is.null(args$size)) args$size = 5
  if(is.null(args$text)) args$text = TRUE
  if(is.null(args$cutoff)) args$cutoff = 0.8
  if(is.null(args$groupnames)) args$groupnames = paste0("N", 1:x$N)
  if(is.null(args$textrotate)) args$textrotate = 90

  # prepare data: explained variance per neighborhood and component
  explvar = matrix(NA, x$N, x$k)
  for(i in 1:x$k){
    for(j in 1:x$N) {
      ind = jth_col(j = j,p = x$p)
      explvar[j,i] =  (t(x$PC[ind, i]) %*%  x$COVS[[j]] %*% x$PC[ind, i])/sum(diag(x$COVS[[j]]))
    }
  }
  colnames(explvar) = paste0(1:x$k)
  rownames(explvar) = args$groupnames

  # cumulative explained variance (CPV)
  meanvar = sapply(1:x$k, function(y) sum(explvar[, 1:y])/x$N)

  # scree plot and global explained variance
  dat_scree = reshape2::melt(explvar)
  g_screeplot = ggplot() +
    geom_hline(aes(yintercept = args$cutoff),
               col ="grey",
               lty = 2) +
    geom_hline(aes(yintercept = 1),
               col ="grey",
               lty = 1)+
    geom_boxplot(data = dat_scree,
                 aes(y = dat_scree$value,
                     x = as.factor(dat_scree$Var2),
                     group = dat_scree$Var2),
                 fill = "grey90",
                 col = "black") +
    geom_line(aes(x = as.factor(1:x$k),
                  y = meanvar,
                  group = 1),
              col = "black") +
    geom_point(aes(x = as.factor(1:x$k),
                   y = meanvar,
                   group = 1),
               col = "black") +
    scale_y_continuous("Explained Variance",
                       limits = c(0, 1.05),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       labels = c("0%", "20%", "40%", "60%", "80%", "100%")) +
    labs(x = "Number of principal components") +
    theme_minimal()


  # cumulative explained variance per neighborhood
  explvar_cumu = t(apply(X = explvar,
                         FUN = cumsum,
                         MARGIN = 1))

  dat1 = explvar_cumu %>% reshape2::melt()
  g_explvar = ggplot() +
    geom_tile(aes(x = dat1$Var1,
                  y = dat1$Var2,
                  fill = dat1$value)) +
    scale_fill_gradient2(midpoint = args$cutoff) +
    geom_text(aes(x = dat1$Var1,
                  y = dat1$Var2,
                  label = 100*round(dat1$value, 2)),
              angle = args$textrotate,
              col = ifelse(args$text, "black", "transparent")) +
    theme(axis.text.x = element_text(angle = 0),
          legend.position = "None") +
    labs(x = "", y = "Number of PCs",
         title = "Cumulative explained variance (%)") +
    theme_minimal()

  return(list(g_screeplot, g_explvar))
}







#' Distance-distance plot of scores of PCA
#'
#' @param X data matrix.
#' @param PC loadings from PCA.
#' @param groups vector containing group assignments.
#' @param ssMRCD ssMRCD object.
#' @param k integer of how many components should be used.
#' @param ... other input arguments, see details.
#'
#' @return Returns distance-distance plot of orthogonal and score distance.
#' @export
#'
#' @details
#' Additional parameters that can be given to the function are: \tabular{ll}{
#'    \code{shape} \tab point shape  \cr
#'    \tab \cr
#'    \code{size} \tab point size \cr
#'    \tab \cr
#'    \code{alpha} \tab transparency  \cr
#'    \tab \cr
#' }
#'
#'
#' @examples
#' # set seed
#' set.seed(236)
#'
#' data = matrix(rnorm(2000), ncol = 4)
#' groups = sample(1:10, 500, replace = TRUE)
#' W = time_weights(N = 10, c(3,2,1))
#'
#' # calculate covariance matrices
#' covs = ssMRCD(data, groups = groups, weights = W, lambda = 0.3)
#'
#' # sparse PCA
#' pca = sparsePCAloc(eta = 0.3, gamma = 0.7, cor = FALSE, COVS = covs$MRCDcov,
#'              n_max = 1000, increase_rho = list(TRUE, 50, 1), trace = FALSE)
#'
#' # plot score distances
#' plot_score_distances(PC = pca$PC,
#'                      groups = groups,
#'                      X = data,
#'                      ssMRCD = covs,
#'                      k = 2,
#'                      alpha = 0.4,
#'                      shape = 16,
#'                      size = 2)
#' @import ggplot2
plot_score_distances = function(X, PC, groups, ssMRCD, k, ...){

  args = list(...)
  if(is.null(args$shape)) args$shape = 43
  if(is.null(args$size)) args$size = 5
  if(is.null(args$alpha)) args$alpha = 0.7

  sd = scores.SD(X = X,
                 PC = PC[, 1:k],
                 groups = groups,
                 ssMRCD = ssMRCD)
  od = scores.OD(X = X,
                 PC = PC[, 1:k],
                 groups = groups,
                 ssMRCD = ssMRCD)
  sc = data.frame(sd = sd,
                  od = od,
                  groups = groups)


  g = ggplot(sc) +
    geom_point(aes(x = sd,
                   y = od,
                   col = as.factor(groups)),
               shape = args$shape,
               size = args$size,
               alpha = args$alpha) +
    theme_minimal() +
    labs(x = "Score distance",
         y = "Orthogonal distance",
         title = "Distance-distance plot") +
    scale_colour_discrete("Group")

  return(g)
}
