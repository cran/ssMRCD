##########################################################################################
#### ALIGN PCs                                                                        ####
##########################################################################################
#' Align Loadings of Principal Components
#'
#' @description
#' Aligns loadings per neighborhood for better visualization and comparison. Different options are available.
#'
#' @param PC Array of loadings of size p x k x N as returned by \code{msPCA$PC}.
#' @param type Character string specifying the sign adjustment method. One of:
#' \tabular{ll}{
#' \code{"largest"} \tab Sets the sign so that the largest (by absolute value) loading is positive, per neighborhood and component. \cr
#' \code{"maxvar"}  \tab Makes the variable with the highest absolute value positive, per neighborhood and component. \cr
#' \code{"nonzero"} \tab Makes the loading with the largest distance from zero over neighborhoods positive, per neighborhood and component. \cr
#' \code{"scalar"}  \tab Adjusts signs so that the scalar product with \code{vec} is positive, per neighborhood and component.
#' \code{vec} can be a p-dimensional vector used across all components and neighborhoods, a (p x k) matrix for individual vectors for components,
#' or a (pN x k) matrix with individual vectors for each neighborhood and component \cr
#' \code{"mean"}    \tab Like \code{"scalar"}, but \code{vec} is the mean of loadings across neighborhoods for each components. \cr
#' \code{"none"}    \tab No sign adjustment; returns the raw \code{PC}. \cr
#' }
#' @param vec \code{NULL} or vector containing vectors for type \code{"scalar"}.
#' If \code{vec} is of dimension \code{p × k}, the same \code{vec} is used for all neighborhoods.

#'
#' @return Returns an array of loadings of size \code{p} times \code{k} times \code{N}.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' # Note, that in this example the vectors are not feasible loadings.
#' COVS = list("a"= matrix(runif(16, -1, 1), 4), "b" = matrix(runif(16, -1, 1), 4))
#' COVS = lapply(COVS, function(x) x %*% t(x))
#'
#' pca = msPCA(eta = 1, gamma = 0.5, COVS = COVS, k = 2)
#' x = pca$PC
#'
#' align(PC = x, type = "largest")
#' align(PC = x, type = "maxvar")
#' align(PC = x, type = "mean")
#' align(PC = x, type = "nonzero")
#' align(PC = x, type = "scalar", vec = c(1,1,1,1))

align = function(PC, type = "largest", vec = NULL){

  N = dim(PC)[3]
  k = dim(PC)[2]
  p = dim(PC)[1]

  type_orig = type
  if(type == "scalar" & is.null(vec)) stop("Necessary 'vec' for option 'type = scalar' is missing.")
  if(!is.null(vec)) {
    vec = as.matrix(vec)
    if(dim(vec)[2] == 1){
      vec = matrix(rep(vec, times = k), ncol = k)
    }
    if(dim(vec)[1] == p) {
      vec = matrix(rep(t(vec), times = N), nrow = p*N, ncol = k, byrow = TRUE)
    }
  }

  if(type == "none") return(PC)

  for(l in 1:k){
    x = PC[, l, ]

    # setup
    if(type == "scalar")  veci = vec[, l]

    if(type_orig == "maxvar"){
      m = t(abs(x))
      var_ind = which.max(colSums(m))
    }
    if(type_orig == "nonzero"){
      m = t(abs(x))
      var_ind = which.max(sapply(as.data.frame(t(m)), min))  # take variable where abs values are furthest away from zero
    }
    if(type_orig == "mean"){
      type ="scalar"
      veci = rep(rowMeans(x), times = N)
    }

    # align
    for(i in 1:N){
      x_i = x[, i]
      if(type == "maxvar" | type == "nonzero"){
        mult = sign(x_i[var_ind])
      }
      if(type == "largest"){
        # largest absolute value will be positive
        mult = sign(x_i[which.max(abs(x_i))])
      }
      if(type == "scalar"){
        # scalar product is positive
        a = matrix(veci[jth_col(j=i, p= p)], ncol = 1)
        b = matrix(x_i, nrow = 1)
        mult = c(sign(b %*% a))
      }
      if(mult == 0) {
        # if orthogonal: default to largest
        mult = unlist(sign(x_i[which.max(abs(x_i))]))
      }
      if(mult == 0) mult = 1

      # redirect
      x[, i] = x[, i]*mult
    }
    PC[, l, ] = x
  }
  return(PC)
}

##########################################################################################
#### SCORES AND DISTANCES                                                             ####
##########################################################################################

#' Calculate Scores and Distances for Multi-Source PCA
#'
#' Computes principal component scores, score distances (SD), and orthogonal distances (OD)
#' for observations grouped into multiple sources using the multi-source PCA model. The function supports
#' either an `ssMRCD` object for robust local centering and scaling or manually provided group-wise
#' means (`mu`) and covariances (`Sigma`).
#'
#' @param PC A 3D array representing the principal component loading matrices for each group
#' (dimensions: variables × components × groups).
#' @param ssMRCD An optional `ssMRCD` object used to robustly center and scale `X`. If `NULL`,
#' then `X`, `groups`, `mu` and `Sigma`, must be provided.
#' @param X An optional matrix of observations (rows are samples, columns are variables), required if `ssMRCD` is not provided.
#' @param groups An optional numeric vector specifying group/source membership for each observation in `X`, required if `ssMRCD` is not provided.
#' @param mu Optional list of group-wise means, required if `ssMRCD` is not provided.
#' @param Sigma Optional list of group-wise covariance matrices, required if `ssMRCD` is not provided.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{scores}{Matrix of principal component scores for each observation.}
#'   \item{SD}{Numeric vector of score distances, i.e., Mahalanobis distances in the PCA space.}
#'   \item{OD}{Numeric vector of orthogonal distances (reconstruction error orthogonal to PCA space).}
#'   \item{X_centered}{Locally centered input data.}
#' }
#'
#' @export
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}}, \code{\link[ssMRCD]{scale.ssMRCD}},  \code{\link[ssMRCD]{msPCA}}
#'
#' @examples
#' # create data set
#' x1 = matrix(runif(200), ncol = 2)
#' x2 = matrix(rnorm(200), ncol = 2)
#' x = list(x1, x2)
#'
#' # create weighting matrix
#' W = matrix(c(0, 1, 1, 0), ncol = 2)
#'
#' # calculate ssMRCD
#' loccovs = ssMRCD(x, weights = W, lambda = 0.5)
#'
#' # calculate PCA
#' pca = msPCA(eta = 1, gamma = 0.5,COVS = loccovs$MRCDcov)
#'
#' # calculate scores
#' scores_all = scores(PC = pca$PC, ssMRCD = loccovs)
#' str(scores_all)
#'
#' scores_all = scores(PC = pca$PC,
#'                     X = rbind(x1, x2),
#'                     groups = rep(c(1,2), each = 100),
#'                     mu = loccovs$MRCDmu,
#'                     Sigma = loccovs$MRCDcov)
#' str(scores_all)


scores = function( PC, ssMRCD = NULL, X = NULL, groups = NULL, mu = NULL, Sigma = NULL){

  N = dim(PC)[3]
  k = dim(PC)[2]
  p = dim(PC)[1]

  ssmrcd_scale = list()
  if(is.null(ssMRCD)){
    if(is.null(mu) | is.null(Sigma) | is.null(groups) | is.null(X)) stop("Either ssMRCD object or list of means 'mu', a list of covariances 'Sigma', data 'X' and vector of 'groups' need to be given.")
    ssmrcd_scale$MRCDmu = mu
    ssmrcd_scale$MRCDcov = Sigma
    ssmrcd_scale$N = N
    X = as.matrix(X)
    ssmrcd_scale$mX = restructure_as_list(X, groups)
    groups = rep(1:N, times= sapply(ssmrcd_scale$mX, nrow))
  } else if("ssMRCD" %in% class(ssMRCD)) {
    ssmrcd_scale = ssMRCD
    X = do.call(rbind, ssMRCD$mX)
    groups = rep(1:N, times= sapply(ssMRCD$mX, nrow))
  }
  class(ssmrcd_scale) = "ssMRCD"

  # center data
  X_centered = scale.ssMRCD(x = ssmrcd_scale, center_only = TRUE)
  n = nrow(X_centered)

  # naming
  TT = matrix(NA, ncol = k, nrow = nrow(X), dimnames = list(NULL, paste0("Score",1:k)))

  # calculate scores
  for(l in 1:N){
    ind_obs = which(groups == l)
    centered = X_centered[ind_obs, ]
    TT[ind_obs, ] = centered %*% PC[, , l ]
  }

  # calculate orthogonal distances
  OD = rep(NA, n)
  for(l in 1:N){
    ind_obs = which(groups == l)
    tmp = X_centered[ind_obs, ] -  TT[ind_obs, ] %*% t(PC[, ,l ])
    OD[ind_obs] = sqrt(diag(tmp %*% t(tmp))) # norm per observation
  }

  # calculate score distances
  COVS = ssmrcd_scale$MRCDcov

  SD = rep(NA, n)
  for(j in 1:N){
    ind_obs = which(groups == j)
    lambdainv = diag(t(PC[, ,j ]) %*% COVS[[j]] %*% PC[, , j ])^(-1)
    SD[ind_obs] = diag(TT[ind_obs, ] %*% diag(lambdainv, ncol = k) %*% t(TT[ind_obs, ]))
  }

  # return scores
  return(list(scores = TT,
              SD = SD,
              OD = OD,
              X_centered = X_centered))  # X_centered also scaled if ssMRCD is used
}
