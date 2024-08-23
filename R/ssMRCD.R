#' Spatially Smoothed MRCD Estimator
#'
#' The ssMRCD function calculates the spatially smoothed MRCD estimator from Puchhammer and Filzmoser (2023).
#'
#' @param x a list of matrices containing the observations per neighborhood sorted which can be obtained by the function \code{\link[ssMRCD]{restructure_as_list}}, or matrix or data frame containing data.
#'          If matrix or data.frame, group vector has to be given.
#' @param groups vector of neighborhood assignments
#' @param weights weighting matrix, symmetrical, rows sum up to one and diagonals need to be zero (see also \code{\link[ssMRCD]{geo_weights}} or \code{\link[ssMRCD]{rescale_weights}} .
#' @param lambda numeric between 0 and 1.
#' @param TM target matrix (optional), default value is the covMcd from robustbase.
#' @param alpha numeric, proportion of values included, between 0.5 and 1.
#' @param maxcond optional, maximal condition number used for rho-estimation.
#' @param maxcsteps maximal number of c-steps before algorithm stops.
#' @param n_initialhsets number of initial h-sets, default is 6 times number of neighborhoods.
#'
#' @return An object of class \code{"ssMRCD"} containing the following elements:\tabular{ll}{
#'    \code{MRCDcov} \tab List of ssMRCD-covariance matrices sorted by neighborhood. \cr
#'    \tab \cr
#'    \code{MRCDicov} \tab List of inverse ssMRCD-covariance matrices sorted by neighborhood. \cr
#'    \tab \cr
#'    \code{MRCDmu} \tab List of ssMRCD-mean vectors sorted by neighborhood. \cr
#'    \tab \cr
#'    \code{mX} \tab List of data matrices sorted by neighborhood.\cr
#'    \tab \cr
#'    \code{N} \tab Number of neighborhoods. \cr
#'    \tab \cr
#'    \code{mT} \tab Target matrix. \cr
#'    \tab \cr
#'    \code{rho} \tab Vector of regularization values sorted by neighborhood. \cr
#'    \tab \cr
#'    \code{alpha} \tab Scalar what percentage of observations should be used. \cr
#'    \tab \cr
#'    \code{h} \tab Vector of how many observations are used per neighborhood, sorted. \cr
#'    \tab \cr
#'    \code{numiter} \tab The number of iterations for the best initial h-set combination. \cr
#'    \tab \cr
#'    \code{c_alpha} \tab Consistency factor for normality. \cr
#'    \tab \cr
#'    \code{weights} \tab The weighting matrix. \cr
#'    \tab \cr
#'    \code{lambda} \tab Smoothing factor. \cr
#'    \tab \cr
#'    \code{obj_fun_values} \tab A matrix with objective function values for all
#' initial h-set combinations (rows) and iterations (columns). \cr
#'    \tab \cr
#'    \code{best6pack} \tab initial h-set combinations with best objective function value
#' after c-step iterations. \cr
#'    \code{Kcov} \tab returns MRCD-estimates without smoothing. \cr
#' }
#'
#'
#'
#' @seealso \code{\link[ssMRCD]{plot.ssMRCD}, \link[ssMRCD]{summary.ssMRCD}, \link[ssMRCD]{restructure_as_list}}
#'
#' @references Puchhammer P. and Filzmoser P. (2023): Spatially smoothed robust covariance estimation for local outlier detection. \doi{10.48550/arXiv.2305.05371}
#'
#' @export
#' @importFrom robustbase covMcd
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
#' ssMRCD(x, weights = W, lambda = 0.5)

ssMRCD = function(x,
                  groups = NULL,
                  weights,
                  lambda,
                  TM = NULL,
                  alpha = 0.75,
                  maxcond = 50,
                  maxcsteps = 200,
                  n_initialhsets = NULL){

  # input checks
  check_input(alpha, "alpha")
  check_input(weights, "W")
  check_input(lambda, "lambda")

  if( (is.data.frame(x) | is.matrix(x)) & !is.null(groups)) {
    x = as.matrix(x)
    x = restructure_as_list(data = x, groups = groups)
  } else if (!is.list(x)) {
    warning("x should either be list of data matrices or data matrix with grouping given by groups input.")
  }

  N <- length(x)
  p <- dim(x[[1]])[2]
  if(is.null(n_initialhsets)) n_initialhsets <- N*6

  # target matrix
  Xsum <- do.call(rbind, x)

  if(is.null(TM)){
    TM <- tryCatch({
      robustbase::covMcd(Xsum, alpha = alpha)$cov
      }, error = function(cond) {
        stop("MCD cannot be calculated automatically.")
      })
  }

  # transpose matrix X
  mXsum <- t(Xsum)

  # setup and allocations
  n <- dim(mXsum)[2]
  p <- dim(mXsum)[1]

  # restore target matrix
  mT <- TM

  # apply SVD and transform u to z, save as mX
  mTeigen <- eigen(mT)
  mQ <- mTeigen$vectors
  mL <- diag(mTeigen$values)
  msqL <- diag(sqrt(mTeigen$values))
  misqL <- diag(sqrt(mTeigen$values)^(-1))
  sdtrafo <- mQ %*% msqL

  # initial step transformation and rho selection per neighborhood
  init <- list()
  for( i in 1:N) {
    temp <- initial_step(x = x[[i]],
                        alpha = alpha,
                        maxcond = maxcond,
                        mQ = mQ,
                        misqL = misqL)
    temp$Vselection <- c(temp$initV, temp$setsV)
    init <- append(init, list(temp))
  }


  # get starting h-set combinations
  V_selection_N <- matrix(NA, n_initialhsets, N)
  for( i in 1:N) {
    V_selection_N[, i] <- sample(x = init[[i]]$Vselection,
                                size = n_initialhsets,
                                replace = TRUE)
  }
  V_selection_N <- V_selection_N[!duplicated(V_selection_N),]
  if(dim(V_selection_N)[1] < 6){
    message("Less than 6 starting values - increase n_initialhsets.")
  }

  # C-steps: set up
  obj_fun_values <- c()

  # C-steps: first initial value for comparison
  k <- V_selection_N[1,]
  ret <- cstep(init,
               maxcsteps = maxcsteps,
               which_indices = k,
               weights = weights,
               lambda = lambda,
               mT =mT)

  objret <- objective_init(ret$out, lambda = lambda, weights= weights)
  best6pack <- V_selection_N[1,i]
  obj_fun_values <- rbind(obj_fun_values, ret$obj_value)


  # C-steps: all starting values
  if(dim(V_selection_N)[1] <= 1) stop("Not enough initial starting values.
                                         You might want to increase n_initialhsets.")
  for(i in 2:dim(V_selection_N)[1]) {
    k <- V_selection_N[i,]
    tmp <- cstep(init,
                 maxcsteps = maxcsteps,
                 which_indices = k,
                 weights = weights,
                 lambda = lambda,
                 mT = mT)
    objtmp <- objective_init(tmp$out, lambda = lambda, weights = weights)
    obj_fun_values <- rbind(obj_fun_values, tmp$obj_value)

    # take minimum objective function values
    if (objtmp < objret){
      ret <- tmp
      objret <- objtmp
      best6pack <- k
    }
    else if (objtmp == objret) {
      best6pack <- c(best6pack, k)
    }
  }

  # back-transform
  temp <- back_transformation(best_ret = ret,
                              TM = TM,
                              alpha = alpha,
                              mQ = mQ,
                              msqL = msqL)
  temp <- append(temp, list(weights = weights,
                            lambda = lambda,
                            obj_fun_values = obj_fun_values,
                            best6pack = best6pack) )

  temp <- linear_combination(temp)

  # check output
  class(temp) <- c("ssMRCD", "list")

  return(temp)
}
