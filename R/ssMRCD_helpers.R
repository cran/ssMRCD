# INITIAL STEP (initial_step) ##################################################
# Initial Steps (1 to 3.4) for one Neighborhood
#
# Applies transformations to data , calculates the neighborhood specific
# regularization paremeter rho and the six initial h-sets for the neighborhood.
#
# @param x matrix with observations
# @param alpha robustness factor
# @param maxcond maximal condition number
# @param mQ eigenvectors of target matrix
# @param misqL diagonal matrix with inverse square root of eigenvalues of target matrix
#
# @return Returns one entry of an init-object.
#' @importFrom robustbase Qn
initial_step <- function(x, alpha, maxcond, mQ, misqL){

  # transpose matrix X for more easy coding to mX
  mX <- t(x)

  # setup and allocations
  n <- dim(mX)[2]
  p <- dim(mX)[1]
  h <- as.integer(ceiling(alpha * n))


  # apply SVD and transform u to w, save as mX
  mW <- t(mX) %*% mQ %*% misqL
  mX <- t(mW)
  mT = diag(p)

  # if no initial h observations are given, follow Hubert 2012 to get six good estimates
  hsets.init <- r6pack(x = t(mX), h = h, full.h = FALSE,
                       adjust.eignevalues = FALSE, scaled = FALSE,
                       scalefn = robustbase::Qn)

  # select h observations (r6pack should make h many observations, so just to check for consistency)
  hsets.init <- hsets.init[1:h, ]

  # get consistency factor for neighborhood cov
  scfac <- MDcons_robustbase(p, h/n)

  # calculate rho
  nsets <- ncol(hsets.init)
  rho_out = rho_estimation(mX = mX, mT = mT, h = h, hsets.init = hsets.init,
                           maxcond = maxcond, scfac = scfac)
  rho = rho_out$rho
  setsV = rho_out$setsV
  initV = rho_out$initV

  return(list(nsets = nsets, mX = mX, rho = rho,
              h = h, scfac = scfac, hsets.init = hsets.init ,
              initV = initV, setsV = setsV,
              n = n, p = p))

}


# GET INITIAL ESTIMATES (r6pack) ###############################################
# Function to get the SIX INITIAL STARTING VALUES
# taken from covMrcd

# Robust Distance based observations orderings based on robust "Six pack" (see robustbase, covMrcd)
#
# @param x a numeric matrix constituting the observation matrix
# @param h number of observations used to calculate robust covariance estimations
# @param full.h logical
# @param adjust.eignevalues logical
# @param scaled if the data is already scaled
# @param scalefn scaling function, used to scale data if scaled = FALSE
#
# @return matrix with 6 starting H-Sets for the MRCD
#' @importFrom robustbase doScale colMedians mahalanobisD classPC
r6pack <- function(x, h, full.h, adjust.eignevalues = TRUE,
                   scaled = TRUE, scalefn = robustbase::Qn) {
  # x ... original matrix with observations (not transposed)
  # h ... number of observations used to calculate
  # full.h .... only if adjust.eignevalues == T
  # adjust.eignevalues ... adjustements for eigenvalues (not clear?)
  # scaled ... if False, data will be scaled with median and scalefn
  # scalefn ... scale function, Qn ... robust
  # RETURN: matrix with 6 columns containing 6 initial subsets of size h, according to Hubert 2012
  #         contains indizes

  initset <- function(data, scalefn, P, h) {
    stopifnot(length(d <- dim(data)) == 2, length(h) ==
                1, h >= 1)
    n <- d[1]
    stopifnot(h <= n)
    lambda <- robustbase::doScale(data %*% P, center = stats::median, scale = scalefn)$scale
    sqrtcov <- P %*% (lambda * t(P))
    sqrtinvcov <- P %*% (t(P)/lambda)
    estloc <- robustbase::colMedians(data %*% sqrtinvcov) %*% sqrtcov
    centeredx <- (data - rep(estloc, each = n)) %*% P
    sort.list(robustbase::mahalanobisD(centeredx, FALSE, lambda))[1:h]
  }

  ogkscatter <- function(Y, scalefn, only.P = TRUE) {
    stopifnot(length(p <- ncol(Y)) == 1, p >= 1)
    U <- diag(p)
    for (i in seq_len(p)[-1L]) {
      sYi <- Y[, i]
      ii <- seq_len(i - 1L)
      for (j in ii) {
        sYj <- Y[, j]
        U[i, j] <- (scalefn(sYi + sYj)^2 - scalefn(sYi -
                                                     sYj)^2)/4
      }
      U[ii, i] <- U[i, ii]
    }
    P <- eigen(U, symmetric = TRUE)$vectors
    if (only.P)
      return(P)
    Z <- Y %*% t(P)
    sigz <- apply(Z, 2, scalefn)
    lambda <- diag(sigz^2)
    list(P = P, lambda = lambda)
  }

  stopifnot(length(dx <- dim(x)) == 2)
  n <- dx[1]
  p <- dx[2]
  if (!scaled) {
    x <- robustbase::doScale(x, center = stats::median, scale = scalefn)$x
  }
  nsets <- 6
  hsets <- matrix(integer(), h, nsets)
  y1 <- tanh(x)
  R1 <- stats::cor(y1)
  P <- eigen(R1, symmetric = TRUE)$vectors
  hsets[, 1] <- initset(x, scalefn = scalefn, P = P, h = h)
  R2 <- stats::cor(x, method = "spearman")
  P <- eigen(R2, symmetric = TRUE)$vectors
  hsets[, 2] <- initset(x, scalefn = scalefn, P = P, h = h)
  y3 <- stats::qnorm((apply(x, 2L, rank) - 1/3)/(n + 1/3))
  R3 <- stats::cor(y3, use = "complete.obs")
  P <- eigen(R3, symmetric = TRUE)$vectors
  hsets[, 3] <- initset(x, scalefn = scalefn, P = P, h = h)
  znorm <- sqrt(rowSums(x^2))
  ii <- znorm > .Machine$double.eps
  x.nrmd <- x
  x.nrmd[ii, ] <- x[ii, ]/znorm[ii]
  SCM <- crossprod(x.nrmd)
  P <- eigen(SCM, symmetric = TRUE)$vectors
  hsets[, 4] <- initset(x, scalefn = scalefn, P = P, h = h)
  ind5 <- order(znorm)
  half <- ceiling(n/2)
  Hinit <- ind5[1:half]
  covx <- stats::cov(x[Hinit, , drop = FALSE])
  P <- eigen(covx, symmetric = TRUE)$vectors
  hsets[, 5] <- initset(x, scalefn = scalefn, P = P, h = h)
  P <- ogkscatter(x, scalefn, only.P = TRUE)
  hsets[, 6] <- initset(x, scalefn = scalefn, P = P, h = h)
  if (!adjust.eignevalues)
    return(hsets)
  if (full.h)
    hsetsN <- matrix(integer(), n, nsets)
  for (k in 1:nsets) {
    xk <- x[hsets[, k], , drop = FALSE]
    svd <- robustbase::classPC(xk, signflip = FALSE)
    score <- (x - rep(svd$center, each = n)) %*% svd$loadings
    ord <- order(robustbase::mahalanobisD(score, FALSE, sqrt(abs(svd$eigenvalues))))
    if (full.h)
      hsetsN[, k] <- ord
    else hsets[, k] <- ord[1:h]
  }
  if (full.h)
    hsetsN
  else hsets
}

# RHO ESTIMATION (rho_estimation) ##############################################
# Rho-Estimation (from CovMrcd)
#
# This function is basically taken from rrcov and calculates the best regularization parameter rho.
#
# @param mX matrix with all observations, with observations as columns.
# @param mT target covariance matrix.
# @param hsets.init matrix, six initial h-sets gained through r6pack function, one column per initial h-set.
# @param maxcond maximal condition number for regularized matrix.
# @param scfac consistency factor.
# @param h number of observations included in covariance calculations.
#
# @return Returns a list with the value for rho and h-sets that are used further.

rho_estimation = function(mX, mT, hsets.init, maxcond, scfac, h){

  # setup
  n = dim(mX)[2]
  p = dim(mX)[1]
  nsets <- ncol(hsets.init)
  rho6pack <- condnr <- c()

  # for all different columns of subsets
  for (k in 1:nsets) {
    mXsubset <- mX[, hsets.init[, k]]
    vMusubset <- rowMeans(mXsubset)
    mE <- mXsubset - vMusubset
    mS <- mE %*% t(mE)/(h - 1)

    # depending on T different condition number calculations (faster)
    if (all(mT == diag(p))) {
      veigen <- eigen(scfac * mS)$values
      e1 <- min(veigen)
      ep <- max(veigen)
      fncond <- function(rho) {
        condnr <- (rho + (1 - rho) * ep)/(rho + (1 -
                                                   rho) * e1)
        return(condnr - maxcond)
      }
    }
    else {
      fncond <- function(rho) {
        rcov <- rho * mT + (1 - rho) * scfac * mS
        temp <- eigen(rcov)$values
        condnr <- max(temp)/min(temp)
        return(condnr - maxcond)
      }
    }

    # find zero value depending on rho for difference
    out <- try(stats::uniroot(f = fncond, lower = 1e-05, upper = 0.99),
               silent = TRUE)
    if (!inherits(out, "try-error")) {
      rho6pack[k] <- out$root
    }
    # try grid search if algorithm is not converging
    else {
      grid <- c(1e-06, seq(0.001, 0.99, by = 0.001),
                0.999999)
      if (all(mT == diag(p))) {
        objgrid <- abs(fncond(grid))
        irho <- min(grid[objgrid == min(objgrid)])
      }
      else {
        objgrid <- abs(apply(as.matrix(grid), 1, "fncond"))
        irho <- min(grid[objgrid == min(objgrid)])
      }
      rho6pack[k] <- irho
    }
  }

  # define cutoff, rho and regular sets of initial values
  cutoffrho <- max(c(0.1, stats::median(rho6pack)))
  rho <- max(rho6pack[rho6pack <= cutoffrho])
  Vselection <- seq(1, nsets)
  Vselection[rho6pack > cutoffrho] = NA
  if (sum(!is.na(Vselection)) == 0) {
    stop("None of the initial subsets is well-conditioned")
  }
  initV <- min(Vselection, na.rm = TRUE)
  setsV <- Vselection[!is.na(Vselection)]
  setsV <- setsV[-1]

  return( list(rho = rho, initV = initV, setsV = setsV))
}

# C-STEPS (cstep) ##############################################################
# C-Steps for ssMRCD
#
# Computes the C-Steps based on the theorem.
#
# @param init list object including all informations like regularization parameter, target matrix, observation values etc.
# @param maxcsteps maximum of c-steps before stopping calculations
# @param which_indices vector, which of the six initial h-sets should be used per neighborhood
# @param lambda scalar, smoothing factor
# @param weights weighting matrix
# @param mT target matrix
#
# @return Returns the optimal h-sets for the initial combination of h-sets.
cstep = function(init, maxcsteps, which_indices, lambda, weights, mT){

  N = length(init)
  p = init[[1]]$p # should be equal for all neighborhoods

  obj_value = rep(NA, maxcsteps)

  # for each neighborhood initial set up of first calculations for cov-matrix
  for( i in 1:N ){

    # get h set for neighborhood i
    init[[i]]$index = init[[i]]$hsets.init[, which_indices[i]]

    # select corresponding observations from X and calculate mean and regularized covariance
    XX <- init[[i]]$mX[, init[[i]]$index]
    init[[i]]$vMu = rowMeans(XX)
    ret = RCOV(XX = XX, vMu = init[[i]]$vMu, rho = init[[i]]$rho, mT = mT,
               scfac = init[[i]]$scfac)
    init[[i]]$rho = ret$rho
    init[[i]]$mS = ret$rcov
  }
  obj_value[1] =  objective_init(init_object = init, lambda = lambda, weights= weights)

  # calculate observations with minimal distances
  for(i in 1:N){
    init[[i]]$vdst = dist_cstep(init = init, i = i, lambda = lambda, weights = weights)
    init[[i]]$index = sort(sort.int( init[[i]]$vdst, index.return = TRUE)$ix[1:init[[i]]$h])
  }

  # start iteration
  iter = 1
  decreasing_obj = NULL

  while (iter < maxcsteps) {

    # calculate mu and Cov
    for (i in  1:N) {
      # same procedure with new h set "index"
      XX <- init[[i]]$mX[, init[[i]]$index]
      init[[i]]$vMu <- rowMeans(XX)
      init[[i]]$ret <- RCOV(XX = XX, vMu = init[[i]]$vMu, rho = init[[i]]$rho, mT = mT,
                            scfac = init[[i]]$scfac)
      init[[i]]$mS <- init[[i]]$ret$rcov
    }

    # select new minimal observations
    for (i in 1:N){
      init[[i]]$vdst <- dist_cstep(init = init, i = i, lambda = lambda, weights = weights)
      nndex <- sort(sort.int(init[[i]]$vdst, index.return = TRUE)$ix[1:init[[i]]$h])

      # if observations stay the same, break and use last calculated values
      if (all(nndex ==  init[[i]]$index)) {
        init[[i]]$stop = 1
      }
      else {
        init[[i]]$stop = 0
      }
      init[[i]]$index <- nndex
    }

    # save information
    obj_value[iter + 1] =  objective_init(init, lambda = lambda, weights = weights)

    # check break condition
    break_sum = 0
    for(i in 1:N) break_sum = sum(break_sum, init[[i]]$stop)
    if(break_sum == N){
      break
    }
    iter <- iter + 1
  }


  return(list(numit = iter, out = init, obj_value = obj_value))
}

# DISTANCE FOR C-STEP (dist_cstep) #############################################
# Distances calculated for the C-Step
#
# Computes the distances used to get the next observations for the ssMRCD.
#
# @param init list, contains all relevant information per neighborhood \eqn{a_i} in i-th entry.
# @param i index, for observations of which neighborhood the distances should be calculated
# @param lambda scalar, smoothing parameter
# @param weights matrix with properties of weighting matrix
#                (symmetric, diagonal equals zero, sum of the rows is equal to one)
#
# @return a vector containing the distances for each observation in neighborhood i
dist_cstep <- function(init, i, lambda, weights){
  # returns vector of distances for each observation in mX

  N = length(init)
  p = init[[i]]$p

  x_centered = init[[i]]$mX - init[[i]]$vMu
  Cov_neighb = matrix(0, p, p)
  for(j in 1:N) {
    Cov_neighb = Cov_neighb + weights[i,j]*init[[j]]$mS  #diagonal entries are 0!
  }
  Cov_matrix = (1-lambda) * init[[i]]$mS + lambda * Cov_neighb
  return(diag(t(x_centered) %*% ( chol2inv(chol(Cov_matrix))) %*% x_centered))
}


# OBJECTIVE FUNCTION VALUE INTERNALLY (objective_init) #########################
# Objective function for init-Object
#
# Calculates the objective function based on the matrices stored in the init object.
#
# @param init_object list object
# @param lambda scalar for spatial smoothing
# @param weights matrix with weights
#
# @return Returns value of objective function.
objective_init = function(init_object, lambda , weights){

  N = length(init_object)

  matrix_list = list()
  for(i in 1:N){
    rcov = init_object[[i]]$mS
    matrix_list = append(matrix_list, list(rcov))
  }

  return(objective_matrix(matrix_list = matrix_list, lambda = lambda, weights = weights))
}


# OBJECTIVE FUNCTION VALUE FOR MATRIX LIST (objective_matrix) ##################
#' Calculation of Objective Function
#'
#' Calculation of the value of the objective function for the \code{\link[ssMRCD]{ssMRCD}} for a given list of matrices,
#' lambda and a weighting matrix according to formula (3) in Puchhammer and Filzmoser (2023).
#'
#' @param matrix_list a list of matrices \eqn{K_i}
#' @param lambda scalar smoothing parameter
#' @param weights matrix of weights
#'
#' @return Returns the value of the objective function using matrices \eqn{K_i}.
#' @references Puchhammer P. and Filzmoser P. (2023): Spatially smoothed robust covariance estimation for local outlier detection. \doi{10.48550/arXiv.2305.05371}
#' @export
#'
#' @examples
#' # construct matrices
#' k1 = matrix(c(1,2,3,4), nrow = 2)
#' k2 = matrix(c(1,3,5,7), nrow = 2)
#'
#' # construct weighting matrix
#' W = matrix(c(0, 1, 1, 0), nrow = 2)
#'
#' objective_matrix(list(k1, k2), 0.5, W)
objective_matrix = function(matrix_list, lambda, weights){

  N = length(matrix_list)
  p = dim(matrix_list[[1]])[1]

  f = 0
  for(i in 1:N){
    P = matrix(0, p, p)
    for(j in 1:N){
      P = P + matrix_list[[j]] * weights[i,j]
    }
    f = f + det( (1-lambda) * matrix_list[[i]] + lambda * P)
  }

  return(f)
}



# BACK TRANSFORMATION (back_transformation) ####################################
# Restructure of Results
#
# Restructures and backtransforms the results by the best c-step and initial h-set combination.
#
# @param best_ret list object, produced by cstep function
# @param TM matrix, target matrix
# @param alpha robustness factor
# @param mQ eigenvectors of target matrix
# @param msqL square root of eigenvalues of target matrix
#
# @return Returns back-transformed observations and covariance matrices.
back_transformation = function(best_ret, TM, alpha, mQ, msqL){

  N = length(best_ret$out)
  p = length(best_ret$out[[1]]$vMu)
  MRCDcov = list()
  MRCDicov = list()
  MRCDmu = list()
  mX = list()
  mt = diag(1, p)
  dist_tmp = list()
  rho = rep(NA, N)
  h = rep(NA, N)
  hset = list()
  c_alpha = rep(NA, N)


  for(i in 1:N){
    c_alpha[i] = best_ret$out[[i]]$scfac
    h[i] = best_ret$out[[i]]$h
    rho[i] = best_ret$out[[i]]$rho

    hindex = best_ret$out[[i]]$index
    hset = append(hset, list(hindex))

    mx = best_ret$out[[i]]$mX

    mE = mx[, hindex] - best_ret$out[[i]]$vMu
    weightedScov <- mE %*% t(mE)/(h[i] - 1)

    # back transformation to original space
    mu = mQ %*% msqL %*% best_ret$out[[i]]$vMu
    mx =  t(t(mx) %*% msqL %*% t(mQ))
    cov = (1 - rho[i]) * c_alpha[i] * weightedScov + rho[i] * mt
    cov =  mQ %*% msqL %*% cov %*% msqL %*% t(mQ)
    icov <- chol2inv(chol(cov))

    dist_tmp <- stats::mahalanobis(t(mx), center = mu, cov = icov,
                                   inverted = TRUE)


    MRCDmu = append(MRCDmu, list(mu))
    MRCDcov = append(MRCDcov, list(cov))
    MRCDicov = append(MRCDicov, list(icov))
    mX = append(mX, list(t(mx)))
    dist = append(dist, dist_tmp)
  }

  out = list(MRCDcov = MRCDcov, MRCDicov = MRCDicov, MRCDmu = MRCDmu, mX = mX,
             N = N, mT = TM, rho = rho, alpha = alpha, h = h, numiter = best_ret$numit, hset = hset,
             c_alpha = c_alpha)
  return(out)

}


# HELPER FUNCTIONS #############################################################
# (RCOV, MDcons_robustbase, monotonic_decreasing, linear_combination)

# Regularized Covariance Matrix
#
# Computes the regularized covariance matrix for a given subset of observations,
# the target matrix, a vector of means, the regularity parameter and the consistency factor for
# consistency of the estimate. See also
# \href{https://link.springer.com/article/10.1007/s11222-019-09869-x}{Boudt et al. (2020)}.
#
# @param XX matrix of selected observations
# @param vMu vector of means, transposed
# @param rho regularity parameter
# @param mT target matrix
# @param scfac scalar, consistency factor

RCOV <- function(XX, vMu, rho, mT, scfac) {
  # XX ... h observations, transposed
  # vMu ... vector of means, transposed
  # rho ... given
  # mT ... target matrix T
  # scfac ... consistency factor (c_alpha)
  # target ... not necessary
  # invert ... if inverse should also be calculated
  # RETURNS: list with rho, Ki, inv(Ki), T, Si

  mE <- XX - vMu
  n <- dim(mE)[2]
  p <- dim(mE)[1]
  mS <- mE %*% t(mE)/n
  rcov <- rho * mT + (1 - rho) * scfac * mS

  return(list(rho = rho, mT = mT, cov = mS, rcov = rcov))
}



MDcons_robustbase = function (p, alpha) {
  qalpha <- stats::qchisq(alpha, p)
  caI <- stats::pgamma(qalpha/2, p/2 + 1)/alpha
  1/caI
}


monotonic_decreasing = function(x){

  x = stats::na.omit(x)
  n = length(x)

  monotonic = TRUE
  for(i in 1:(n-1)){
    if(x[i] < x[i+1]){
      monotonic = FALSE
      break
    }
  }

  return(monotonic)
}


linear_combination = function(x){

  N = x$N
  p = dim(x$MRCDcov[[1]])
  weight = x$weights
  covs = list()
  icovs = list()
  for(i in 1:N){
    sums = matrix(0, p, p)
    for(j in 1:N){
      sums = sums + weight[i, j]* x$MRCDcov[[j]]
    }
    cov = x$MRCDcov[[i]] * (1-x$lambda) + x$lambda * sums
    covs = c(covs, list(cov))
    icovs = c(icovs, list(solve(cov)))
  }
  x$Kcov = x$MRCDcov

  x$MRCDcov = covs
  x$MRCDicov = icovs

  return(x)
}
