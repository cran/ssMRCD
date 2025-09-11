# MG-GMM: MAIN AND HELPER FUNCTIONS

# Regularization according to MRCD approach with either fix rho or if a
# maximal condition number is given then the rho providing the condition is calculated.
# The function needs to be applied to each Sigma separately.
regularize = function(Sigma, Sigma_reg,
                      rho = NULL, rho_fixed = FALSE,
                      maxcond = 100){

  fncond <- function(rho) {
    rcov <- rho * diag(diag(Sigma_reg)) + (1 - rho) * Sigma
    temp <- eigen(rcov, symmetric = TRUE)$values
    condnr <- abs(max(temp)/min(temp))
    return(condnr - maxcond)
  }

  p = dim(Sigma)

  if(is.null(rho) & rho_fixed)  {
    stop("You must either provide a fixed 'rho', or allow the function to estimate it by setting 'rho_fixed = FALSE'.")
  }
  if(cond(Sigma_reg) > maxcond & !rho_fixed)  {
    stop("The condition number of Sigma_reg is bigger than the selected maxcond.")
  }

  if(rho_fixed & !is.null(rho)) rho_out = rho
  if(is.null(rho) | !rho_fixed){
    rho_out = 1e-6

    if(fncond(rho_out) > 0){
      out <- try(stats::uniroot(f = fncond, lower = 1e-6, upper = 0.99), silent = TRUE)
      if (!methods::is(out, "try-error")) {
        rho_out <- out$root
      } else {
        print(out[1])
        stop(paste("No rho found to achieve condition number", maxcond, "in regularization."))
      }
    }
  }
  Sigma_reg = rho_out * Sigma_reg + (1 - rho_out) * Sigma

  return(list(Sigma = Sigma_reg, rho = rho_out))
}



# Calculates the EM-step as described including the regularization of the covariance matrices.
em_step = function(X, mu, Sigma, Sigmai, pi_groups, W, groups, alpha, Sigma_reg, rho){

  # X: data matrix.
  # mu: list of mean vectors.
  # Sigma: list of covariances.
  # Sigmai: list of inverse covariances.
  # pi_groups: matrix of mixture probabilities.
  # W: binary matrix.
  # groups: vector of grouping assignment.
  # alpha: scalar between 0.5 and 1.
  # Sigma_reg: list of regularization matrices.
  # rho: vector of regularization factors

  n = dim(X)[1]
  p = dim(X)[2]
  N = length(Sigma)

  # calculate class probabilities t_k
  probs = probabs(X = X, mu = mu, Sigma = Sigma, pi_groups = pi_groups, W = W, groups = groups)
  Nk = colSums(probs)

  # initialize imputed data matrix
  Ximp_out = matrix(0, nrow = n, ncol = p)

  # repeat the parameter estimation for all mixture components
  for(k in 1:N){
    bias <- matrix(0, p, p)
    Ximp <- X
    Ximp_p = matrix(NA, n, p)

    # iterate over all unique row structures in W
    W_pattern = unique(W)
    wn = dim(W_pattern)[1]

    for(l in 1:wn){
      # get missingness pattern
      obs = which(W_pattern[l, ] == 1)
      mis = which(W_pattern[l, ] == 0)

      # check if there are missing observations and subset matrix (speed)
      if (length(mis) != 0){
        Sigmai_tmp = chol2inv(chol(Sigmai[[k]][mis, mis, drop = F]))
      }

      # get all observation indices that have the same missingness pattern
      ind_w =  which(colSums(abs(t(W) - W_pattern[l, ])) == 0)

      # iterate over all observation with the same pattern
      for(i in ind_w){
        ximp = X[i, ]
        g = groups[i]

        # if no observed variables insert mean
        if(length(obs) == 0){
          ximp = mu[[k]]
          bias = bias + probs[i, k] * Sigma[[k]]
        }

        # if some variables are observed use conditional expectation and variance bias
        if(length(mis) != 0 & length(obs) != 0){
          ximp[mis] = mu[[k]][mis] - Sigmai_tmp %*% Sigmai[[k]][mis, obs, drop = F] %*% (c(ximp[obs])-c(mu[[k]][obs]))
          bias[mis, mis] =  bias[mis, mis] + probs[i, k] * Sigmai_tmp
        }

        # collect imputed values and weighted imputed values
        Ximp[i, ] = ximp
        Ximp_p[i, ] = probs[i, k] * Ximp[i, ]
      }
    }

    # update weighted imputed value over all mixture components
    Ximp_out = Ximp_out + Ximp_p

    # update mean estimation
    mu[[k]] = colSums(Ximp_p)/Nk[k]

    # update covariance estimation including regularization and bias
    weighted_X <- sqrt(probs[, k]) * sweep(Ximp, 2, mu[[k]])  # center Ximp by mu[[k]]
    tmp <- t(weighted_X) %*% weighted_X
    Sigma[[k]] = (1-rho[k])*(tmp + bias)/Nk[k] + rho[k]*Sigma_reg[[k]]
  }

  # calculate probabilities per group
  pi_groups = pis(probs = probs, groups = groups, alpha = alpha)

  return(list(mu = mu,
              Sigma = Sigma,
              pi_groups = pi_groups,
              probs = probs,
              Ximp = Ximp_out))
}


# This function iterates the W-step and the EM-step. It starts with initial values provided and returns
# the estimates parameters.
iter_cstep = function(X,
                      mu, Sigma, Sigma_reg, pi_groups, W,
                      groups,
                      Q,
                      alpha = 0.5, hperc = 0.75,  # number of unflagged cells
                      nsteps = 100, crit = 1e-04, silent = TRUE,
                      rho = NULL,
                      plot_conv = FALSE){

  # X: data matrix.
  # mu: list of initial mu values.
  # Sigma: list of initial covariance matrices.
  # Sigma_reg: list of regularization matrices.
  # pi_groups: matrix of initial values for mixture probabilities.
  # W: initial binary matrix.
  # groups: vector of group assignements.
  # Q: penalty weights for more efficiency based on initialization.
  # alpha: scalar between 0.5 and 1.
  # hperc: minimal percentage of observation to use per variable per group.
  # nsteps: maximal number of iterative steps.
  # crit: tolerance for convergence (connected to change in Sigma).
  # silent: boolean, whether progress should be reported.
  # rho: vector of fixed regularization factor
  # maxcond: maximal condition number of covariance matrices.
  # plot_conv: boolean, whether values of objective function should be plotted over iteration steps


  # initialize variables
  X = as.matrix(X)
  Ximp = X

  N = length(Sigma)
  p = dim(X)[2]
  n = dim(X)[1]

  Sigmai = lapply(Sigma, function(x) chol2inv(chol(x)))
  probs = matrix(NA, nrow = n, ncol = p)
  prevSigma = Sigma

  nn = unname(table(groups))
  h = ceiling(hperc * nn)

  # set variables for loop iterations
  convcrit <- 1
  nusteps <- 1

  # initialize objective function and get value for the initial values
  objvals <- rep(NA, nsteps)
  objvals[nusteps] = objective(X, mu, Sigma, pi_groups, W, groups, Q)

  if(!silent)  cat(paste0("Objective at step ", nusteps, " = ", round(objvals[nusteps], 5), "\n"))

  # calculate as long as there is no convergence and the maximal number of steps are not reached
  while (convcrit > crit && nusteps < nsteps) {

    # update number of steps
    nusteps = nusteps + 1

    # W-step (cpp)
    W = w_step(X, mu, Sigma, pi_groups, W, groups, Q, h)

    # EM-step
    EM = em_step(X = X, mu = mu, Sigma = Sigma, Sigmai = Sigmai, pi_groups = pi_groups, W = W,
                 groups = groups, alpha = alpha, Sigma_reg = Sigma_reg, rho = rho)

    # get new estimates
    Sigma = EM$Sigma
    Sigmai = lapply(EM$Sigma, function(x) chol2inv(chol(x)))
    pi_groups = EM$pi_groups
    probs = EM$probs
    mu = EM$mu
    Ximp = EM$Ximp

    # check convergence criteria
    convcrit = max(sapply(1:N, function(x) abs(Sigma[[x]] - prevSigma[[x]])))
    prevSigma = Sigma

    # calculate objective function value
    objvals[nusteps] = objective(X, mu, Sigma, pi_groups, W, groups, Q)
    if(!silent) cat(paste0("Objective at step ", nusteps, " = ", round(objvals[nusteps], 4), "\n"))

    if(objvals[nusteps] > objvals[nusteps-1] & abs(objvals[nusteps] - objvals[nusteps-1])> 0.0001) {
      message(paste0("Iteration finished due to computational inaccuracies (increase of objective function in last step of size: ",
                     abs(objvals[nusteps -1] - objvals[nusteps]), "). Check convergence plot."))
      break
    }

    # check convergence criteria and break loop
    if(nusteps == nsteps) {
      message(paste0("Iteration finished - maximal step number reached. Increase number of maximal steps."))
      break
    }
    if(convcrit <= crit) {
      message(paste0("Iteration finished - algorithm converged."))
      break
    }
  }

  # plot objective function value over iterations
  if(plot_conv | (objvals[nusteps] > objvals[nusteps-1] & abs(objvals[nusteps] - objvals[nusteps-1])> 0.0001)) {
    plot(objvals,
         type = "l",
         xlab = "Iteration",
         ylab = "Objective Function Value"
    )
  }

  # calculate most likely class
  class = apply(X = probs, MARGIN = 1, FUN = which.max)

  return(list(X = X,
              Ximp = Ximp,
              groups = groups,
              class = class,
              mu = mu,
              Sigma = Sigma,
              Sigmai = Sigmai,
              probs = probs,
              pi_groups = pi_groups,
              W = W,
              Q = Q,
              Sigma_reg = Sigma_reg,
              rho = rho,
              alpha = alpha,
              hperc = hperc,
              nsteps = nusteps,
              objvals = objvals))
}


######################################################################
# calculates the condition number of a covariance matrix S
cond = function(S){

  # S: symmtric matrix

  p = dim(S)[1]
  eig = eigen(S, symmetric = TRUE)
  return(eig$values[1]/eig$values[p])
}



# This function estimates the initial estimates for each group separately based on the
# description in the paper. Regularization with a target matrix is used.

#' @importFrom cellWise estLocScale
#' @importFrom cellWise DDC
#' @import stats
initial_est = function(X, hperc = 0.75, maxcond = 100){

  # X: data matrix
  # hperc: scalar, percentage of cells per variable used
  # mxcond: scalar, condition number if possible

  # helper function to flag outlying rows
  RR <- function(Z, cov, b = 2, quant = 0.99) {
    # helper function to remove cases in step 4
    p <- dim(Z)[2]
    MDs <- stats::mahalanobis(pmin(pmax(Z, -b), b), rep(0, p), cov)
    rowinds <- which(MDs/stats::median(MDs) * stats::qchisq(0.5, p) > stats::qchisq(quant, p))
    return(rowinds)
  }

  n <- dim(X)[1]
  p <- dim(X)[2]

  # get scales of data and calculate target matrix
  locsca <- cellWise::estLocScale(X)
  rscales <- locsca$scale
  Target = diag(rscales^2)

  # get the maximal condition number for regularization based on target matrix
  maxcond = max(maxcond, cond(Target)*1.1)

  # scale data, use identity as target matrix
  Xs <- scale(X, center = locsca$loc, scale = rscales)

  # options for DDC
  maxCol = 1-hperc
  tolProbCell = 0.9
  lmin = 1e-04

  # run DDC
  DDCout <- cellWise::DDC(Xs, list(silent = TRUE, tolProbCell = 0.9, standType = "wrap"))

  # W matrix for flagged cells to be replaced with imputed values (1: replace, 0:keep)
  Wna <- matrix(0, n, p)
  Wna[DDCout$indcells] <- 1

  # variables where more than 25% of the cells are flagged and replaced
  overflag <- which(colSums(Wna) > maxCol * n)
  if (length(overflag) > 0) {
    for (i in seq_len(length(overflag))) {
      ind <- overflag[i]  # which variable

      # replace only the 25% most outlying cells
      ord <- order(abs(DDCout$stdResid[, ind]), decreasing = TRUE)
      replacement <- rep(0, n)
      replacement[ord[seq_len(floor(maxCol * n))]] <- 1
      Wna[, ind] <- replacement
    }

    # impute outlying cells
    DDCout$indcells <- which(Wna == 1)
    DDCout$Ximp <- X
    DDCout$Ximp[DDCout$indcells] <- DDCout$Xest[DDCout$indcells]
  }

  # re-scale imputed and original data according to DDC loc and scale
  locScale <- list(loc = DDCout$locX, scale = DDCout$scaleX)
  Z <- scale(X, locScale$loc, locScale$scale)
  Zimp <- scale(DDCout$Ximp, locScale$loc, locScale$scale)
  Zimporig <- Zimp
  Zorig <- Z

  # remove rows that are flagged as outliers by DDC for covariance estimation
  if (length(DDCout$indrows) > 0) {
    Z <- Z[-DDCout$indrows, ]
    Zimp <- Zimp[-DDCout$indrows, ]
  }


  # STEP 3 on imputed values

  # eigen decomposition and keep eigenvectors with high enough variance
  eig = eigen(cov(Zimp), symmetric = TRUE)
  keep = which(eig$values >= lmin) # only project on eigenvectors with enough variance
  eigenvectors <- eig$vectors[, keep]

  # project data on selected eigenvectors and get scores
  Zimp_proj <- Zimp %*% eigenvectors

  # estimate location and scale (at least lmin)
  locscale_proj <- cellWise::estLocScale(Zimp_proj, silent = TRUE)
  locscale_proj$scale = pmax(locscale_proj$scale, lmin)

  # wrap imputed data scores
  Zimp_proj_w <- cellWise::wrap(X = Zimp_proj,
                                locX = locscale_proj$loc,
                                scaleX=locscale_proj$scale)$Xw


  # STEP 4

  # get covariance in wrapped score space and transfer to original Z imputed space
  cov <- eigenvectors %*% cov(Zimp_proj_w) %*% t(eigenvectors)
  # regularize covariance in z-space with transformed target matrix
  cov = regularize(Sigma = cov,
                   Sigma_reg = diag(locScale$scale^2),
                   maxcond = maxcond)$Sigma
  # transfer covariance in original data space (scaled)
  cov <- t(stats::cov2cor(cov) * locScale$scale) * locScale$scale

  Zimp <- Zimporig
  Z <- Zorig

  # remove cases based on cov estimates
  rowinds <- RR(Z, cov, 2, 0.99)
  if (length(rowinds) > 0) {
    Z <- Z[-rowinds, ]
    Zimp <- Zimp[-rowinds, ]
  }

  # STEP 5: redo everything with DDC output with removed rows

  # project on eigenvectors
  eig = eigen(cov(Zimp), symmetric = TRUE)
  keep = which(eig$values >= lmin)
  eigenvectors <- eig$vectors[, keep]
  Zimp_proj <- Zimp %*% eigenvectors

  # scale and wrap the scores
  locscale_proj <- estLocScale(Zimp_proj, silent = TRUE)
  locscale_proj$scale = pmax(locscale_proj$scale, lmin)
  Zimp_proj_w <- cellWise::wrap(X = Zimp_proj,
                                locX = locscale_proj$loc,
                                scaleX = locscale_proj$scale)$Xw

  # calculate covariance and transfer it to original space
  Zcov.raw <- eigenvectors %*% cov(Zimp_proj_w) %*% t(eigenvectors)

  # regularize covariances
  Zcov <- stats::cov2cor(Zcov.raw)
  Zcov = regularize(Sigma = Zcov,
                    Sigma_reg = diag(locScale$scale^2),
                    maxcond = maxcond)$Sigma

  # transfer to original space
  cov <- t(t(Zcov) * locScale$scale) * locScale$scale
  S <- diag(rscales) %*% stats::cov2cor(cov) %*% diag(rscales)

  # return
  return(list(mu = locScale$loc * rscales + locsca$loc, Sigma = S, Target = Target))
}



# calculates initial probabilities given initial parameter estimates
initial_probs = function(X, Sigma_reg, mu, groups, alpha = 0.5){

  # X: data matrix
  # Sigma_reg: list of regularized covariance matrices
  # mu: list of location vectors
  # groups: vectors of grouping assignement
  # alpha: scalar, flexibility parameter

  n = dim(X)[1]
  p = dim(X)[2]
  N = length(Sigma_reg)

  # get initialization for mixture probabilities pi
  tt = matrix(1/(N-1), N, N)
  diag(tt) = 0
  pi_groups = diag(alpha, N) + (1-alpha)*tt

  # get probabilities t_k for each observation
  probs = probabs(X = X,
                  mu = mu,
                  Sigma = Sigma_reg,
                  pi_groups = pi_groups,
                  W = matrix(1, n, p),
                  groups = groups)

  # return initial values
  return(list(probs = probs,
              pi_groups = pi_groups,
              Sigma = Sigma_reg,
              W = matrix(1, n, p)))
}



# estimate penalty weights based on initialization parameters
initial_Q = function(Sigma_init, probs, quantq = 0.99){

  # Sigma_init: list of matrices
  # probs: matrix of probabilities for each observation
  # quantq: quantile to use for chi-square distribution

  N = length(Sigma_init)
  p = dim(Sigma_init[[1]])[1]
  n = dim(probs)[1]

  # initialize
  qhelp = matrix(NA, nrow = N, ncol = p)
  Q = matrix(NA, nrow = n, ncol = p)

  # for all variables, classes and observations
  for(j in 1:p){
    for(k in 1:N){
      qhelp[k, j] = stats::qchisq(p = quantq, df = 1) + log(2*pi) - log(chol2inv(chol(Sigma_init[[k]]))[j,j])
    }
  }
  for(i in 1:n){
    for(j in 1:p){
      Q[i, j] = sum(probs[i, ] * qhelp[ , j])
    }
  }

  # return matrix Q
  return(Q)
}



# general initialization function
initialization = function(X, groups, alpha = 0.5, maxcond = 100){

  # X: data matrix
  # groups: vector of group assignements
  # alpha: scalar, flexibility parameter
  # maxcond: maximal condition number suggested

  # get N and p
  N = length(unique(groups))
  p = dim(X)[2]

  # initialize variables
  Sigma = list()
  Sigma_reg = list()
  mu = list()
  rho = rep(NA, N)
  maxcond_all = 0

  # repeat for each groups
  for(g in 1:N){

    # get initial mean, covariance and target matrix estimate for the procedure
    tmp = initial_est(X[groups == g, ], hperc = 0.75, maxcond = maxcond)

    # get the real condition number used and maximimum over all groups
    condmax = max(maxcond, cond(tmp$Target)*1.1)
    maxcond_all = max(maxcond_all, condmax)

    # regularize and get the regularization factor
    tmp2 = regularize(tmp$Sigma, Sigma_reg = tmp$Target, maxcond = condmax)

    # save results
    rho[g] = tmp2$rho
    Sigma_reg = c(Sigma_reg, list(tmp$Target))
    Sigma = c(Sigma, list(tmp2$Sigma))
    mu = c(mu, list(tmp$mu))
  }

  # calculate initial values for mixing probabilities pi and t_k
  probs_pi = initial_probs(X = X,
                           Sigma_reg = Sigma_reg,
                           mu = mu,
                           groups = groups,
                           alpha = alpha)

  # calculate penalty weights
  Q = initial_Q(Sigma_init = probs_pi$Sigma, probs = probs_pi$probs)

  return(list(Sigma = probs_pi$Sigma,
              Sigma_reg = Sigma_reg,
              mu = mu,
              pi_groups = probs_pi$pi_groups,
              probs = probs_pi$probs,
              Q = Q,
              W = probs_pi$W,
              rho = rho,
              maxcond = maxcond_all))
}




####################################################################################
#' Calculation of the cellwise robust multi-group Gaussian mixture model
#'
#' Performs robust estimation of multivariate location and scatter within predefined groups
#' using an iiterative EM-based algorithm.
#'
#' @param X A numeric matrix or data frame with observations in rows and variables in columns.
#' @param groups A vector indicating the group membership for each observation (length must match `nrow(X)`).
#' @param alpha A non-negative numeric value between `0.5` and `1` controlling the flexibility degree. Default is `0.5`.
#' @param hperc A numeric value in `[0,1]` controlling robustness of the estimation. Default is `0.75`.
#' @param nsteps Number of main iteration steps in the algorithm. Default is `100`.
#' @param crit Convergence criterion for iterative updates. Default is `1e-4`.
#' @param silent Logical; if `TRUE`, suppresses progress output. Default is `FALSE`.
#' @param maxcond Maximum allowed condition number for covariance matrices. Default is `100`.
#'
#' @return A list containing: \describe{
#'   \item{\code{X}}{The original data matrix.}
#'   \item{\code{Ximp}}{The imputed and/or scaled data matrix.}
#'   \item{\code{groups}}{Vector specifying group assignments from the input.}
#'   \item{\code{class}}{Vector indicating the most likely group membership for each observation, as inferred by the model.}
#'   \item{\code{mu}}{A list of estimated location (mean) vectors for each group.}
#'   \item{\code{Sigma}}{A list of estimated covariance matrices for each group.}
#'   \item{\code{Sigmai}}{A list of estimated inverse covariance matrices for each group.}
#'   \item{\code{probs}}{A matrix of class probabilities for each observation (rows = observations, columns = groups).}
#'   \item{\code{pi_groups}}{A matrix of estimated mixture probabilities, where rows correspond to groups and columns to distributions.}
#'   \item{\code{W}}{A binary matrix indicating outlying cells (0 = outlier, 1 = no outlier).}
#'   \item{\code{Q}}{A matrix of penalty weights.}
#'   \item{\code{Sigma_reg}}{A list of estimated target (regularization) matrices.}
#'   \item{\code{rho}}{A vector of regularization factors used in the estimation.}
#'   \item{\code{alpha}}{Flexibility parameter, as provided in the function input.}
#'   \item{\code{hperc}}{A matrix or vector indicating the percentage of outlying cells per variable and group, based on input.}
#'   \item{\code{nsteps}}{The number of iteration steps taken until convergence.}
#'   \item{\code{objvals}}{The values of the objective function across the iteration steps.}
#' }
#'
#'
#' @seealso \code{\link[ssMRCD]{residuals_mggmm}}
#'
#' @references Puchhammer, P., Wilms, I., & Filzmoser, P. (2025). A smooth multi-group Gaussian Mixture Model for cellwise robust covariance estimation. \emph{ArXiv preprint} \doi{10.48550/arXiv.2504.02547}.
#'
#'
#' @importFrom cellWise estLocScale
#' @export
#'
#' @examples
#' data("weatherAUT2021")
#' cut_lon = c(min(weatherAUT2021$lon)-0.2, 12, 16, max(weatherAUT2021$lon) + 0.2)
#' cut_lat = c(min(weatherAUT2021$lat)-0.2, 48, max(weatherAUT2021$lat) + 0.2)
#' groups = groups_gridbased(weatherAUT2021$lon, weatherAUT2021$lat, cut_lon, cut_lat)
#' N = length(unique(groups))
#' model = cellMGGMM(X = weatherAUT2021[, c("p", "s", "vv", "t", "rsum", "rel")],
#'                  groups = groups,
#'                  alpha = 0.5)

cellMGGMM = function(X,
                    groups,
                    alpha = 0.5,
                    hperc = 0.75,
                    nsteps = 100,
                    crit = 1e-04,
                    silent = TRUE,
                    maxcond = 100)  {


  # convert X to matrix if data.frame
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) stop("'X' must be a matrix or data frame that can be converted to a numeric matrix.")
  storage.mode(X) <- "numeric"

  # get variables
  n <- nrow(X)
  p <- ncol(X)
  N <- length(unique(groups))

  # check groups
  if (length(groups) != n) stop("Length of 'groups' must equal number of rows in 'X'.")

  # check scalars
  if(alpha < 0.5 | alpha > 1) stop("The flexibility parameter alpha needs to be within [0.5, 1].")
  if(hperc < 0.5 | hperc > 1) stop("The percentage of cells used per variable hperc needs to be within [0.5, 1].")

  # collect names
  cnames = colnames(X)
  rnames = rownames(X)

  # collect group names
  groups_isnum = is.numeric(groups)

  groups_char = as.character(factor(groups))
  groups = as.numeric(factor(groups))
  gnames = unique(groups_char)[sort.int(unique(groups), index.return = TRUE)$ix]

  # scale everything
  locsca <- cellWise::estLocScale(X, type = "mcd")
  Xs <- scale(X, center = locsca$loc, scale = locsca$scale)

  # initialization
  initEst = initialization(X = Xs,
                        groups = groups,
                        alpha = alpha,
                        maxcond = maxcond)


  rho = initEst$rho
  Sigma_reg = initEst$Sigma_reg

  # run model
  out = iter_cstep (X = Xs,
                    mu = initEst$mu,
                    Sigma = initEst$Sigma,
                    Sigma_reg = Sigma_reg,
                    pi_groups = initEst$pi_groups,
                    W = matrix(1, n, p),
                    groups = groups,
                    Q = initEst$Q,
                    alpha = alpha,
                    hperc = hperc,
                    nsteps = nsteps,
                    crit = crit,
                    silent = silent,
                    rho = rho)

  # rescale to original X space (should be included in the main function)
  rscales = locsca$scale
  out$Ximp = out$Ximp %*% diag(rscales) + matrix(locsca$loc, ncol = p, nrow = n,  byrow = T)
  out$mu <- lapply(out$mu, function(x) locsca$loc + x * rscales)
  out$Sigma = lapply(out$Sigma, function(x)  diag(rscales) %*% x %*% diag(rscales))

  # rename variables
  if(!is.null(cnames)){
    out$mu = lapply(out$mu, function(x) {names(x) = cnames;x})
    out$Sigma = lapply(out$Sigma, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
    out$Sigma_reg = lapply(out$Sigma_reg, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
    out$Sigmai = lapply(out$Sigmai, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
    colnames(out$W) = cnames
    colnames(out$Ximp) = cnames
  }

  # rename observations
  if(!is.null(rnames)){
    rownames(out$W) = rnames
    rownames(out$Ximp) = rnames
    rownames(out$probs) = rnames
    names(out$class) = rnames
  }

  # rename groups
  if(!is.null(gnames)){
    names(out$Sigma) = names(out$Sigmai) = names(out$Sigma_reg) = names(out$mu) = gnames
    colnames(out$pi_groups) = rownames(out$pi_groups) = gnames
    colnames(out$probs) = gnames
    out$groups = ifelse(groups_isnum, as.numeric(groups_char), groups_char)
    out$class =  gnames[out$class]
  }

  return(out)
}




#' Calculation of Residuals for the Multi-Group GMM
#'
#' This function calculates the cell-wise residuals for each observation based on the fitted parameters
#' of a multi-group Gaussian Mixture Model (GMM) and the cellwise outlyingness pattern in matrix `W`.
#'
#' @param X A numeric data matrix or data frame with observations in rows and variables in columns.
#' @param groups A vector indicating pre-defined group membership for each observation (length must match `nrow(X)`).
#' @param Sigma A list of estimated covariance matrices.
#' @param mu A list of estimated mean vectors.
#' @param probs A matrix of posterior probabilities for each observation (rows) and group (columns).
#' @param W A binary matrix indicating which entries are considered non-outlying (1 = clean, 0 = outlying). Same dimensions as `X`.
#' @param set_to_zero A boolean indicating whether residuals of non-outlying cells should be set to zero.
#'
#' @details
#' Positive values of residuals mean that the observed value of the outlying variable is higher than would have been expected based on the other observed variables,
#' negative values mean that the observed value is lower than expected.
#' For non-outlying cells (i.e. where `W[i, j] == 1`), the residual is set to zero.
#'
#' @return A numeric matrix of residuals of the same dimension as `X`, where each cell represents the standardized deviation
#'         from the model-based conditional expectation, or zero if the cell was not flagged as outlying in `W`.
#'
#' @seealso \code{\link[ssMRCD]{cellMGGMM}}
#'
#' @references Puchhammer, P., Wilms, I., & Filzmoser, P. (2025). A smooth multi-group Gaussian Mixture Model for cellwise robust covariance estimation. \emph{ArXiv preprint} \doi{10.48550/arXiv.2504.02547}.
#'
#' @export
#'
#' @examples
#' data("weatherAUT2021")
#' cut_lon = c(min(weatherAUT2021$lon)-0.2, 12, 16, max(weatherAUT2021$lon) + 0.2)
#' cut_lat = c(min(weatherAUT2021$lat)-0.2, 48, max(weatherAUT2021$lat) + 0.2)
#' groups = ssMRCD::groups_gridbased(weatherAUT2021$lon, weatherAUT2021$lat, cut_lon, cut_lat)
#' N = length(unique(groups))
#' model = cellMGGMM(X = weatherAUT2021[, c("p", "s", "vv", "t", "rsum", "rel")],
#'                  groups = groups,
#'                  alpha = 0.5)
#' res = residuals_mggmm(X =  weatherAUT2021[, c("p", "s", "vv", "t", "rsum", "rel")],
#'                 groups = groups,
#'                 Sigma = model$Sigma,
#'                 mu = model$mu,
#'                 probs = model$probs,
#'                 W = model$W)

residuals_mggmm = function(X, groups, Sigma, mu, probs, W, set_to_zero = TRUE){

  # convert X to matrix if data.frame
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) stop("'X' must be a matrix or data frame that can be converted to a numeric matrix.")
  storage.mode(X) <- "numeric"

  # get variables
  n <- nrow(X)
  p <- ncol(X)
  N <- length(Sigma)

  # check groups
  if (length(groups) != n) stop("Length of 'groups' must equal number of rows in 'X'.")

  # check Sigma and mu
  if (!is.list(Sigma) || length(Sigma) != N) stop("'Sigma' must be a list of length equal to the number of groups.")
  if (!is.list(mu) || length(mu) != N) stop("'mu' must be a list of length equal to the number of groups.")
  if (!all(sapply(Sigma, function(s) all(dim(s) == c(p, p))))) stop("Each matrix in 'Sigma' must be p x p.")
  if (!all(sapply(mu, function(m) length(m) == p))) stop("Each mean vector in 'mu' must be of length p.")

  # check probs
  if (!is.matrix(probs) || nrow(probs) != n || ncol(probs) != N) {
    stop("'probs' must be an n x N matrix of posterior probabilities.")
  }

  # check W
  if (!is.matrix(W)) stop("'W' must be a numeric matrix.")
  if (!all(dim(W) == dim(X))) stop("'W' must have the same dimensions as 'X'.")
  if (!all(W %in% c(0, 1))) stop("'W' must only contain 0 (outlying) or 1 (non-outlying) entries.")

  # initialize residuals
  residuals = matrix(0, ncol = p, nrow = n)
  colnames(residuals) = colnames(X)
  rownames(residuals) = rownames(X)

  # get unique row patterns in W
  W_pattern = unique(W)
  wn = dim(W_pattern)[1]

  # iterate over unique row patterns
  for(l in 1:wn){
    # get missingness pattern
    obs = which(W_pattern[l, ] == 1)
    mis = which(W_pattern[l, ] == 0)

    # get observations with the same missingness pattern
    ind_w =  which(colSums(abs(t(W) - W_pattern[l, ])) == 0)

    # iterate over all distributions and variables
    for(k in 1:N){
      for(j in 1:p){

        # calculate help variable C and S
        obsj = obs[obs != j]
        S = Sigma[[k]][j, obsj, drop = FALSE] %*% chol2inv(chol(Sigma[[k]][obsj, obsj, drop = FALSE]))
        C = Sigma[[k]][j,j] - Sigma[[k]][j, obsj, drop = FALSE] %*% chol2inv(chol(Sigma[[k]][obsj, obsj, drop = FALSE])) %*% Sigma[[k]][obsj, j, drop = FALSE]
        C = as.numeric(C)

        # calculate conditional expectation and residual for each observation
        for(i in ind_w){
          xhat = mu[[k]][j] + S %*% c(X[i, obsj]-mu[[k]][obsj])
          residuals[i,j] = residuals[i,j] + probs[i,k] * (X[i, j] - xhat)/sqrt(C)
        }
      }
    }
  }

  # set non-outlying residuals to zero
  if(set_to_zero) {
    residuals[W == 1] = 0
  }

  return(residuals)
}

