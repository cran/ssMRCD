#########################################
#### MAIN ssMRCD FUNCTION - EXTERNAL ####
#########################################


#' Spatially Smoothed MRCD Estimator
#'
#' The ssMRCD function calculates the spatially smoothed MRCD estimator from Puchhammer and Filzmoser (2023).
#'
#' @param X a list of matrices containing the observations per neighborhood sorted, or matrix or data frame containing data.
#'          If matrix or data.frame, group vector has to be given.
#' @param groups vector of neighborhood assignments
#' @param weights weighting matrix, symmetrical, rows sum up to one and diagonals need to be zero (see also \code{\link[ssMRCD]{geo_weights}} or \code{\link[ssMRCD]{time_weights}} .
#' @param lambda numeric between 0 and 1.
#' @param tuning default NULL. List of tuning specifications if lambda contains more than one value. See Details.
#' @param TM target matrix (optional), default value is the covMcd from robustbase.
#' @param alpha numeric, proportion of values included, between 0.5 and 1.
#' @param maxcond optional, maximal condition number used for rho-estimation.
#' @param maxcsteps maximal number of c-steps before algorithm stops.
#' @param n_initialhsets number of initial h-sets, default is 6 times number of neighborhoods.
#'
#'
#' @details
#' The necessary list elements for the parameter \code{tuning} depend on the method specified.
#' For both tuning approaches (residual-based or contamination-based) the element \code{method} needs to be specified to
#' \code{"residuals"} and \code{"local contamination"}, respectively. The boolean list element \code{plot} is available for both methods and
#' specifies if a plot should be constructed after tuning.
#'
#' For \code{tuning$method = "local contamination"}, additional information needs to be passed.
#' The number of nearest neighbors \code{tuning$k} used for the local outlier detection method
#' is \code{10} by default. The percentage of exchanged/contaminated observations is specified
#' by \code{tuning$cont} and is set to \code{0.05} by default. Also the coordinates must be given in  \code{tuning$coords}
#' and the number of repetitions for the switching procedure, \code{tuning$repetitions}.
#'
#' For \code{tuning$method = "local contamination"} no optimal value is returned but the choice has to
#' be made by the user. Be aware that the FNR does not take into account that there are also natural outliers
#' included in the data set that might or might not be found. The best parameter selection depends on the goal of the analysis
#' and whether false negatives should be avoided or whether the number of flagged outliers should be low.
#'
#'
#' @return The output depends on whether parameters are tuned.
#' If there is no tuning the output is an object of class \code{"ssMRCD"} containing the following elements:\tabular{ll}{
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
#' If parameters are tuned, the output consists of:\tabular{ll}{
#'    \code{ssMRCD} \tab Object of class ssMRCD with optimally selected parameter \code{lambda}. \cr
#'    \tab \cr
#'    \code{tuning_grid} \tab Vector of lambda to tune over given by the input. \cr
#'    \tab \cr
#'    \code{tuning_values} \tab If \code{tuning$method = "residuals"} then a vector returning
#'    the values of the residual criteria for the corresponding values of lambda in \code{tuning_grid}.\cr
#'    \tab If \code{tuning$method = "local contamination"}, then matrix with false negative rates
#'     and the total number of flagged outliers. \cr
#'    \tab \cr
#'    \code{plot} \tab If \code{tuning$plot = TRUE}, then a plot for parameter tuning is added. \cr
#'
#' }
#'
#'
#'
#' @seealso \code{\link[ssMRCD]{plot.ssMRCD}}
#'
#' @references Puchhammer P. and Filzmoser P. (2023). Spatially Smoothed Robust Covariance Estimation for Local Outlier Detection. \emph{Journal of Computational and Graphical Statistics, 33(3), 928–940}.  \doi{10.1080/10618600.2023.2277875}
#'
#' @export
#' @import ggplot2
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
#' out = ssMRCD(X = x, weights = W, lambda = 0.5)
#' str(out)


ssMRCD = function(X,
                  groups = NULL,
                  weights,
                  lambda = 0.5,
                  tuning = list(method = NULL, plot = FALSE,
                                k = 10, repetitions = 5, cont = 0.05),
                  TM = NULL,
                  alpha = 0.75,
                  maxcond = 50,
                  maxcsteps = 200,
                  n_initialhsets = NULL){

  # input checks
  if(alpha > 1 | alpha  < 0.5) stop("The robustness parameter alpha needs to be inside [0.5, 1].")
  if(all(lambda > 1) | all(lambda < 0)) stop("The smoothing parameter lambda needs to be inside [0, 1].")

  # collect names
  cnames = NULL
  rnames = NULL

  if(is.list(X) | is.null(groups)){
    cnames = colnames(X[[1]])
    rnames = do.call(rbind, sapply(X, rownames))
  }

  if( (is.data.frame(X) | is.matrix(X)) & !is.null(groups)) {
    X = as.matrix(X)
    cnames = colnames(X)
    rnames = rownames(X)
    if(nrow(X)!= length(groups)) stop("The length of groups must match the number of observations in matrix X.")
    X = restructure_as_list(data = X, groups = groups)
  } else if (!is.list(X)) {
    stop("X should either be list of data matrices or a data matrix with grouping given by groups input.")
  }



  weights = as.matrix(weights)
  if(!sum(diag(weights)) == 0){
    warning("Diagonal elements of weights are not all equal to zero and are set to zero for further calculations.")
    diag(weights)= 0
  }
  if(!all.equal(weights/rowSums(weights), weights)) {
    warning("The rows of weights are not scaled to sum up to one. They are rescaled by the algorithm.")
    weights = weights/rowSums(weights)
    if(any(!is.finite(weights))){
      stop("No rescaling of weight matrix possible. There should be at least one positive value per row.")
    }
  }
  if(any(weights < 0)) stop("There are negative values in weights. Please choose non-negative weights.")


  if(!is.null(groups)){
    groups_isnum = is.numeric(groups)
    groups_char = as.character(factor(groups))
    groups = as.numeric(factor(groups))
  }
  if(is.null(groups)) {
    groups_isnum = TRUE
    groups = rep(1:length(X), times = sapply(X, nrow))
    groups_char = as.character(factor(groups))
  }

  # collect group names
  gnames = unique(groups_char)[sort.int(unique(groups), index.return = TRUE)$ix]


  # no tuning
  if(length(lambda) == 1){
    out = ssMRCD_internal(X = X,
                          weights = weights,
                          lambda = lambda,
                          TM = TM,
                          alpha = alpha,
                          maxcond = maxcond,
                          maxcsteps = maxcsteps,
                          n_initialhsets = n_initialhsets)

    # rename variables
    if(!is.null(cnames)){
      out$MRCDmu = lapply(out$MRCDmu, function(x) {names(x) = cnames;x})
      out$MRCDcov = lapply(out$MRCDcov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
      out$MRCDicov = lapply(out$MRCDicov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
      out$Kcov = lapply(out$Kcov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
      colnames(out$mT) = rownames(out$mT) = cnames
    }

    # rename groups
    if(!is.null(gnames)){
      names(out$MRCDmu) = names(out$MRCDcov) = names(out$MRCDicov) = names(out$Kcov)= gnames
      colnames(out$weights) = rownames(out$weights) = gnames
      #out$groups = ifelse(groups_isnum, as.numeric(groups_char), groups_char)
    }
  }


  # with tuning of lambda
  if(length(lambda) > 1 & !is.null(tuning$method)){

    # check input
    if(!tuning$method %in% c("residuals", "local contamination")) stop("Tuning method must be either 'residuals' or 'local contamination'.")
    message(paste0("Hyperparameter tuning of lambda based on method '", tuning$method, "'."))

    # make target matrices
    if(is.null(TM)) {
      # get the same target matrix
      Xsum <- do.call(rbind, X)
      TM = tryCatch({
        robustbase::covMcd(Xsum, alpha = alpha)$cov
      }, error = function(cond){
        message("Singularity issues on the global level. Use MRCD estimator with default values as target matrix.")
        rrcov::CovMrcd(Xsum, alpha = alpha)$cov
      })
    }



    # tuning residual based
    if(tuning$method == "residuals"){
      out = tuning_residuals(X = X,
                             weights = weights,
                             lambda = lambda,
                             TM = TM,
                             alpha = alpha,
                             maxcond = maxcond,
                             maxcsteps = maxcsteps,
                             n_initialhsets = n_initialhsets)

      if(tuning$plot){
        g = ggplot2::ggplot() +
          ggplot2::geom_line(aes(x = out$tuning_grid,
                                 y = out$tuning_values)) +
          ggplot2::geom_point(aes(x = out$tuning_grid,
                                  y = out$tuning_values)) +
          ggplot2::geom_point(aes(x = out$lambda_optimal,
                                  y = out$tuning_values[out$tuning_grid == out$lambda_optimal]),
                              col = "red") +
          ggplot2::theme_bw() +
          ggplot2::labs(x = expression(lambda),
                        y = "Residual Value",
                        title = "Optimal Smoothing Based on Residuals")
        print(g)
        out$plot = g

        # rename variables
        if(!is.null(cnames)){
          out$ssMRCD_optimal$MRCDmu = lapply(out$ssMRCD_optimal$MRCDmu, function(x) {names(x) = cnames;x})
          out$ssMRCD_optimal$MRCDcov = lapply(out$ssMRCD_optimal$MRCDcov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
          out$ssMRCD_optimal$MRCDicov = lapply(out$ssMRCD_optimal$MRCDicov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
          out$ssMRCD_optimal$Kcov = lapply(out$ssMRCD_optimal$Kcov, function(x) {colnames(x) = cnames; rownames(x) = cnames;x})
          colnames(out$ssMRCD_optimal$mT) = rownames(out$ssMRCD_optimal$mT) = cnames
        }

        # rename groups
        if(!is.null(gnames)){
          names(out$ssMRCD_optimal$MRCDmu) = names(out$ssMRCD_optimal$MRCDcov) = names(out$ssMRCD_optimal$MRCDicov) = names(out$ssMRCD_optimal$Kcov) = names(out$ssMRCD_optimal$mT) = gnames
          colnames(out$ssMRCD_optimal$weights) = rownames(out$ssMRCD_optimal$weights) = gnames
        }
      }
    }


    # tuning contamination based
    if(tuning$method == "local contamination"){
      if(is.null(tuning$k) |
         is.null(tuning$cont) |
         is.null(tuning$coords) |
         is.null(tuning$repetitions)){
        stop("Some specifications for 'local contamination' tuning are missing in the list tuning.")
      }
      out = tuning_contamination(data = Xsum,
                                 coords = tuning$coords,
                                 groups = groups,
                                 lambda = lambda,
                                 weights = weights,
                                 k = tuning$k,
                                 cont = tuning$cont,
                                 repetitions = tuning$repetitions)

      if(tuning$plot){
        g1 = ggplot2::ggplot() +
          ggplot2::geom_boxplot(aes(x = out$tuning_values$lambda,
                                    y = out$tuning_values$FNR,
                                    group = out$tuning_values$lambda)) +
          ggplot2::theme_bw() +
          ggplot2::labs(x = expression(lambda),
                        y = "False Negative Rate",
                        title = "False Negative Rate Based on Additonal Contamination")
        print(g1)

        g2 = ggplot2::ggplot() +
          ggplot2::geom_boxplot(aes(x = out$tuning_values$lambda,
                                    y = out$tuning_values$n_outliers,
                                    group = out$tuning_values$lambda)) +
          ggplot2::theme_bw() +
          ggplot2::labs(x = expression(lambda),
                        y = "Number of Flagged Local Outliers",
                        title = "Number of Flagged Outliers Based on Additonal Contamination")
        print(g2)

        out$plot_fnr = g1
        out$plot_nout = g2
      }

      return(out)
    }
  }

  out$rnames = rnames
  out$cnames = cnames
  out$gnames = gnames

  return(out)
}


#########################################
#### MAIN ssMRCD FUNCTION - INTERNAL ####
#########################################


#' @importFrom robustbase covMcd
#' @importFrom rrcov CovMrcd
ssMRCD_internal = function(X,
                  weights,
                  lambda = 0.5,
                  TM = NULL,
                  alpha = 0.75,
                  maxcond = 50,
                  maxcsteps = 200,
                  n_initialhsets = NULL){

  N <- length(X)
  p <- dim(X[[1]])[2]
  if(is.null(n_initialhsets)) n_initialhsets <- N*6

  # bind data matrix
  Xsum <- do.call(rbind, X)

  # get target matrix
  if(is.null(TM)) {
    TM = tryCatch({
      robustbase::covMcd(Xsum, alpha = alpha)$cov
    }, error = function(cond){
      message("Singularity issues on the global level. Use MRCD estimator with default values as target matrix.")
      rrcov::CovMrcd(Xsum, alpha = alpha)$cov
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
  mTeigen <- eigen(mT, symmetric = TRUE)
  mQ <- mTeigen$vectors
  mL <- diag(mTeigen$values)
  msqL <- diag(sqrt(mTeigen$values))
  misqL <- diag(sqrt(mTeigen$values)^(-1))
  sdtrafo <- mQ %*% msqL

  # initial step transformation and rho selection per neighborhood
  init <- list()
  for( i in 1:N) {
    temp <- initial_step(x = X[[i]],
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
  obj_fun_values <- rbind(obj_fun_values, c(ret$obj_value))


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
    obj_fun_values <- rbind(obj_fun_values, c(tmp$obj_value))

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





###########################
#### TUNING FUNCTIONS #####
###########################
# Optimal Smoothing Parameter for ssMRCD based on Local Outliers
#
# This function provides insight into the effects of different parameter settings.
#
# @param data matrix with observations.
# @param coords matrix of coordinates of these observations.
# @param groups numeric vector, the neighborhood structure that should be used for \code{\link[ssMRCD]{ssMRCD}}.
# @param lambda scalar, the smoothing parameter.
# @param weights weighting matrix used in \code{\link[ssMRCD]{ssMRCD}}.
# @param k vector of possible k-values to evaluate.
# @param dist vector of possible dist-values to evaluate.
# @param cont level of contamination, between 0 and 1.
# @param repetitions number of repetitions wanted to have a good picture of the best parameter combination.

tuning_contamination = function(data,
                            coords,
                            groups,
                            lambda = c(0, 0.25, 0.5, 0.75, 0.9),
                            weights = NULL,
                            k = NULL,
                            cont = 0.05,
                            repetitions = 5){

  # function to introduce outliers
  contamination_random = function(cont, data){

    n = dim(data)[1]
    n_cont = 2*round(n*cont/2)

    # select outliers
    outliers = sample(1:n, size = n_cont , replace = FALSE)

    # switch first outlier observation with last outlier and so forth
    for (i in 1:(n_cont/2)){
      tmp = data[outliers[i], ]
      data[outliers[i],] =  data[outliers[n_cont - (i-1)],]
      data[outliers[n_cont - (i-1)],] = tmp
    }
    data_switched = data

    return(list(data = data_switched, out = outliers))

  }

  coords = as.matrix(coords)
  parameter_combinations = expand.grid(lambda, 1:repetitions)
  colnames(parameter_combinations) = c("lambda", "repetitions")

  # parameter settings
  n_par = dim(parameter_combinations)[1]
  n_outliers = rep(NA, n_par)
  FNR = rep(NA, n_par)

  for(i in 1:n_par){

    reps = parameter_combinations[i, 2]

    # contamination
    set.seed(reps)
    contaminated = contamination_random(cont = cont, data = data)
    data_contam = contaminated$data
    outliers = contaminated$out

    # performance
    lambda_i = parameter_combinations[i, 1]
    tmp = locOuts(data_contam,
                  coords,
                  groups,
                  lambda = lambda_i,
                  k = k,
                  dist = NULL,
                  weights = weights)

    outs_method = tmp$outliers
    n_outliers[i] = length(outs_method)
    FNR[i] = 1 - sum(outliers %in% outs_method)/length(outliers)
  }

  param_results = cbind(parameter_combinations, "FNR" = FNR, "n_outliers" = n_outliers)

  return(list(tuning_values = param_results,
              tuning_grid = lambda))
}




# Optimal Smoothing Parameter for ssMRCD based on Residuals
#
# The optimal smoothing value for the ssMRCD estimator is based on the residuals and the
# trimmed mean of the norm.
#
# @param X list of data matrix containing observations.
# @param weights weight matrix for groups, see \code{\link[ssMRCD]{rescale_weights}}, and \code{\link[ssMRCD]{geo_weights}}.
# @param lambda vector of parameter values for smoothing, between 0 and 1.
# @param TM target matrix.
# @param alpha percentage of outliers to be expected.
tuning_residuals = function(X,
                            weights,
                            lambda = seq(0, 1, 0.1),
                            TM,
                            alpha = 0.75,
                            maxcond = 50,
                            maxcsteps = 100,
                            n_initialhsets = NULL){


  COV_list = list()
  for(i in 1:length(lambda)){
    COVS = ssMRCD_internal( X = X,
                            weights = weights,
                            maxcond = maxcond,
                            TM = TM,
                            lambda = lambda[i],
                            alpha = alpha,
                            maxcsteps = maxcsteps,
                            n_initialhsets = n_initialhsets)
    COV_list = c(COV_list, list(COVS))

  }

  names(COV_list) = paste0("lambda", lambda)
  resid = sapply(X = COV_list,
                 function(x) residuals.ssMRCD(object = x,
                                              type = "trimmed_mean"))
  names(resid) = paste0("lambda", lambda)

  lambda_opt = lambda[which.min(resid)]

  return(list(ssMRCD_optimal = COV_list[[which(lambda == lambda_opt)]],
              lambda_optimal = lambda_opt,
              tuning_grid = lambda,
              tuning_values = resid))
}



####################################
#### INTERNAL HELPER FUNCTIONS #####
####################################


# Applies transformations to data, calculates the neighborhood-specific regularization parameter (rho), and generates six initial h-sets for the neighborhood.
# @param x A matrix of observations.
# @param alpha Robustness factor.
# @param maxcond Maximum condition number.
# @param mQ Eigenvectors of the target matrix.
# @param misqL Diagonal matrix with the inverse square root of eigenvalues of the target matrix.
# @return A list containing elements such as `nsets`, `mX`, `rho`, `h`, `scfac`, `hsets.init`, `initV`, `setsV`, `n`, and `p`.
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


# Calculates six robust initial h-sets based on robust covariance estimations from package covMrcd.
# @param x A numeric matrix constituting the observation matrix.
# @param h Number of observations used for robust covariance estimations.
# @param full.h Logical, only if adjust.eignevalues is TRUE
# @param adjust.eignevalues Logical, to adjust eigenvalues.
# @param scaled Logical, indicating whether data is already scaled.
# @param scalefn Scaling function (default: robustbase::Qn).
# @return A matrix containing six h-sets used in MRCD initialization.
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


# Calculates the best regularization parameter rho for neighborhood covariance estimation.
# @param mX Matrix with all observations (observations as columns).
# @param mT Target covariance matrix.
# @param hsets.init Matrix with six initial h-sets, one column per h-set.
# @param maxcond Maximum condition number for the regularized matrix.
# @param scfac Consistency factor.
# @param h Number of observations for covariance calculations.
# @return A list with values for `rho`, `initV`, and `setsV`.
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


# Back-transforms results from the c-step procedure using the best combination of h-set and initial estimates.
# @param best_ret List object produced by the cstep function.
# @param TM Target matrix.
# @param alpha Robustness factor.
# @param mQ Eigenvectors of the target matrix.
# @param msqL Square root of the eigenvalues of the target matrix.
# @return A list containing back-transformed observations and covariance matrices.
back_transformation = function(best_ret, TM, alpha, mQ, msqL){

  #browser()
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

    hindex = sort(c(best_ret$out[[i]]$index))
    hset = append(hset, list(hindex))

    mx = best_ret$out[[i]]$mX
    best_ret$out[[i]]$vMu = c(best_ret$out[[i]]$vMu)

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




# Computes the consistency factor for robust MCD estimator.
# @param p Dimension of the data.
# @param alpha Percentage of observations used for estimation.
# @return Consistency factor.
MDcons_robustbase = function (p, alpha) {
  qalpha <- stats::qchisq(alpha, p)
  caI <- stats::pgamma(qalpha/2, p/2 + 1)/alpha
  1/caI
}


# Check if Values are Monotonically Decreasing
# @param x Numeric vector.
# @return Logical, TRUE if the values are monotonically decreasing, FALSE otherwise.
monotonic_decreasing = function(x, tol = 1e-08){

  x = stats::na.omit(x)
  n = length(x)

  monotonic = TRUE
  for(i in 1:(n-1)){
    if(x[i] + tol < x[i+1]){
      monotonic = FALSE
      break
    }
  }

  return(monotonic)
}

# Linear Combination of Covariances
# @param x Object containing multiple MRCD covariance matrices as list.
# @return An updated object with combined covariance matrices.
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
