##########################################################################################
# OBJECTIVE FUNCTION: EXPLAINED VARIANCE, SPARSITY, PENALTY                           ####
##########################################################################################

#' Objective function value for local sparse PCA
#'
#' @param PC vectorised component to evaluate.
#' @param eta degree of sparsity.
#' @param gamma distribution of sparsity between groupwise (\eqn{\gamma = 1}) and entrywise (\eqn{\gamma = 0}) sparsity.
#' @param COVS list of covariance matrices used for PCA
#'
#' @return Returns value of the objective function for given \code{v}.
#' @export
#'
#' @examples
#' S1 = matrix(c(1, 0.9, 0.8, 0.5,
#'               0.9, 1.1, 0.7, 0.4,
#'               0.8, 0.7, 1.5, 0.2,
#'               0.5, 0.4, 0.2, 1), ncol = 4)
#' S2 = t(S1)%*% S1
#' S2 = S2/2
#'
#' eval_objective(PC = c(1,0,0,0,sqrt(2),0,0,-sqrt(2)),
#'                eta = 1, gamma = 0.5,
#'                COVS = list(S1, S2))

eval_objective = function(PC, eta, gamma, COVS){

  N = length(COVS)
  p = dim(COVS[[1]])[1]
  S = (eta*gamma)*sum(abs(PC))  #L1- penalty
  for(i in 1:N){
    ind = jth_col(j = i, p = p)
    S = S - t(PC[ind])%*% COVS[[i]] %*% PC[ind]  #Variance
  }
  for (i in 1:p){ # groupwise penalty
    ind = jth_row(j = i, p = p, N = N)
    S = S + eta* (1-gamma)* sqrt(t(PC[ind])%*%PC[ind])
  }
  return(S)
}



#' Entry-wise Sparsity in the Loadings
#'
#' @param PC matrix-like object of PCs.
#' @param N integer, number of groups.
#' @param p integer, number of variables.
#' @param tolerance tolerance for sparsity.
#' @param k integer or integer vector of which component should be used.
#' @param scaled logical, if total number or percentage of possible sparse entries should be returned.
#'
#' @return Returns either a percentage (\code{scaled = TRUE}) or the amount of zero-values entries (\code{scaled = FALSE}).
#' @export
#'
#' @examples
#' PC = matrix(c(1,0,2,3,0,7,0,1,0,1,0.001,0), ncol = 2)
#' sparsity_entries(PC, N = 2, p = 3, tolerance = 0, k = 1, scaled = FALSE)
#' sparsity_entries(PC, N = 2, p = 3, tolerance = 0.001, k = 2, scaled = TRUE)
#'
sparsity_entries = function(PC,
                            N,
                            p,
                            tolerance = 0,
                            k = 1,
                            scaled = TRUE){

  PC = as.matrix(PC)
  sparse_entries = sum( abs(PC[, k]) <= tolerance)
  if(scaled){
    return( sparse_entries/(length(k)*N*(p-1))    )
  } else {
    return(sparse_entries)
  }
}




#' Group-wise Sparsity in the Loadings
#'
#' @param PC matrix-like object of PCs.
#' @param N integer, number of groups.
#' @param p integer, number of variables.
#' @param tolerance tolerance for sparsity.
#' @param k integer, which components should be used. Does not work for multiple PCs simultaneously.
#' @param scaled logical, if total number or percentage of possible sparse entries should be returned.
#'
#' @return Returns either a matrix of percentages (\code{scaled = TRUE}) or the amounts
#'  of zero-values entries (\code{scaled = FALSE}) for each group/neighborhood.
#' @export
#'
#' @examples
#' PC = matrix(c(1,0,2,3,0,7,0,1,0,1,0.001,0), ncol = 2)
#' sparsity_group(PC, N = 2, p = 3, tolerance = 0, k = 1, scaled = FALSE)
#' sparsity_group(PC, N = 2, p = 3, tolerance = 0.001, k = 2, scaled = TRUE)
sparsity_group = function(PC,
                          N,
                          p,
                          tolerance = 0,
                          k = 1,
                          scaled = TRUE){

  PC = as.matrix(PC)
  n_sparse_groups = 0
  for(j in 1:p){
    ind = jth_row(j = j, p = p, N = N)
    sparse_entries = sum(abs(matrix(PC[ind,k], ncol = 1)) <= tolerance)
    n_sparse_groups = n_sparse_groups + as.numeric(sparse_entries == N) # per variable there are N many entries
  }

  if (scaled ) {
    return(n_sparse_groups/((p-1)*length(k)))
  }  else{
    return(n_sparse_groups)
  }
}






#' Mixed Sparsity of the Loadings
#'
#' @param PC matrix-like object of PCs.
#' @param N integer, number of groups.
#' @param p integer, number of variables.
#' @param k integer, which components should be used. Does not work for multiple PCs simultaneously.
#' @param tolerance tolerance for sparsity.
#' @param mean if \code{"arithmetic"} or \code{"geometric"} mean should be used.
#'
#' @return Returns the geometric mean of the percentage of entry-wise and group-wise sparsity.
#' @export
#'
#' @examples
#' PC = matrix(c(1,0,2,3,0,7,0,1,0,1,0.001,0), ncol = 2)
#' sparsity_mixed(PC, N = 2, p = 3, tolerance = 0, k = 1)
#' sparsity_mixed(PC, N = 2, p = 3, tolerance = 0.001, k = 2, mean = "geometric")
sparsity_mixed = function(PC,
                          p,
                          N,
                          k = 1,
                          tolerance = 0.001,
                          mean = "arithmetic"){

  PC = as.matrix(PC)

  sg = sparsity_group(PC = PC,
                      N = N,
                      p = p,
                      tolerance = tolerance,
                      k = k,
                      scaled = TRUE)
  se = sparsity_entries(PC = PC,
                        N = N,
                        p = p,
                        tolerance = tolerance,
                        k = k,
                        scaled = TRUE)

  if(mean == "geometric") return(sqrt(se*sg))
  if(mean == "arithmetic") return((se+sg)/2)
}




#' Entry-wise Sparsity in the Loadings per Group
#'
#' @param PC matrix-like object of PCs.
#' @param N integer, number of groups.
#' @param p integer, number of variables.
#' @param tolerance tolerance for sparsity.
#' @param k integer or integer vector of which component should be used.
#' @param scaled logical, if total number or percentage of possible sparse entries should be returned.
#'
#' @return Returns either a matrix of percentages (\code{scaled = TRUE}) or the amounts
#'  of zero-values entries (\code{scaled = FALSE}) for each group/neighborhood.
#' @export
#'
#' @examples
#' PC = matrix(c(1,0,2,3,0,7,0,1,0,1,0.001,0), ncol = 2)
#' sparsity_summary(PC, N = 2, p = 3, tolerance = 0, k = 1, scaled = FALSE)
#' sparsity_summary(PC, N = 2, p = 3, tolerance = 0.001, k = 2, scaled = TRUE)

sparsity_summary = function(PC, N, p, tolerance = 0, k = 1, scaled = FALSE){

  PC = as.matrix(PC)
  sparse = matrix(NA, nrow = N, ncol = 1)
  for(j in 1:N){
    ind = jth_col(j = j, p = p)
    sparse[j, ] =  sparsity_entries(PC = PC[ind, ],
                                    N = 1,
                                    p = p,
                                    tolerance = tolerance,
                                    k = k,
                                    scaled = scaled)
  }
  return(sparse)
}



#' Explained Variance summarized over Groups
#'
#' @param COVS list of covariance matrices
#' @param PC matrix-like object holding the loadings of length np
#' @param k which component should be evaluated
#' @param type character, either \code{"scaled"} for scaling using the extremes solutions or \code{"percent"} as percentage of overall variance.
#' @param cor logical, if \code{COVS} is a correlation matrix or not
#' @param gamma scalar between 0 and 1 indicatig distribution of sparsity.
#'
#' @return Returns scalar
#' @export
#'
#' @importFrom Matrix bdiag
#'
#' @examples
#' S1 = matrix(c(1, 0.9, 0.8, 0.5,
#'               0.9, 1.1, 0.7, 0.4,
#'               0.8, 0.7, 1.5, 0.2,
#'               0.5, 0.4, 0.2, 1), ncol = 4)
#' S2 = t(S1)%*% S1
#' S2 = S2/2
#'
#' explained_var(COVS = list(S1, S2),
#'               PC = c(1,0,0,0,sqrt(2),0,0,-sqrt(2)),
#'               k = 1,
#'               cor = FALSE,
#'               gamma = 0.5)
#'
#' explained_var(COVS = list(cov2cor(S1), cov2cor(S2)),
#'               PC = c(1,0,0,0,sqrt(2),0,0,-sqrt(2)),
#'               k = 1,
#'               cor = TRUE,
#'               gamma = 0.5)

explained_var = function(COVS,
                         PC,
                         k,
                         type = "scaled",
                         cor = FALSE,
                         gamma = 0.5){

  PC = as.matrix(PC)
  SIGMA = as.matrix(Matrix::bdiag(COVS))
  var_full = sum(diag(SIGMA))
  var_expl = sum(diag(t(PC[, k]) %*% SIGMA %*% PC[, k]))

  if(type == "percent") {
    return(c("percent" = var_expl/var_full,
             "explained_tot" = var_expl,
             "total_tot" = var_full))
  }

  if(type == "scaled"){
    starts = starting_value_ADMM (Sigma = COVS,
                                  lambda = 1,
                                  alpha = gamma,
                                  k = k,
                                  return_all = TRUE,
                                  Xi = NULL,
                                  cor = cor)

    var_min = sum(diag(t(starts$penalty_solution) %*% SIGMA %*% starts$penalty_solution))
    var_max = sum(diag(t(starts$eigenvector) %*% SIGMA %*% starts$eigenvector))

    var_scaled = (var_expl - var_min)/(var_max - var_min)
    return(var_scaled)
  }
}



# for summary_diagnostic
explained_var_singular= function(x, COVS){
  # COVS: covariance matrix
  # x: p vector
  var_all = sum(diag(COVS))
  var_x = sum(diag(t(x) %*% COVS %*% x))
  return(c(var_x/var_all, var_x, var_all))
}

# for summary_diagnostic
explained_var_N = function(x, COVS){
  # COVS: list of covariance matrices
  # x: N*p vector
  N = length(COVS)
  p = dim(COVS[[1]])[1]
  vars = matrix(NA, nrow = N, ncol = 3)
  for(i in 1:N){
    ind = jth_col(j = i, p = p)
    vars[i, ] = explained_var_singular(x[ind], COVS[[i]])
  }
  colnames(vars) = c("varPC_perc", "varPC", "varbig")
  return(vars)
}


##########################################################################################
# STRUCTURAL                                                                          ####
##########################################################################################
# GET INDEXES
jth_col = function(j, p){
  # get p-dim index vector for j-th neighborhood
  seq( (j-1)*p + 1, j*p, by = 1)
}

jth_row = function(j, p, N){
  # get N-dim index vector for j-th variable
  j + p*(0:(N-1))
}

renorm = function(x, p, N){
  # renorm per neighborhoods for all neighborhoods
  for (j in 1:N){
    ind = jth_col(j = j, p = p)
    x[ind] = x[ind] / norm(matrix(x[ind]), "F")
  }
  return(x)
}

project_to_orthogonal = function(PC, Xi, N, p, renorm = FALSE){
  # project PC onto orthogonal space of Xi and renorm per one neighborhood
  # PC: vector (N*p)
  # Xi: matrix or vector

  k = 1

  if(is.null(dim(Xi))) Xi = matrix(Xi, ncol = 1)
  proj = rep(0, N*p)

  for (j in 1:N){
    # get neighborhood indices
    ind = jth_col(j = j, p = p)
    A_tmp = PC[ind]

    # get neighborhood specific Xi
    X_tmp = as.matrix(Xi, nrow = N*p, ncol = k)[ind, ]

    # calculate projection vector
    zz = matrix(0, ncol = p, nrow = dim(Xi)[2])
    diag(zz) = c(t(X_tmp) %*% A_tmp)
    sumi = rowSums(X_tmp %*% zz)

    # remove projection vector to get in orthogonal space
    proj[ind] = A_tmp - sumi
  }

  # renorm if necessary
  if(renorm) proj = renorm(x = proj, p = p, N = N)

  return(proj)
}



#' Align Loadings of Principal Components
#'
#' @description
#' Aligns loadings per neighborhood for better visualization and comparison. Different options are available.
#'
#' @param PC matrix of loadings of size Np x k
#' @param N integer, number of groups/neighborhoods
#' @param p integer, number of variables
#' @param type character indicating how loadings are aligned (see details),
#'             options are \code{"largest", "maxvar","nonzero","mean", "scalar", "none"}.
#' @param vec \code{NULL} or vector containing vectors for type \code{"scalar"}
#'
#' @return Returns a matrix of loadings of size \code{Np} times \code{k}.
#'
#' @details
#' For input \code{type} possible values are \code{"largest", "maxvar","nonzero","mean","scalar"}.
#' For option \code{"maxvar"} the variable with the highest absolute value in the loading
#' is scaled to be positive (per neighborhood, per loading).
#' For option \code{"nonzero"} the variable with largest distance to zero in the entries is
#' scaled to be positive (per neighborhood, per loading).
#' For option \code{"scalar"} the variable is scaled in a way, that the scalar product
#' between the loading and the respective part of \code{vec} is positive (per neighborhood, per loading).
#' If \code{vec} is of size \code{p} times \code{k}, the same vector is used for all neighborhoods.
#' Option \code{"mean"} is option \code{"scalar"} with \code{vec} being the mean of the loadings per variable across neighborhoods.
#' Option \code{"largest"} scales the largest absolute value to be positive per neighborhood and per PC.
#' Option \code{"none"} does nothing and returns \code{PC}.
#'
#' @export
#'
#' @examples
#' x = matrix(c(1, 0, 0, 0, sqrt(0.5), -sqrt(0.5), 0, 0,
#'              0, sqrt(1/3), -sqrt(1/3), sqrt(1/3), sqrt(0.5), sqrt(0.5), 0, 0),
#'            ncol = 2)
#' align_PC(PC = x, N = 2, p = 4, type = "largest")
#' align_PC(PC = x, N = 2, p = 4, type = "mean")

align_PC = function(PC,
                    N,
                    p,
                    type = "largest",
                    vec = NULL){

  type_orig = type
  if(type == "scalar" & is.null(vec)) stop("Necessary vec for option scalar is missing.")

  PC = as.matrix(PC)
  k = dim(PC)[2]
  if(!is.null(vec)) {
    vec = as.matrix(vec)
    if(dim(vec)[1] == p) {
      vec = matrix(rep(t(vec), times = N), nrow = p*N, ncol = k, byrow = TRUE)
    }
  }

  if(type == "none") return(PC)

  for(l in 1:k){
    x = PC[, l]
    veci = vec[, l]

    # setup
    if(type_orig == "maxvar"){
      m = matrix(abs(x), nrow = p, ncol = N, byrow = FALSE)
      var_ind = which.max(colSums(m))
    }
    if(type_orig == "nonzero"){
      m = matrix(abs(x), nrow = p, ncol = N, byrow = FALSE)
      var_ind = which.max(sapply(as.data.frame(t(m)), min))  # take variable where abs values are furthest away from zero
    }
    if(type_orig == "mean"){
      type ="scalar"
      veci = rep(rowMeans(matrix(x, ncol = N, byrow = FALSE)), times = N)
    }

    # align
    for(i in 1:N){
      ind = jth_col(j = i, p = p)
      x_i = x[ind]
      if(type == "maxvar" | type == "nonzero"){
        mult = sign(x_i[var_ind])
      }
      if(type == "largest"){
        # largest absolute value will be positive
        mult = sign(x_i[which.max(abs(x_i))])
      }
      if(type == "scalar"){
        # scalar product is positive
        a = matrix(veci[ind], ncol = 1)
        b = matrix(x_i, nrow = 1)
        mult = c(sign(b %*% a))
      }
      if(mult == 0) {
        # if orthogonal: default to largest
        mult = unlist(sign(x_i[which.max(abs(x_i))]))
      }
      if(mult == 0) mult = 1

      # redirect
      x[ind] = x[ind]*mult
    }
    PC[, l] = x
  }
  return(PC)
}


##########################################################################################
# ADJUSTED LAMBDA                                                                     ####
##########################################################################################
adjusted_eta = function(COVS, Xi, k, eta){

  # Xi: matrix(p*N x k) of PCs calculated
  # k: which PC is calculated next
  # eta: given eta by user
  # COVS: matrix list

  N = length(COVS)
  p = dim(COVS[[1]])[1]
  if(k == 1) {
    v = sum(sapply(COVS, function(x) eigen(x)$values[1]))
  } else {
    Xi = as.matrix(Xi)
    v = 0
    for(i in 1:N){
      Xii = Xi[jth_col(j = i, p = p), 1:(k-1)]
      P = diag(1, nrow = p) - Xii %*% t(Xii)
      S = Re(eigen(P %*% COVS[[i]] %*% t(P))$values[1])
      v = v + S
    }
  }
  return(eta * v/sum(sapply(COVS, function(x) eigen(x)$values[1])) )
}


##########################################################################################
# SCORES                                                                              ####
##########################################################################################

#' Calculate Scores for local sparse PCA
#'
#' @param X data set as matrix.
#' @param PC loading matrix.
#' @param groups vector of grouping structure (numeric).
#' @param ssMRCD ssMRCD object used for scaling \code{X}. If \code{NULL} no scaling and centering is performed.
#'
#' @return Returns a list with scores and univariately and locally centered and scaled observations.
#' @export
#'
#' @seealso  \code{\link[ssMRCD]{ssMRCD}}, \code{\link[ssMRCD]{scale_ssMRCD}}
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
#' pca = sparsePCAloc(eta = 1, gamma = 0.5, cor = FALSE,
#'                    COVS = loccovs$MRCDcov,
#'                    increase_rho = list(FALSE, 20, 1))
#'
#' # calculate scores
#' scores(X = rbind(x1, x2), PC = pca$PC,
#'        groups = rep(c(1,2), each = 100), ssMRCD = loccovs)

scores = function(X, PC, groups, ssMRCD = NULL){
  # calculates scores for all nieghborhoods

  # PC: PC matrix from PCA, already scaled with corresponding "eigenvalues"
  # ssMRCD: ssMRCD object

  PC = as.matrix(PC)
  X = as.matrix(X)
  N = length(unique(groups))
  p = dim(X)[2]
  n = dim(X)[1]
  k = dim(PC)[2]

  TT = matrix(NA, ncol = k, nrow = n)
  if(is.null(ssMRCD)) {
    centered_all = X
  } else {
    centered_all = scale_ssMRCD(ssMRCD = ssMRCD,
                                X = X,
                                groups = groups,
                                multivariate = TRUE,
                                center_only = TRUE)
  }

  for(l in 1:N){
    ind = jth_col(j = l, p = p)
    ind_obs = which(groups == l)
    centered = centered_all[ind_obs, ]
    TT[ind_obs, ] = centered %*% PC[ind, ]
  }

  return(list(scores = matrix(TT, ncol = k, nrow = n),
              X_centered = centered_all))  # X_centered also scaled if ssMRCD is used
}

#' Orthogonal Distances for PCAloc
#'
#' @param X data matrix of observations.
#' @param PC loadings of sparse local PCA.
#' @param groups grouping vector for locality.
#' @param ssMRCD ssMRCD object used for PCA calculation.
#'
#' @return Returns vector of orthogonal distances of observations.
#' @export
#'
#' @seealso  \code{\link[ssMRCD]{scores}}, \code{\link[ssMRCD]{scores.SD}}, \code{\link[ssMRCD]{sparsePCAloc}}, \code{\link[ssMRCD]{scale_ssMRCD}}
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
#' pca = sparsePCAloc(eta = 1, gamma = 0.5, cor = FALSE,
#'                    COVS = loccovs$MRCDcov,
#'                    increase_rho = list(FALSE, 20, 1))
#'
#' # calculate scores
#' scores.OD(X = rbind(x1, x2), PC = pca$PC,
#'           groups = rep(c(1,2), each = 100), ssMRCD = loccovs)
#'

scores.OD = function(X, PC, groups, ssMRCD){

  PC = as.matrix(PC)
  X = as.matrix(X)
  TT_all = scores(X = X, PC = PC, groups = groups, ssMRCD = ssMRCD)

  TT = as.matrix(TT_all$scores)
  XTT = as.matrix(TT_all$X_centered)

  N = length(unique(groups))
  p = dim(X)[2]
  n = dim(X)[1]
  k = dim(PC)[2]

  OD = rep(NA, n)
  for(l in 1:N){
    ind = jth_col(j = l, p = p)
    ind_obs = which(groups == l)
    tmp = XTT[ind_obs, ] -  TT[ind_obs, ] %*% t(PC[ind, ])
    OD[ind_obs] = sqrt(diag(tmp %*% t(tmp))) # norm per observation
  }
  return(OD)
}

#' Score Distances for PCAloc
#'
#' @param X data matrix of observations.
#' @param PC loadings of sparse local PCA.
#' @param groups grouping vector for locality.
#' @param ssMRCD ssMRCD object used for PCA calculation.
#'
#' @return Returns vector of score distances of observations.
#' @export
#'
#' @seealso  \code{\link[ssMRCD]{scores}}, \code{\link[ssMRCD]{scores.OD}}, \code{\link[ssMRCD]{sparsePCAloc}}, \code{\link[ssMRCD]{scale_ssMRCD}}
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
#' pca = sparsePCAloc(eta = 1, gamma = 0.5, cor = FALSE,
#'                    COVS = loccovs$MRCDcov,
#'                    increase_rho = list(FALSE, 20, 1))
#'
#' # calculate scores
#' scores.SD(X = rbind(x1, x2), PC = pca$PC,
#'           groups = rep(c(1,2), each = 100), ssMRCD = loccovs)
#'
scores.SD = function(X, PC, groups, ssMRCD){

  PC = as.matrix(PC)
  X = as.matrix(X)

  TT_all = scores(X = X, PC = PC, groups = groups, ssMRCD = ssMRCD)
  TT = as.matrix(TT_all$scores)

  N = length(unique(groups))
  p = dim(X)[2]
  n = dim(X)[1]
  k = dim(PC)[2]
  COVS = ssMRCD$MRCDcov

  SD = rep(NA, n)
  for(j in 1:N){
    ind = jth_col(j = j, p = p)
    ind_obs = which(groups == j)
    lambdainv = diag(t(PC[ind, ]) %*% COVS[[j]] %*% PC[ind, ])^(-1)
    SD[ind_obs] = diag(TT[ind_obs, ] %*% diag(lambdainv, ncol = k) %*% t(TT[ind_obs, ]))
  }

  return(SD)
}
