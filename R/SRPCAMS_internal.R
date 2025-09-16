##########################################################################################
#### STARTING VALUE                                                                   ####
##########################################################################################

solveL1_cov = function(Sigma, k){

  # Calculate Extreme Solutions for L1 Penalty and Covariance Matrix
  # @param Sigma A list of covariance matrices, where each matrix corresponds to a neighborhood.
  # @param k An integer specifying which eigenvector component to use for the solution.
  # @return A numeric vector of length \code{N*p}, where \code{N} is the number of covariance matrices in \code{Sigma} and \code{p} is the dimension of each covariance matrix.

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]

  sol_k = rep(0, N*p)
  for(i in 1:N){
    # get variable with highest variance per neighborhood
    ind = sort.int(diag(Sigma[[i]]),
                   index.return = T,
                   decreasing = T)$ix[k]
    # get value of eigenvector connected to variable with highest variance
    tmp = sign(eigen(Sigma[[i]])$vectors[ind,k])
    if(tmp == 0) tmp = 1
    sol_k[jth_col(j = i, p = p)][ind] = tmp
  }
  return(sol_k)
}

solveL1_cor = function(Sigma, k){

  # Calculate Extreme Solutions for L1 Penalty and Correlation Matrix
  # @param Sigma A list of correlation matrices.
  # @param k An integer specifying which eigenvector to use.
  # @return A numeric vector of length \code{N*p}, where \code{N} is the number of matrices and \code{p} is the matrix dimension.

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]

  # take variable with largest entry in eigenvector as entry-wise variable
  sol_k = rep(0, N*p)
  for(i in 1:N){
    EV = eigen(Sigma[[i]])$vectors[,k]
    ind = which.max(abs(EV))
    sign = sign(EV[ind])  # since absolute max, always non zero
    sol_k[jth_col(j = i, p = p)][ind] = sign
  }
  return(sol_k)
}

solveL2_cov = function(Sigma, k){

  # Calculate Extreme Solutions for L2 Penalty (Group) and Covariance Matrix
  # @param Sigma A list of covariance matrices.
  # @param k An integer specifying which variable to select based on variance.
  # @return A numeric vector of length \code{N*p}, where \code{N} is the number of matrices
  # and \code{p} is the matrix dimension.

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]

  var_k = rep(NA, N*p)
  sol_k = rep(0, p)
  sign_vec =  rep(0, N)

  # get correct rotation for sparse solution and calculate variance
  for(i in 1:N){
    var_k[jth_col(j = i, p = p)] = diag(Sigma[[i]])
    ind = sort.int(diag(Sigma[[i]]),
                   index.return = T,
                   decreasing = T)$ix[k]   # get variable with highest variance per neighborhood
    sign_vec[i] = sign(eigen(Sigma[[i]])$vectors[ind,k])  # get value of eigenvector connected to variable with highest variance
    if(sign_vec[i] == 0) sign_vec[i] = 1  # if orthogonal to original eigenvector!
  }

  # group-wise sparsity
  var_variables = rep(NA, p)
  for(i in 1:p){
    ind = jth_row(j = i, p = p, N = N)
    var_variables[i] = sum(var_k[ind])
  }

  # construct sparsity vector with right rotation
  sorted_ind = sort.int(var_variables, index.return = T, decreasing = T)$ix
  sol_k[sorted_ind[k]] = 1
  groupwise_k = rep(sol_k, times = N)* rep(sign_vec, each = p)

  return(groupwise_k)
}

solveL2_cor = function(Sigma, k){

  # Calculate Extreme Solutions for L2 Penalty (Group) and Correlation Matrix
  # @param Sigma A list of correlation matrices.
  # @param k An integer specifying which eigenvector component to use.
  # @return A numeric vector of length \code{N*p}, where \code{N} is the number of matrices and \code{p} is the matrix dimension.

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]

  vec = rep(0, p)
  sign_vec = rep(0, N)

  # get loadings
  EVEC_orig = lapply(Sigma, function(x) eigen(x)$vectors[,k] * eigen(x)$values[k])

  # orient to the same direction/halfspace
  EVEC_orient = EVEC_orig
  for(i in 1:N){
    sign = t(EVEC_orient[[1]]) %*% EVEC_orient[[i]]
    if(sign < 0) EVEC_orient[[i]] = EVEC_orient[[i]] * (-1)
  }

  # sum over neighborhoods
  EVEC_mean = sapply(EVEC_orient, function(x) x)
  EVEC_mean = rowMeans(EVEC_mean)

  # take largest entry as groupwise variable
  which_var = which.max(abs(EVEC_mean))

  # get correct rotation for sparse vector
  for(i in 1:N){
    sign_vec[i] = sign(eigen(Sigma[[i]])$vectors[which_var,k])  # get value of eigenvector connected to variable with highest variance
    if(sign_vec[i] == 0) sign_vec[i] = 1  # if orthogonal to original eigenvector!
  }
  vec[which_var] = 1
  groupwise_k = rep(vec, times = N)* rep(sign_vec, each = p)

  return(groupwise_k)
}

starting_value_ADMM = function(Sigma, eta, gamma, k = 1, Xi = NULL, cor = FALSE, return_all = FALSE){

  # Get Starting Value for ADMM
  # @param Sigma A list of covariance or correlation matrices.
  # @param eta A numeric value representing the regularization parameter.
  # @param gamma A numeric value controlling the mixture of penalty types (1 for entrywise, otherwise groupwise).
  # @param k An integer specifying which eigenvector component to use (default is 1).
  # @param return_all A logical indicating whether to return all computed components (default is FALSE).
  # @param Xi Optional matrix used for orthogonal projection (only relevant if \code{k > 1}).
  # @param cor A logical indicating whether to use correlation matrices (TRUE) or covariance matrices (FALSE).
  # @return A numeric matrix of size \code{N*p} where \code{N} is the number of matrices and \code{p} is the matrix dimension.
  # If \code{return_all} is TRUE, a list with additional components is returned:
  # \itemize{
  #   \item \code{starting_value}: The computed starting value for ADMM.
  #   \item \code{penalty_solution}: The penalty solution vector.
  #   \item \code{eigenvector}: The eigenvector estimates.
  # }

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]
  var_k = rep(NA, N*p)

  # eigenvalue solution
  eigenvec_start_k = rep(NA, N*p)
  for(i in 1:N){
    tmp = eigen(Sigma[[i]], symmetric = TRUE)$vectors[,k]
    eigenvec_start_k[jth_col(j = i, p = p)] = tmp
    var_k[jth_col(j = i, p = p)] = diag(Sigma[[i]])
  }

  # extreme sparse solutions
  if(!cor) {
    groupwise_start_k = solveL2_cov(Sigma = Sigma, k = k)
    entrywise_start_k = solveL1_cov(Sigma = Sigma, k = k)
  }
  if(cor) {
    groupwise_start_k = solveL2_cor(Sigma = Sigma, k = k)
    entrywise_start_k = solveL1_cor(Sigma = Sigma, k = k)
  }

  # best full-sparse extreme value
  if(gamma == 1)  penalty_start = entrywise_start_k
  if(gamma != 1) penalty_start = groupwise_start_k

  # take average
  starting_value = 0.5*eigenvec_start_k + 0.5*penalty_start

  # project to correct space
  if(k > 1 & !is.null(Xi)) {
    starting_value = project_to_orthogonal(PC = starting_value,
                                           p = p,
                                           N = N,
                                           renorm = TRUE,
                                           Xi = Xi)
  }

  # special case
  if(eta == 0) starting_value = eigenvec_start_k

  # return
  if(!return_all) return(starting_value)

  return(list(starting_value = starting_value,
              penalty_solution = penalty_start,
              eigenvector = eigenvec_start_k))
}


##########################################################################################
#### SOFT THRESHOLDING                                                                ####
##########################################################################################

soft_thresholding_scalar = function(a, kappa){

  # Soft Thresholding for Scalar Values
  # @param a A numeric scalar that represents the value to be thresholded.
  # @param kappa A numeric scalar that represents the threshold parameter.

  if(a > kappa)
    return(a-kappa)
  if(a < -kappa)
    return(a+kappa)
  if(abs(a) <= kappa)
    return(0)
}

solve_minimization_scalar_softthreshold = function(eta, gamma, rho, U2, A_mean){

  # Solve Minimization with Scalar Soft Thresholding
  # @param eta A numeric value representing the sparsity parameter.
  # @param gamma A numeric value that scales the sparsity parameter.
  # @param rho A numeric value for the ADMM parameter.
  # @param U2 A numeric vector representing the dual variable.
  # @param A_mean A numeric vector representing the mean of matrices.

  p = length(A_mean)
  A_res = rep(NA, p)
  kappa = (eta*gamma)/rho
  a = A_mean - U2/rho

  for(i in 1:p) {
    A_res[i] = soft_thresholding_scalar(a = a[i],
                                        kappa = kappa)
  }

  return(A_res)
}

soft_thresholding_group = function(a, kappa){

  # Soft Thresholding for Group Variables
  # @param a A numeric vector representing observations in one group.
  # @param kappa A numeric scalar representing the thresholding parameter.
  # @return A numeric vector where each element is scaled based on the soft thresholding rule.

  a = as.matrix(a)
  Nor = norm(a, type = "F")
  if(Nor == 0) {
    return (a)
  } else{
    mult = max(1 - kappa/Nor, 0)
    return(a * mult)
  }
}

solve_minimization_group_softthreshold = function(eta, gamma, rho, U3, A_mean, N){

  # Soft Thresholding for Group Variables
  # @param eta A numeric scalar representing the sparsity parameter.
  # @param gamma A numeric scalar between 0 and 1 that determines the balance between L1 and L2 penalties.
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param U3 A numeric vector representing the dual variable for the group minimization.
  # @param A_mean A numeric vector representing the mean of the group estimates to be adjusted.
  # @param N An integer representing the number of groups.
  # @return A numeric vector of the same length as \code{A_mean}, with group-wise soft thresholding applied.

  eta = eta * sqrt(N)
  kappa = eta*(1-gamma)/rho
  a = A_mean - U3/rho
  p = length(A_mean)/N

  A_res = rep(NA, p*N)
  for(j in 1:p) {
    ind = jth_row(j = j, p = p, N = N)
    A_res[ind] = soft_thresholding_group(a = a[ind],
                                         kappa = kappa)
  }

  return(A_res)
}

##########################################################################################
#### ADMM HELPER FUNCTION                                                             ####
##########################################################################################

f_find_root_of = function(x_mu_lambda, Sigma, U, rho, A, k, w, Xi = NULL){

  # Function to Find the Root of
  # @param x_mu_lambda A numeric vector combining the variables \(x\), \(mu\), and \(lambda\). The vector should be structured as follows:
  #        - The first \(p\) elements represent the variable \(x\),
  #        - The next element represents \(mu\),
  #        - The remaining elements represent \(lambda\) (only if \(k > 1\)).
  # @param Sigma A numeric matrix representing the covariance matrix for one neighborhood.
  # @param U A numeric vector representing the dual variable \(U\).
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param A A numeric vector representing the mean vector \(A\).
  # @param k An integer representing the number of principal components considered.
  # @param w A numeric vector representing the projection vector for uniqueness.
  # @param Xi A numeric matrix with prior principal components (dimensions \(p x (k-1)\)). It can be `NULL` if \(k = 1\).
  # @return A list containing:
  #   \item{res}{A numeric vector of the gradient of the Lagrangian function with respect to \(x\), \(mu\), and \(lambda\) (if \(k > 1\)).}
  #   \item{mu}{A numeric scalar representing \(mu * x'w\).}

  p = dim(Sigma)[1]
  x = x_mu_lambda[1:p]
  mu = x_mu_lambda[p+1]

  z = U*(rho^(-1)) - A
  lambda0 = c(t(x) %*% Sigma %*% x) - (rho/2)*(c(t(x) %*% z) + 1)
  mu_1 = mu
  sum_i = rep(0, p)
  if(k > 1){
    lambda = x_mu_lambda[(p+2) : length(x_mu_lambda)]
    xxj = rep(0, k-1)
    for(i in 1:(k-1)){
      sum_i = sum_i + lambda[i] * Xi[, i]
      xxj[i] = t(x) %*%  Xi[, i]
    }
  }

  # Lagrange Gradient
  res =  -2 * (Sigma %*% x) + rho * (z+x) - mu_1 * w + 2* lambda0 * x + sum_i
  stopifnot(length(res) == p)

  mux = mu_1 * (t(x) %*% w)
  res =  unname(rbind(res, mux))
  if(k > 1){
    res =  c(res, xxj)
  }

  return(list(res = res, mu = mux))
}

#' @importFrom rootSolve multiroot
find_root = function(Sigma, U, rho, A, k, Xi, w,
                     starting_value = NULL, mu_start = 20, eps_root = 1e-2, maxiter_root = 100){

  # Find the Root for an Optimization Problem
  # @param Sigma A numeric matrix representing the covariance matrix.
  # @param U A numeric vector representing the dual variables.
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param A A numeric vector representing the mean of all As from the ADMM iteration.
  # @param k An integer specifying which of principal components is calculated.
  # @param Xi A numeric matrix of available principal components (k-1 components), or NULL if not applicable.
  # @param w A numeric vector of projection vectors for uniqueness.
  # @param starting_value A numeric vector representing the starting value for the root-finding algorithm, or NULL to use default.
  # @param mu_start A numeric scalar specifying the initial guess for the mu parameter (default is 20).
  # @param eps_root A numeric scalar specifying the tolerance for convergence in the root-finding algorithm (default is 1e-2).
  # @param maxiter_root An integer specifying the maximum number of iterations for the root-finding algorithm (default is 100).
  # @return A list with a component:
  # \itemize{
  #   \item \code{roots}: A numeric matrix of the roots found, where each row represents a solution.
  # }

  p = dim(Sigma)[1]
  f = function(x_mu_lambda) f_find_root_of(x_mu_lambda = x_mu_lambda,
                                           Sigma = Sigma,
                                           U = U,
                                           rho = rho,
                                           A = A,
                                           k = k,
                                           Xi = Xi,
                                           w= w)$res

  # update starting value for toot finder
  starting_value = c(starting_value, mu_start)
  if(k > 1){
    n = length(starting_value)
    proj = project_to_orthogonal(PC = starting_value[1: (n-1)],
                                 Xi = Xi,
                                 N = 1 , # only one neighborhood
                                 p = p,
                                 renorm = TRUE)
    starting_value[1: (n-1)] = proj

    # calculate lambda_j for starting value
    lambdaj = c()
    for(j in 1:(k-1)){
      tmp =   t(Xi[, j]) %*% ( 2* Sigma %*% starting_value[1:(n-1)] - rho * (U/rho - A) + starting_value[n] * w)
      lambdaj = c(lambdaj, tmp)
    }
    starting_value = c(starting_value, lambdaj)
  }

  # calculate root of function
  tmp = rootSolve::multiroot(f = f,
                             start = c(starting_value),
                             useFortran = FALSE,
                             atol = eps_root*1e-1,
                             rtol = eps_root*1e-1,
                             maxiter = maxiter_root)

  root = tmp$root
  return(list(roots = matrix(root[1:(p+1)], nrow = 1)))
}

f_to_minimize = function(x, Sigma, U, rho, A) {

  # @param x A numeric vector representing the current solution or variable.
  # @param Sigma A numeric matrix representing the covariance matrix.
  # @param U A numeric vector representing the dual variable.
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param A A numeric vector representing the mean of the variables.

  -t(x) %*% Sigma %*% x + (rho/2) * t(x+U/rho-A) %*% (x+U/rho-A)
}

check_constraints = function(x, Xi, k, mu_1, w, eps = 1e-4, return_vals = FALSE){

  # Check Constraints for Optimization Solution
  # @param x A numeric vector representing the solution to be checked.
  # @param Xi A numeric matrix containing the principal components, or NULL if not applicable.
  # @param k An integer specifying the number of principal components.
  # @param mu_1 A numeric scalar representing the mu parameter.
  # @param w A numeric vector of projection vectors.
  # @param eps A numeric scalar specifying the tolerance for constraint checks (default is 1e-4).
  # @param return_vals A logical indicating whether to return detailed values of constraints (default is FALSE).
  # @return A logical value indicating whether all constraints are satisfied if `return_vals` is FALSE.
  #         If `return_vals` is TRUE, returns a list with:
  #         \itemize{
  #           \item \code{check}: A logical indicating whether all constraints are satisfied.
  #           \item \code{values}: A numeric vector of the constraint values.
  #         }

  check_all = TRUE

  const = c(t(x) %*% x, t(x)%*% w,  mu_1* (t(x)%*% w))

  if( abs(const[1] - 1) > eps) {
    check_all = FALSE
  }
  if( const[2] < -eps) {
    check_all = FALSE
  }
  if( mu_1 < -eps) {
    check_all = FALSE
  }
  if( abs(const[3]) > eps) {
    check_all = FALSE
  }
  if(k > 1){
    for( i in 1:(k-1)){
      const = c(const, t(x) %*% Xi[, i])
      if(abs(t(x) %*% Xi[, i]) > eps){
        check_all = FALSE
      }
    }
  }

  if(return_vals){
    return(list(check = check_all,
                values = const))
  }
  return(check_all)
}

find_minimum = function(Sigma, U, rho, A, k, Xi = NULL, w,
                        eps_root = 1e-2, maxiter_root = 50){


  # Find Minimum for Optimization Problem
  # @param Sigma A numeric matrix representing the covariance matrix for the current neighborhood.
  # @param U A numeric vector representing the dual variables.
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param A A numeric vector representing the mean of all As in the current ADMM iteration.
  # @param k An integer specifying the number of principal components or a similar parameter.
  # @param Xi A numeric matrix of available principal components (k-1 components), or NULL if not applicable.
  # @param w A numeric vector of projection vectors for uniqueness.
  # @param eps_root A numeric scalar specifying the tolerance for convergence in the root-finding algorithm (default is 1e-2).
  # @param maxiter_root An integer specifying the maximum number of iterations for the root-finding algorithm (default is 50).
  # @return A list with components:
  # \itemize{
  #   \item \code{min}: A numeric vector of the optimal solution.
  #   \item \code{roots}: A numeric vector of roots found.
  #   \item \code{mu}: The value of the mu parameter.
  #   \item \code{c_check}: A logical indicating whether the constraints are satisfied.
  #   \item \code{val}: The value of the objective function to minimize.
  # }

  p = dim(Sigma)[1]

  # find root
  roots_res = find_root(Sigma = Sigma,
                        U = U,
                        rho = rho,
                        A = A,
                        k = k,
                        Xi = Xi,
                        w = w,
                        starting_value = A,
                        mu_start = 20,
                        eps_root = eps_root,
                        maxiter_root = maxiter_root)


  # check for feasibility
  roots = roots_res$roots[1:p]
  mu = roots_res$roots[(p+1)]

  c_check_values = check_constraints(x = roots,
                                     k = k,
                                     Xi = Xi,
                                     eps = eps_root*1e1,
                                     mu_1 = mu,
                                     w = w,
                                     return_vals = TRUE)
  c_check = c_check_values$check

  # if not feasible print message and stop
  if(!c_check) {
    vec_out = c(roots_res$roots, c_check_values$values)
    if(k == 1) {
      names(vec_out) = c(paste0("A", 1:p), "mu", "norm(x)", "x'w","mu*(x'w)")
    } else {
      names(vec_out) = c(paste0("A", 1:p), "mu", "norm(x)", "x'w",  "mu*(x'w)", paste0("x'x_", 1:(k-1)))
    }
    print(vec_out)
    stop("Found root is not feasible regarding inequality constraint. Try to increase rho.")
  }

  f_vals = f_to_minimize(x = roots,
                         Sigma = Sigma,
                         U = U,
                         rho = rho,
                         A = A)

  return( list(min = roots,
               roots = roots,
               mu = mu,
               c_check = c_check,
               val = f_vals))
}

solve_minimization_PCA = function(Sigma, rho, U1, A_mean, k, Xi, w,
                                  eps_root = 1e-3, maxiter_root = 50){

  # Solve Minimization for PCA in ADMM
  # @param Sigma A list of covariance matrices, each corresponding to a neighborhood.
  # @param rho A numeric scalar representing the ADMM parameter.
  # @param U1 A numeric vector of dual variables used in ADMM.
  # @param A_mean A numeric vector of means from the ADMM iteration.
  # @param k An integer specifying which principal component to find.
  # @param Xi A matrix of available principal components (k-1 components), or NULL if k=1.
  # @param w A numeric vector of projection vectors for uniqueness.
  # @param eps_root A numeric scalar specifying the tolerance for convergence (default is 1e-3).
  # @param maxiter_root An integer specifying the maximum number of iterations for the root-finding algorithm (default is 50).

  N = length(Sigma)
  p = dim(Sigma[[1]])[1]
  if(!is.null(Xi)) Xi = as.matrix(Xi)

  optimal = rep(NA, N*p)
  val = 0

  # separable across neighborhoods
  for(j in 1:N){
    ind = jth_col(j = j, p = p)

    if(k == 1) {
      tmp = find_minimum(Sigma = Sigma[[j]],
                         U = U1[ind],
                         rho = rho,
                         A = A_mean[ind],
                         k = k,
                         Xi = NULL,
                         w = w[ind],
                         eps_root = eps_root,
                         maxiter_root = maxiter_root)
    } else {
      tmp = find_minimum(Sigma = Sigma[[j]],
                         U = U1[ind],
                         rho = rho,
                         A = A_mean[ind],
                         k = k,
                         Xi = as.matrix(Xi[ind, ], ncol = 1),
                         w = w[ind],
                         eps_root = eps_root,
                         maxiter_root = maxiter_root)
    }

    optimal[ind] = tmp$min
    val = val + tmp$val
  }

  return(list(A_optim = optimal,
              val = val))

}

check_constraints_all = function(x, Xi, k, p, w, eps = 1e-4){

  # Check Constraints for All Neighborhoods
  # @param x A numeric vector containing the concatenated results for all neighborhoods, where each neighborhood's vector is \(p\)-dimensional.
  # @param Xi A numeric matrix of prior principal components, where each column represents a principal component. The dimensions are \(p x (k-1)\).
  # @param k An integer representing the number of principal components considered.
  # @param p An integer representing the dimensionality of each neighborhood.
  # @param w A numeric vector representing the projection vectors for uniqueness.
  # @param eps A numeric scalar representing the tolerance for checking constraints. Defaults to \(10^{-4}\).
  # @return A list containing:
  #   \item{checks}{A logical value indicating whether all constraints are satisfied.}
  #   \item{values}{A matrix where each row corresponds to a neighborhood, and columns represent constraint values including normalization, orthogonality checks, and others. Column names vary based on the value of \(k\).}

  N = length(x)/p
  check_res = TRUE
  values = matrix(NA, nrow = N, ncol = (2 + k))

  for (i in 1:N){
    ind = jth_col(j = i, p = p)
    c = check_constraints(x = x[ind],
                          Xi = as.matrix(Xi[ind, ]),
                          k = k,
                          eps = eps,
                          mu_1 = 0,
                          w = w[ind],
                          return_vals = TRUE)
    check_res = check_res & c$check
    values[i, ] = c$values
  }

  coln = c("xx", "xw", "muxw")
  if(k > 1) coln = c(coln, paste0("xxi_", 1:(k-1)))
  colnames(values) = coln

  if(!check_res) {
  }
  return(list(checks =check_res,
              values = values) )
}

eval_objective = function(PC, eta, gamma, COVS){

  # Objective function value for local sparse PCA
  # @param PC vectorised component to evaluate.
  # @param eta degree of sparsity.
  # @param gamma distribution of sparsity.
  # @param COVS list of covariance matrices used for PCA

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



##########################################################################################
#### ADMM MAIN FUNCTION                                                               ####
##########################################################################################

#' @importFrom graphics abline legend text points
solve_ADMM = function(eta,
                      gamma,
                      Sigma,
                      k = 1,
                      Xi = NULL,
                      rho = NULL,
                      n_max = 100,
                      cor = FALSE,
                      eps_root = 1e-2,
                      eps_ADMM = 1e-4,
                      eps_threshold = NULL,
                      maxiter_root = 100,
                      convergence_plot = TRUE,
                      show_progress = TRUE){

  # ADMM Main function
  # @param eta scalar, non-negative
  # @param gamma scalar, between 0 and 1
  # @param Sigma list of covariance or correlation matrices
  # @param k integer bigger equal 1, number of component to calculate
  # @param Xi array or NULL, if k bigger than 1 contains prior components
  # @param rho positive scalar
  # @param n_max integer, number of maximal ADMM steps
  # @param cor logical, if correlation specific starting value should be used
  # @param eps_root positive small error for root finder
  # @param eps_ADMM positive small error for convergence of ADMM
  # @param eps_threshold positive small error for final threshold
  # @param maxiter_root integer, maximal number of steps for root finder
  # @param convergence_plot logical, if convergence plot should be plotted

  # extract information
  N = length(Sigma)
  p = dim(Sigma[[1]])[1]
  eps_threshold_given = eps_threshold


  # set parameters
  if(is.null(rho)) rho = sum(sapply(Sigma, diag))/N + eta
  A_mean = starting_value_ADMM(Sigma = Sigma,
                               eta = eta,
                               gamma = gamma,
                               k = k,
                               Xi = Xi,
                               cor = cor)
  w = A_mean


  # setup  convergence plots
  residual_plot = FALSE
  if((p > 10 | N > 10) & convergence_plot){
    convergence_plot = FALSE
    residual_plot = TRUE
  }
  if(convergence_plot){
    plot(x = 0:n_max,
         y = rep(NA, n_max+1),
         ylim=c(-1.25, 1.25),
         main = "Convergence (Entries)",
         xlab = "Iterations",
         ylab = "Values PC")
    text(0, 1.1,
         paste("eta:", round(eta,3),
               ", gamma:", round(gamma,3),
               "; rho:", round(rho, 3),
               "; k:", k),
         pos = 4)
    abline(0,0)
    abline(1,0, col = "grey")
    abline(-1,0, col = "grey")
    colors = c("blue", "orange", "red", "cyan", "green",
               "darkgreen", "magenta", "magenta4", "brown", "yellow4")
    shapes = c(0, 1, 2, 3, 4, 5, 6, 20, 16, 18)
    for(i in 1:(N*p)){
      ind2 = (i-1) %% p +1
      ind1 = 1 +  floor((i-1)/p)
      points(x = 0,
             y = A_mean[i],
             col = colors[ind1],
             pch = shapes[ind2])
    }
  }

  if(residual_plot){
    plot(x = 0:n_max,
         y = rep(NA, n_max+1),
         ylim = c(0, 10.25),
         main = "Convergence (Residuals)",
         xlab = "Iterations",
         ylab = "Residual values")
    text(0, 10.1,
         paste("eta:", round(eta, 3),
               ", gamma:", round(gamma, 3),
               "; rho:", round(rho, 3),
               "; k:", k),
         pos = 4)
  }

  # show progress only when interactive
  if(!interactive()) show_progress = FALSE



  # initialize vectors
  U1 = rep(0, N*p)
  U2 = rep(0, N*p)
  U3 = rep(0, N*p)

  A1 = rep(0, N*p)
  A2 = rep(0, N*p)
  A3 = rep(0, N*p)

  # setup iteration
  A_mean_old = A_mean
  norm_residual_dual = c()
  norm_residual_prime = c()
  value_objective = rep(NA, n_max)


  # make progress bar
  progress_bar <- function(i, total = n_max, prefix = paste("PC", k), bar_length = 50) {
    filled <- floor(i / total * bar_length)
    bar <- paste0(paste0(rep("=", filled), collapse = ""),  paste0(rep(" ", bar_length - filled), collapse = ""))
    cat(sprintf("\r%s |%s| %d iterations", prefix, bar, i))
    flush.console()
  }


  # iterate between problems
  for( i in 1:n_max){
    if(show_progress) progress_bar(i)

    # PCA subproblem
    A1 = solve_minimization_PCA(Sigma = Sigma,
                                rho = rho,
                                U1 = U1,
                                A_mean = A_mean,
                                k = k,
                                Xi = Xi,
                                w = w,
                                eps_root = eps_root,
                                maxiter_root = maxiter_root)$A_optim


    # soft thresholding - scalar
    A2 = solve_minimization_scalar_softthreshold(eta = eta,
                                                 gamma = gamma,
                                                 rho = rho,
                                                 U2 = U2,
                                                 A_mean = A_mean)


    # soft thresholding- group
    A3 = solve_minimization_group_softthreshold(eta = eta,
                                                gamma = gamma,
                                                rho = rho,
                                                U3 = U3,
                                                A_mean = A_mean,
                                                N = N)


    A_mean = (1/3) * (A1 + A2 + A3)  + (1/(3*rho)) * (U1 + U2 + U3)
    if(gamma == 1)  {
      A_mean = (1/2) * (A1 + A2)  + (1/(2*rho)) * (U1 + U2)
    }


    # project to feasible space
    if( k > 1){
      A_mean = project_to_orthogonal(PC = A_mean, p = p, N = N, Xi = Xi, renorm = TRUE)
    } else {
      A_mean = renorm(A_mean, p = p, N = N)
    }


    # update dual variables
    U1 = U1 + rho*(A1 - A_mean)
    U2 = U2 + rho*(A2 - A_mean)
    U3 = U3 + rho*(A3 - A_mean)


    # calculate metrics for convergence (page 51 Boyd)
    if(gamma < 1){
      residual_prime = norm(as.matrix(A1 - A_mean), "F")^2 +
        norm(as.matrix(A2 - A_mean), "F")^2 +
        norm(as.matrix(A3 - A_mean), "F")^2
      residual_dual = rho^2 * 3 * norm(as.matrix(A_mean_old - A_mean), "F")^2
    }
    if(gamma == 1){
      residual_prime = norm(as.matrix(A1 - A_mean), "F")^2 + norm(as.matrix(A2 - A_mean), "F")^2
      residual_dual = rho^2 * 2 * norm(as.matrix(A_mean_old - A_mean), "F")^2
    }
    eps_prime = sqrt(N*p)*eps_ADMM + eps_ADMM * max(c(norm(as.matrix(A1), "F"),
                                                      norm(as.matrix(A2), "F"),
                                                      norm(as.matrix(A3), "F"),
                                                      norm(as.matrix(A_mean), "F")))
    eps_dual = sqrt(N*p)*eps_ADMM + eps_ADMM * max(c(norm(as.matrix(U1), "F"),
                                                     norm(as.matrix(U2), "F"),
                                                     norm(as.matrix(U3), "F")))


    # save residuals and threshold
    norm_residual_dual = c(norm_residual_dual, residual_dual)
    norm_residual_prime = c(norm_residual_prime, residual_prime)
    if(is.null(eps_threshold_given)) {
      eps_threshold = max(abs(A_mean_old - A_mean))
    } else {eps_threshold = eps_threshold_given}

    # calculate value of objective function
    value_objective[i] = eval_objective(PC = A_mean,
                                        eta = eta,
                                        gamma = gamma,
                                        COVS = Sigma)


    # update plots
    if(convergence_plot){
      for(j in 1:(N*p)){
        ind2 = (j-1) %% p +1
        ind1 = 1 + floor((j-1)/p)
        points(x = i,
               y = A_mean[j],
               col = colors[ind1],
               pch = shapes[ind2])
      }
    }
    if(residual_plot){
      points(x = i,
             y = norm_residual_dual[i],
             col = "lightblue",
             pch = 3)
      points(x = i,
             y = norm_residual_prime[i],
             col = "darkblue",
             pch = 4)
      legend("topright",
             legend = c("Primal residual", "Dual residual"),
             col = c("darkblue", "lightblue"),
             lty = c(1,1))
    }


    # check for convergence
    if(residual_dual < eps_dual  & residual_prime < eps_prime) {
      break
    } else {
      A_mean_old = A_mean
    }
  }
  if(show_progress) cat("\n")

  # thresholding and projection to feasible space
  if(k > 1){
    A_mean = project_to_orthogonal(PC = A_mean, Xi = Xi, p = p, N = N, renorm = TRUE)
  }
  A_mean[abs(A_mean) < eps_threshold] = 0
  A_mean = renorm(x = A_mean, p = p, N = N)


  # message for convergence
  if (i == n_max){
    warning(paste0("\nAlgorithm did not converge - maximal number of iterations (n_max = ", n_max ,") reached! ",
                   "Primal residual: ", round(residual_prime, 7), " (adaptive threshold: ", round(eps_prime, 7), "), ",
                   "dual residual :", round(residual_dual, 7), " (adaptive threshold: ", round(eps_dual, 7), ").\n"
    ))
  }

  return(list(PC = A_mean,
              converged = (i != n_max),
              n_steps = i,
              starting_value = w,
              value_objective = value_objective[1:i],
              residuals = cbind(norm_residual_prime[i], norm_residual_dual[i]),
              eps_threshold = eps_threshold))
}




##########################################################################################
# ADJUSTED ETA                                                                     ####
##########################################################################################

adjusted_eta = function(COVS, Xi, k, eta){

  # Adjust Eta Value Based on Principal Components and Covariance Matrices
  # @param COVS A list of covariance matrices (one for each neighborhood).
  # @param Xi A matrix of principal components (PCs) calculated. Dimensions are \( p x N x k \) where \( p \) is the dimension, \( N \) is the number of neighborhoods, and \( k \) is the number of principal components.
  # @param k An integer indicating the index of the principal component being calculated.
  # @param eta A numeric value provided by the user to be adjusted.
  # @return A numeric value representing the adjusted `eta`.

  N = length(COVS)
  p = dim(COVS[[1]])[1]
  if(k == 1) {
    return(eta)
  } else {
    Xi = as.matrix(Xi, ncol = k-1)
    v = 0
    for(i in 1:N){
      Xii = Xi[jth_col(j = i, p = p), 1:(k-1), drop = FALSE]
      P = diag(1, nrow = p) - Xii %*% t(Xii)
      S = Re(eigen(P %*% COVS[[i]] %*% t(P))$values[1])
      v = v + S
    }
  }
  return(eta * v/sum(sapply(COVS, function(x) eigen(x)$values[1])) )
}

##########################################################################################
# STRUCTURAL                                                                          ####
##########################################################################################
jth_col = function(j, p){

  # This function generates a vector of column indices for the \( j \)-th neighborhood in a matrix of dimension \( p \).
  # @param j An integer representing the neighborhood index.
  # @param p An integer representing the dimension of each neighborhood.
  # @return A numeric vector containing the column indices for the \( j \)-th neighborhood.

  seq( (j-1)*p + 1, j*p, by = 1)
}

jth_row = function(j, p, N){

  # This function generates a vector of row indices for the \( j \)-th variable across all neighborhoods.
  # @param j An integer representing the variable index.
  # @param p An integer representing the dimension of each neighborhood.
  # @param N An integer representing the total number of neighborhoods.
  # @return A numeric vector containing the row indices for the \( j \)-th variable.

  j + p*(0:(N-1))
}

renorm = function(x, p, N){

  # This function renormalizes vectors within each neighborhood based on the Frobenius norm.
  # @param x A numeric vector to be renormalized.
  # @param p An integer representing the dimension of each neighborhood.
  # @param N An integer representing the total number of neighborhoods.
  # @return A numeric vector where each neighborhood is renormalized based on its Frobenius norm.

  for (j in 1:N){
    ind = jth_col(j = j, p = p)
    x[ind] = x[ind] / norm(matrix(x[ind]), "F")
  }
  return(x)
}

project_to_orthogonal = function(PC, Xi, N, p, renorm = FALSE){

  # This function projects a vector `PC` onto the orthogonal space of a given matrix `Xi` and optionally renormalizes the result (per neighborhood).
  # @param PC A numeric vector of length \( N x p \) representing the principal components.
  # @param Xi A numeric matrix or vector representing the basis for orthogonal projection.
  # @param N An integer representing the number of neighborhoods.
  # @param p An integer representing the dimension of each neighborhood.
  # @param renorm A logical value indicating whether to renormalize the projected vector (default is `FALSE`).
  # @return A numeric vector of length \( N x p \) representing the projection of `PC` onto the orthogonal space of `Xi`.

  if(is.null(dim(Xi))) Xi = matrix(Xi, ncol = 1)
  k = dim(Xi)[2]
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
