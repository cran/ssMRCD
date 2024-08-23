# INCLUDES:
#   - ADMM main function
#   - starting value calculation
#   - subroutines

##########################################################################################
# ADMM Main function

# @param lambda scalar, non-negative
# @param alpha scalar, between 0 and 1
# @param Sigma list of covariance or correlation matrices
# @param k integer bigger equal 1, number of component to calculate
# @param Xi matrix or NULL, if k bigger than 1 contains prior components
# @param rho positive scalar
# @param n_max integer, number of maximal ADMM steps
# @param cor logical, if correlation specific starting value should be used
# @param starting_value vector, if non-default starting value should be used
# @param eps_root positive small error for root finder
# @param eps_ADMM positive small error for convergence of ADMM
# @param eps_threshold positive small error for final threshold
# @param maxiter_root integer, maximal number of steps for root finder
# @param trace logical, if additional information should be printed
# @param convergence_plot logical, if convergence plot should be plotted

#' @importFrom graphics abline legend text
#' @importFrom utils setTxtProgressBar txtProgressBar
solve_ADMM = function(lambda,
                      alpha,
                      Sigma,
                      k = 1,
                      Xi = NULL,
                      rho = NULL,
                      n_max = 100,
                      cor = FALSE,
                      starting_value = NULL,
                      eps_root = 1e-2,
                      eps_ADMM = 1e-4,
                      eps_threshold = NULL,
                      maxiter_root = 100,
                      trace = TRUE,
                      convergence_plot = TRUE){

  # extract information
  N = length(Sigma)
  p = dim(Sigma[[1]])[1]
  eps_threshold_given = eps_threshold

  if(trace) coloured_print(paste0("Start ADMM for k-th = ", k, " principal component."),
                           colour = "message",
                           trace = trace)


  # set parameters
  if(!is.null(Xi)) Xi = as.matrix(Xi)
  if(is.null(rho)) rho = sum(sapply(Sigma, diag))/N + lambda

  if(is.null(starting_value)) {
    coloured_print("Starting value not given, resort to default.",
                   colour = "message",
                   trace = trace)
    A_mean = starting_value_ADMM(Sigma = Sigma,
                                 lambda = lambda,
                                 alpha = alpha,
                                 k = k,
                                 Xi = Xi,
                                 cor = cor)
  } else if(is.numeric(starting_value)){
    A_mean = c(starting_value)
  } else warning("Given starting value not valid. Choose numeric vector of size N*p or NULL.")

  w = A_mean

  # setup  convergence plots
  residual_plot = FALSE
  if((p > 10 | N > 10) & convergence_plot){
    convergence_plot = FALSE
    residual_plot = TRUE
    coloured_print("Entrywise convergence plot not possible for p > 10 or N > 10. Instead, norm of residuals will be plotted for convergence.",
                   colour = "message",
                   trace = trace)
  }
  if(convergence_plot){
    plot(x = 0:n_max,
         y = rep(NA, n_max+1),
         ylim=c(-1.25, 1.25),
         main = "Convergence (Entries)",
         xlab = "Iterations",
         ylab = "Values PC")
    text(0, 1.1,
         paste("lambda:", round(lambda,3),
               ", alpha:", round(alpha,3),
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
         paste("lambda:", round(lambda, 3),
               ", alpha:", round(alpha, 3),
               "; rho:", round(rho, 3),
               "; k:", k),
         pos = 4)
  }

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

  # iterate between problems
  if(trace) pb = txtProgressBar(min = 0, max = n_max, initial = 0)
  for( i in 1:n_max){
    if(trace) setTxtProgressBar(pb,i)

    # PCA subproblem
    A1 = solve_minimization_PCA(Sigma = Sigma,
                                rho = rho,
                                U1 = U1,
                                A_mean = A_mean,
                                k = k,
                                Xi = Xi,
                                w = w,
                                eps_root = eps_root,
                                maxiter_root = maxiter_root,
                                trace = trace)$A_optim

    # soft thresholding - scalar
    A2 = solve_minimization_scalar_softthreshold(lambda = lambda,
                                                 alpha = alpha,
                                                 rho = rho,
                                                 U2 = U2,
                                                 A_mean = A_mean)

    # soft thresholding - group
    A3 = solve_minimization_group_softthreshold(lambda = lambda,
                                                alpha = alpha,
                                                rho = rho,
                                                U3 = U3,
                                                A_mean = A_mean,
                                                N = N)


    A_mean = (1/3) * (A1 + A2 + A3)  + (1/(3*rho)) * (U1 + U2 + U3)
    if(alpha == 1)  {
      A_mean = (1/2) * (A1 + A2)  + (1/(2*rho)) * (U1 + U2)
      #A3 = rep(0, N*p)
      #U3 = rep(0, N*p)
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
    if(alpha < 1){
      residual_prime = norm(as.matrix(A1 - A_mean), "F")^2 +
        norm(as.matrix(A2 - A_mean), "F")^2 +
        norm(as.matrix(A3 - A_mean), "F")^2
      residual_dual = rho^2 * 3 * norm(as.matrix(A_mean_old - A_mean), "F")^2
    }
    if(alpha == 1){
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
                                        eta = lambda,
                                        gamma = alpha,
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
      coloured_print(paste(i, " many iteration steps"),
                     colour = "message",
                     trace = trace)
      break
    } else {
      A_mean_old = A_mean
    }
  }

  if(trace) {
    close(pb)
    coloured_print(paste0("Thresholding applied to entries: ", round(eps_threshold, 8)),
                   colour = "message",
                   trace = trace)
  }


  # thresholding and projection to feasible space
  if(k > 1){
    A_mean = project_to_orthogonal(PC = A_mean,
                                   Xi = Xi,
                                   p = p,
                                   N = N,
                                   renorm = TRUE)
  }
  A_mean[abs(A_mean) < eps_threshold] = 0
  A_mean = renorm(x = A_mean, p = p, N = N)


  # message for convergence
  if (i == n_max){
    coloured_print(paste0("\nAlgorithm did not converge - maximal iteration number (n_max = ", n_max ,") reached!\n",
                          "Primal residual is ", round(residual_prime, 8), " (adaptive threshold: ", eps_prime, "). \n",
                          "Dual residual is ", round(residual_dual, 8), " (adaptive threshold: ", eps_dual, "). \n"
    ),
    colour = "problem",
    trace = trace)
  }

  # check_constraints
  summary = summary_diagnostic(A_mean,
                               Xi = Xi,
                               k = k,
                               p = p,
                               w = w,
                               COVS = Sigma,
                               trace = trace,
                               eps = eps_ADMM)
  if(trace) print(summary)

  return(list(PC = A_mean,
              converged = (i != n_max),
              n_steps = i,
              summary = summary,
              starting_value = w,
              value_objective = value_objective[1:i],
              residuals = cbind(norm_residual_prime[i],
                                norm_residual_dual[i]),
              eps_threshold = eps_threshold))
}




##########################################################################################
# STARTING VALUE
solveL1_cov = function(Sigma, k){
  # calculate extreme solutions for L1 penalty and covariance matrix

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
  # calculate extreme solutions for L1 penalty and correlation matrix

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
  # calculate extreme solutions for L2 penalty (group) and covariance matrix

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
  # calculate extreme solutions for L2 penalty (group) and correlation matrix

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

starting_value_ADMM = function(Sigma,
                               lambda,
                               alpha,
                               k = 1,
                               return_all = FALSE,
                               Xi = NULL,
                               cor = FALSE){

  # get starting value for ADMM
  N = length(Sigma)
  p = dim(Sigma[[1]])[1]
  var_k = rep(NA, N*p)

  # eigenvalues
  eigenvec_start_k = rep(NA, N*p)
  for(i in 1:N){
    tmp = eigen(Sigma[[i]])$vectors[,k]
    eigenvec_start_k[jth_col(j = i, p = p)] = tmp
    var_k[jth_col(j = i, p = p)] = diag(Sigma[[i]])
  }

  if(!cor) {
    groupwise_start_k = solveL2_cov(Sigma = Sigma, k = k)
    entrywise_start_k = solveL1_cov(Sigma = Sigma, k = k)
  }
  if(cor) {
    groupwise_start_k = solveL2_cor(Sigma = Sigma, k = k)
    entrywise_start_k = solveL1_cor(Sigma = Sigma, k = k)
  }

  # average starting values
  if(alpha == 1)  penalty_start = entrywise_start_k
  if(alpha != 1) penalty_start = groupwise_start_k

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
  if(lambda == 0) starting_value = eigenvec_start_k

  # return
  if(return_all){
    return(list(starting_value = starting_value,
                penalty_solution = penalty_start,
                eigenvector = eigenvec_start_k))
  }

  return(starting_value)
}


##########################################################################################
### solve L1-norm
# see Appendix B4 in Tag Lasso of Ines: solution is soft thresholding

# solve minimization problem 2 (L1 Norm)
solve_minimization_scalar_softthreshold = function(lambda,
                                                   alpha,
                                                   rho,
                                                   U2,
                                                   A_mean){

  # lambda, alpha: sparsity parameters
  # rho: ADMM parameter
  # U2: dual variable
  # A_mean: mean of all three As, where we want to find minimum

  p = length(A_mean)
  A_res = rep(NA, p)
  kappa = (lambda*alpha)/rho
  a = A_mean - U2/rho

  for(i in 1:p) {
    A_res[i] = soft_thresholding_scalar(a = a[i],
                                        kappa = kappa)
  }

  return(A_res)
}

# soft thresholding for single values
soft_thresholding_scalar = function(a, kappa){

  if(a > kappa)
    return(a-kappa)
  if(a < -kappa)
    return(a+kappa)
  if(abs(a) <= kappa)
    return(0)
}


##########################################################################################
### solve group-wise norm
solve_minimization_group_softthreshold = function(lambda,
                                                  alpha,
                                                  rho,
                                                  U3,
                                                  A_mean,
                                                  N){

  # lambda, alpha: sparsity parameters
  # rho: ADMM parameter
  # U3: dual variable
  # A_mean: mean of all three As, where we want to find minimum

  lambda = lambda * sqrt(N)
  kappa = lambda*(1-alpha)/rho
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

# soft thresholding for singular group
soft_thresholding_group = function(a, kappa){
  # a: vector of observations of one group
  # kappa: thresholding parameter

  a = as.matrix(a)
  Nor = norm(a, type = "F")
  if(Nor == 0) {
    return (a)
  } else{
    mult = max(1 - kappa/Nor, 0)
    return(a * mult)
  }
}



##########################################################################################
### solve PCA

# solve for all neighborhoods
solve_minimization_PCA = function(Sigma,
                                  rho,
                                  U1,
                                  A_mean,
                                  k,
                                  Xi,
                                  w,
                                  eps_root = 1e-3,
                                  maxiter_root = 50,
                                  trace = TRUE){

  # Sigma: list covariance matrices
  # rho: parameter ADMM
  # U1: dual variable
  # A_mean: from ADMM
  # k: number of principal components
  # Xi: available principal components (k-1 many)
  # w: projection vector for uniqueness

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
                         maxiter_root = maxiter_root,
                         trace = trace)
    } else {
      tmp = find_minimum(Sigma = Sigma[[j]],
                         U = U1[ind],
                         rho = rho,
                         A = A_mean[ind],
                         k = k,
                         Xi = as.matrix(Xi[ind, ], ncol = 1),
                         w = w[ind],
                         eps_root = eps_root,
                         maxiter_root = maxiter_root,
                         trace = trace)
    }

    optimal[ind] = tmp$min
    val = val + tmp$val
  }

  return(list(A_optim = optimal,
              val = val))

}



# solve minimization for single neighborhood
find_minimum = function(Sigma,
                        U,
                        rho,
                        A,
                        k,
                        Xi = NULL,
                        w,
                        eps_root = 1e-2,
                        trace = FALSE,
                        maxiter_root = 50){



  p = dim(Sigma)[1]

  # find root
  roots_res = find_root(Sigma = Sigma,
                        U = U,
                        rho = rho,
                        A = A,
                        k = k,
                        Xi = Xi,
                        w= w,
                        starting_value = A,
                        mu_start = 20,
                        eps_root = eps_root,
                        maxiter_root = maxiter_root,
                        trace = trace)


  # check for feasibility
  roots = roots_res$roots[1:p]
  mu = roots_res$roots[(p+1)]

  c_check_values = check_constraints(x = roots,
                                     k = k,
                                     Xi = Xi,
                                     eps = eps_root*1e1,
                                     trace = trace,
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

    stop("Found root is not feasible regarding inequality constraint. Try higher rho.")
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

#' @importFrom rootSolve multiroot
find_root = function(Sigma,
                     U,
                     rho,
                     A,
                     k,
                     Xi,
                     w,
                     starting_value = NULL,
                     mu_start = 20,
                     eps_root = 1e-2,
                     maxiter_root = 100,
                     trace = TRUE){

  # lambda_j_calc: not really faster
  # starting value will be re-projected and -normed if k > 1

  p = dim(Sigma)[1]
  f = function(x_mu_lambda) f_find_root_of(x_mu_lambda = x_mu_lambda,
                                           Sigma = Sigma,
                                           U = U,
                                           rho = rho,
                                           A = A,
                                           k = k,
                                           Xi = Xi,
                                           w= w)$res

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


  tmp = rootSolve::multiroot(f = f,
                             start = c(starting_value),
                             useFortran = FALSE,
                             atol = eps_root*1e-1,
                             rtol = eps_root*1e-1,
                             maxiter = maxiter_root)

  root = tmp$root
  if(max(abs(tmp$f.root)) > 0.1){
    coloured_print(paste0("Root finder stopped at non-root: max(abs(f(x*)))=", max(abs(tmp$f.root)), "!"),
                   colour = "problem",
                   trace = trace)
  }
  return(list(roots = matrix(root[1:(p+1)], nrow = 1)))
}


check_constraints = function(x,
                             Xi,
                             k,
                             mu_1,
                             w,
                             eps = 1e-4,
                             trace = TRUE,
                             return_vals = FALSE){

  check_all = TRUE

  const = c(t(x) %*% x, t(x)%*% w,  mu_1* (t(x)%*% w))

  if( abs(const[1] - 1) > eps) {
    check_all = FALSE
    coloured_print(paste("Constraint Issue: x not normalized. Norm:", const[1]),
                   colour = "problem",
                   trace = trace)
  }
  if( const[2] < -eps) {
    check_all = FALSE
    coloured_print(paste("Constraint Issue: Sign of x'w: ", const[2]),
                   colour = "problem",
                   trace = trace)
  }
  if( mu_1 < -eps) {
    check_all = FALSE
    coloured_print(paste("Constraint Issue: Sign of mu: ", mu_1),
                   colour = "problem",
                   trace = trace)
  }
  if( abs(const[3]) > eps) {
    check_all = FALSE
    coloured_print(paste("Constraint Issue: Value of mu*x'w: ", abs(const[3])),
                   colour = "problem",
                   trace = trace)
  }
  if(k > 1){
    for( i in 1:(k-1)){
      const = c(const, t(x) %*% Xi[, i])
      if(abs(t(x) %*% Xi[, i]) > eps){
        check_all = FALSE
        coloured_print(paste("Constraint Issue: x not orthogonal to", i,
                             "th principal component. Scalar Product:", t(x) %*% Xi[, i]),
                       colour = "problem",
                       trace = trace)
      }
    }
  }

  if(return_vals){
    return(list(check = check_all,
                values = const))
  }
  return(check_all)
}

# objective function
f_to_minimize = function(x, Sigma, U, rho, A) {
  -t(x) %*% Sigma %*% x + (rho/2) * t(x+U/rho-A) %*% (x+U/rho-A)
}

f_find_root_of = function(x_mu_lambda,
                          Sigma,
                          U,
                          rho,
                          A,
                          k,
                          w,
                          Xi = NULL){

  # Sigma: cov matrix for one neighborhood
  # U: U3 vector
  # A: A_k vector
  # k: number of principle component looked for
  # Xi: matrix with prior PC (dim = p x (k-1))

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



# for summary
check_constraints_all = function(x,
                                 Xi,
                                 k,
                                 p,
                                 w,
                                 eps = 1e-4,
                                 trace = TRUE){

  N = length(x)/p
  check_res = TRUE
  values = matrix(NA, nrow = N, ncol = (2 + k))

  for (i in 1:N){
    ind = jth_col(j = i, p = p)
    c = check_constraints(x = x[ind],
                          Xi = as.matrix(Xi[ind, ]),
                          k = k,
                          eps = eps,
                          trace = trace,
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
    coloured_print(paste0("Result does not fulfill constraints of minimization problem with tolerance ", eps, " !"),
                   colour = "problem",
                   trace = trace)
  }
  return(list(checks =check_res,
              values = values) )
}


# summary for one set of PCs
summary_diagnostic = function(x,
                              Xi,
                              k,
                              p,
                              w,
                              COVS,
                              trace = FALSE,
                              eps = 0){

  N = length(COVS)

  # constraints
  constraints = check_constraints_all(x = x,
                                      Xi = Xi,
                                      k = k,
                                      p = p,
                                      w = w,
                                      eps = eps,
                                      trace = trace)$values
  constraints = constraints[, colnames(constraints) != "muxw"]

  # variance
  expl_var = explained_var_N(x, COVS)
  glob_var = colSums(expl_var)
  glob_var[1] = glob_var[2] / glob_var[3]

  # sparsity
  sparsity = sparsity_summary(PC = x,
                              N = N,
                              p =p,
                              scaled = FALSE)
  colnames(sparsity) = "sparsity (#0)"
  glob_sparsity = sum(sparsity)

  # summary
  diagnostic_matrix = cbind(constraints, expl_var, sparsity)
  diagnostic_matrix = rbind(diagnostic_matrix,
                            c(rep(NA, dim(constraints)[2]), glob_var, glob_sparsity ))
  rownames(diagnostic_matrix) = c(paste0("N", 1:N), "global")

  return(diagnostic_matrix)
}


