##########################################################################################
#### MAIN USER FUNCTION                                                               ####
##########################################################################################
#' Compute Sparse Multi-Source Principal Components
#'
#' Estimates sparse principal components from multiple covariance or correlation matrices
#' using an ADMM-based optimization routine (see Puchhammer, Wilms and Filzmoser, 2024).
#'
#' @param eta Numeric or numeric vector. Controls the overall sparsity level. If a single value is provided, it will be used directly.
#' If a vector is given, the optimal value will be selected via internal model selection.
#' @param gamma Numeric or numeric vector. Controls the distribution of sparsity across components.
#' If a single value is provided, the optimal \code{eta} is selected automatically.
#' @param COVS A list of covariance or correlation matrices (one per data source or group).
#' @param k Integer. Number of principal components to compute. If not specified, all components are estimated.
#' @param adjust_eta Logical. If \code{TRUE} (default), the sparsity parameter \code{eta} is adjusted based on the variance structure.
#' @param convergence_plot Logical. If \code{TRUE}, a convergence diagnostic plot is displayed of either the residuals or the loading entries (default: \code{FALSE}).
#' @param rho List of parameters controlling the ADMM penalty parameter \code{rho}, with the following elements:
#' \describe{
#'   \item{1}{Initial value for \code{rho} (default: \code{NA}).}
#'   \item{2}{Logical; whether to increase \code{rho} if convergence is not reached (default: \code{TRUE}).}
#'   \item{3}{Maximum value for \code{rho} (default: 100). You may need to increase this for high-dimensional problems.}
#'   \item{4}{Step size for increasing \code{rho} (default: 1).}
#' }
#' @param eps Numeric vector of tolerance parameters used in optimization. Includes:
#' \describe{
#'   \item{1}{Tolerance for soft-thresholding (default: \code{1e-5}).}
#'   \item{2}{Tolerance for ADMM convergence (default: \code{1e-4}).}
#'   \item{3}{Tolerance for convergence of the internal root-finding step (default: \code{1e-1}).}
#'   \item{4}{Maximum number of iterations for the root finder (default: 50).}
#' }
#' @param n_max Integer. Maximum number of ADMM iterations (default: 200).
#' @param show_progress Logical. Indicates whether progress bars should be displayed.
#'
#' @return An object of class \code{"msPCA"} containing the following elements:\tabular{ll}{
#'    \code{PC} \tab Array of dimension p x k x N of loading vectors.  \cr
#'    \tab \cr
#'    \code{p} \tab Number of variables. \cr
#'    \tab \cr
#'    \code{N} \tab Number of neighborhoods. \cr
#'    \tab \cr
#'    \code{k} \tab Number of components. \cr
#'    \tab \cr
#'    \code{COVS} \tab List of covariance matrices sorted by neighborhood. \cr
#'    \tab \cr
#'    \code{gamma} \tab Sparsity distribution. \cr
#'    \tab \cr
#'    \code{eta} \tab Amount of sparsity. \cr
#'    \tab \cr
#'    \code{converged} \tab Logical, if ADMM converged with given specifications. \cr
#'    \tab \cr
#'    \code{n_steps} \tab Number of steps used. \cr
#'    \tab \cr
#'    \code{residuals} \tab Primary and secondary residuals. \cr
#'    \tab \cr
#' }
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}},  \code{\link[ssMRCD]{plot.msPCA}},
#' \code{\link[ssMRCD]{summary.msPCA}}, \code{\link[ssMRCD]{biplot.msPCA}},
#' \code{\link[ssMRCD]{screeplot.msPCA}}, \code{\link[ssMRCD]{scores}}, \code{\link[ssMRCD]{align}}
#'
#' @references Puchhammer, P., Wilms, I., & Filzmoser, P. (2024). Sparse outlier-robust PCA for multi-source data. *ArXiv Preprint*.  \doi{10.48550/arXiv.2407.16299}
#'
#'
#' @export
#'
#' @examples
#'
#' C1 = diag(c(1.1, 0.9, 0.6, 0.5, 2))
#' C2 = matrix(runif(25, -1, 1), 5, 5)
#' C2 = t(C2) %*% C2
#' C3 = matrix(runif(25, -1, 1), 5, 5)
#' C3 = t(C3) %*% C3
#'
#' pca1 = msPCA(eta = 1, gamma = 0.5, COVS = list(C1, C2, C3), k = 3,
#'              n_max = 100, rho = list(NA, TRUE, 100, 1), show_progress = FALSE)
#' summary(pca1)
#'
#' pca2 = msPCA(eta = seq(0, 3, 0.25), gamma = 1, COVS = list(C1, C2, C3), k = 3,
#'              n_max = 100, rho = list(NA, TRUE, 100, 1), show_progress = FALSE)
#' summary(pca2)

msPCA = function(eta,
                 gamma,
                 COVS,
                 k = ncol(COVS[[1]]),
                 adjust_eta = TRUE,
                 convergence_plot = FALSE,
                 n_max = 200,
                 rho = list(NA, TRUE, 100, 1),
                 eps = c(1e-5, 1e-4, 1e-1, 50),
                 show_progress = FALSE){


  # input checks
  ssmrcd_input = FALSE
  tuning = FALSE

  # eta, gamma: numeric
  if (!is.numeric(eta)) {
    stop("'eta' must be numeric.")
  }
  if (!is.numeric(gamma)) {
    stop("'gamma' must be numeric value.")
  }

  # COVS: list of matrices
  gnames = NULL
  if("ssMRCD" %in% class(COVS)){
    ssmrcd_input = TRUE
    ssmrcd = COVS
    COVS = ssmrcd$MRCDcov
    p = ncol(COVS[[1]])
    gnames = ssmrcd$gnames
  } else if (!is.list(COVS) | !all(sapply(COVS, is.matrix))) {
    stop("'COVS' must be a list of matrices or an ssMRCD object.")
  }
  if(!"ssMRCD" %in% class(COVS)){
    p = ncol(COVS[[1]])
    gnames = 1:length(COVS)
    if(!is.null(names(COVS)))  gnames = names(COVS)
  }

  cnames = 1:p
  if(!is.null(colnames(COVS[[1]]))) cnames = colnames(COVS[[1]])



  # rho: list
  if (!is.list(rho)) {
    stop("'increase_rho' must be a list.")
  }
  # k: positive integer
  if (!is.numeric(k) | length(k) != 1 | k <= 0 | k != floor(k)) {
    stop("'k' must be a positive integer.")
  }
  # n_max: positive integer
  if (!is.numeric(n_max) | length(n_max) != 1 | n_max <= 0 | n_max != floor(n_max)) {
    stop("'n_max' must be a positive integer.")
  }

  # convergence_plot: logical scalar
  if (!is.logical(convergence_plot) | length(convergence_plot) != 1) {
    stop("'convergence_plot' must be a logical (TRUE/FALSE).")
  }
  # adjust_eta: logical scalar
  if (!is.logical(adjust_eta) | length(adjust_eta) != 1) {
    stop("'adjust_eta' must be a logical (TRUE/FALSE).")
  }


  # get additional information
  cor = all(sapply(COVS, function(x) diag(x) == 1))
  eps_threshold = eps[1]
  eps_ADMM = eps[2]
  eps_root = eps[3]
  maxiter_root = eps[4]
  rho_init = rho[[1]]
  increase_rho = rho[-1]
  if(length(eta) > 1) tuning = TRUE


  # parameter tuning
  tuning_add = NULL
  if(tuning){
    if(length(gamma) != 1 & length(eta) == 1) stop("No parameter tuning can be achieved for single values of 'eta'.")

    out_tuning = select_sparsity (COVS = COVS,
                               rho = rho_init,
                               cor = cor,
                               eta = eta,
                               gamma = gamma,
                               eps_threshold = eps[1],
                               eps_root = eps[3],
                               eps_ADMM = eps[2],
                               n_max = n_max,
                               increase_rho = increase_rho,
                               convergence_plot = FALSE,
                               plot = TRUE)

    eta = out_tuning$eta
    gamma = out_tuning$gamma
    tuning_add = out_tuning[c("tpo", "auc", "eval", "plot_auc", "plot_tpo")]
  }

  p = ncol(COVS[[1]])
  N = length(COVS)

  # calculate pca to be returned
  out = sparsePCAloc(eta = eta,
                     adjust_eta = adjust_eta,
                     gamma = gamma,
                     COVS = COVS,
                     cor = cor,
                     k = k,
                     n_max = n_max,
                     convergence_plot = convergence_plot,
                     eps_threshold = eps[1],
                     eps_ADMM = eps[2],
                     eps_root = eps[3],
                     maxiter_root = eps[4],
                     rho = rho_init,
                     increase_rho = increase_rho,
                     show_progress = show_progress)

  # add tuning info
  if(tuning) out = c(out, tuning_add)

  # transform PCs to array
  PC = array(out$PC, dim = c(p, N, k), dimnames = list(cnames, gnames, paste0("PC", 1:k)))
  PC = aperm(PC, c(1, 3, 2))

  # name everything
  out$PC = PC
  names(out$converged) = names(out$n_steps) = paste0("PC", 1:k)
  out$residuals = do.call(rbind,out$residuals)
  colnames(out$residuals) = c("primal", "dual")
  rownames(out$residuals) =  paste0("PC", 1:k)
  out$COVS = COVS
  class(out) <- c("msPCA", "list")

  return(out)

}



##########################################################################################
#### INTERNAL HELPER FUNCTIONS                                                        ####
##########################################################################################

sparsity_entries = function(PC, N, p, tolerance = 0, k = 1, scaled = TRUE){

  # Entry-wise Sparsity in the Loadings
  # @param PC matrix-like object of PCs.
  # @param N integer, number of groups.
  # @param p integer, number of variables.
  # @param tolerance tolerance for sparsity.
  # @param k integer or integer vector of which component should be used.
  # @param scaled logical, if total number or percentage of possible sparse entries should be returned.
  # @return Returns either a percentage (\code{scaled = TRUE}) or the amount of zero-values entries (\code{scaled = FALSE}).

  sparse_entries = sum( abs(PC[, k ]) <= tolerance)
  if(scaled){
    return(sparse_entries/(length(k)*N*(p-1))    )
  } else {
    return(sparse_entries)
  }
}

sparsity_group = function(PC, N, p, tolerance = 0, k = 1, scaled = TRUE){

  # Group-wise Sparsity in the Loadings
  # @param PC matrix-like object of PCs.
  # @param N integer, number of groups.
  # @param p integer, number of variables.
  # @param tolerance tolerance for sparsity.
  # @param k integer, which components should be used. Does not work for multiple PCs simultaneously.
  # @param scaled logical, if total number or percentage of possible sparse entries should be returned.
  # @return Returns either a matrix of percentages (\code{scaled = TRUE}) or the amounts
  #  of zero-values entries (\code{scaled = FALSE}) for each group/neighborhood.

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

explained_var_N = function(x, COVS){

  # This function calculates the explained variance for a vector \( x \) with respect to a list of covariance matrices.
  # It computes the explained variance for each covariance matrix individually and returns the results in a matrix.
  # @param x A numeric vector of length \( N x p \), where \( N \) is the number of covariance
  # matrices and \( p \) is the dimension of each covariance matrix.
  # @param COVS A list of numeric matrices, where each matrix has dimensions \( p x p \),
  # representing the covariance matrices.
  # @return A numeric matrix with \( N \) rows and 3 columns:
  #   \itemize{
  #     \item Percentage of explained variance for each covariance matrix.
  #     \item Explained variance for each covariance matrix.
  #     \item Total variance of each covariance matrix.
  #   }
  # @details
  # The function computes the variance explained by the vector \( x \) for each covariance matrix in the list `COVS`. It provides results in a matrix where each row corresponds to one covariance matrix.
  # used in summary and to adjust rho

  explained_var_singular= function(x, COVS){
    # COVS: covariance matrix
    # x: p vector
    var_all = sum(diag(COVS))
    var_x = sum(diag(t(x) %*% COVS %*% x))
    return(c(var_x/var_all, var_x, var_all))
  }

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
#### INTERNAL WRAPPER FUNCTIONS                                                       ####
##########################################################################################

sparsePCAloc = function(eta,
                      gamma,
                      COVS,
                      cor = FALSE,
                      rho = NA,
                      k = NULL,
                      eps_threshold = NULL,
                      eps_ADMM = 1e-4,
                      n_max = 200,
                      eps_root = 1e-1,
                      maxiter_root = 50,
                      increase_rho = list(TRUE, 100, 1),
                      convergence_plot = TRUE,
                      adjust_eta = TRUE,
                      show_progress = TRUE){

  # Calculate Sparse Principle Components
  # @param eta numeric, degree of sparsity.
  # @param gamma numeric, distribution of sparsity.
  # @param COVS list of covariance or correlation matrices.
  # @param cor logical, if starting value for correlation or covariance matrices should be used.
  # @param rho numeric bigger than zero, penalty for ADMM.
  # @param k number of components to calculate.
  # @param eps_threshold tolerance for thresholding.
  # @param eps_ADMM tolerance for ADMM convergence.
  # @param n_max number of maximal iterations.
  # @param eps_root tolerance for root finder.
  # @param maxiter_root maximal number of iterations for root finder.
  # @param increase_rho list with entries for stable convergence. See Details.
  # @param convergence_plot logical, if convergence plot should be displayed.
  # @param starting_value optional given starting value.
  # @param adjust_eta logical, if eta should be adjusted by the variance.

  p = dim(COVS[[1]])[1]
  N = length(COVS)

  if(is.null(k)) k = p
  if(is.na(maxiter_root)) maxiter_root = 50
  if(is.na(rho)) rho = sum(sapply(COVS, diag))/N
  changing_rho = increase_rho[[1]]
  changing_rho_max = max(rho, increase_rho[[2]])

  PC_flat = matrix(NA, nrow = N*p, ncol = 0)
  converged_k = rep(NA, k)
  n_steps = rep(NA, k)
  residuals = list()
  eta_adjusted = rep(NA, k)
  names(eta_adjusted) = paste0("PC", 1:k)
  eta_j = eta

  for(l in 1:k){

    converged = FALSE
    rho_while = rho + eta_j
    if(adjust_eta) {
      eta_j = adjusted_eta(COVS = COVS, Xi = PC_flat, k = l,  eta = eta)
    }

    while(!converged & rho_while <= changing_rho_max) {
      tmp = tryCatch({
        solve_ADMM(rho = rho_while,
                   eta = eta_j,
                   gamma = gamma,
                   Sigma = COVS,
                   k = l,
                   Xi = PC_flat,
                   n_max = n_max,
                   convergence_plot = convergence_plot,
                   cor = cor,
                   maxiter_root = maxiter_root,
                   eps_root = eps_root,
                   eps_ADMM = eps_ADMM,
                   eps_threshold = eps_threshold,
                   show_progress = show_progress)
      }, error = function(cond){
        message(cond)
        return(NULL)
      })

      if(!is.null(tmp)) converged = T
      if(!changing_rho) break
      if(is.null(tmp) & changing_rho) {
        rho_while = rho_while + increase_rho[[3]]
      }
    }

    PC_flat = cbind(PC_flat, as.matrix(tmp$PC, ncol = 1))
    n_steps[l] = tmp$n_steps
    converged_k[l] = tmp$converged
    eta_adjusted[l] = eta_j
    residuals = c(residuals, list(tmp$residuals))

    # adjust rho
    exvar = sum(explained_var_N(tmp$PC, COVS)[, 2])/  sum(explained_var_N(tmp$PC, COVS)[, 3])
    rho = rho * (1-exvar)
  }

  msPCA_out = list(PC = PC_flat,
                converged = converged_k,
                n_steps = n_steps,
                residuals = residuals,
                gamma = gamma,
                eta = eta,
                N = N,
                p = p,
                k = k,
                eta_adjusted = eta_adjusted)

  class(msPCA_out) <- c("msPCA", "list")

  return(msPCA_out)
}







#' @importFrom DescTools AUC
#' @import ggplot2
#' @importFrom Matrix bdiag
#' @importFrom utils flush.console
select_sparsity = function(COVS,
                           rho = NULL,
                           cor = FALSE,
                           eta = seq(0, 5, by = 0.2),
                           gamma = seq(0, 1, 0.05),
                           eps_threshold = 1e-3,
                           eps_root = 1e-1,
                           eps_ADMM = 1e-4,
                           n_max = 300,
                           increase_rho = list(TRUE, 100, 1),
                           convergence_plot = FALSE,
                           plot = TRUE){

  # Optimal Sparsity Parameter Selection for PCA
  # @param COVS list of covariance or correlation matrices.
  # @param k number of components to be returned.
  # @param rho penalty parameter for ADMM.
  # @param cor logical, if starting values for covariances or correlation matrices should be used.
  # @param eta vector of possible values for degree of sparsity.
  # @param gamma vector of possible values for distribution of sparsity. If only one value is provided, the optimal eta is calculated.
  # @param eps_threshold tolerance for thresholding.
  # @param eps_root tolerance for root finder.
  # @param eps_ADMM tolerance for ADMM iterations.
  # @param n_max maximal number of ADMM iterations.
  # @param adjust_eta if eta should be adjusted for further components.
  # @param increase_rho list of settings for improved automated calculation and convergence. See Details.
  # @param convergence_plot logical, if convergence plot should be plotted. Not applicable for \code{cores > 1}.
  # @param stop.sparse calculate if AUC should be calculated for PCAs until full sparsity is reached (\code{TRUE})
  # or over the whole eta range (\code{FALSE}). Set to \code{TRUE}.
  # @return Returns list with \tabular{ll}{
  #    \code{PCA} \tab object of type PCAloc.\cr
  #    \tab \cr
  #    \code{PC} \tab local loadings of PCA\cr
  #    \tab \cr
  #    \code{gamma} \tab optimal value for gamma. \cr
  #    \tab \cr
  #    \code{eta} \tab optimal value for eta. \cr
  #    \tab \cr
  #    \code{eta_tpo} \tab values of Trade-Off-Product for eta from optimization process. \cr
  #    \tab \cr
  #    \code{auc} \tab area under the curve for varying gamma values. \cr
  #    \tab \cr
  #    \code{pars} \tab parameters and respective sparsity entrywise and mixed and explained variance. \cr
  #    \tab \cr
  #    \code{plot} \tab ggplot object for optimal parameter selection. \cr
  #    \tab \cr
  #    \code{plot_info} \tab additional data for plotting functions. \cr
  # }


  sparsity = function(PC, p, N, k = 1, tolerance = 0, mean = "arithmetic"){

    sg = sparsity_group(PC = PC, N = N, p = p, tolerance = tolerance, k = 1, scaled = TRUE)
    se = sparsity_entries(PC = PC, N = N, p = p, tolerance = tolerance, k = 1, scaled = TRUE)

    return(c((se+sg)/2, se, sg))
  }

  explained_var = function(COVS, PC, k, type = "scaled", cor = FALSE, gamma = 0.5){

    # used for TPO and AUC

    # Explained Variance summarized over Groups
    # @param COVS list of covariance matrices
    # @param PC matrix-like object holding the loadings of length np
    # @param k which component should be evaluated
    # @param type character, either \code{"scaled"} for scaling using the extremes solutions or \code{"percent"} as percentage of overall variance.
    # @param cor logical, if \code{COVS} is a correlation matrix or not
    # @param gamma scalar between 0 and 1 indicatig distribution of sparsity.
    # @return Returns scalar

    SIGMA = as.matrix(Matrix::bdiag(COVS))
    var_full = sum(diag(SIGMA))
    var_expl = sum(diag(t(PC[, k]) %*% SIGMA %*% PC[, k]))

    starts = starting_value_ADMM (Sigma = COVS,
                                  eta = 1,
                                  gamma = gamma,
                                  k = k,
                                  return_all = TRUE,
                                  Xi = NULL,
                                  cor = cor)

    var_min = sum(diag(t(starts$penalty_solution) %*% SIGMA %*% starts$penalty_solution))
    var_max = sum(diag(t(starts$eigenvector) %*% SIGMA %*% starts$eigenvector))

    var_scaled = (var_expl - var_min)/(var_max - var_min)
    return(var_scaled)
  }

  # matrix of lists
  PCA_list <- matrix( list(), ncol = length(gamma), nrow = length(eta))
  colnames(PCA_list) = paste0("gamma", gamma)
  rownames(PCA_list) = paste0("eta", eta)

  p = dim(COVS[[1]])[1]
  N = length(COVS)

  if(is.null(rho)) rho = sum(sapply(COVS, diag))/N + max(eta)

  # run all PCAs
  j = 0
  PCA_list = list()
  pars = expand.grid("gamma" = gamma, "eta" = eta)
  sparsitym = rep(NA, j)
  sparsitye = rep(NA, j)
  expvar = rep(NA, j)

  # make progress bar
  progress_bar <- function(j, total = dim(pars)[1], prefix = "", bar_length = 50) {
    filled <- floor(j / total * bar_length)
    bar <- paste0(paste0(rep("=", filled), collapse = ""),  paste0(rep(" ", bar_length - filled), collapse = ""))
    cat(sprintf("\rTuning |%s| %d parameter combinations of %d", bar, j, dim(pars)[1]))
    utils::flush.console()
  }


  for(j in 1:dim(pars)[1]){
    progress_bar(j)
    PCA1 = sparsePCAloc(rho = rho,
                      eta = pars$eta[j],
                      gamma = pars$gamma[j],
                      COVS = COVS,
                      k = 1,
                      n_max = n_max,
                      eps_root = eps_root,
                      eps_ADMM = eps_ADMM,
                      increase_rho = increase_rho,
                      convergence_plot = convergence_plot,
                      adjust_eta = TRUE,
                      eps_threshold = eps_threshold,
                      show_progress = FALSE)
    PCA_list = c(PCA_list, list(PCA1))

    sparsitym[j] = sparsity(PC = PCA1$PC,p = p, N = N)[1]
    sparsitye[j] = sparsity(PC = PCA1$PC,p = p, N = N)[2]
    expvar[j] = explained_var(COVS = COVS, PC = PCA1$PC, k = 1, cor = cor, gamma = PCA1$gamma)
  }
  cat("\n")


  # GET OPTIMAL GAMMA
  auc = rep(NA, length(gamma))
  for(g in gamma){
    ind_gamma = which(pars$gamma == g)
    max_eta_ind = min(which(sparsitye[ind_gamma] == max(sparsitye[ind_gamma])))
    auc[which(gamma == g)] = DescTools::AUC(x = sparsitym[ind_gamma][1:max_eta_ind],
                                            y = expvar[ind_gamma][1:max_eta_ind])

    if(any(is.na(auc[which(gamma == g)]))){
      message("NA value for AUC. Check sparsity and explained variance.")
    }
  }
  names(auc) = gamma

  optimal_gamma = gamma[ which.max(auc)]
  optimal_gamma_etaindex = which(pars$gamma == optimal_gamma)

  g_auc = ggplot() +
    geom_line(aes( x = gamma, y = auc)) +
    geom_point(aes( x = gamma, y = auc)) +
    geom_point(aes( x = optimal_gamma,
                    y = auc[ which.max(auc)]),
               col = "red") +
    labs(title = "Optimal Distribution of Sparsity",
         x = expression(gamma),
         y = "AUC") +
    theme_bw()

  # GET OPTIMAL ETA
  eta_set = pars$eta[optimal_gamma_etaindex]
  sort_ind = sort.int(eta_set, index.return = TRUE)$ix

  eta_set = eta_set[sort_ind]
  spe = sparsitye[optimal_gamma_etaindex][sort_ind]
  exv = expvar[optimal_gamma_etaindex][sort_ind]

  TPO = spe * exv
  optimal_eta = eta_set[which.max(TPO)]
  optimal_eta_index = which.max(TPO)
  names(TPO) = eta

  g_tpo = ggplot() +
    geom_line(aes( x = eta_set, y = TPO)) +
    geom_point(aes( x = eta_set, y = TPO)) +
    geom_point(aes( x = optimal_eta,
                    y = TPO[ which.max(TPO)]),
               col = "red") +
    labs(title = "Optimal Sparsity",
         x = expression(eta),
         y = "TPO") +
    theme_bw()


  # RETURN OPTIMAL PCA
  PC = PCA_list[[(optimal_gamma_etaindex[optimal_eta_index])]]

  return(c(PC,
           list(
              tpo = TPO,
              auc = auc,
              eval = cbind(pars, "s_mixed" = sparsitym, "explained_var" = expvar, "s_entry" = sparsitye),
              plot_auc = g_auc,
              plot_tpo = g_tpo)))
}

