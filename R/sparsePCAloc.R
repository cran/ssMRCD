# USER FUNCTIONS
##########################################################################################
#' Calculate Sparse Principle Components
#'
#' @param eta numeric, degree of sparsity.
#' @param gamma numeric, distribution of sparsity.
#' @param COVS list of covariance or correlation matrices.
#' @param cor logical, if starting value for correlation or covariance matrices should be used.
#' @param rho numeric bigger than zero, penalty for ADMM.
#' @param k number of components to calculate.
#' @param eps_threshold tolerance for thresholding.
#' @param eps_ADMM tolerance for ADMM convergence.
#' @param n_max number of maximal iterations.
#' @param eps_root tolerance for root finder.
#' @param maxiter_root maximal number of iterations for root finder.
#' @param increase_rho list with entries for stable convergence. See Details.
#' @param convergence_plot logical, if convergence plot should be displayed.
#' @param starting_value optional given starting value.
#' @param adjust_eta logical, if eta should be adjusted by the variance.
#' @param trace logical, if messages should be displayed.
#'
#' @details The input \code{increase_rho} consists of a logical indicating if rho should be adjusted
#' if algorithm did not converged within the given maximal number of iterations. Two integers specify the
#' maximal \code{rho} that is allowed and the step size.
#'
#' @return An object of class \code{"PCAloc"} containing the following elements:\tabular{ll}{
#'    \code{PC} \tab Matrix of dimension Np x k of stacked loading vectors.  \cr
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
#'    \code{summary} \tab Description of result per component. \cr
#'    \tab \cr
#'    \code{residuals} \tab Primary and secondary residuals. \cr
#'    \tab \cr
#' }
#'
#' @export
#'
#' @examples
#' C1 = diag(c(1.1, 0.9, 0.6))
#' C2 = matrix(c(1.1, 0.1, -0.1,
#'               0.1, 1.0, -0.2,
#'              -0.1, -0.2, 0.7), ncol = 3)
#' C3 = (C1 + C2)/2
#'
#' sparsePCAloc(eta = 1, gamma = 0.5, cor = FALSE, COVS = list(C1, C2, C3),
#'              n_max = 100, increase_rho = list(FALSE, 100, 1))
#'
sparsePCAloc = function(eta,
                      gamma,
                      COVS,
                      cor = FALSE,
                      rho = NULL,
                      k = NULL,
                      eps_threshold = NULL,
                      eps_ADMM = 1e-4,
                      n_max = 200,
                      eps_root = 1e-1,
                      maxiter_root = 50,
                      increase_rho = list(TRUE, 100, 1),
                      convergence_plot = TRUE,
                      starting_value = NULL,
                      adjust_eta = TRUE,
                      trace = TRUE){

  p = dim(COVS[[1]])[1]
  N = length(COVS)
  if(is.null(k)) k = p

  if(is.null(rho)) rho = sum(sapply(COVS, diag))/N
  changing_rho = increase_rho[[1]]
  changing_rho_max = max(rho, increase_rho[[2]])

  PC = matrix(NA, ncol = 0, nrow = N*p)
  converged_k = rep(NA, k)
  n_steps = rep(NA, k)
  summary = list()
  residuals = list()
  eta_j = eta

  # Conditions to check:
  # stopifnot(all.equal(COVS, t(COVS)))
  # stopifnot(length(U) == length(A) & dim(COVS)[1] == length(A))
  # stopifnot(rho > 0, k > 0)
  # if(k > 1) stopifnot(dim(Xi)[1] == dim(COVS)[1])
  # stopifnot(eta >= 0 & gamma >= 0 & gamma <= 1 & rho > 0 & k > 0 & is.list(COVS)) -- > add to user functions

  for(l in 1:k){
    converged = FALSE
    rho_while = rho + eta_j

    if(adjust_eta) {
      eta_j = adjusted_eta(COVS = COVS, Xi  = PC, k = l,  eta = eta)
      coloured_print(paste("eta_j:", eta_j),
                     colour = "message",
                     trace = trace)
    }

    while(!converged & rho_while <= changing_rho_max) {
      tmp = tryCatch({
        solve_ADMM(rho = rho_while,
                   lambda = eta_j,
                   alpha = gamma,
                   Sigma = COVS,
                   k = l,
                   Xi = PC,
                   starting_value = starting_value,
                   trace = trace,
                   n_max = n_max,
                   convergence_plot = convergence_plot,
                   cor = cor,
                   maxiter_root = maxiter_root,
                   eps_root = eps_root,
                   eps_ADMM = eps_ADMM,
                   eps_threshold = eps_threshold)
      }, error = function(cond){
        print(cond)
        if(changing_rho & trace) coloured_print(paste("Increase rho to:", rho_while + increase_rho[[3]]),
                                                colour = "message",
                                                trace = trace)
        return(NULL)
      })

      if(!is.null(tmp)) converged = T
      if(!changing_rho) break
      if(is.null(tmp) & changing_rho) {
        rho_while = rho_while + increase_rho[[3]]
      }
      if(!is.null(tmp)){
        if(changing_rho & max(tmp$residuals) > 1 & trace) {
          coloured_print(paste("Residuals to high, increase rho to:", rho_while + increase_rho[[3]]),
                         colour = "message",
                         trace = trace)
          rho_while = rho_while + increase_rho[[3]]
        }
      }
    }

    PC = cbind(PC, as.matrix(tmp$PC, ncol = 1))
    n_steps[l] = tmp$n_steps
    converged_k[l] = tmp$converged
    summary = c(summary, list(tmp$summary))
    residuals = c(residuals, list(tmp$residuals))

    # adjust rho
    exvar = sum(explained_var_N(tmp$PC, COVS)[, 2])/  sum(explained_var_N(tmp$PC, COVS)[, 3])
    #rho = rho * (1-exvar)
  }

  colnames(PC) = paste0("PC", 1:k)
  PCAloc = list(PC = PC,
                converged = converged_k,
                n_steps = n_steps,
                summary = summary,
                residuals = residuals,
                gamma = gamma,
                eta = eta,
                N = N,
                p = p,
                k = k,
                eps_threshold = eps_threshold,
                COVS = COVS)
  class(PCAloc) <- c("PCAloc", "list")

  return(PCAloc)
}


##########################################################################################
#' Optimal Sparsity Parameter Selection for PCA
#'
#' @param COVS list of covariance or correlation matrices.
#' @param k number of components to be returned.
#' @param rho penalty parameter for ADMM.
#' @param cor logical, if starting values for covariances or correlation matrices should be used.
#' @param eta vector of possible values for degree of sparsity.
#' @param gamma vector of possible values for distribution of sparsity. If only one value is provided, the optimal eta is calculated.
#' @param eps_threshold tolerance for thresholding.
#' @param eps_root tolerance for root finder.
#' @param eps_ADMM tolerance for ADMM iterations.
#' @param n_max maximal number of ADMM iterations.
#' @param adjust_eta if eta should be adjusted for further components.
#' @param cores number of cores for parallel computing.
#' @param increase_rho list of settings for improved automated calculation and convergence. See Details.
#' @param convergence_plot logical, if convergence plot should be plotted. Not applicable for \code{cores > 1}.
#' @param trace logical, if messages should be displayed. Not applicable for \code{cores > 1}.
#' @param stop.sparse calculate if AUC should be calculated for PCAs until full sparsity is reached (\code{TRUE})
#' or over the whole eta range (\code{FALSE}). Set to \code{TRUE}.
#'
#' @return Returns list with \tabular{ll}{
#'    \code{PCA} \tab object of type PCAloc.\cr
#'    \tab \cr
#'    \code{PC} \tab local loadings of PCA\cr
#'    \tab \cr
#'    \code{gamma} \tab optimal value for gamma. \cr
#'    \tab \cr
#'    \code{eta} \tab optimal value for eta. \cr
#'    \tab \cr
#'    \code{eta_tpo} \tab values of Trade-Off-Product for eta from optimization process. \cr
#'    \tab \cr
#'    \code{auc} \tab area under the curve for varying gamma values. \cr
#'    \tab \cr
#'    \code{pars} \tab parameters and respective sparsity entrywise and mixed and explained variance. \cr
#'    \tab \cr
#'    \code{plot} \tab ggplot object for optimal parameter selection. \cr
#'    \tab \cr
#'    \code{plot_info} \tab additional data for plotting functions. \cr
#' }
#' @export
#'
#' @details The input \code{increase_rho} consists of a logical indicating if rho should be adjusted
#' if algorithm did not converged within the given maximal number of iterations. Two integers specify the
#' maximal \code{rho} that is allowed and the step size.
#'
#' @importFrom DescTools AUC
#' @import parallel
#' @import foreach
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' C1 = matrix(c(1,0,0,0.9), ncol = 2)
#' C2 = matrix(c(1.1, 0.1, 0.1, 1), ncol = 2)
#' C3 = matrix(c(1.2, 0.2, 0.2, 1), ncol = 2)
#'
#' select_sparsity(COVS = list(C1, C2, C3),
#'                 k = 1,
#'                 rho = 5,
#'                 eta = c(0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5 ,0.75,  1),
#'                 gamma =c(0, 0.25, 0.5, 0.75, 1),
#'                 eps_threshold = 0.005,
#'                 increase_rho = list(FALSE, 20, 5))
#' }

select_sparsity = function(COVS,
                           k = 1,
                           rho = NULL,
                           cor = FALSE,
                           eta = seq(0, 5, by = 0.2),
                           gamma = seq(0, 1, 0.05),
                           eps_threshold = 1e-3,
                           eps_root = 1e-1,
                           eps_ADMM = 1e-4,
                           n_max = 300,
                           adjust_eta = FALSE,
                           cores = 1,
                           increase_rho = list(TRUE, 100, 1),
                           convergence_plot = FALSE,
                           trace = FALSE,
                           stop.sparse = TRUE){

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

  if(cores !=  1) {
    cl_select = parallel::makeCluster(cores)
    PCA_list = foreach::foreach(j= 1:dim(pars)[1],
                                .combine = "c",
                                .packages = c("ssMRCD")) %dopar% {

                                  gamma = pars$gamma[j]
                                  eta = pars$eta[j]

                                  PCA1 = sparsePCAloc(rho = rho,
                                                    eta = eta,
                                                    gamma = gamma,
                                                    COVS = COVS,
                                                    k = k,
                                                    n_max = n_max,
                                                    eps_root = eps_root,
                                                    eps_ADMM = eps_ADMM,
                                                    increase_rho = increase_rho,
                                                    convergence_plot = FALSE,
                                                    trace = FALSE,
                                                    adjust_eta = adjust_eta,
                                                    eps_threshold = eps_threshold,
                                                    cor = cor)
                                  list(PCA1)
                                }
    stopCluster(cl_select)
  }
  if(cores == 1){

    for(j in 1:dim(pars)[1]){
      PCA1 = sparsePCAloc(rho = rho,
                        eta = pars$eta[j],
                        gamma = pars$gamma[j],
                        COVS = COVS,
                        k = k,
                        n_max = n_max,
                        eps_root = eps_root,
                        eps_ADMM = eps_ADMM,
                        increase_rho = increase_rho,
                        convergence_plot = convergence_plot,
                        trace = trace,
                        adjust_eta = adjust_eta,
                        eps_threshold = eps_threshold,
                        cor = cor)
      PCA_list = c(PCA_list, list(PCA1))
    }
  }

  sparsitym = sapply(X = PCA_list, function(x) sparsity_mixed(PC = x$PC,
                                                              p = p,
                                                              N = N,
                                                              k = 1,
                                                              tolerance = 0.0,
                                                              mean = "arithmetic"))
  expvarm = sapply(X = PCA_list, function(x) explained_var(COVS = COVS,
                                                           PC = x$PC,
                                                           k = 1,
                                                           type = "scaled",
                                                           cor = cor,
                                                           gamma = x$gamma)[1])
  sparsitye = sapply(X = PCA_list, function(x) sparsity_entries(PC = x$PC,
                                                                k = 1,
                                                                N = N,
                                                                p = p,
                                                                tolerance = 0,
                                                                scaled = TRUE)
  )

  # GET OPTIMAL GAMMA
  auc = rep(NA, length(gamma))
  for(g in gamma){
    ind_gamma = which(pars$gamma == g)

    if(stop.sparse) {
      max_eta_ind = min(which(sparsitye[ind_gamma] == max(sparsitye[ind_gamma])))
      auc[which(gamma == g)] = DescTools::AUC(x = sparsitym[ind_gamma][1:max_eta_ind],
                                              y = expvarm[ind_gamma][1:max_eta_ind])
    } else {
      auc[which(gamma == g)] = DescTools::AUC(x = sparsitym[ind_gamma],
                                              y = expvarm[ind_gamma])
    }

    if(is.na(auc[which(gamma == g)])){
      coloured_print("NA values for AUC. Check sparsity and explained variance.",
                     colour = "message")
    }
  }

  optimal_gamma_index = which.max(auc)
  optimal_gamma = gamma[optimal_gamma_index]
  optimal_gamma_etaindex = which(pars$gamma == optimal_gamma)

  g = ggplot() +
    geom_line(aes( x = gamma, y = auc)) +
    geom_point(aes( x = gamma, y = auc)) +
    geom_point(aes( x = optimal_gamma,
                    y = auc[optimal_gamma_index]),
               col = "red") +
    labs(title = "Optimal Distribution of Sparsity",
         x = "gamma",
         y = "AUC") +
    theme_bw()


  # GET OPTIMAL ETA
  eta_set = pars$eta[optimal_gamma_etaindex]
  sort_ind = sort.int(eta_set, index.return = TRUE)$ix

  eta_set = eta_set[sort_ind]
  spe = sparsitye[optimal_gamma_etaindex][sort_ind]
  exv = expvarm[optimal_gamma_etaindex][sort_ind]

  TPO = spe * exv
  optimal_eta = eta_set[which.max(TPO)]
  optimal_eta_index = which.max(TPO)


  # RETURN OPTIMAL PCA
  PC = PCA_list[[(optimal_gamma_etaindex[optimal_eta_index])]]

  return(list(PCA = PC,
              PC = PC$PC,
              gamma = optimal_gamma,
              eta = optimal_eta,
              eta_tpo = TPO,
              auc = auc,
              pars = cbind(pars,
                           "sparsem" = sparsitym,
                           "explained_var" = expvarm,
                           "sparsee" = sparsitye),
              plot = g,
              plot_info = list(PCA_list = PCA_list,
                               parameters = pars)))
}

