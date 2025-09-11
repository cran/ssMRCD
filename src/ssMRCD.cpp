#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @title Regularized Covariance Matrix Calculation
//' @description Computes a regularized covariance matrix using a convex combination of a target matrix and the sample covariance matrix.
//' @name RCOV
//'
//' @param XX A numeric matrix containing the observations, where each column represents a different observation.
//' @param vMu A numeric vector representing the mean vector to center the columns of \code{XX}.
//' @param rho A numeric scalar representing the regularization parameter (between 0 and 1) that balances between the target matrix and the sample covariance.
//' @param mT A numeric matrix representing the target covariance matrix for regularization.
//' @param scfac A numeric scalar representing a scaling factor applied to the sample covariance matrix.
//'
//' @return A named list containing:
//'   \item{rho}{The regularization parameter used in the calculation.}
//'   \item{mT}{The target covariance matrix provided as input.}
//'   \item{cov}{The sample covariance matrix calculated from the centered observations.}
//'   \item{rcov}{The regularized covariance matrix, which is a convex combination of the target matrix and scaled sample covariance matrix.}
//'
//' @details The function calculates the sample covariance matrix of the centered data \code{XX} and combines it with the target covariance matrix \code{mT}, scaled by the factor \code{scfac}. The regularization is controlled by the parameter \code{rho}, where \code{rho = 1} results in using only the target matrix, and \code{rho = 0} uses only the sample covariance.
//'
//' @importFrom Rcpp sourceCpp
//' @keywords internal

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List RCOV(const arma::mat& XX, const arma::vec& vMu, double rho, const arma::mat& mT, double scfac) {
  arma::mat mE = XX.each_col() - vMu; // Center the matrix
  int n = mE.n_cols; // Number of observations

  arma::mat mS = (mE * mE.t()) / n; // Sample covariance
  arma::mat rcov = rho * mT + (1 - rho) * scfac * mS; // Regularized covariance

  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("mT") = mT,
    Rcpp::Named("cov") = mS,
    Rcpp::Named("rcov") = rcov
  );
}



//' @title Calculation of Objective Function
//' @description Calculation of the value of the objective function for a given list of matrices,
//' lambda, and a weighting matrix.
//' @param matrix_list A list of matrices \eqn{K_i}.
//' @param lambda Scalar smoothing parameter.
//' @param weights Matrix of weights.
//' @return Returns the value of the objective function.
//' @name objective_matrix
//' @keywords internal
// [[Rcpp::export]]
 double objective_matrix(const List& matrix_list, double lambda, const arma::mat& weights) {

   // Convert List to vector of matrices
   int N = matrix_list.size();
   int p = Rcpp::as<arma::mat>(matrix_list[0]).n_rows;

   double f = 0.0;
   std::vector<arma::mat> matrices(N);

   for (int i = 0; i < N; ++i) {
     matrices[i] = Rcpp::as<arma::mat>(matrix_list[i]);
   }



   for (int i = 0; i < N; ++i) {
     arma::mat P = arma::zeros(p, p);

     for (int j = 0; j < N; ++j) {
       P += matrices[j] * weights(i, j);
     }

     arma::mat A = (1 - lambda) * matrices[i] + lambda * P;
     f += arma::det(A);
   }

   return f;
 }



//' @title Objective Function for Init Object
//' @description Calculates the objective function based on the matrices stored in the init object.
//' @param init_object List object containing matrices
//' @param lambda Scalar for spatial smoothing
//' @param weights Matrix with weights
//' @return Returns the value of the objective function.
//' @keywords internal
// [[Rcpp::export]]
double objective_init(const List& init_object, double lambda, const arma::mat& weights) {
 int N = init_object.size();

 // Create a List to hold the matrices
 List matrix_list(N);

 for (int i = 0; i < N; ++i) {
   List item = init_object[i];
   // Convert item["mS"] to arma::mat
   try {
     matrix_list[i] = Rcpp::as<arma::mat>(item["mS"]);
   } catch (const Rcpp::not_compatible& e) {
     Rcpp::stop("Error converting item['mS'] to arma::mat: %s", e.what());
   }
 }

 return objective_matrix(matrix_list, lambda, weights);
}



//' Computes Mahalanobis Distances for a Given Set of H-Subsets
//'
//' This function calculates the Mahalanobis distances for a set of observations
//' by centering them with the mean vector and using a covariance matrix computed
//' as a weighted combination of the covariance matrix of the current item and
//' the covariance matrices of its neighbors.
//'
//' @param init A list of items where each item contains the following elements:
//'   \itemize{
//'     \item \code{mX} A matrix of observations (one column per observation).
//'     \item \code{vMu} A vector of means.
//'     \item \code{mS} A covariance matrix of the observations in \code{mX}.
//'   }
//' @param i An integer index specifying which item from the \code{init} list to use.
//' @param lambda A numeric value representing the weight for the covariance matrix of the current item.
//' @param weights A matrix of weights where each element \code{weights(i, j)} specifies the weight of the \code{j}-th item for the \code{i}-th item.
//'
//' @return A numeric vector of distances for each observation in the centered matrix.
//'
//' @details The Mahalanobis distances are computed using the covariance matrix,
//' which is a weighted combination of the current item's covariance matrix and
//' those of its neighbors. The covariance matrix is smoothed using the parameter
//' \code{lambda} and the distances are computed as \code{(x_centered^T * Cov_matrix_chol_inv * x_centered)} for each observation.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::vec dist_cstep(const List& init, int i, double lambda, const arma::mat& weights) {
  int N = init.size();

  // Extract the matrix and vector from the init list
  List item = init[i];

  arma::mat X = as<arma::mat>(item["mX"]);
  arma::vec vMu = as<arma::vec>(item["vMu"]);


  // Ensure vMu is a column vector
  if (vMu.n_rows != X.n_rows) {
    stop("Dimension mismatch: vMu length does not match number of rows in mX.");
  }

  int p = vMu.n_rows;

  // Initialize covariance matrix for neighbors
  arma::mat Cov_neighb(p, p, arma::fill::zeros);
  // Center the matrix X by subtracting vMu from each column
  arma::mat x_centered = X.each_col() - vMu;

  for (int j = 0; j < N; ++j) {
    List tmpj = init[j];
    Cov_neighb += weights(i, j) * as<arma::mat>(tmpj["mS"]);;
  }
  arma::mat Cov_matrix = (1 - lambda) * as<arma::mat>(item["mS"]) + lambda * Cov_neighb;

  // Compute the Cholesky decomposition and its inverse
  arma::mat Cov_matrix_chol_inv = arma::inv_sympd(Cov_matrix);

  // Initialize vector to store distances
  arma::vec distances(x_centered.n_cols);

  // Compute distances for each observation
  for (arma::uword j = 0; j < x_centered.n_cols; ++j) {
    arma::vec x_col = x_centered.col(j);
    distances(j) = as_scalar(x_col.t() * Cov_matrix_chol_inv * x_col);
  }

  return distances;
}


//' Perform Concentration Step
//'
//' This function performs concentration steps by iteratively updating
//' the neighborhoods of items based on Mahalanobis distances. The function computes
//' the covariance matrix and updates the neighborhood list until convergence or
//' a maximum number of iterations is reached.
//'
//' @param init A list of items where each item contains:
//'   \itemize{
//'     \item \code{mX} A matrix of observations (one column per observation).
//'     \item \code{hsets.init} A matrix where each column specifies initial indices for neighbors.
//'     \item \code{rho} A regularization parameter for covariance estimation.
//'     \item \code{scfac} A scaling factor for consistency.
//'     \item \code{mS} A covariance matrix of the observations in \code{mX}.
//'     \item \code{index} (Initially not present; will be updated with indices of neighbors).
//'     \item \code{vdst} (Initially not present; will be updated with Mahalanobis distances).
//'     \item \code{stop} (Initially not present; used to indicate convergence).
//'     \item \code{ret} (Initially not present; will be updated with results from \code{RCOV}).
//'   }
//' @param maxcsteps An integer specifying the maximum number of iterations for the optimization.
//' @param which_indices An integer vector specifying which initial indices to use for each item.
//' @param lambda A numeric value representing the weight for the covariance matrix in the optimization.
//' @param weights A matrix of weights where each element \code{weights(i, j)} specifies the weight of the \code{j}-th item for the \code{i}-th item.
//' @param mT A matrix used for regularization in the covariance matrix calculation.
//'
//' @return A list containing:
//'   \itemize{
//'     \item \code{numit} An integer representing the number of iterations performed.
//'     \item \code{out} The updated list of items with updated neighborhoods and additional information.
//'     \item \code{obj_value} A numeric vector of objective values at each iteration, including the initial value.
//'   }
//'
//' @details The function updates the neighborhoods of each item based on Mahalanobis distances, recalculates means and covariances,
//' and checks for convergence. If the neighborhoods do not change between iterations, the optimization stops early.
//'
//' @keywords internal
// [[Rcpp::export]]
List cstep(const List& init, int maxcsteps, const arma::uvec& which_indices,
           double lambda, const arma::mat& weights, const arma::mat& mT) {

  int N = init.size();

  // Initialize the obj_value vector and fill with NaN
  arma::vec obj_value(maxcsteps +1);
  obj_value.fill(arma::datum::nan);  // Fill the vector with NaN values
  List updated_init = clone(init);  // Clone the list to avoid modifying the input directly

  // Initialize calculations for covariance matrix for each neighborhood
  for (int i = 0; i < N; ++i) {
    List neighborhood = as<List>(updated_init[i]);
    arma::mat hset_indices = neighborhood["hsets.init"];
    arma::vec sel_col = hset_indices.col(which_indices[i]-1); // no R to C++ index shift in storage index
    neighborhood["index"] = sel_col;


    // Select corresponding observations and calculate mean and regularized covariance
    arma::mat XX = as<arma::mat>(neighborhood["mX"]).cols(as<arma::uvec>(neighborhood["index"]) - 1); // R to C++ index shift when accessing
    arma::vec vMu = arma::mean(XX, 1);
    neighborhood["vMu"] = vMu;
    List ret = RCOV(XX, vMu, as<double>(neighborhood["rho"]), mT, as<double>(neighborhood["scfac"]));
    neighborhood["rho"] = ret["rho"];
    neighborhood["mS"] = ret["rcov"];
    updated_init[i] = neighborhood;  // Update the list with modified neighborhood
  }

  obj_value[0] = objective_init(updated_init, lambda, weights);

  // Calculate observations with minimal distances
  for (int i = 0; i < N; ++i) {
    // Call dist_cstep to calculate distances
    arma::vec vdst = dist_cstep(updated_init, i, lambda, weights);

    // Get the current neighborhood list
    List neighborhood = as<List>(updated_init[i]);

    // Sort the indices based on the distances
    arma::uvec sorted_indices = arma::sort_index(vdst);  // C++ indexing
    arma::uvec nndex = sorted_indices.subvec(0, as<int>(neighborhood["h"]) - 1) + 1;  // h-many of the smallest indices, plus one for R storage


    // Update the neighborhood with the new indices
    neighborhood["vdst"] = vdst; // Add 'vdst' to the neighborhood list
    neighborhood["index"] = nndex;

    // Save the updated neighborhood list back to 'init'
    updated_init[i] = neighborhood;
  }

  // Start iteration
  int iter = 0;
  while (iter < (maxcsteps-1)) {
    // Calculate mu and Cov
    for (int i = 0; i < N; ++i) {
      List neighborhood = as<List>(updated_init[i]);
      arma::mat XX = as<arma::mat>(neighborhood["mX"]).cols(as<arma::uvec>(neighborhood["index"]) - 1); // R to C++ index shift when accessing
      arma::vec vMu = arma::mean(XX, 1);
      neighborhood["vMu"] = vMu;

      List ret = RCOV(XX, vMu, as<double>(neighborhood["rho"]), mT, as<double>(neighborhood["scfac"]));
      neighborhood["mS"] = ret["rcov"];
      neighborhood["ret"] = ret;
      updated_init[i] = neighborhood;
    }

    // Select new minimal observations
    for (int i = 0; i < N; ++i) {
      List neighborhood = as<List>(updated_init[i]);

      // Convert the proxy to an arma::vec
      arma::vec vdst = dist_cstep(updated_init, i, lambda, weights);
      neighborhood["vdst"] = vdst; // Add 'vdst' to the neighborhood list

      // Sort the indices based on the distances and get the top 'h' smallest indices
      arma::uvec sorted_indices = arma::sort_index(vdst);
      arma::uvec nndex = sorted_indices.subvec(0, as<int>(neighborhood["h"]) - 1) + 1;  // h-many of the smallest indices, plus one for R storage

      // If observations stay the same, break and use last calculated values
      if (arma::all(nndex == as<arma::uvec>(neighborhood["index"]))) {
        neighborhood["stop"] = 1;
      } else {
        neighborhood["stop"] = 0;
      }

      neighborhood["index"] = nndex;
      updated_init[i] = neighborhood;
    }

    // Save information
    obj_value[iter+1]= objective_init(updated_init, lambda, weights);

    // Check break condition
    int break_sum = 0;
    for (int i = 0; i < N; ++i) {
      break_sum += as<int>(as<List>(updated_init[i])["stop"]);
    }
    if (break_sum == N) {
      break;
    }
    iter++;
  }

  return List::create(
    Named("numit") = iter,
    Named("out") = updated_init,
    Named("obj_value") = obj_value
  );
}

