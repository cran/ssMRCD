#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]
double density_mnorm(arma::vec x, arma::vec mean, arma::mat sigmai) {

  // Ensure `x` and `mean` are column vectors (1-column matrices)
  arma::mat x_mat = arma::reshape(x, x.n_elem, 1);
  arma::mat mean_mat = arma::reshape(mean, mean.n_elem, 1);

  //Rcpp::Rcout << "Herea " << std::endl;
  // Dimension of the mean vector
  int p = mean.n_elem;

  // Mahalanobis distance term
  arma::mat diff = x_mat - mean_mat;
  double exponent = -0.5 * arma::as_scalar(diff.t() * sigmai * diff);

  //Rcpp::Rcout << "Hereb " << std::endl;
  // Normalization constant
  double norm_const = std::pow(2 * M_PI, -0.5 * p) * std::sqrt(arma::det(sigmai));

  //Rcpp::Rcout << "Herec " << std::endl;
  // Density calculation
  double density = std::exp(exponent) * norm_const;

  //Rcpp::Rcout << "Hered " << std::endl;
  return density;
}





// Function to check if a row already exists in the result
bool row_exists(const arma::mat& result, const arma::rowvec& row) {
  for (size_t i = 0; i < result.n_rows; ++i) {
    if (arma::all(result.row(i) == row)) {
      return true;
    }
  }
  return false;
}

// Function to find unique rows manually
arma::mat find_unique_rows(const arma::mat& A) {
  // Create an empty matrix to store unique rows
  //Rcpp::Rcout << "A " << A << std::endl;
  arma::mat unique_rows(0, A.n_cols);

  // Iterate through each row of matrix A
  for (size_t i = 0; i < A.n_rows; ++i) {
    arma::rowvec current_row = A.row(i);
    //Rcpp::Rcout << "current_row " << current_row << std::endl;

    //Rcpp::Rcout << "unique_rows " << unique_rows << std::endl;
    // Check if the current row already exists in unique_rows
    if (!row_exists(unique_rows, current_row)) {
      // If it doesn't exist, add the current row to the result
      unique_rows.insert_rows(unique_rows.n_rows, current_row);
    }
  }

  return unique_rows;
}

// [[Rcpp::export]]
double objective(const arma::mat& X,
                 const Rcpp::List& mu,
                 const Rcpp::List& Sigma,
                 const arma::mat& pi_groups,
                 const arma::mat& W,
                 const arma::uvec& groups,
                 const arma::mat& Q) {


  //Rcpp::Rcout << "Input W: " << W << std::endl;
  int N = pi_groups.n_rows;
  int p = X.n_cols;

  double obj = 0;

  // Find unique rows in W (equivalent to unique(W) in R)
  arma::mat W_pattern = find_unique_rows(W);
  //Rcpp::Rcout << "Input W: " << W_pattern << std::endl;
  int wn = W_pattern.n_rows;

  for (int l = 0; l < wn; l++) {
    arma::uvec obs = arma::find(W_pattern.row(l) == 1);
    //Rcpp::Rcout << "obs " << obs << std::endl;
    arma::uvec ind_i;

    // Finding indices where colSums(abs(t(W) - W_pattern[l,])) == 0
    if (wn == 1) {
      // If there is only one unique pattern, use all rows
      ind_i = arma::regspace<arma::uvec>(0, W.n_rows - 1);
    } else {

      //Rcpp::Rcout << "Input W: " << W.row(i) << std::endl;
      //Rcpp::Rcout << "Input W: " << W_pattern.row(l) << std::endl;

      // Find indices where the row in W matches the pattern (for multiple patterns)
      for (arma::uword j = 0; j < W.n_rows; j++) {
        //Rcpp::Rcout << "sum " << arma::accu(arma::abs(W.row(j) - W_pattern.row(l)))  << std::endl;

        if (arma::accu(arma::abs(W.row(j) - W_pattern.row(l))) == 0) {
          ind_i.insert_rows(ind_i.n_rows, arma::uvec{(unsigned int)j});
        }
      }
      //Rcpp::Rcout << "ind_i " << ind_i << std::endl;
    }

    //Rcpp::Rcout << "ind_i: " << ind_i << std::endl;
    //Rcpp::Rcout << "Here1 " << std::endl;
    if (!obs.empty()) {
      // If there are observed variables
      std::vector<arma::mat> Sigmai_tmp(N);

      //Rcpp::Rcout << "Here3 " << std::endl;

      // Compute the inverse of Sigma for each group
      for (int k = 0; k < N; k++) {
        // Convert the list element to an arma::mat
        arma::mat Sigma_k = Rcpp::as<arma::mat>(Sigma[k]); // Convert list element Sigma[k] to matrix
        Sigmai_tmp[k] = arma::inv(Sigma_k.submat(obs, obs));  // Inverse of Sigma for observed indices
      }

      //Rcpp::Rcout << "Here4 " << std::endl;

      for (auto i : ind_i) {
        int g = groups[i] - 1;  // Adjust for 0-based indexing
        double dens = 0.0;
        for (int k = 0; k < N; k++) {
          if (pi_groups(g, k) > 0) {
            //Rcpp::Rcout << "Here5 " << std::endl;
            arma::rowvec Xi = X.row(i);  // Extract the i-th row from X
            arma::rowvec Xi_obs = Xi.cols(obs);  // Extract the observed columns of Xi

            arma::rowvec mu_k = mu[k];  // Extract the mean vector for group k
            arma::rowvec mu_k_obs = mu_k.cols(obs);  // Extract the observed columns of mu_k

            //Rcpp::Rcout << "Xi_obs " << Xi_obs << std::endl;
            //Rcpp::Rcout << "mu_k_obs " << mu_k_obs << std::endl;
            //Rcpp::Rcout << "Sigmai_tmp " << Sigmai_tmp[k] << std::endl;
            double tmp = density_mnorm(Xi_obs.t(), mu_k_obs.t(), Sigmai_tmp[k]);  // Calculate density
            dens += tmp * pi_groups(g, k);
            //Rcpp::Rcout << "Here7 " << std::endl;
          }
        }
        obj += std::log(dens);
      }
    }
  }
  //Rcpp::Rcout << "Here2 " << std::endl;
  obj = -2 * obj;
  // Add penalty term
  for (int j = 0; j < p; j++) {
    obj += arma::accu(Q.col(j) % (1 - W.col(j)));
  }
  //Rcpp::Rcout << "End " << std::endl;
  return obj;
}




// [[Rcpp::export]]
arma::mat w_step(arma::mat& X,
                 Rcpp::List& mu,
                 Rcpp::List& Sigma,
                 arma::mat& pi_groups,
                 arma::mat& W,
                 arma::uvec& groups,
                 arma::mat& Q,
                 arma::vec& h) {

  int p = X.n_cols;
  int n = X.n_rows;
  int N = Sigma.size();

  for (int j = 0; j < p; j++) {  // over each variable/column
    arma::mat W_pattern =  find_unique_rows(W);  // Get unique rows
    int wn = W_pattern.n_rows;

    arma::vec delta = arma::zeros<arma::vec>(n);  // Initialize delta vector

    for (int l = 0; l < wn; l++) {  // Over different W patterns
      arma::uvec obs = arma::find(W_pattern.row(l) == 1);
      arma::uvec ind_w;

      // Finding the indices where the row in W matches the pattern W_pattern[l]
      for (arma::uword i = 0; i < W.n_rows; i++) {
        if (arma::accu(arma::abs(W.row(i) - W_pattern.row(l))) == 0) {
          ind_w.insert_rows(ind_w.n_rows, arma::uvec{(unsigned int)i});
        }
      }

      arma::mat W1 = W;
      arma::mat W1ind = W1.rows(ind_w);  // Extract the i-th row from X
      W1ind.col(j).fill(1);  // Extract the observed columns of Xi
     // W1indj.fill(1);  // Set W1[j] = 1 for these rows
      arma::uvec obs1 = arma::find(W1ind.row(0) == 1);

      arma::mat W0 = W;
      arma::mat W0ind = W0.rows(ind_w);  // Extract the i-th row from X
      W0ind.col(j).fill(0);  // Set W0[j] = 0 for these rows
      arma::uvec obs0 = arma::find(W0ind.row(0) == 1);

      // Temporary sigma inverses
      std::vector<arma::mat> Sigmai_tmp1(N);
      for (int k = 0; k < N; k++) {
        arma::mat Sigma_k = Rcpp::as<arma::mat>(Sigma[k]); // Convert list element Sigma[k] to matrix
        Sigmai_tmp1[k] = arma::inv(Sigma_k.submat(obs1, obs1));
      }

      std::vector<arma::mat> Sigmai_tmp0(N);
      if (obs0.n_elem > 0) {
        for (int k = 0; k < N; k++) {
          arma::mat Sigma_k = Rcpp::as<arma::mat>(Sigma[k]); // Convert list element Sigma[k] to matrix
          Sigmai_tmp0[k] = arma::inv(Sigma_k.submat(obs0, obs0));
        }
      }

      for (auto i : ind_w) {  // For each observation
        int g = groups[i] - 1;  // Adjust for 0-based indexing

        // If cell is observed
        double dens = 0.0;

        for (int k = 0; k < N; k++) {
          if (pi_groups(g, k) > 0) {
            arma::rowvec Xi = X.row(i);  // Extract the i-th row from X
            arma::rowvec Xi_obs = Xi.cols(obs1);  // Extract the observed columns of Xi

            arma::rowvec mu_k = mu[k];  // Extract the mean vector for group k
            arma::rowvec mu_k_obs = mu_k.cols(obs1);  // Extract the observed columns of mu_k

            double tmp = density_mnorm(Xi_obs.t(), mu_k_obs.t(), Sigmai_tmp1[k]);
            dens += tmp * pi_groups(g, k);
          }
        }
        double o1 = -2 * std::log(dens);

        // If cell is unobserved
        double o0 = 0.0;
        if (obs0.n_elem > 0) {
          dens = 0.0;
          for (int k = 0; k < N; k++) {
            if (pi_groups(g, k) > 0) {
              arma::rowvec Xi = X.row(i);  // Extract the i-th row from X
              arma::rowvec Xi_obs = Xi.cols(obs0);  // Extract the observed columns of Xi

              arma::rowvec mu_k = mu[k];  // Extract the mean vector for group k
              arma::rowvec mu_k_obs = mu_k.cols(obs0);  // Extract the observed columns of mu_k
              double tmp = density_mnorm(Xi_obs.t(), mu_k_obs.t(), Sigmai_tmp0[k]);
              dens += tmp * pi_groups(g, k);
            }
          }
          o0 = -2 * std::log(dens);
        }

        delta[i] = o1 - o0 - Q(i, j);
        if (std::isinf(o1)) {
          delta[i] = std::numeric_limits<double>::infinity();
        }
      }
    }

    // For all groups, set outlier cells to 0
    for (int g = 0; g < N; g++) {
      arma::uvec ind = arma::find(groups == (g+1));
      arma::vec sorted_delta = arma::sort(delta(ind)); // Sort the values
      double cutoff = std::max(0.0, sorted_delta(h(g) - 1)); // Index into the sorted vector
      for (auto i : ind) {
        if (delta[i] > cutoff) {
          W(i, j) = 0;
        } else {
          W(i, j) = 1;
        }
      }
    }
  }

  return W;  // Return the updated W matrix
}


// Function to calculate probabilities
// [[Rcpp::export]]
arma::mat probabs(const arma::mat& X,
                  const Rcpp::List& mu,
                  const Rcpp::List& Sigma,
                  const arma::mat& pi_groups,
                  const arma::mat& W,
                  const arma::uvec& groups) {

  int n = X.n_rows;
  int N = Sigma.size();

  arma::mat probs(n, N, arma::fill::zeros);

  // Get unique rows of W
  arma::mat W_pattern = find_unique_rows(W);  // Get unique rows
  int wn = W_pattern.n_rows;

  for (int l = 0; l < wn; l++) {
    arma::uvec obs = find(W_pattern.row(l) == 1);
    arma::uvec mis = find(W_pattern.row(l) == 0);

    // Find indices of observations with the same missing pattern
    arma::uvec ind_w;
    for (arma::uword i = 0; i < W.n_rows; i++) {
      if (accu(abs(W.row(i) - W_pattern.row(l))) == 0) {
        ind_w.insert_rows(ind_w.n_rows, arma::uvec{(unsigned int)i});
      }
    }
    // Precompute inverses of Sigma submatrices if there are observed values
    std::vector<arma::mat> Sigmai_tmp(N);
    if (!obs.empty()) {
      // Compute the inverse of Sigma for each group
      for (int k = 0; k < N; k++) {
        // Convert the list element to an arma::mat
        arma::mat Sigma_k = Rcpp::as<arma::mat>(Sigma[k]); // Convert list element Sigma[k] to matrix
        Sigmai_tmp[k] = arma::inv(Sigma_k.submat(obs, obs));  // Inverse of Sigma for observed indices
      }
    }

    for (auto i : ind_w) {
      int g = groups[i] - 1;  // Adjust for 0-based indexing

      if (obs.empty()) {
        // No observed variables, use group probability
        probs.row(i) = pi_groups.row(g);
        Rcpp::Rcout << "No observed variables in probabs function" << std::endl;
      } else {
        arma::vec count(N, arma::fill::zeros);
        for (int k = 0; k < N; k++) {
          arma::rowvec Xi = X.row(i);  // Extract the i-th row from X
          arma::rowvec Xi_obs = Xi.cols(obs);  // Extract the observed columns of Xi
          arma::rowvec mu_k = mu[k];  // Extract the mean vector for group k
          arma::rowvec mu_k_obs = mu_k.cols(obs);  // Extract the observed columns of mu_k
          count[k] = pi_groups(g, k) * density_mnorm(Xi_obs.t(), mu_k_obs.t(), Sigmai_tmp[k]);
        }
        double devi = sum(count);
        if (devi == 0) {
          probs.row(i) = pi_groups.row(g);
          Rcpp::Rcout << "Class-probability p_ik^g of observation i = " << i << " cannot be calculated exactly (Division by 0)." << std::endl;
        } else {
          probs.row(i) = trans(count / devi);
        }
      }
    }
  }

  return probs;
}




// Function to compute pi_groups
// [[Rcpp::export]]
arma::mat pis(const arma::mat& probs,
              const arma::uvec& groups,
              double alpha = 0.5) {

  int N = probs.n_cols;
  arma::mat pi_groups(N, N, arma::fill::none); // Initialize with NaN values

  for (int g = 0; g < N; g++) {
    arma::uvec indices = find(groups == (g + 1));

    if (!indices.empty()) {
      arma::mat probs_r = probs.rows(indices);  // Extract the i-th row from X
      arma::colvec probs_rc = probs_r.col(g);  // Extract the observed columns of Xi

      double Nk = sum(probs_rc);
      pi_groups(g, g) = std::max(alpha, Nk / indices.n_elem);
    }

    for (int g = 0; g < N; g++) {
      for (int l = 0; l < N; l++) {
        if(!(l==g)){
          arma::uvec indices = find(groups == (g + 1));

          if (!indices.empty()) {
            arma::mat probs_r = probs.rows(indices);  // Extract the i-th row from X
            arma::colvec probs_rc = probs_r.col(g);  // Extract the observed columns of Xi
            double Ng = sum(probs_rc);

            arma::mat probs_rl = probs.rows(indices);  // Extract the i-th row from X
            arma::colvec probs_rcl = probs_rl.col(l);  // Extract the observed columns of Xi
            double Nl = sum(probs_rcl);

            if (1 - pi_groups(g, g) != 0) {
              pi_groups(g, l) = (1 - pi_groups(g, g)) * (Nl / indices.n_elem) / (1 - Ng / indices.n_elem);
            } else {
              pi_groups(g, l) = 0.0;
            }
          }
        }
      }
    }
  }
  return pi_groups;
}

