# HELPING FUNCTIONS FOR APPLICATIONS


#' Restructure Data Matrix as List
#'
#' This function restructures neighborhood information given by a data matrix
#' containing all information and one neighborhood assignment vector. It returns a list
#' of data matrices used in \code{\link[ssMRCD]{ssMRCD}}.
#'
#' @param data data matrix with all observations.
#' @param neighborhood_vec numeric neighborhood assignment vector.
#'                         Should contain numbers from \code{1} to \code{N} and not leave integers out.
#'
#' @return Returns a list containing the observations per neighborhood assignment.
#'
#' @examples
#'
#' # data matrix
#' data = matrix(rnorm(n = 3000), ncol = 3)
#' N_assign = sample(x = 1:10, size = 1000, replace = TRUE)
#'
#' restructure_as_list(data, N_assign)
#'
#' @export
restructure_as_list = function(data, neighborhood_vec){

  data = as.matrix(data)
  check_input(data, "matrix")
  check_input(neighborhood_vec, "vector")
  check_input(neighborhood_vec, "valid_N_vec")

  p = dim(data)[2]
  N = length(unique(neighborhood_vec))

  x = list()
  for(i in 1:N){
    x = append(x, list(data[neighborhood_vec == i, 1:p]))
  }

  return(x)
}


#' Rescale Weight Matrix
#'
#' Given a matrix with values for neighborhood influences the function rescales
#' the matrix in order to get an appropriate weight matrix used for the function \code{\link[ssMRCD]{ssMRCD}}.
#'
#' @param W weight matrix with diagonals equal to zero and at least one positive entry per row.
#'
#' @return An appropriately scaled weight matrix.
#'
#' @seealso \code{\link[ssMRCD]{ssMRCD}, \link[ssMRCD]{local_outliers_ssMRCD}, \link[ssMRCD]{geo_weights}}
#'
#' @examples
#'
#' W = matrix(c(0, 1, 2,
#'              1, 0, 1,
#'              2, 1, 0), nrow = 3)
#' rescale_weights(W)
#'
#' @export
rescale_weights = function(W){

  W = as.matrix(W)
  stopifnot(all.equal(unname(diag(W)), rep(0, dim(W)[1])))

  W = W * (rowSums(W)^(-1))
  if(any(!is.finite(W))){
    stop("No finite weight matrix possible. There should be at least one positive value per row.")
  }
  check_input(W, "W")

  return(W)
}



#' Inverse Geographic Weight Matrix
#'
#' Calculates a inverse-distance based weight matrix for the function \code{\link[ssMRCD]{ssMRCD}} (see details).
#'
#' @param coordinates matrix of coordinates of observations.
#' @param N_assignments vector of neighborhood assignments.
#'
#' @details
#' First, the centers (means of the coordinates given) \eqn{c_i} of each neighborhood is calculated.
#' Then, the Euclidean distance between the centers is calculated and the weight is based on
#' the inverse distance between two neighborhoods, \deqn{w_{ij} = \frac{1}{dist(c_i, c_j)}. }
#' It is scaled according to a weight matrix.
#'
#'@return Returns a weighting matrix \code{W} and the coordinates of the centers per neighborhood \code{centersN}.
#'
#' @seealso \code{\link[ssMRCD]{rescale_weights}}
#'
#' @examples
#' coordinates = matrix(rnorm(1000), ncol = 2, nrow = 500)
#' N_ass = sample(1:5, 500, replace = TRUE)
#'
#' geo_weights(coordinates, N_ass)
#'
#' @export
geo_weights = function(coordinates, N_assignments){

  coordinates = as.matrix(coordinates)
  check_input(N_assignments, "vector")

  N = length(unique(N_assignments))
  p_coord = dim(coordinates)[2]

  centersN = matrix(NA, N, p_coord)
  for(i in 1:N) {
    tmp = coordinates[N_assignments == i,]
    centersN[i, ] = colMeans(tmp)
  }

  W = as.matrix(stats::dist(centersN)^(-1))
  W = rescale_weights(W)

  return(list(W = W, centersN = centersN))
}



#' Creates Grid-Based Neighborhood Structure
#'
#' This function creates a grid-based neighborhood structure for the \code{\link[ssMRCD]{ssMRCD}} function using cut-off values for two coordinate axis.
#'
#' @param x vector of first coordinate of data set.
#' @param y vector of second coordinate of data set.
#' @param cutx cut-offs for first coordinate.
#' @param cuty cut-offs for second coordinate.
#'
#' @return Returns a neighborhood assignment vector for the coordinates \code{x} and \code{y}.
#' @export
#'
#' @examples
#' # get data
#' data(weatherAUT2021)
#'
#' # set cut-off values
#' cut_lon = c(9:16, 18)
#' cut_lat = c(46, 47, 47.5, 48, 49)
#'
#' # create neighborhood assignments
#' N_structure_gridbased(weatherAUT2021$lon,
#'                       weatherAUT2021$lat,
#'                       cut_lon,
#'                       cut_lat)

N_structure_gridbased = function(x, y, cutx, cuty){

  check_input(x, "vector")
  check_input(y, "vector")
  check_input(cutx, "vector")
  check_input(cuty, "vector")

  N = c()
  Nvec = seq(1,(length(cutx)-1)*(length(cuty)-1))
  N_matrix = matrix(Nvec, nrow = length(cuty)-1, ncol = length(cutx)-1)
  for(i in 1:length(x)){
    xi = sum(x[i] >= cutx)
    yi = sum(y[i] >= cuty)
    N[i] = N_matrix[yi, xi]
  }

  N = as.numeric(as.factor(N))
  return(N)
}


