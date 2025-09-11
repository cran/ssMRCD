# HELPING FUNCTIONS FOR APPLICATIONS

restructure_as_list = function(data, groups){

  # Restructure Data Matrix as List
  # This function restructures neighborhood information given by a data matrix
  # containing all information and one neighborhood assignment vector. It returns a list
  # of data matrices used in \code{\link[ssMRCD]{ssMRCD}}.
  # @param data data matrix with all observations.
  # @param groups numeric neighborhood assignment vector.
  # @return Returns a list containing the observations per neighborhood assignment.
  # The list is sorted according to the order of the first appearance in the groups vector.

  data = as.matrix(data)
  groups = as.numeric(as.factor(groups))

  p = dim(data)[2]
  N = length(unique(groups))

  x = list()
  for(i in 1:N){
    x = append(x, list(data[groups == i, 1:p]))
  }

  return(x)
}



#' Inverse Geographic Weight Matrix
#'
#' Calculates a inverse-distance based weight matrix for the function \code{\link[ssMRCD]{ssMRCD}} (see details).
#'
#' @param coordinates matrix of coordinates of observations.
#' @param groups vector of neighborhood groups.
#'
#' @details
#' First, the centers (means of the coordinates given) \eqn{c_i} of each neighborhood is calculated.
#' Then, the Euclidean distance between the centers is calculated and the weight is based on
#' the inverse distance between two neighborhoods, \deqn{w_{ij} = \frac{1}{dist(c_i, c_j)}. }
#' It is scaled according to a weight matrix.
#'
#' @return Returns a weighting matrix \code{W} and the coordinates of the centers per neighborhood \code{centersN}.
#'
#' @seealso \code{\link[ssMRCD]{time_weights}}
#'
#' @examples
#' coordinates = matrix(rnorm(1000), ncol = 2, nrow = 500)
#' groups = sample(1:5, 500, replace = TRUE)
#'
#' geo_weights(coordinates, groups)
#'
#' @export
geo_weights = function(coordinates, groups){

  coordinates = as.matrix(coordinates)

  N = length(unique(groups))
  p_coord = dim(coordinates)[2]
  groups = as.numeric(as.factor(groups))

  centersN = matrix(NA, N, p_coord)
  for(i in 1:N) {
    tmp = coordinates[groups == i,]
    centersN[i, ] = colMeans(tmp)
  }

  W = as.matrix(stats::dist(centersN)^(-1))
  W = W * (rowSums(W)^(-1))

  return(list(W = W, centersN = centersN))
}



#' Band weight matrix for time series groupings
#'
#' @param N number of groups.
#' @param off_diag vector for off-diagonal values unequal to zero.
#'
#' @return Returns weight matrix for time series groups appropriate for \code{\link[ssMRCD]{ssMRCD}}.
#' @export
#'
#' @seealso \code{\link[ssMRCD]{geo_weights}}
#'
#' @examples
#' time_weights(N = 10, off_diag = c(2,1))
#'
time_weights = function(N, off_diag) {
  w = diag(0, N)
  for( i in 1:length(off_diag)){
    diag(w[-c(1:i), -c((N+1-i):N) ]) = off_diag[i]
    diag(w[-c((N+1-i):N), -c(1:i)]) = off_diag[i]
  }
  w =  w * (rowSums(w)^(-1))
  return(w)
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
#' groups_gridbased(weatherAUT2021$lon,
#'                       weatherAUT2021$lat,
#'                       cut_lon,
#'                       cut_lat)

groups_gridbased = function(x, y, cutx, cuty){

  groups = c()
  Nvec = seq(1,(length(cutx)-1)*(length(cuty)-1))
  N_matrix = matrix(Nvec, nrow = length(cuty)-1, ncol = length(cutx)-1)
  for(i in 1:length(x)){
    xi = sum(x[i] >= cutx)
    yi = sum(y[i] >= cuty)
    groups[i] = N_matrix[yi, xi]
  }

  groups = as.numeric(as.factor(groups))
  return(groups)
}
