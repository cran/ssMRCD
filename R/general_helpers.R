# HELPING FUNCTIONS FOR APPLICATIONS




#' Restructure Data Matrix as List
#'
#' This function restructures neighborhood information given by a data matrix
#' containing all information and one neighborhood assignment vector. It returns a list
#' of data matrices used in \code{\link[ssMRCD]{ssMRCD}}.
#'
#' @param data data matrix with all observations.
#' @param groups numeric neighborhood assignment vector.
#'
#' @return Returns a list containing the observations per neighborhood assignment.
#' The list is sorted according to the order of the first appearance in the groups vector.
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
restructure_as_list = function(data, groups){

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

  W = W * (rowSums(W)^(-1))
  if(any(!is.finite(W))){
    stop("No finite weight matrix possible. There should be at least one positive value per row.")
  }

  return(W)
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
#'@return Returns a weighting matrix \code{W} and the coordinates of the centers per neighborhood \code{centersN}.
#'
#' @seealso \code{\link[ssMRCD]{rescale_weights}}
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
  W = rescale_weights(W)

  return(list(W = W, centersN = centersN))
}




# Smoothing structure
#' Band weight matrix for time series groupings
#'
#' @param N number of groups.
#' @param off_diag vector for off-diagonal values unequal to zero.
#'
#' @return Returns weight matrix for time series groups appropriate for \code{\link[ssMRCD]{ssMRCD}}.
#' @export
#'
#' @seealso \code{\link[ssMRCD]{geo_weights}}, \code{\link[ssMRCD]{rescale_weights}}
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
  return(rescale_weights(w))
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


##########################################################################################
colour_to_ansi <- function(colour) {
  # Note ANSI colour codes
  colour_codes <- list(
    "black" = 30,
    "red" = 31,
    "green" = 32,
    "yellow" = 33,
    "blue" = 34,
    "magenta" = 35,
    "cyan" = 36,
    "white" = 37
  )

  # Check colour provided in codes above
  if ((colour %in% names(colour_codes) == FALSE & is.character(colour)) |
      (is.numeric(colour) & (colour > 255 | colour < 0) ) ) {
    stop(
      paste0(
        "Colour provided (", colour, ") can't be converted to ANSI. ",
        "Must be one of: \n", paste(names(colour_codes), collapse = ",")
      )
    )
  }

  # Create ANSI version of colour
  if(is.character(colour)){
    ansi_colour <- paste0("\033[", colour_codes[[colour]], "m")
  }
  if(is.numeric(colour)){
    ansi_colour <- paste0("\033[38;5;", colour, "m")
  }

  return(ansi_colour)
}




coloured_print <- function(text, colour = "green", trace = T) {

  # algo messages
  if(colour == "parsetting_message") colour = "blue"
  if(colour == "general_message") colour = "blue"
  if(colour == "message") colour = "blue"
  # algo problems
  if(colour == "root_finder_issue") colour = 9
  if(colour == "constraint_issue") colour = 9
  if(colour == "convergence_issue") colour = 9
  if(colour == "problem") colour = 9
  # additional information
  if(colour == "trace") colour = 117
  if(colour == "info") colour = 117

  if(trace) cat(colour_to_ansi(colour), text, "\033[0m\n")
}
