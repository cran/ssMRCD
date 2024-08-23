# CHECK FUNCTION INPUTS AND OUTPUTS


check_input = function(object, category){

  if(category == "vector"){
    stopifnot(is.vector(object), is.numeric(object))
  }
  if(category == "matrix"){
    if(is.data.frame(object)) {
      object = as.matrix.data.frame(object)
    }
    stopifnot(is.matrix(object))
  }
  if(category == "numeric"){
    stopifnot(is.numeric(object))
  }
  if(category == "list"){
    stopifnot(is.list(object))
  }
  if(category == "scalar"){
    stopifnot(is.numeric(object), length(object) == 1)
  }
  if(category == "alpha"){
    stopifnot(is.numeric(object), length(object) == 1, object >= 0.5, object <= 1)
  }
  if(category == "lambda"){
    stopifnot(is.numeric(object), length(object) == 1, object >= 0, object <= 1)
  }
  if(category == "valid_N_vec"){
    stopifnot(is.vector(object), as.numeric(as.factor(object)) == object)
  }
  if(category == "W"){
    stopifnot(is.matrix(object),
              is.numeric(object),
              sum(diag(object)) == 0,
              all.equal(object/rowSums(object), object),
              all(object >= 0))
  }
}


check_locOut = function(object){
  # check local outlier object

  stopifnot("locOuts" %in% class(object))
  stopifnot(names(object) == c("outliers", "next_distance", "cutoff", "coords",
                               "data", "groups", "k", "dist", "centersN",
                               "matneighbor", "ssMRCD"))
  stopifnot(is.numeric(object$outliers) & is.vector(object$outliers))
  stopifnot(is.numeric(object$next_distance) & is.vector(object$next_distance))
  stopifnot(is.numeric(object$cutoff) & length(object$cutoff))
  stopifnot(is.numeric(object$coords) & is.matrix(object$coords))
  stopifnot(is.numeric(object$data) & is.matrix(object$data))
  stopifnot(dim(object$data)[1] == dim(object$coords)[1] & dim(object$data)[1] ==length(object$next_distance))
  stopifnot(dim(object$matneighbor)[1] == dim(object$data)[1] & dim(object$matneighbor)[2] == dim(object$data)[1])
}
