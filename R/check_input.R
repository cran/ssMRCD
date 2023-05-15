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
  if(category == "ssMRCD"){
    check_ssMRCD(object)

  }

}

check_ssMRCD = function(object){

  stopifnot(is.list(object),
            "ssMRCD" %in% class(object))
  stopifnot(names(object) == c("MRCDcov", "MRCDicov", "MRCDmu", "mX", "N", "mT", "rho", "alpha", "h",
                               "numiter", "hset", "c_alpha", "weights", "lambda", "obj_fun_values", "best6pack", "Kcov"))
  # lambda
  stopifnot(is.numeric(object$lambda),
            length(object$lambda) == 1,
            object$lambda <= 1,
            object$lambda >= 0)
  # N
  stopifnot(is.numeric(object$N),
            length(object$N) == 1,
            object$N > 0)
  # TM
  stopifnot(is.matrix(object$mT),
            dim(object$mT)[1] == dim(object$mT)[2],
            object$mT == t(object$mT))
  # cov [[1]]
  stopifnot(is.matrix(object$MRCDcov[[1]]),
            dim(object$MRCDcov[[1]])[1] == dim(object$MRCDcov[[1]])[2],
            all.equal(object$MRCDcov[[1]], t(object$MRCDcov[[1]])))
  # icov[[1]]
  stopifnot(is.matrix(object$MRCDicov[[1]]),
            dim(object$MRCDicov[[1]])[1] == dim(object$MRCDicov[[1]])[2],
            all.equal(object$MRCDicov[[1]], t(object$MRCDicov[[1]])))
  # icov and cov [[1]]
  stopifnot(is.matrix(object$MRCDcov[[1]]),
            is.matrix(object$MRCDicov[[1]]),
            dim(object$MRCDicov[[1]])[1] == dim(object$MRCDcov[[1]])[2],
            all.equal(object$MRCDicov[[1]] %*% object$MRCDcov[[1]], diag(1, dim(object$MRCDcov[[1]]))))
  # dims cov, mu, x, tm
  stopifnot(dim(object$MRCDcov[[1]])[1] == dim(object$MRCDmu[[1]])[1],
            dim(object$MRCDcov[[1]])[1] == dim(object$mX[[2]])[2],
            dim(object$MRCDcov[[1]])[1] == dim(object$mT[[1]])[1])
  # dim weights, N, rho
  stopifnot(dim(object$weights) == object$N,
            object$N == length(object$c_alpha),
            object$N == length(object$rho),
            object$N == length(object$hset),
            object$N == length(object$MRCDmu),
            object$N == length(object$MRCDcov),
            object$N == length(object$MRCDicov))
}


check_locOut = function(object){
  # check local outlier object

  stopifnot("locOuts" %in% class(object))
  stopifnot(names(object) == c("outliers", "next_distance", "cutoff", "coords",
                               "data", "N_assignments", "k", "dist", "centersN",
                               "matneighbor", "ssMRCD"))
  check_ssMRCD(object$ssMRCD)
  stopifnot(is.numeric(object$outliers) & is.vector(object$outliers))
  stopifnot(is.numeric(object$next_distance) & is.vector(object$next_distance))
  stopifnot(is.numeric(object$cutoff) & length(object$cutoff))
  stopifnot(is.numeric(object$coords) & is.matrix(object$coords))
  stopifnot(is.numeric(object$data) & is.matrix(object$data))
  stopifnot(dim(object$data)[1] == dim(object$coords)[1] & dim(object$data)[1] ==length(object$next_distance))
  stopifnot(dim(object$matneighbor)[1] == dim(object$data)[1] & dim(object$matneighbor)[2] == dim(object$data)[1])
}
