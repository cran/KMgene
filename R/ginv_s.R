ginv_s <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xe <- eigen(X,symmetric=TRUE)
  Positive <- Xe$values > max(tol * Xe$values[1L], 0)
  if (all(Positive)) 
    #       Xe$v %*% (1/Xe$d * t(Xe$u))
    Xe$vectors %*% (1/Xe$values * t(Xe$vectors)) # 2014-10-07 - HC
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xe$vectors[, Positive, drop = FALSE] %*% ((1/Xe$values[Positive]) * 
                                                   t(Xe$vectors[, Positive, drop = FALSE]))
}
