# get the cluster result from the estimated interaction matrices from HGMND
getCluster <- function(est.HGMND, method = "F", tol = 1e-5) {

  if (class(est.HGMND) != "est.HGMND") stop("theta must be an object of class \"est.HGMND\" from function HGMND.")

  theta       <- est.HGMND[["Theta"]]
  M           <- dim(theta)[3]
  mat.compare <- matrix(0, nrow = M, ncol = M)

  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      mat.1             <- theta[,,i] - diag(diag(theta[,,i]))
      mat.2             <- theta[,,j] - diag(diag(theta[,,j]))
      mat.compare[i, j] <- norm(mat.1 - mat.2, type = method) <= tol
    }
  }
  tri.comapre       <- mat.compare
  diag(tri.comapre) <- 1
  mat.compare       <- mat.compare + t(mat.compare)
  diag(mat.compare) <- 1


  candidate   <- 1:M
  est.cluster <- rep(0, M)
  t           <- 1
  while (length(candidate) > 0) {

    temp.cluster <- which(tri.comapre[min(candidate),] == 1)
    candidate    <- setdiff(candidate, temp.cluster)

    est.cluster[temp.cluster] <- t

    t <- t + 1
  }


  cluster.HGMND <- list(mat.compare = mat.compare, est.cluster = est.cluster)
  class(cluster.HGMND) <- "cluster.HGMND"
  return(cluster.HGMND)

}


