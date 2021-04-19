# turn result of group fused lasso into origin matrix form of Z
.constructTensor <- function(A, P){
  M      <- dim(A)[2]
  result <- array(0, dim = c(P, P, M))
  for (m in 1:M) {
    temp                  <- matrix(0, nrow = P, ncol = P)
    temp[lower.tri(temp)] <- A[, m]
    result[,, m]          <- temp + t(temp)
  }
  return(result)
}


# form a matrix from the upper-triangular elements of Z
.constructMatrix <- function(A, P){
  M   <- dim(A)[3]
  A.m <- matrix(0, nrow = M, ncol = P*(P - 1)/2)
  for (t in 1:M) {
    A.m[t,] <- A[,, t][lower.tri(A[,, t])]
  }
  return(A.m)
}


# group lars
.block.optimize <- function(beta, AS, lambda, C, n, X, tol = 1e-8, maxit = 1e5){

  a <- dim(beta)[1]
  p <- dim(beta)[2]

  norm.beta <- sqrt(rowSums(beta^2))

  tol  <- tol*p
  gain <- 2*tol*rep(1, a)

  mean.X <- apply(X, 2, function(x) scale(x, center = TRUE, scale = FALSE))

  # main loop
  if(a > 0){

    # optimize each block in turn
    itc <- 0 # iteration count

    while((any(gain > tol)) && (itc < maxit)){

      i    <- itc %% a + 1 # index of block to update in AS
      AS.i <- AS[i]        # block to update

      XitX <- (matrix(mean.X[, AS.i], nrow = 1) %*% mean.X)[AS]
      gamma.i <- XitX[i]

      indices.without.i <- c(seq_len(i-1), seq_len((a-i-1>0)*(a-i-1))+i)

      # compute the vector S
      S <- matrix(C, ncol = p)[i,] - matrix(XitX[indices.without.i], nrow = 1) %*% beta[indices.without.i,]

      # check the norm of S to decide where the optimum is
      norm.S <- svd(S)[['d']][1]
      if(norm.S<lambda){
        beta.new <-  matrix(0, nrow = 1, ncol = p)
      } else {
        beta.new <- matrix(S * (1-lambda/norm.S)/gamma.i, nrow = 1)
      }

      norm.beta.new <- norm(beta.new, type = 'f')
      gain[i]       <- (gamma.i*(norm.beta[i] + norm.beta.new)/2 + lambda) * (norm.beta[i]-norm.beta.new) +
        S %*% (t(beta.new) - matrix(beta[i,], ncol = 1))

      # update beta
      beta[i,]     <- beta.new
      norm.beta[i] <- norm.beta.new

      itc <- itc + 1

    } # end of while
  } # end of if
  return(beta)
}


# group fused lasso solver
.gflasso <- function(Y, X, lambda){
  maxit  <- 1e2  # for convergence
  tol    <- 1e-3
  Y.mean <- colMeans(Y)
  n      <- dim(Y)[1]
  p      <- dim(Y)[2]

  R   <- apply(X, 2, function(x) scale(x, center = T, scale = F))
  XtY <- t(R) %*% Y

  beta <- matrix(nrow = 0, ncol = p)
  AS   <- NULL

  globalSol <- 0

  while (!globalSol) {

    beta <- .block.optimize(beta, AS, lambda, XtY[AS,], n, X, tol, maxit)
    beta <- matrix(beta, ncol = p)

    # remove from active set the zero coefficients
    nonzero.coef <- which(rowSums(beta^2) != 0)
    AS           <- AS[nonzero.coef]
    beta         <- matrix(beta[nonzero.coef,], ncol = p)

    # check optimality
    temp.beta <- matrix(0, nrow = n - 1, ncol = p)
    temp.beta[AS,] <- beta

    S      <- XtY - t(R) %*% R %*% temp.beta
    norm.S <- rowSums(S^2)

    if(is.null(AS)){

      max.norm.S <- max(norm.S)
      i.max      <- which.max(norm.S)

      if(max.norm.S < lambda^2 + tol){
        # the optimal solution is the null matrix
        globalSol <- 1
      } else {
        # this should only occur at the first iteration
        AS   <- i.max
        beta <- matrix(0, nrow = 1, ncol = p)
      }
    } else {

      # At optimality we must have normS(i)=lambda^2 for i in AS and normS(i)<lambda^2 for i not in AS
      lagr  <- max(norm.S[AS])
      lagr  <- min(c(lagr, lambda^2))
      nonas <- setdiff(1:(n-1), AS)

      # check optimality conditions for blocks not in the active set
      y <- max(norm.S[nonas])
      i <- which.max(norm.S[nonas])

      if(is.null(nonas) || (y < lagr + tol)){
        # optimality conditions are fulfilled, we have found the global solution
        globalSol = 1
      } else {

        # otherwise we add the block that violates most the optimality condition
        AS   <- c(AS, nonas[i])
        beta <- rbind(beta, rep(0, p))
      }
    }
  } # end of while

  result <- list(active = AS, beta = beta, mean = Y.mean, lambda = lambda)
  return(result)
}





HGMND <- function(x, setting, h, centered, mat.adj, lambda1, lambda2, gamma = 1, maxit = 200, tol = 1e-5, silent = TRUE) {
  # check input
  if (!is.function(h)) stop("h must be a function returning the value of h and its derivative with data.")
  if (!is.list(x)) stop("x must be a list containing multiple matrices of datasets.")
  if (!length(unique(unlist(lapply(x, ncol)))) == 1) stop("each dataset in x must have same number of columns.")
  if (!setting %in% c("gaussian", "gamma", "exp")) stop("dist must be chosen from \"gaussian\", \"gamma\", and \"exp\".")

  # preparations --------------------------------------------------
  maxsub   <- 100  # Dykstra maximum number of iterations
  tolDual  <- tol  # tolerance for convergence of ADMM
  tolPrime <- tol
  tolDyk   <- 1e-3

  M <- length(x)
  P <- ncol(x[[1]])

  # reconstruction matrix in group fused lasso part
  if (nrow(mat.adj) != M | ncol(mat.adj) != M) stop("mat.adj must be a square matrix with size the same as the number of datasets.")
  if (!all(mat.adj %in% c(0, 1))) stop("mat.adj must contain only 0s and 1s.")
  if (0 %in% rowSums(mat.adj)) stop("the network among the multiple datasets must be connected.")

  adjMat                    <- matrix(0, nrow = M, ncol = M)
  adjMat[upper.tri(adjMat)] <- mat.adj[upper.tri(mat.adj)]
  pos.edge <- which(adjMat == 1, arr.ind = TRUE)
  mat.D    <- matrix(0, nrow = nrow(pos.edge), ncol = M)
  mat.D[cbind(1:nrow(pos.edge), pos.edge[, 1])] <- -1
  mat.D[cbind(1:nrow(pos.edge), pos.edge[, 2])] <-  1
  mat.D    <- rbind(c(1, rep(0, M-1)), mat.D)
  mat.R    <- solve(mat.D)[, -1]

  # calculate Gamma and g in the loss function
  domain        <- make_domain("R+", p = P)
  mat.g         <- array(0, dim = c(P  , P, M))
  mat.Gamma     <- array(0, dim = c(P^2, P, M))
  mat.Gamma.inv <- array(0, dim = c(P^2, P, M))

  for (m in 1:M) {
    elts <- get_elts(h_hp                  = h,
                     x                     = x[[m]],
                     setting               = setting,
                     domain                = domain,
                     centered              = centered,
                     profiled_if_noncenter = TRUE,
                     scale                 = "",
                     diagonal_multiplier   = 1)

    mat.g[,,m]          <- matrix(elts$g_K, ncol = P)
    mat.Gamma[,,m]      <- t(elts$Gamma_K)
    mat.Gamma.inv[,, m] <- do.call(rbind, lapply(1:P, function(j) solve(mat.Gamma[((j-1)*P + 1):(j*P),, m] + gamma*diag(P))))
  }

  # initialize auxiliary and dual variables of ADMM
  Z     <- array(0, dim = c(P, P, M))
  U     <- Z
  Theta <- Z

  # initialize active set and solutions for fused group lasso
  active <- NULL
  beta   <- matrix(nrow = 0, ncol = P)

  aditer   <- 1            # iteration count of ADMM
  difPrime <- tolPrime + 1 # difference of Theta
  difDual  <- tolDual  + 1 # difference of Z


  # ADMM --------------------------------------------------
  while ((difPrime > tolPrime) && (difDual > tolDual) && (aditer < maxit)) {
    # * Step 1: Update Theta ===========================
    tic <- Sys.time()
    for (m in 1:M) {
      Theta[,,m] <- sapply(1:P,
                           function(j) mat.Gamma.inv[((j-1)*P + 1):(j*P),, m] %*% (mat.g[, j, m] + gamma * (Z[,,m] - U[,,m])[,j]))
      Theta[,,m] <- 1/2*(Theta[,,m] + t(Theta[,,m]))
    }
    time.step1 <- Sys.time() - tic


    # * Step 2: Update Z ===========================
    Z.old <- Z
    A     <- .constructMatrix(U + Theta, P) # reconstruct matrix for step 2

    if (lambda2 > 0) {
      # initialize Dykstra Algorithm
      X.k    <- A
      difDyk <- 0
      P.k    <- matrix(0, nrow = M, ncol = P*(P-1)/2)
      Q.k    <- P.k
      n.iter <- 1 # iteration count of Dykstra

      while ((difDyk>tolDyk) && (n.iter<maxsub) || (n.iter==1)) {

        # ** group fused step ###################
        result <- .gflasso(X.k + P.k, mat.R, lambda2/gamma)
        beta   <- result[['beta']]
        active <- result[['active']]
        Beta   <- matrix(0, nrow = M-1, ncol = dim(beta)[2])
        if(!is.null(active)) Beta[active,] <- beta[1:length(active),]
        offSet <- matrix(1, nrow = 1, ncol = M) %*% (X.k + P.k - mat.R %*% Beta)/M

        # reconstruct signal
        Y.k <- matrix(1, nrow = M, ncol = 1) %*% offSet + mat.R %*% Beta

        # update the 1st auxiliary variable
        P.k <- P.k + X.k - Y.k

        # ** lasso part ###################
        X.old <- X.k

        # soft threshold
        X.k   <- ifelse(abs(Y.k + Q.k) - lambda1/gamma > 0,
                        sign(Y.k + Q.k)*(abs(Y.k + Q.k) - lambda1/gamma),
                        0)

        # update th 2nd auxiliary variable
        Q.k <- Q.k + Y.k - X.k

        # update iteration count
        n.iter <- n.iter + 1

        difDyk <- svd(X.k - X.old)[['d']][1]
      } # end of while

      Z.m <- X.k

    } else {
      # update of auxiliary without the group fused lasso penalty
      Z.m <- ifelse(abs(A)-lambda1/gamma > 0,
                    sign(A)*(abs(A)-lambda1/gamma),
                    0)
    }

    # update off-diagonal elements
    Z <- .constructTensor(t(Z.m), P)

    # update diagonal elements
    for (m in 1:M) {
      diag(Z[,,m]) <- diag(Theta[,,m] + U[,,m])
    }

    time.step2 <- Sys.time() - tic

    # * Step 3: Update U ===========================
    U.old <- U
    U     <- U.old + (Theta - Z)

    # * Step 4: dual and primal feasibility ===========================
    difPrime <- sum(sapply(1:M, function(m) norm(Theta[,,m] - Z[,,m], 'f')^2))     # prime feasibility
    difDual  <- sum(sapply(1:M, function(m) norm(Z[,,m]     - Z.old[,,m], 'f')^2)) # dual feasibility
    aditer   <- aditer + 1 # update iteration count of ADMM

    if (!silent) {
      print(paste0("HGMND: "    ,
                   "diffPrime: ", round(difPrime, round(-log10(tol))), "; ",
                   "diffDual:  ", round(difDual , round(-log10(tol))), "; ",
                   "ADMM loop: ", aditer))
      print(paste0("Theta update: ", round(time.step1, 3), ";  ", "Z update: ", round(time.step2, 3)))
    }
  } # end of ADMM

  est.HGMND <- list(Theta = Z, M = M, P = P, dataset.edge = mat.D, mat.adj = mat.adj, h = h, centered = centered, lambda1 = lambda1, lambda2 = lambda2)
  class(est.HGMND) <- "est.HGMND"
  return(est.HGMND)
}







