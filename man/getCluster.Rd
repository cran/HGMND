\name{getCluster}
\alias{getCluster}
\title{Get the cluster structure of the HGMND estimate}
\description{
After estimating the conditional dependence matrices of the multiple datasets using the HGMND method, the cluster structure can be revealed by comparison of these matrices.
}
\usage{
  getCluster(est.HGMND, method = "F", tol = 1e-5)
}
\arguments{
  \item{est.HGMND}{a list, the result of the function \code{HGMND} with class "est.HGMND".}

  \item{method}{the method of evaluating the difference of two conditional dependence matrices. The function \code{norm} from the \code{base} package is used to calculate the matrix norm of the element-wise difference of two matrices. It must be chosen from \code{"O"}, \code{"I"}, \code{"F"}, \code{"M"}, \code{"2"}, corresponding to the same settings in the function \code{norm}. Default to \code{"F"}, the Frobenius norm.}

  \item{tol}{tolerance in evaluating the difference of two conditional dependence matrices. If the calculated difference is no larger than \code{tol}, they are regarded as in one cluster. Default to \code{1e-5}.}

}
\value{
  the function \code{getCluster} returns the clustering structure of the multiple conditional dependence matrices.

  \item{mat.comapre}{a matrix of 0 or 1. If the element on the \eqn{i}th row and \eqn{j}th column of the matrix is 1, the \eqn{i}th and the \eqn{j}th conditional dependence matrices are in the same cluster, 0 otherwise.}

  \item{est.cluster}{a vector with length same as the number of conditional dependence matrices indicating the cluster label of each matrix.}

}
\examples{
# This is an example of HGMND with simulated data
data(HGMND_SimuData)
h              <- genscore::get_h_hp("mcp", 1, 5)
HGMND_SimuData <- lapply(HGMND_SimuData, function(x) scale(x, center = FALSE))
mat.chain      <- diag(length(HGMND_SimuData))
diag(mat.chain[-nrow(mat.chain), -1]) <- 1

result <- HGMND(x        = HGMND_SimuData,
                setting  = "gaussian",
                h        = h,
                centered = FALSE,
                mat.adj  = mat.chain,
                lambda1  = 0.086,
                lambda2  = 3.6,
                gamma    = 1,
                tol      = 1e-3,
                silent  = TRUE)
Theta       <- result[["Theta"]]
res.cluster <- getCluster(result)
}
