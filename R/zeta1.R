#' Multivariate dependence measure
#'
#' @param X a numeric matrix of dimension d indicating the conditioning (predictor) variables
#' @param Y a numeric vector indicating the response (uni-variate)
#' @param ties.correction logical indicating if the dependence measure uses a correction term w.r.t. ties mentioned in the paper. Default = FALSE.
#' @param approx logical indicating if an approximated version of the dependence measure is computed (fast version). It should be only applied on data with no (or only a few) ties. Default = FALSE.
#' @param resolution an integer indicating the resolution N of the checkerboard aggregation. Default = NULL, which is recommended (N(n) = floor(n^(1/(d+1)))).
#'
#' @description Given a d-dimensional random vector X and a uni-variate random variable Y, the multivariate copula-based dependence measure
#' \eqn{\zeta}1 (zeta 1) quantifies the influence of the random vector X on the random variable Y.
#' More precisely, \eqn{\zeta}1 fulfills the following properties:
#' \itemize{
#'    \item[N] \eqn{\zeta}1(X,Y) attains values in [0,1] (normalization).
#'    \item[I] \eqn{\zeta}1(X,Y) = 0 if and only if X and Y are independent (independence).
#'    \item[C] \eqn{\zeta}1(X,Y) = 1 if and only if Y is a function of X (complete dependence).
#'    \item[S] \eqn{\zeta}1(X,Y) is scale-invariant.
#'    \item[IG] \eqn{\zeta}1(X,Y) fulfills the information gain inequality, i.e., adding variables to X can not decrease the dependence measure.
#' }
#' For more information, see Griessenberger, Junker, Trutschnig (2021), On a multivariate copula-based dependence measure and its estimation, <arXiv:2109.12883> .
#'
#' The derived checkerboard estimator is proved to be strongly consistent and implemented in the function zeta1. Note, that the estimator only attains positive values (within [0,1]).
#' Therefore, interpretation of low values have to be done with care and always under consideration of the sample size.
#'
#' @return A numeric value indicating the dependence between X and Y (or, equivalently, the influence of X on Y).
#'
#' @references Griessenberger, Junker, Trutschnig (2021), On a multivariate copula-based dependence measure and its estimation, <arXiv:2109.12883>.
#'
#' @examples
#' #(complete dependence for dimension 4)
#' n <- 300
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- x1 + x2 + rnorm(n)
#' y <- x1 + x2 + x3
#' zeta1(X = cbind(x1,x2,x3), Y = y)
#' zeta1(X = cbind(x1,x2,x3), Y = y, approx = TRUE)
#'
#' #(independence for dimension 4)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- x1 + x2 + rnorm(n)
#' y <- runif(n)
#' zeta1(X = cbind(x1,x2,x3), Y = y)
#' zeta1(X = cbind(x1,x2,x3), Y = y, approx = TRUE)
#'
#' #(binary output for dimension 3)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- ifelse(x1 + x2 < 1, 0, 1)
#' zeta1(X = cbind(x1,x2), Y = y, ties.correction = FALSE)
#' zeta1(X = cbind(x1,x2), Y = y, ties.correction = TRUE)
#'

zeta1 <- function(X, Y, ties.correction = FALSE, approx = FALSE, resolution = NULL){
  if(approx){
    res <- .zeta1_approx(X,Y,N = resolution)
    return(c(zeta1 = res))
  }else{
    res <- .zeta1(X, Y, resolution = resolution, ties.correction = ties.correction)
    return(c(zeta1 = res))
  }
}
