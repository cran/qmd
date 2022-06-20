#' Multivariate dependence measure
#'
#' @param X a numeric matrix or data.frame of dimension d containing the explanatory variables
#' @param y a numeric vector containing the uni-variate response variable
#' @param ties.correction logical indicating if the measure of dependence should be calculated with ties-correction (experimental version). Default = FALSE.
#' @param bin.size either "fixed", "adaptive" or "sparse.adaptive", indicating whether the checkerboard copula may vary its bin sizes (defaults to "fixed"). Setting this to "adaptive" might affect the results but will be faster if the sample has many ties.
#' @param resolution an integer indicating the resolution N of the checkerboard aggregation. We recommend to use the default configuration (resolution = NULL), which uses the resolution N(n) = floor(n^(1/(d+1))), where d denotes the number of explanatory variables.
#'
#' @description Function for estimating the non-parametric copula-based multivariate measure of dependence \eqn{\zeta}1.
#' This measure quantifies the extent of dependence between a d-dimensional random vector X and a uni-variate random variable y (i.e.,
#' it measures the influence of d explanatory variables X1,...,Xd on a univariate variable y).
#'
#' @return A numeric value indicating the extent of dependence between the vector X and the variable y (or, equivalently, the influence of X on y).
#'
#'@details see function qmd(...).
#'
#' @references
#' Griessenberger, F., Junker, R.R. and Trutschnig, W. (2022). On a multivariate copula-based dependence measure and its estimation, Electronic Journal of Statistics, 16, 2206-2251.
#'
#' @examples
#' #(complete dependence for dimension 4)
#' n <- 300
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- x1 + x2 + rnorm(n)
#' y <- x1 + x2 + x3
#' zeta1(X = cbind(x1,x2,x3), y = y)
#'
#' #(independence for dimension 4)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- x1 + x2 + rnorm(n)
#' y <- runif(n)
#' zeta1(X = cbind(x1,x2,x3), y = y)
#'
#' #(binary output for dimension 3)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- ifelse(x1 + x2 < 1, 0, 1)
#' zeta1(X = cbind(x1,x2), y = y)


zeta1 <- function(X, y, ties.correction = FALSE, bin.size = "fixed", resolution = NULL){
  res <- .zeta1(X, y, resolution = resolution, ties.correction = ties.correction, bin.size = bin.size)
  return(c(zeta1 = res))
}

