#' Compute empirical checkerboard copula in arbitrary dimension
#'
#' @param X a numeric matrix of dimension rho indicating a sample of rho variables
#' @param resolution an integer indicating the resolution N of the checkerboard copula
#' @param bin.size either "fixed" or "adaptive", indicating whether the checkerboard copula may vary its bin sizes (defaults to "fixed")
#'
#' @description The function ECBC computes the mass distribution of the empirical (checkerboard) copula,
#' given a rho-dimensional sample X. If resolution equals sample size, the bi-linearly extended empirical copula is returned.
#' Note, if there are ties in the sample an adjusted empirical copula is calculated.
#' If bin.size is set to "adaptive" the sizes of the bins will be adjusted to fit the data without overspilling into neighboring bins.
#' This might affects the result, but is more efficient with samples having many ties as no adjustment is needed.
#'
#' @return array of dimension resolution^rho.
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- x1 + x2 + rnorm(n)
#' M <- ECBC(X = cbind(x1,x2,y), resolution = 8)
#'

ECBC <- function(X, resolution, bin.size = "fixed") {
    U <- apply(X, 2, qmdrank)
    if (bin.size == "adaptive") return(.EACBC(U, resolution))
    else return(.ECBC(U, resolution))
}
