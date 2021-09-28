#' Compute empirical checkerboard copula in arbitrary dimension
#'
#' @param X a numeric matrix of dimension d indicating the conditioning (predictor) variables
#' @param Y a numeric vector indicating the response variable
#' @param resolution an integer indicating the resolution N of the checkerboard copula. 
#' 
#' @description The function ECBC computes the mass distribution of the empirical (checkerboard) copula, 
#' given a (d+1)-dimensional sample (X,Y). If resolution equals sample size, the bi-linearly extended empirical copula is returned. 
#' Note, if there are ties in the sample an adjusted empirical copula is calculated. 
#' 
#' @return array of dimension dim(X,Y).
#' 
#' @examples
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- x1 + x2 + rnorm(n)
#' M <- pure_r_ECBC(X = cbind(x1,x2), Y = y, resolution = 8)
#' 

pure_r_ECBC <- function(X, Y, resolution){
  N <- resolution
  n <- NROW(X)
  #Define the overall vector with predictor and response
  Z <- cbind(X,Y)
  #PIT (ranks)
  U <- apply(Z, 2, function(x) rank(x, ties.method = "max"))
  U_range <- apply(U, 2, .pure_r_range)
  #dimension of vector
  rho <- NCOL(U)
  #define results array
  results <- array(0, dim = rep(N,rho))
  for(i in 1:n){
    lambda <- low_upp <- z_range <- list()
    for(j in 1:rho){
      rZ <- U_range[,j][U[i,j]]
      z_range[[j]] <- rZ
      upper = unname(ceiling(U[i,j]/n * N))
      lower = max(ceiling((U[i,j] - rZ)/n * N), 1)
      lambda[[j]] <- as.vector(rep(0,length(lower:upper)))
      low_upp[[j]] <- lower:upper
      k <- 1
      for(z in lower:upper){
        lambda[[j]][k] <- min(U[i,j], z/N * n) - max(U[i,j] - rZ, (z - 1)/N * n)
        k <- k + 1
      }
    }
    
    LAMBDA <- .expand.grid.new(lambda)
    LOW_UPP <- .expand.grid.new(low_upp)
    z_RANGE <- .expand.grid.new(z_range)
    #LAMBDA <- as.matrix(expand.grid(lambda))
    #LOW_UPP <- as.matrix(expand.grid(low_upp))
    #z_RANGE <- as.matrix(expand.grid(z_range))
    results[LOW_UPP] <- results[LOW_UPP] + apply(LAMBDA, 1, prod)/(n * apply(z_RANGE, 1, prod))
  }
  return(results)
}
