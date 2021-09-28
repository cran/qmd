#dependence measures

.zeta1 <- function(X, Y, resolution = NULL, ties.correction = FALSE){

  N <- resolution

  #dimension of conditioning
  d <- NCOL(X)
  #dimension of vector
  rho <- d+1
  #sample size
  n <- NROW(X)
  #resolution
  if(is.null(N)){
    N <- max(c(floor(n^(1/rho)),2))
  }else{
    N <- N
  }


  #Compute checkerboard weights
  M <- pure_r_ECBC(X,Y, N)

  #Compute marginal checkerboard weights
  M_margin <- apply(M, MARGIN = 1:d, sum)

  #Compute the integral and hence zeta1
  M_new <- matrix(M, ncol = N, byrow = F)
  index <- which(M_margin > 0)
  k_sum <- 0
  for(i in 1:length(index)){
    k1 <- M_new[index[i],]/M_margin[index[i]]
    k_indep <- rep(1/N,N)
    k_int <- .r_local_kernel_integral(k1,k_indep)
    k_sum <- k_sum + k_int*M_margin[index[i]]
  }

  if(ties.correction){
    return(3*k_sum/qad::zeta1(Y,Y,resolution = n))
  }else{
    return(k_sum*3)
  }
}




.zeta1_approx <- function(X,Y,N=NULL){

  X <- as.matrix(X)
  #dimension of conditioning
  d <- NCOL(X)
  #dimension of vector
  rho <- d+1
  #sample size
  n <- NROW(X)
  #resolution
  if(is.null(N)){
    N <- max(c(floor(n^(1/rho)),2))
  }else{
    N <- N
  }

  U <- apply(X, 2, function(x)ecdf(x)(x))
  V <- ecdf(Y)(Y)
  W <- cbind(U,V)
  b <- seq(0,1,1/N)

  test <- apply(W, 2, function(x) cut(x,b))
  l <- list()
  for(i in 1:NCOL(test)){
    l[[i]] <- factor(test[,i], levels = levels(cut(W[,1],b)))
  }

  M <- table(l)
  M <- as.array(M)/n



  #Compute marginal checkerboard weights
  M_margin <- apply(M, MARGIN = 1:d, sum)

  #Compute the integral and hence zeta1
  M_new <- matrix(M, ncol = N, byrow = F)
  index <- which(M_margin > 0)
  k_sum <- 0
  for(i in 1:length(index)){
    k1 <- M_new[index[i],]/M_margin[index[i]]
    k_indep <- rep(1/N,N)
    k_int <- .r_local_kernel_integral(k1,k_indep)
    k_sum <- k_sum + k_int*M_margin[index[i]]
  }
  return(unname(k_sum*3))
}




.zeta1_pvalue <- function(X, Y, ties.correction = FALSE, approx = FALSE, resolution = NULL,p.value = TRUE, R = 1000){
  if(approx){
    res <- .zeta1_approx(X,Y,N = resolution)
    if(p.value){
      res_R <- rep(0,R)
      for(i in 1:R){
        y <- runif(length(Y))
        res_R[i] <- .zeta1_approx(X,y,N = resolution)
      }
      pvalue <- 1 - ecdf(res_R)(res)
      return(c(zeta1 = res, p.value = pvalue))
    }else{
      return(c(zeta1 = res))
    }
  }else{
    res <- .zeta1(X, Y, resolution = resolution, ties.correction = ties.correction)
    if(p.value){
      res_R <- rep(0,R)
      for(i in 1:R){
        y <- runif(length(Y))
        res_R[i] <- .zeta1(X,y, resolution = resolution, ties.correction = ties.correction)
      }
      pvalue <- 1 - ecdf(res_R)(res)
      return(c(zeta1 = res, p.value = pvalue))
    }else{
      return(c(zeta1 = res))
    }
  }
}
