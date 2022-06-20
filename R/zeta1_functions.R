#dependence measures

.zeta1 <- function(X, Y, resolution = NULL, ties.correction = FALSE, bin.size = "fixed"){

  #dimension of conditioning
  d <- NCOL(X)
  #dimension of vector
  rho <- d+1

  W <- apply(cbind(X, Y), 2, qmdrank)
  U <- W[,1:(rho-1)]
  V <- W[,rho]

  return(.zeta1_ranks(U, V, resolution, ties.correction, bin.size))
}

.zeta1_ranks <- function(U, V, resolution = NULL, ties.correction = FALSE, bin.size = "fixed"){

  N <- resolution

  #dimension of conditioning
  d <- NCOL(U)
  #dimension of vector
  rho <- d+1
  #sample size
  n <- NROW(U)
  #resolution
  if(is.null(N)){
    N <- floor(n^(1/rho))
    if (N == 1) {
      warning("the automatically calculated resolution of the checkerboard approximation is 1; caused by to many dimensions for this sample size")
    }
  }else{
    N <- N
  }

  #Compute checkerboard weights
  if (bin.size == "adaptive")
    M <- .EACBC(cbind(U, V), N)
  else if (bin.size == "sparse.adaptive")
    M <- .EACBC_nonzero(cbind(U, V), N)
  else M <- .ECBC(cbind(U, V), N)

  mY = rep(1/N,N)
  if (bin.size == "adaptive" || bin.size == "sparse.adaptive") {
    #Compute Y-masses
    mY <- .adaptive_masses(V, N)
  }

  #Compute marginal checkerboard weights
  if (bin.size == "sparse.adaptive")
    M_margin <- lapply(M, sum)
  else
    M_margin <- apply(M, MARGIN = 1:d, sum)

  #Compute the integral and hence zeta1
  if (bin.size != "sparse.adaptive") {
    M_new <- matrix(M, ncol = N, byrow = F)
    index <- which(M_margin > 0)
    k_sum <- 0
    for(i in 1:length(index)){
        k1 <- M_new[index[i],]/M_margin[index[i]]
        k_indep <- mY                                   # USE THE MASSES OF THE TRANSFORMED Y-BLOCKS FOR PI
        k_int <- .local_kernel_integral(k1,k_indep,mY)
        k_sum <- k_sum + k_int*M_margin[index[i]]
    }
  } else {
    k_sum <- 0
    for(i in 1:length(M)){
        k1 <- M[[i]]/M_margin[[i]]
        k_indep <- mY                                   # USE THE MASSES OF THE TRANSFORMED Y-BLOCKS FOR PI
        k_int <- .local_kernel_integral(k1,k_indep,mY)
        k_sum <- k_sum + k_int*M_margin[[i]]
    }
  }

  if(ties.correction){
    return(3*k_sum/qad::zeta1(V,V,resolution = n))
  }else{
    return(k_sum*3)
  }
}

.zeta1CB <- function(CB){

  D <- dim(CB)
  N <- D[1]

  #dimension of vector
  rho <- length(D)
  d <- rho-1

  M=CB

  #Compute marginal checkerboard weights
  M_margin <- apply(M, MARGIN = 1:d, sum)

  #Compute the integral and hence zeta1
  M_new <- matrix(M, ncol = N, byrow = F)
  index <- which(M_margin > 0)
  k_sum <- 0
  for(i in 1:length(index)){
    k1 <- M_new[index[i],]/M_margin[index[i]]
    k_indep <- rep(1/N,N)
    k_int <- .local_kernel_integral(k1,k_indep,k_indep)
    k_sum <- k_sum + k_int*M_margin[index[i]]
  }


  return(k_sum*3)

}

.zeta1_pvalue <- function(X, Y, ties.correction = FALSE, resolution = NULL, p.value = TRUE, R = 1000){
  rho <- NCOL(X) + 1
  W <- apply(cbind(X, Y), 2, qmdrank)
  U <- W[,1:(rho-1)]
  V <- W[,rho]
  res <- .zeta1_ranks(U, V, resolution = resolution, ties.correction = ties.correction)
  if(p.value){
    res_R <- rep(0,R)
    for(i in 1:R){
      y <- sample(length(Y), length(Y))
      res_R[i] <- .zeta1_ranks(U, y, resolution = resolution, ties.correction = ties.correction)
    }
    pvalue <- 1 - mean(res_R <= res)
    return(c(zeta1 = res, p.value = pvalue))
  }else{
    return(c(zeta1 = res))
  }
}
