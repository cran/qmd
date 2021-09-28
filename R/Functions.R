#Expand grid of matrices (fast version of original expand.grid)
.expand.grid.new <- function(list_seq){
  rep.fac <- 1
  res <- c()
  orep <- prod(sapply(list_seq, length))

  for(i in 1:length(list_seq)){
    nx <- length(list_seq[[i]])
    x <- list_seq[[i]]
    orep <- orep/nx
    x <- x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)]
    rep.fac <- rep.fac * nx
    res <- cbind(res,x)
  }
  return(res)
}


# Compute difference of integrals exact

.r_local_kernel_integral <- function(k1, k2){
  #k1...checkerboard weights of first kernel (local)
  #k2...checkerboard weights of second kernel (local)

  #get resolution
  N <- length(k1)

  #append 0
  k1 <- cumsum(c(0,k1))
  k2 <- cumsum(c(0,k2))

  #Compute the integral for each segment
  ysum <- 0
  for(i in 1:N){
    #case 1: there is no intersection
    if((k1[i+1]-k2[i+1])*(k1[i]-k2[i]) >= 0){
      I1 <- (k1[i+1] - k1[i])/(2 * N) + k1[i]/N
      I2 <- (k2[i+1] - k2[i])/(2 * N) + k2[i]/N
      I <- abs(I1 - I2)
    }else{
      x_intersect <- (k1[i] - k2[i])/(((k2[i+1] - k2[i]) - (k1[i+1] - k1[i])) * N)
      y_intersect <- k1[i] + N * (k1[i+1] - k1[i]) * x_intersect

      I11 <- ((y_intersect - k1[i]) * x_intersect)/2 + k1[i] * x_intersect
      I12 <- ((y_intersect - k2[i]) * x_intersect)/2 + k2[i] * x_intersect
      I21 <- ((k1[i+1] - y_intersect) * (1/N - x_intersect))/2 + y_intersect * (1/N - x_intersect)
      I22 <- ((k2[i+1] - y_intersect) * (1/N - x_intersect))/2 + y_intersect * (1/N - x_intersect)
      I <- abs(I11 - I12) + abs(I21 - I22)
    }
    ysum <- ysum + I
  }
  return(ysum)
}


#get number of ranges
.pure_r_range = function(X) {
  r = rep(0, length(X))
  for (x in X) {
    r[x] = r[x] + 1
  }
  return(r)
}


