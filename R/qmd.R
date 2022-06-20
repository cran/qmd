#' Quantification of Multivariate Dependence
#'
#' @param X a numeric matrix or data.frame of dimension d containing the explanatory variables
#' @param y a numeric vector containing the uni-variate response variable
#' @param p.value logical indicating if a p-value is returned using permutations of Y
#' @param R integer indicating the number of repetitions for the calculation of the p-value (default = 1000)
#' @param ties.correction logical indicating if the measure of dependence should be calculated with ties-correction (experimental version). Default = FALSE.
#' @param resolution an integer indicating the resolution N of the checkerboard aggregation. We recommend to use the default configuration (resolution = NULL), which uses the resolution N(n) = floor(n^(1/(d+1))), where d denotes the number of explanatory variables.
#' @param print logical indicating whether the results of the function are printed
#' @param na.exclude logical if all rows containing NAs should be removed.
#'
#' @description Function for estimating the non-parametric copula-based multivariate measure of dependence \eqn{\zeta}1.
#' This measure quantifies the extent of dependence between a d-dimensional random vector X and a uni-variate random variable y (i.e.,
#' it measures the influence of d explanatory variables X1,...,Xd on a univariate variable y). Further details can be found in the section Details and the corresponding references.
#'
#'
#' @return qmd returns a list object containing the following components:
#' \itemize{
#'     \item input: data containing the explanatory variables (X)
#'     \item output: data containing the response (y)
#'     \item q(X,y): dependence measure indicating the extent of dependence between X and y
#'     \item results: data.frame containing the dependence measure and the corresponding p-value
#'     \item resolution: an integer indicating the resolution of the aggregated checkerboard copula
#'     \item Sample size
#' }
#'
#'
#' @details In the following we will simply write q for the dependence measure \eqn{\zeta}1. Furthermore, X denotes a random vector consisting of d random variables and y denotes a univariate random variable.
#' Then the theoretical dependence measure q fulfills the following essential properties of a dependence measure:
#' \itemize{
#'    \item[N]  q(X,y) attains values in [0,1] (normalization).
#'    \item[I]  q(X,y) = 0 if and only if X and y are independent (independence).
#'    \item[C]  q(X,y) = 1 if and only if y is a function of X (complete dependence).
#' }
#' Further properties of q and the exact mathematical definition can be found in Griessenberger et al. (2022). This function qmd() contains
#' the empirical checkerboard-estimator (ECB-estimator), which is strongly consistent and attains always positive values between 0 and 1.
#' Note, that interpretation of low values has to be done with care and always under consideration of the sample size. For instance, values of 0.2 can point towards independence in small sample settings.
#' An additional p-value (testing for independence and being based on permutations of y) helps in order to correctly understand the dependence values.
#' Since independence constitutes the null hypothesis a p-value above the significance level (e.g., 0.05) indicates independence between X and y.
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
#' qmd(X = cbind(x1,x2,x3), y = y, p.value = TRUE)
#'
#' #(independence for dimension 4)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- x1 + x2 + rnorm(n)
#' y <- runif(n)
#' qmd(X = cbind(x1,x2,x3), y = y, p.value = TRUE)
#'
#' #(binary output (classification) for dimension 3)
#' n <- 500
#' x1 <- runif(n)
#' x2 <- runif(n)
#' y <- ifelse(x1 + x2 < 1, 0, 1)
#' qmd(X = cbind(x1,x2), y = y, p.value = TRUE)
#' #(independence)
#' y <- runif(n)
#' qmd(X = cbind(x1,x2), y = y, p.value = TRUE)



qmd <- function(X, y,
                ties.correction = FALSE,
                resolution = NULL,
                p.value = FALSE,
                R = 1000,
                print = TRUE,
                na.exclude = FALSE){

  X <- as.matrix(X)

  if(!all(apply(X, 2, is.numeric))){
    stop("Not every input variable of X is numeric!")
  }
  if(!is.numeric(y)){
    stop("The response y is not numeric!")
  }
  if(NROW(X) != NROW(y)){
    stop("Number of observations differ!")
  }
  if(!na.exclude){
    if(any(apply(cbind(X,y), 2, is.na))){
      stop("Remove NAs or use na.exclude = TRUE")
    }
  }

  Z <- cbind(X,y)
  Z <- na.omit(Z)
  X_new <- Z[,1:NCOL(X)]
  Y_new <- Z[,NCOL(X) + 1]

  #The next lines contain the calculation of zeta1

  res <- .zeta1(X_new, Y_new, resolution = resolution, ties.correction = ties.correction)
  pvalue <- as.numeric(NA)
  if(p.value){
    res_R <- rep(0,R)
    for(i in 1:R){
      y_perm <- sample(Y_new, length(Y_new))
      res_R[i] <- .zeta1(X_new,y_perm, resolution = resolution, ties.correction = ties.correction)
    }
    #pvalue <- 1 - ecdf(res_R)(res)
    pvalue <- mean(ifelse(res_R >= res, 1, 0))
  }
  #---end calculation of zeta1

  #Get names of variables
  if(is.null(names(y))){
    name_y <- deparse(substitute(y))
  }else{
    name_y <- names(y)
  }
  if(is.null(names(X))){
    names_X <- paste0(names(apply(X, 2, function(x) deparse(substitute(x)))), collapse = ", ")
  }else{
    names_X <- paste0(names(X), collapse = ", ")
  }

  #Get resolution
  d <- NCOL(X)
  #dimension of vector
  rho <- d+1
  if(is.null(resolution)){
    reso <- floor(NROW(X_new)^(1/rho))
  }else{
    reso <- resolution
  }

  #output
  names <- c('q(X,y)')
  q.values <- c(res)
  p.values <- c(pvalue)
  output <- data.frame(names, q.values, p.values)
  names(output) <- c('','coef','p.value')

  output_q <- output[1,]
  names(output_q) <- c('','qmd','p.value')
  output_q[,2:3] <- round(output_q[,2:3],3)

  if(print){
    cat("\n")
    cat("Quantification of Multivariate Dependence (qmd):", "\n")
    cat("\nNumber of input variables (random vector X):", NCOL(X))
    cat("\nVariable names: X :=", names_X)
    cat("\n                y :=", name_y)
    cat("\n")
    cat(paste("\nOriginal sample size:"),NROW(X))
    cat(paste("\nNumber of rows (with NAs) removed:", NROW(X)-NROW(X_new)))
    cat(paste("\nResolution of checkerboard:",reso))
    cat("\n\nDependence measure:")
    cat("\n")
    if(is.na(pvalue)){
      print.data.frame(format(output_q[,1:2], justify='left', digits=3), row.names = FALSE)
    }else{
      print.data.frame(format(output_q, justify='left', digits=3), row.names = FALSE)
    }

  output <- list(input = X,
                     output = y,
                     `q(X,y)` = res,
                     results = output,
                     resolution = reso,
                     n = NROW(X))
  invisible(output)
  }

}
