#' Variable selection using the qmd-dependence values
#'
#' @param X a numeric matrix or data.frame of dimension d containing the explanatory variables
#' @param y a numeric vector containing the uni-variate response variable
#' @param method possible options are c("combVar", "addVar"), see Details.
#' @param bin.size either "fixed", "adaptive" or "sparse.adaptive", indicating whether the checkerboard copula may vary its bin sizes (defaults to "fixed"). Setting this to "adaptive" might affect the results but will be faster if the sample has many ties.
#' @param plot logical indicating whether the feature selection plot is printed
#' @param na.exclude logical if all rows containing NAs should be removed.
#' @param max_num_features maximal number of explanatory variables to be selected
#' @param plot.title a label for the title
#' @param plot.color a colour for the selected variables
#'
#' @description Given a d-dimensional random vector X containing the explanatory variables and a uni-variate response variable y,
#' this function uses the qmd-dependence values to select the most relevant (influential) explanatory variables.
#' Two different methods are available and are explained in the section Details.
#'
#' @details
#' method 1 (default) - "combVar": computes all qmd-dependence scores, i.e., calculates the dependence of every combination of explanatory variables to the response variable y and selects for each number of explanatory variables the combination with the greatest dependence score. This procedure is computational expensive and is only available up to 15 explanatory variables.
#'
#' method 2 - "addVar": stepwise procedure which calculates all bi-variate dependence values q(X_i,Y) and selects the variable X_j exhibiting the greatest dependence value. In the next step all three-dimensional combinations q((X_j, X_i), Y) (for every i =1,..., d and i not j) are computed and the variable exhibiting again the greatest dependence score is added. In this manner the procedure works up to dimension d.
#'
#'
#' @return a list containing a data.frame (result) and the corresponding plots.
#' The data.frame result contains the number of explanatory variables (nVars),
#' the combination of selected variables (selVars),
#' the dependence measure zeta1 (qmd) of the selected variables to the response y
#' and the resolution of the empirical checkerboard copula (ECBC_resolution).
#' For the method "combVar" the dependence value zeta1 (qmd) is returned for all combinations of explanatory variables and is sorted in decreasing order according to zeta1.
#'
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n)
#' x2 <- rexp(n)
#' x3 <- x1 + log(x2) + rnorm(n)
#' x4 <- rnorm(n)
#' x5 <- x4^2
#' x6 <- x1 + x5 + rnorm(n)
#' x7 <- 1:n
#' y <- x2 + x4*x7 + runif(n)
#' X <- data.frame(x1,x2,x3,x4,x5,x6,x7)
#' fit <- feature_selection(X, y, method = "combVar", plot = TRUE)
#' fit <- feature_selection(X, y, method = "addVar", plot = TRUE)


feature_selection <- function(X, y,
                              method = "combVar",
                              bin.size = "fixed",
                              plot = TRUE,
                              na.exclude = FALSE,
                              max_num_features = NULL,
                              plot.title = NULL,
                              plot.color = "hotpink"){
  .feature_selection(X, y,
                     method = method,
                     bin.size = bin.size,
                     plot = plot,
                     na.exclude = na.exclude,
                     max_Var = max_num_features,
                     plot.title = plot.title,
                     plot.color = plot.color)
}


.feature_selection <- function(X, y,
                               method = "addVar",
                               bin.size = "fixed",
                               plot = TRUE,
                               na.exclude = FALSE,
                               max_Var = NULL,
                               plot.title = NULL,
                               plot.color = "hotpink"){

  X <- as.matrix(X)
  if(!all(apply(X, 2, is.numeric))){
    stop("Not every input variable of X is numeric!")
  }
  if(!is.numeric(y)){
    stop("The response y is not numeric!")
  }
  if(NROW(X) != NROW(y)){
    stop("Number of observations between explanatory variables and response variable differ!")
  }
  if(!na.exclude){
    if(any(apply(cbind(X,y), 2, is.na))){
      stop("Remove NAs or use na.exclude = TRUE")
    }
  }

  Z <- cbind(X,y)
  Z <- apply(na.omit(Z), 2, qmdrank)
  X_new <- Z[,1:NCOL(X)]
  Y_new <- Z[,NCOL(X) + 1]


  if(method == "addVar"){
    d <- NCOL(X_new)
    k <- data.frame(step0 = 1:d)
    add <- add_values <- res <-  rep(0,d)

    max_useful_d <- floor(NROW(X_new)^(1/(1:d+1)))
    max_useful_d <- max(which(max_useful_d > 1))
    if(!is.null(max_Var)){
      max_useful_d <- min(max_useful_d, max_Var)
    }
    for(i in 1:max_useful_d){

      #resolution
      N <- floor(NROW(X_new)^(1/(i + 1)))
      res[i] <- N
      q_values <- apply(k, 1, function(x) .zeta1_ranks(U = X_new[,x], V = Y_new, resolution = N, bin.size = bin.size))
      if(i == 1){
        add[i] <- which.max(q_values)
      }else{
        add[i] <- as.numeric(names(which.max(q_values)))
      }
      add_values[i] <- max(q_values)

      k[[paste0("step",i)]] <- rep(add[i], d-i+1)
      k <- k[-which(k$step0 == add[i]),]

    }
    #sapply(add, function(x) paste(add[1:which(add == x)], collapse = ","))
    add <- add[1:max_useful_d]
    add_values <- add_values[1:max_useful_d]
    res <- res[1:max_useful_d]

    selVars <- rep(NA, max_useful_d)
    for(i in 1:max_useful_d){
      selVars[i] <- paste0(colnames(X)[add[1:i]], collapse = ",")
    }

    result <- data.frame(nVars = 1:length(add),
                         selVars = selVars,
                         zeta1 = round(add_values,4),
                         #zeta1Delta = round(c(add_values[1],add_values[-1] - add_values[-length(add_values)]), 4),
                         ECBC_resolution = res)
    row.names(result) <- NULL

    if(max_useful_d < d & is.null(max_Var)){
      message(paste0("Note: The feature selection algorithm only considers maximal ", max_useful_d, " variables! Otherwise the checkerboard resolution would be 1 (due to the small sample size or the high dimension)!"))
    }
#____________________________________________________________________________________________________________________________________
  }else if(method == "biVar"){
    # method 3 - "biVar": computes all bivariate dependence values (i.e., q(X_1,Y), q(X_2,Y),...) and returns the values in decreasing order. Note, that this procedure ignores multivariate dependence structures.

    #here: all bivariate zeta 1 values are calculated
    d <- NCOL(X_new)
    k <- data.frame(step0 = 1:d)
    add <- add_values <- rep(0,d)

    if(is.null(max_Var)){
      max_Var <- d
    }

    #resolution
    N <- floor(NROW(X_new)^(1/(1 + 1)))
    res <- rep(N,d)
    q_values <- apply(k, 1, function(x) .zeta1_ranks(U = X_new[,x], V = Y_new, resolution = N, bin.size = bin.size))
    names(q_values) <- colnames(X)
    q_values <- sort(q_values, decreasing = TRUE)

    selVars <- rep(NA, max_Var)
    for(i in 1:max_Var){
      selVars[i] <- paste0(names(q_values)[1:i], collapse = ",")
    }

    result <- data.frame(nVars = 1:max_Var,
                         selVars = selVars,
                         zeta1 = round(q_values[1:max_Var],4),
                         ECBC_resolution = res[1:max_Var])
    row.names(result) <- NULL
    #____________________________________________________________________________________________________________________________________

  }else if(method == "combVar"){
    d <- NCOL(X)

    if(d > 16){
      stop("Method 'combVar' only works for at most 16 explanatory variables. Please reduce the number of explanatory variables or use another method. Alternative methods are method = c('addVar', 'biVar')")
    }

    Z <- cbind(X,y)
    Z <- apply(na.omit(Z), 2, qmdrank)
    X_new <- Z[,1:NCOL(X)]
    Y_new <- Z[,NCOL(X) + 1]

    sample_size <- NROW(Z)

    dimension_limit <- min(d, floor(log2(sample_size)) - 1)
    if(!is.null(max_Var)){
      dimension_limit <- min(dimension_limit, max_Var)
    }
    if (dimension_limit < d && is.null(max_Var)) {
        message(paste0("Note: The feature selection algorithm only considers maximal ", dimension_limit, " variables! Otherwise the checkerboard resolution would be 1 (due to the small sample size or the high dimension)!"))
    }

    if (!is.null(max_Var)){
      if (max_Var > d) stop("too many dimensions selected")
      if (max_Var > dimension_limit) warning("the automatically calculated resolution of the checkerboard approximation is 1, this is caused by having a too many dimensions for this sample size")
    }

    if(is.null(colnames(X))){
      names_X <- colnames(apply(X, 2, function(x) deparse(substitute(x))))
    }else{
      names_X <- colnames(X)
    }

    #Get all combinations
    df <- data.frame()
    combn_sizes <- 1:dimension_limit

    if(!is.null(max_Var))
      combn_sizes <- 1:max_Var

    for(i in combn_sizes){
      combinations <- combn(1:NCOL(X_new), i)
      for(j in 1:NCOL(combinations)){
        N <- floor(NROW(X_new)^(1/(NROW(combinations[,j]) + 1)))
        df <- rbind(df, data.frame(nVars = i,selVars = paste(names_X[combinations[,j]], collapse = ","), zeta1 = .zeta1_ranks(U = X_new[,combinations[,j]], V = Y_new, resolution = N, bin.size = bin.size), ECBC_resolution = N))
      }
      if(is.null(max_Var)) print(paste0(i,"/", dimension_limit))
    }

    result <- df[order(-df$zeta1),]
    result$zeta1 <- round(result$zeta1, 4)
    row.names(result) <- NULL
  }else{
    stop("No method selected!")
  }

  #maximal numbers of variables to plot
  max_vars <- 16

  df1 <- data.frame(dplyr::summarise(dplyr::group_by(result, across("nVars")), zeta1 = max(.data$zeta1), ECBC_resolution = max(.data$ECBC_resolution)))
  df1 <- df1[df1$nVars <= max_vars,]
  df1$resolution_trans <- df1$ECBC_resolution/max(df1$ECBC_resolution)

  index <- which(!duplicated(dplyr::select(result, all_of("nVars"))))
  df2 <- result[index,]

  df_new <- data.frame()
  for(i in 1:NROW(df2)){
    df_new <- rbind(df_new, data.frame(nVars = df2[i,1], selVars = unlist(strsplit(df2[i, 2], ","))))
  }
  df_new$fill <- "1"

  tile <- expand.grid(nVars = 1:max(df1$nVars), selVars = sort(df_new$selVars))
  tile <- merge(tile, df_new, by = c("nVars", "selVars"), all.x = TRUE)



  p <- ggplot(df1, aes_string(x = "nVars", y = "zeta1"))
  p <- p + geom_point(size = 2, shape = 15)
  p <- p + geom_line()
  p <- p + geom_label(aes_string(label = "zeta1"), vjust = 1.3)
  p <- p + geom_point(aes_string(x = "nVars", y = "resolution_trans"), color = "blue")
  p <- p + geom_line(aes_string(x = "nVars", y = "resolution_trans"), color = "blue", linetype = "dashed")
  p <- p + theme_bw()
  p <- p + scale_x_continuous(limits = c(0.5, max(df1$nVars)+0.5), breaks = 1:max(df1$nVars))
  p <- p + scale_y_continuous(limits = c(0,1), breaks = function(x) pretty(x, 6),
                              sec.axis = sec_axis(trans=~.*(max(df1$ECBC_resolution)), name="Checkerboard resolution N",
                                                  breaks = function(x) pretty(x, 10)))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.text.x = element_blank(),
                 panel.grid = element_blank())
  p <- p + labs(caption = "Black line depicts the q-values. \nBlue dashed line depicts the resolution N (second axis).")
  if(method == "biVar"){
    p <- p + ylab("bivariate q-value of added variable")
  }else{
    p <- p + ylab("q-value")
  }
  if(is.null(plot.title)){
    p <- p + ggtitle(paste0("Feature selection qmd (method='", method, "')"))
  }else{
    p <- p + ggtitle(plot.title)
  }
  p1 <- p

  p <- ggplot(df_new, aes_string(x = "nVars", y = "selVars"))
  p <- p + geom_tile(data = tile, aes_string(x = "nVars", y = "selVars", fill = "fill"), color = "black")
  p <- p + scale_fill_manual(values = c("1" = plot.color), na.value = "transparent")
  p <- p + theme_bw() + xlab("Number of explanatory variables") + ylab("Variables")
  p <- p + theme(panel.grid = element_blank())
  p <- p + scale_x_continuous(breaks =  df1$nVars)
  p <- p + guides(fill = "none")
  p2 <- p

  if(plot){
    print(cowplot::plot_grid(p1,p2, ncol = 1, align = "v", rel_heights = c(0.7,0.3)))
  }

  print(utils::head(result, n = 10))

  return(invisible(list(result = result,
              p1 = p1,
              p2 = p2,
              plot = cowplot::plot_grid(p1,p2, ncol = 1, align = "v", rel_heights = c(0.7,0.3)))))

}



