#' Equivalent to rank(x, ties.method = "max") but not as stupidly slow
#'
#' @param x A numeric vector
#' @return An integer vector specifying for each value in x the rank
#'         within x. If one value appears multiple time the maximum
#'         is used.
qmdrank <- function(x) {
  order <- order(x)
  res <- rep(0, length(x))
  x <- x[order]

  val <- seq_until_changes(x)

  res[order] <- val

  res
}
