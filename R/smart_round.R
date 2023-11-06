#' @rdname smart_round
#' @title Smart rounding of non integer vector
#'
#' @description Perform rounding of a numeric vector which conservs the total sum 
#'
#' @param x a numeric vector
#' @return A rounded vector which sum is the same as the initial vector  
#' @export

smart_round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}