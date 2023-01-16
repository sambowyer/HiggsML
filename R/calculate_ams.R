#' Calculate Approximate Median Significance (AMS)
#'
#' This formula is from: \url{http://opendata.cern.ch/record/328}.
#'
#' @param y The true class labels (signal = 0, background = 1).
#' @param yhat The predicted class labels (signal = 0, background = 1).
#' @param w Weights of the data points in question (as provided in the original dataset).
#' @param original_sum_b Sum of background sample weights. Default is 125317.9 to match the Kaggle training set.
#' @param original_sum_s Sum of signal sample weights. Default is 212.418 to match the Kaggle training set.
#' @param b_r Regularisation term. Default is 10 as per the HiggsML Challenge document.
#'
#' @return The AMS of the predicted class labels given the true values.
#' @export
#'
#' @examples
#' y = rep(0, 5)
#' w = rep(1, 5)
#' calculate_ams(y, y, w)
#'
#' @tests
#' y = rep(0, 5)
#' w = rep(1, 5)
#' expect_equal(calculate_ams(y, y, w), 30.90362, tolerance=1e-3)
#'
#' yhat = rep(1, 5)
#' expect_equal(calculate_ams_unormalized(y, yhat, w), 0)
calculate_ams <- function(y, yhat, w, original_sum_b = 125317.9, original_sum_s = 212.418, b_r = 10) {
  w_new <- vector(length = length(w))
  sum_b <- original_sum_b/sum(w[y==1])
  sum_s <- original_sum_s/sum(w[y==0])
  w_new[y ==1 ] <- w[y==1]*sum_b
  w_new[y ==0 ] <- w[y==0]*sum_s
  s <- sum(w_new * (1-y) * (1-yhat))
  b <- sum(w_new * y * (1-yhat))
  ams <- sqrt(2 * ((s + b + b_r) * log(1 + s/(b + b_r)) - s))
  ams
}

#' Calculate Approximate Median Significance (AMS) Without Sample Weights
#'
#' This formula is from: \url{http://opendata.cern.ch/record/328}.
#'
#' @param y The true class labels (signal = 0, background = 1).
#' @param yhat The predicted class labels (signal = 0, background = 1).
#' @param w Weights of the data points in question (as provided in the original dataset).
#' @param b_r Regularisation term. Default is 10 as per the HiggsML Challenge document.
#'
#' @return The AMS of the predicted class labels given the true values.
#' @export
#'
#' @examples
#' y = rep(0, 5)
#' w = rep(1, 5)
#' calculate_ams(y, y, w)
#'
#' @tests
#' y = rep(0, 5)
#' w = rep(1, 5)
#' expect_equal(calculate_ams_unormalized(y, y, w), 0)
#'
#' yhat = rep(1, 5)
#' expect_equal(calculate_ams_unormalized(y, yhat, w), 0)
calculate_ams_unormalized <- function(y, yhat, w, b_r = 10) {
  s <- sum(w * y * yhat)
  b <- sum(w * (1-y) * yhat)
  ams <- sqrt(2 * ((s + b + b_r) * log(1 + s/(b + b_r)) - s))
  ams
}
