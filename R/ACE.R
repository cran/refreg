#' Alternating Conditional Expectation (ACE) algorithm
#'
#' This function implements Alternating Conditional Expectation algorithm
#' (ACE, Hastie and Tibshirani, 1990). This is a bivRegr's inner function
#' for estimating variance, and correlation models.
#'
#' @param y A character defining the response variable.
#' @param predictor The regression predictor, following the mgcv package structure (see example below).
#' @param restriction Type of restriction to be imposed to the response variable,
#' "positive" for variance, and "correlation" for the correlation.
#' @param eps A number defining the allowed estimation error, default = 0.01.
#' @param itmax Maximum number of iterations of the algorithm, default = 10.
#' @param data A data frame containing the response, and predictor variables.
#' @param ... Additional mgcv::gam() parameters to be modified by the user.
#' @return This function returns a mgcv::gam() fit for a transformed response.
#' @references Hastie, T. & Tibshirani, J. (1990) Generalized additive models.
#' CRC press. London.
#' @examples
#' n <- 100
#' x <- runif(n, -1, 1)
#' y <- x^2 + rnorm(n, sd = 0.1)
#' df <- data.frame(y, x)
#' plot(df$x, df$y)
#' m1 <- ACE(
#'   y = "y", predictor = "~s(x)", restriction = "positive",
#'   eps = 0.01, itmax = 10, data = df
#' )$fit
#' nw <- data.frame(x = seq(-1, 1, 0.1))
#' abline(h = 0, col = 2)
#' lines(exp(predict(m1, newdata = nw)) ~ nw$x, col = 3, lwd = 2)
#' legend("top", legend = c("ACE fit", "Zero"), lty = 1, lwd = 2, col = c(3, 2), bty = "n")
#' @importFrom "mgcv" "gam" "predict.gam"
#' @export

ACE <- function(y = "y", predictor = "~s(x)", restriction = "positive", eps = 0.01, itmax = 10, data = data, ...) {
  y <- data[, y]

  if (restriction == "positive") {
    H <- function(x) exp(x)
    H_ <- function(x) log(x)
    H1 <- function(x) exp(x)
  }

  if (restriction == "correlation") {
    H <- function(x) {
      x[x > 3] <- 3
      x[x < (-3)] <- -3
      tanh(x)
    }

    H_ <- function(x) atanh(x)

    H1 <- function(x) {
      x[x > 3] <- 3
      x[x < (-3)] <- -3
      1 - tanh(x)**2
    }
  }

  if (restriction == "correlation" & mean(y) >= 1) {
    eta0 <- H_(0.99)
  } else {
    if (restriction == "correlation" & mean(y) <= -1) {
      eta0 <- H_(-0.99)
    } else {
      eta0 <- H_(mean(y))
    }
  }

  ## Estimation loop
  mu0 <- H(eta0)
  iter <- 0
  continuar <- TRUE
  while (continuar) {
    iter <- iter + 1
    res <- log((y - mu0)**2)
    var <- exp(predict(gam(as.formula(paste0("res", parse(text = predictor))), data = data)))
    der <- H1(eta0)
    var[which(var < 0.01)] <- 0.01
    w <- der**2 / var
    z <- eta0 + (y - mu0) / der
    fit <- gam(as.formula(paste0("z", parse(text = predictor))), data = data, weights = w, ...)
    eta <- predict(fit)
    mu <- H(eta)
    MSE0 <- mean(y - mu0)**2
    MSE <- mean(y - mu)**2
    error <- abs(MSE - MSE0) / MSE0
    if (error < eps | iter > itmax) continuar <- FALSE
    eta0 <- eta
    mu0 <- mu
  }
  return(list(
    fit = fit, error = error
  ))
}
