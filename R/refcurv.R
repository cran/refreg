#' Univariate reference curve model
#'
#' This function obtain univariate conditional quantiles as described
#' in Martinez-Silva et. al (2016).
#'
#' @param mu A formula object for the response mean model following the mgcv package structure (see example below).
#' @param sigma a formula object for fitting a model to the response variance (see example below).
#' @param data A data frame containing both the response, and predictor variables.
#' @return This function returns univariate conditional quantiles estimated using a non parametric location scale model.
#' @details In the Martinez Silva et. al (2016) the non linear effects of the continuous
#' covariates are estimating through polynomial kernel smoother, in this package we implement the
#' same methodology but using penalized splines in order to reduce computational cost.
#' @references Martinez--Silva, I., Roca--Pardinas, J., & Ordonez, C. (2016). Forecasting SO2 pollution incidents by means of quantile curves based on additive models. Environmetrics, 27(3), 147--157.
#' @export
#' @importFrom "stats" "get_all_vars"
#' @examples
#' #--- Glycation hemoglobin reference curve depending on age
#' dm_no <- subset(aegis, aegis$dm == "no")
#' fit1 <- refcurv(mu = "hba1c~s(age)", sigma = "~s(age)", data = dm_no)
#' plot(fit1, newdata = data.frame(age = 18:90), tau = c(0.025, 0.05, 0.10, 0.90, 0.95, 0.975))
#'
refcurv <- function(mu = "y~s(x)", sigma = "~s(x)", data = data) {
  modelo_m <- gam(as.formula(mu), data = data, method = "REML")
  data$res <- residuals(modelo_m)^2

  modelo_v <- ACE(y = "res", predictor = sigma, restriction = "positive", eps = 0.01, itmax = 10, data = data)$fit
  mean <- predict(modelo_m)
  sd <- sqrt(exp(predict(modelo_v)))
  data$res <- (data[, all.vars(modelo_m$formula)[1]] - mean) / sd
  ret <- list(modelo_m = modelo_m, modelo_v = modelo_v, res = data$res)
  class(ret) <- "refcurve"
  return(ret)
}
