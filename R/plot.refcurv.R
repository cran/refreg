#' Plot univariate conditional quantile models curves (i.e. reference curves)
#'
#' This function depict the univariate conditional quantile model based
#' on the non parametric location scale model fitted with the refcurv function.
#'
#' @param x A refcurv object.
#' @param newdata A data frame defining a sequence of the predictor variables values.
#' @param tau A number or vector defining desired quantile.
#' @param ... Additional plot options.
#' @return This function returns a plot of the refcurve model.
#' @export

plot.refcurve <- function(x, newdata = data.frame(x = seq(0,1,0.01)), tau = seq(0.1,0.9,0.2),...){
mean <- predict(x$modelo_m, newdata = newdata)
sd <- sqrt(exp(predict(x$modelo_v, newdata = newdata)))
epsilon <- quantile(x$res, probs = tau)

cuants = matrix(0, ncol = length(epsilon), nrow = dim(newdata)[1])
for(i in 1:length(epsilon))  cuants[,i] = mean + sd*epsilon[i]

plot(rev(x$modelo_m$model), col = "grey")
invisible(apply(cuants, 2, function(x) lines(x ~ newdata[,1])))
}
