#' Predict method for bivRegr
#'
#' Obtains predictions for a bivRegr object, that is the bivRegr means,
#' variances, and correlation for given covariate values.
#'
#' @param object A bivRegr fit.
#' @param newdata A data frame defining the covariate values for prediction. Default
#' is NULL and the prediction will be done in the original data.
#' @param ... Additional predict options.
#' @return This function returns prediction of bivRegr mean, variance and correlation models.
#' @export

predict.bivRegr = function(object,newdata=NULL,...){

  vars = unique(c(names(object$mu1$model)[-1],names(object$mu2$model)[-1],
                  names(object$var1$model)[-1],names(object$var2$model)[-1],
                  names(object$rho$model)[-1]))
  vars = vars[-which(vars=="(weights)")]

hat_mu1 = if(is.null(newdata)){predict(object$mu1)}else{predict(object$mu1,newdata = newdata)}
hat_mu2 = if(is.null(newdata)){predict(object$mu2)}else{predict(object$mu2,newdata = newdata)}
hat_sd1 = if(is.null(newdata)){sqrt(exp(predict(object$var1)))}else{sqrt(exp(predict(object$var1,newdata = newdata)))}
hat_sd2 = if(is.null(newdata)){sqrt(exp(predict(object$var2)))}else{sqrt(exp(predict(object$var2,newdata = newdata)))}
if(is.null(newdata)){R = predict(object$rho); R[R>1.8]=1.8; R[R<(-1.8)]=-1.8;R = tanh(R)}else{R = predict(object$rho,newdata=newdata); R[R>1.8]=1.8; R[R<(-1.8)]=-1.8;R = tanh(R)}
hat_rho = R

if(is.null(newdata)){`Predicted Parameters` = as.data.frame(cbind(object$data[,vars],hat_mu1,hat_sd1,hat_mu2,hat_sd2,hat_rho))}else{`Predicted Parameters` = as.data.frame(cbind(newdata,hat_mu1,hat_sd1,hat_mu2,hat_sd2,hat_rho))}
if(is.null(newdata)){l=dim(object$data[,vars])[2]}else{l = dim(newdata)[2]}
names(`Predicted Parameters`)[l+c(1:5)] = c(paste0("Mean ",names(object$mu1$model)[1]),
                                                             paste0("Std.dev. ",names(object$mu1$model)[1]),paste0("Mean ",names(object$mu2$model)[1]),
                                                             paste0("Std.dev. ",names(object$mu2$model)[1]),
                                                             paste0("Corr. (",names(object$mu1$model)[1],"-",names(object$mu2$model)[1],")"))
return(`Predicted Parameters`)
}

