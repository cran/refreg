#' Trivariate reference region estimation
#'
#' This functions estimate a probabilistic/reference region for trivariate data.
#' It is based on a non parametric kernel density estimation. It can only be applied to
#' a trivRegr object, and for one single tau.
#'
#' @param fit A trivRegr object.
#' @param tau A number defining the desired coverage of the trivariate
#' reference region.
#' @return This function return a region containing a given percentage of trivariate data points.
#' @references Duong, T. (2019) ks: Kernel Smoothing. R package version 1.11.6. https://CRAN.R--project.org/package=ks.
#' @export
#' @importFrom "ks" "kde"

trivRegion=function(fit,tau=0.90){
  kernel.3d = kde(fit$trivres,H=diag(c(0.5,0.5,0.5)),
                  approx.cont = TRUE,gridsize = c(49,50,51))
  ks_hat = predict(kernel.3d,x=as.matrix(fit$trivres))


  names(kernel.3d$cont) =as.numeric(substr(names(kernel.3d$cont),1,nchar(names(kernel.3d$cont))-1))/100
  k = as.numeric(which(ks_hat<kernel.3d$cont[as.character(1-tau)]))
  l = list(kernel_fit = kernel.3d, predict_kernel = ks_hat, which_out = k,
           fit=fit,trivres = fit$trivres,k_limit= as.numeric(kernel.3d$cont[as.character(1-tau)]))
  class(l) = "trivRegion"
  return(l)
}
