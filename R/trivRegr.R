#' Trivariate regression model
#'
#' This function estimates the covariates effects on the means vector,
#' and variance covariance matrix of a trivariate variable. Non linear effects
#' might be estimated for continuous covariates using penalized splines.
#'
#' @param f A list of 9 formulas defining the covariates effects
#' in three responses means, in their variances, and in their correlations. The formulas
#' follow the mgcv::gam() structure.
#' @param data A data frame containing the reponses, and predictor variables values.
#' @return This function returns the covariates effect on the means,
#' variances, and correlation for a trivariate response variable.
#' @importFrom "mbend" "bend"
#' @importFrom "matrixcalc" "is.positive.semi.definite"
#' @examples
#'dm_no = subset(aegis,aegis$dm=="no")
#'
#' # Model formulas
#' mu1 = fpg ~ s(age)
#' mu2 = hba1c ~ s(age)
#' mu3 = fru~s(age)
#' var1 = ~ s(age)
#' var2 = ~ s(age)
#' var3 = ~s(age)
#' theta12 =  ~ s(age)
#' theta13 =  ~ s(age)
#' theta23 =  ~ s(age)
#' f = list(mu1,mu2,mu3,var1,var2,var3,theta12,theta13,theta23)
#'
#' # Model fit
#' fit = trivRegr(f,data=dm_no)

#' # Trivariate region estimation
#' region = trivRegion(fit,tau=0.95)
#' plot(region,col=2,planes=TRUE)
#' plot(region,cond=TRUE,newdata=data.frame(age=c(20,80)),
#'      xlab="FPG, mg/dl",ylab="HbA1c, %",zlab="Fru, mg/dL")
#'
#' @export

trivRegr = function(f=f,data=data){
  modelo_m1=gam(f[[1]],data = data)
  modelo_m2=gam(f[[2]],data = data)
  modelo_m3=gam(f[[3]],data = data)

  medias=cbind(predict(modelo_m1),predict(modelo_m2),predict(modelo_m3))

  # Fit variance models
  # 1) Obtain the residuals
  data$res1 = residuals(modelo_m1)^2
  data$res2 = residuals(modelo_m2)^2
  data$res3 = residuals(modelo_m3)^2

  # 2) Fit the models
  modelo_v1 = ACE(y = "res1",predictor = f[[4]],data=data,restriction = "positive")
  modelo_v2 = ACE(y = "res2",predictor = f[[5]],data=data,restriction = "positive")
  modelo_v3 = ACE(y = "res3",predictor = f[[6]],data=data,restriction = "positive")

  sds=cbind(sqrt(exp(predict(modelo_v1$fit))),sqrt(exp(predict(modelo_v2$fit))),sqrt(exp(predict(modelo_v3$fit))))

  # Fit correlation model
  # 1) Obtain the residuals
  data$res1 = (eval(parse(text = paste0("data$",all.vars(f[[1]])[1]))) - medias[,1]) / sds[,1]
  data$res2 = (eval(parse(text = paste0("data$",all.vars(f[[2]])[1]))) - medias[,2]) / sds[,2]
  data$res3 = (eval(parse(text = paste0("data$",all.vars(f[[3]])[1]))) - medias[,3]) / sds[,3]

  data$res12 = data$res1*data$res2
  data$res13 = data$res1*data$res3
  data$res23 = data$res2*data$res3

  # 2) Fit the model
  modelo_rho1 = ACE(y = "res12",predictor = f[[7]],data=data,restriction = "correlation")
  modelo_rho2 = ACE(y = "res13",predictor = f[[8]],data=data,restriction = "correlation")
  modelo_rho3 = ACE(y = "res23",predictor = f[[9]],data=data,restriction = "correlation")

  R = cbind(tanh(predict(modelo_rho1$fit)),tanh(predict(modelo_rho2$fit)),tanh(predict(modelo_rho3$fit)))

  # Obtain the bivariate residuals
  n = dim(data)[1]
  YC=cbind(data$res1,data$res2,data$res3)

  for (i in 1:n) {
    Sigma=matrix(c(1,R[i,1],R[i,2],R[i,1],1,R[i,3],R[i,2],R[i,3],1),ncol=3,nrow=3,byrow=TRUE)

    if(is.positive.semi.definite(Sigma, tol=1e-8)==FALSE){Sigma = bend(Sigma)$bent}


    P=chol(solve(Sigma))
    YC[i,]=P%*%YC[i,]}
  YC = as.data.frame(YC)
  names(YC)=c(all.vars(f[[1]])[1],all.vars(f[[2]])[1],all.vars(f[[3]])[1])
  return(list(trivres = YC,mean1 = modelo_m1,mean2=modelo_m2,mean3=modelo_m3,
              var1=modelo_v1,var2=modelo_v2,var3=modelo_v3,
              rho1=modelo_rho1,rho2=modelo_rho2,rho3=modelo_rho3))
}
