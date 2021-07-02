#' Bivariate regression model
#'
#' This function estimates the covariates effects on the means vector,
#' and variance covariance matrix from a bivariate variable. Non linear effects
#' might be estimated for continuous covariates using penalized spline smoothers.
#'
#' @param f A list of five formulas defining the covariates effects
#' in both responses means, variances, and in their correlation. The formulas
#' follow the same structure as mgcv::gam() (see example below).
#' @param data A data frame containing the reponses, and predictor variables values.
#' @param ace.eps A number defining the error rate in the ACE algorithm.
#' @param ace.itmax A number defining the maximum number of ACE algorithm iterations.
#' @return This function returns the covariates effect on the means,
#' variances, and correlation of a bivariate response variable.
#' @importFrom "stats" "residuals" "pnorm" "rnorm" "sd"
#' @export
#' @examples
#' # Bivariate reference region for fasting plasma glucose (fpg)
#' #and glycated hemoglobin (hba1c) levels depending on age
#'
#' dm_no = subset(aegis,aegis$dm == "no")#select healthy patients
#' # 1.1) Define predictors
#' mu1 = fpg ~ s(age)
#' mu2 = hba1c ~ s(age)
#' var1 = ~ s(age)
#' var2 = ~ s(age)
#' rho = ~ s(age)
#' f = list(mu1,mu2,var1,var2,rho)
#'
#' fit = bivRegr(f,data=dm_no)
#'
#' # 1.2) Depict the estimated covariates effects
#' plot(fit,eq=1)
#' plot(fit,eq=2)
#' plot(fit,eq=3)
#' plot(fit,eq=4)
#' plot(fit,eq=5)
#' # 1.2.1) Depict the estimated covariates effects with CI (Not Run)
#' \donttest{
#' s0 = summary_boot(fit,B=100) #no parallelization
#' #s1 = summary_boot(fit,B=100,parallel=TRUE) #parallelization
#' plot(s0,eq=1)
#' }
#'
#' # 1.3) Obtain the reference region in the standarized residuals
#' region = bivRegion(fit,tau=0.95,shape=2)
#' plot(region)
#'
#' # 1.4) Identify those patients located outside the reference region
#' summary(region)
#'
#' # 1.5) Depict the conditional reference region for two ages
#' plot(region,cond=TRUE,newdata=data.frame(age=c(20,50)),col="grey",pch="*",
#' reg.lwd = 2,reg.lty=2)

bivRegr <- function(f = f, data = data, ace.eps =0.01 ,ace.itmax=25){
data = na.omit(data[,unique(c(all.vars(f[[1]])[1],all.vars(f[[2]])[1],
                              do.call(c,lapply(f,all.vars))))])

### Fit mean models
modelo_m1 <- gam(f[[1]], data = data, method = "REML")
modelo_m2 <- gam(f[[2]], data = data, method = "REML")
medias <- cbind(predict(modelo_m1), predict(modelo_m2))

### Fit variance models
# 1) Obtain the residuals
data$res1 <- residuals(modelo_m1)^2
data$res2 <- residuals(modelo_m2)^2

# 2) Fit the models
modelo_v1 <- ACE(y = "res1", predictor = f[[3]], restriction = "positive",
                 eps = ace.eps, itmax=ace.itmax, data = data)
modelo_v2 <- ACE(y = "res2", predictor = f[[4]], restriction = "positive",
                 eps = ace.eps, itmax=ace.itmax, data = data)


if(ace.itmax >=10&(modelo_v1$error>=ace.eps|modelo_v2$error>=ace.eps)){
  data$res1 <- log(residuals(modelo_m1)^2)
  data$res2 <- log(residuals(modelo_m2)^2)

  g1=as.formula(paste0("res1~",f[[3]])[2])
  g2=as.formula(paste0("res2~",f[[4]])[2])

  modelo_v1 <- gam(g1,data=data)
  modelo_v2 <- gam(g2,data=data)
}else{
  modelo_v1 = modelo_v1$fit
  modelo_v2 = modelo_v2$fit
}

sds <- cbind(sqrt(exp(predict(modelo_v1))), sqrt(exp(predict(modelo_v2))))

### Fit correlation model
# 1) Obtain the residuals
data$res1 <- (eval(parse(text = paste0("data$", all.vars(f[[1]])[1]))) - medias[,1]) / sds[,1]
data$res2 <- (eval(parse(text = paste0("data$", all.vars(f[[2]])[1]))) - medias[,2]) / sds[,2]

if(sd(data$res1)>=1.1|sd(data$res2)>=1.1){
  res_mu=c(mean(data$res1), mean(data$res2))
  res_sd=c(sd(data$res1), sd(data$res2))
  data$res1=(data$res1-res_mu[1])/res_sd[1]
  data$res2=(data$res2-res_mu[2])/res_sd[2]
}


data$res12 <- data$res1*data$res2

# 2) Fit the model
modelo_rho <- ACE(y = "res12", predictor = f[[5]], restriction = "correlation", eps = ace.eps, itmax = ace.itmax,data = data)$fit
R <- predict(modelo_rho);R[R>1.8]=1.8;R[R<(-1.8)]=-1.8
R <- tanh(R)


#### Obtain the bivariate residuals
n <- dim(data)[1]
YC <- cbind(data$res1, data$res2)

for (i in 1:n) {
  Sigma <- matrix(c(1, R[i], R[i], 1), ncol = 2, nrow = 2)
  P <- chol(solve(Sigma))
  YC[i,] <- P %*% YC[i,]
}
YC = as.data.frame(YC)
names(YC) = c(all.vars(f[[1]])[1],all.vars(f[[2]])[1])


data = data[,-which(names(data)%in%c("res1","res2","res12"))]
ret <- list(formula = f, bivres = YC,mu1 = modelo_m1, mu2 = modelo_m2, var1 = modelo_v1, var2 = modelo_v2,
              rho = modelo_rho,means = medias, std.variations = sds, correlation = R,data=data)
class(ret) <- "bivRegr"
return(ret)
}






