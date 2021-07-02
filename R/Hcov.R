#' Kernel bandwidth selection method based on bivariate density
#' contours coverage
#'
#' This function implements a method for estimating bivariate kernel bandwidth
#' based on data covarage. The method starts with the plug in estimate (which
#' usually overfits the data), and then increase this bandwidth value until the
#' desired coverage is obtained. Region coverage is evaluated in an out sample
#' design, using a k fold cross validation scheme.
#'
#' @param Y A matrix containing bivariate data values.
#' @param shape A sequence of values which controls plug in estimator increasing.
#' @param k A number indicating k fold cross validations to be performed.
#' @param tau The desired region coverage
#' @param display_plot A logical indicating if a plot must be displaying, during
#' the function estimation process, summarizing the results.
#' @return This function return a diagonal kernel bandwidth matrix.
#' @export
#' @importFrom "KernSmooth" "dpik" "bkde2D"
#' @importFrom "grDevices" "contourLines"
#' @importFrom "graphics" "abline" "lines" "plot"
#' @importFrom "sp" "SpatialPolygons" "SpatialPoints" "over"
#' @importFrom "pracma" "akimaInterp"
#' @importFrom "grDevices" "dev.off"
#' @importFrom "graphics" "boxplot"
#' @importFrom "stats" "as.formula" "na.omit" "predict" "quantile"

# Incluir funcion Hcv
Hcov = function(Y,shape = seq(1,10,0.5),k=20,tau=0.90,display_plot=TRUE){
data = Y
train_test_coverage = function(data=data,tau = 0.90, shape = 1){
    ## Evaluate the bivariate region coverage in a CV scheme

    # Define a train sample
    train = sample(nrow(data),0.70*nrow(data))

    # Bivariate kernel estimation
    x_train = data[train,1]
    y_train = data[train,2]
    bandwidth <- c(dpik(x_train), dpik(y_train))#plug in bandwidht
    surf <- bkde2D(cbind(x_train,y_train),bandwidth=bandwidth*shape)

    # Estimate f(u,v) such that P(e1,e2) in f(u,v) = tau
    pts = SpatialPoints(cbind(x_train,y_train))
    cl = contourLines(surf[[1]], surf[[2]], surf[[3]],nlevels = 1000)
    cl = cl[which(do.call(rbind,lapply(cl,function(x) length(x$x)))>=10)]
    pin = as.numeric()

    spols = list()
    for(i in 1:length(cl)){
      spol = Polygon(cbind(cl[[i]]$x,cl[[i]]$y))
      spol = Polygons(list(spol),ID = " ")
      spol = SpatialPolygons(list(spol))
      pin[i] = sum(!is.na(over(pts,spol)))
      spols[[i]] = spol
      if(abs(round(pin[i]/length(x_train),2)-tau)==0){break}#esto engadino para que non busque sen falta
    }
    pin =  pin/length(x_train)


    # Pick up the R_tau region
    rexion = cob = NULL
    best = which.min(abs(pin - tau))
    best.spol = spols[[best]]@polygons[[1]]@Polygons[[1]]
    rexion = cbind(xcord = best.spol@coords[,1],ycord = best.spol@coords[,2])

    ## Evaluate the region
    pts_test = SpatialPoints(cbind(data[-train,1],data[-train,2]))
    spol = SpatialPolygons(list(Polygons(list(Polygon(cbind(rexion[,1],rexion[,2]))),ID = " ")))

    return(test_cov = 100*(sum(!is.na(over(pts_test,spol)))/length(data[-train,1])))
  }




  #Perform the a train test scheme k times, for a sequence of shape parameters
  cov_shape=NULL
  for(i in 1:length(shape)){
    cov_rep = replicate(k,train_test_coverage(data = data,tau=tau,shape=shape[i]))

    plug_in = c(dpik(data[,1]), dpik(data[,2]))
    bandwidth = format(round(shape[i] * plug_in,2),nsmall=2)
    bandwidth = paste0("(",bandwidth[1],",",bandwidth[2],")")

    cov_shape = rbind(cov_shape,cbind(rep(bandwidth,k),rep(shape[i],k),cov_rep))

    if(display_plot==TRUE){
      if(i!=1){dev.off()}
      boxplot(as.numeric(as.character(cov_shape[,3]))~cov_shape[,1],
              xlab="Kernel bandwidth",ylab="Test coverage, %",
              main="Best coverage kernel bandwidth",outline=F)
      abline(h=tau*100,col=2,lwd=2)}

    if(median(cov_rep)>=(100*tau)){break}# if it overcome the nominal coverage STOP
  }

  # Interpolation of the coverage for a sequence of shape parameters
  cov_shape = as.data.frame(cov_shape)
  colnames(cov_shape) = c("bandwidth","shape","cov")
  cov_shape$shape = as.numeric(as.character(cov_shape$shape))
  cov_shape$cov = as.numeric(as.character(cov_shape$cov))

  t = tapply(cov_shape$cov,cov_shape$shape,median)

  if(length(t)>1){interpol = akimaInterp(x=as.numeric(names(t)),y=as.numeric(t),
                                     xi=seq(1,max(as.numeric(names(t))),0.01))
  best_shape=seq(1,max(as.numeric(names(t))),0.01)[which.min(abs(interpol-(100*tau)))]}else{best_shape=as.numeric(t)}

  return(list(cov_shape=cov_shape,best_shape=best_shape,
              interpolados = cbind(seq(1,max(as.numeric(names(t))),0.01),
                                   interpol)))
}
