#' Plot method for bivRegr fit
#'
#' This function takes an bivRegr object and plots the estimated effects for the
#' conditional response means, variances or correlation. summary_boot function
#' must be applied by the user in order to get the estimated effects confidence
#' intervals.
#'
#' @param x A bivRegr fit.
#' @param eq A number indicating the model effects to be depicted; 1 = first response mean,
#' 2 = second response mean, 3 = first response variance, 4 = second response variance, and
#' 5 = correlation model.
#' @param newdata An optional data frame containing covariates values.
#' @param ... Additional plot arguments.
#' @return This function returns a plot for bivRegr mean, variance and correlation models.
#' @importFrom "gridExtra" "grid.arrange"
#' @importFrom "graphics" "legend"
#' @importFrom "stringr" "str_to_title"
#' @import "ggplot2"
#' @import "RColorBrewer"
#' @export

plot.bivRegr = function(x,eq=1,newdata=NULL,...){

    if(eq == 1) model = x$mu1
    if(eq == 2) model = x$mu2
    if(eq == 3) model = x$var1
    if(eq == 4) model = x$var2
    if(eq == 5) model = x$rho

    model_pred = predict(model,type="terms")


    results = list()
    for(i in 1:ncol(model_pred)){
      if(str_detect(colnames(model_pred)[i],"s")==TRUE){
        covar_data = x$data[,which(str_detect(colnames(model_pred)[i],names(x$data))),drop=FALSE]
      }else{
        covar_data = x$data[,colnames(model_pred)[i],drop=FALSE]}

      fit_pred = cbind(covar_data,model_pred[,colnames(model_pred)[i],drop=FALSE])
      fit_pred = unique(fit_pred)
      fit_pred = as.data.frame(fit_pred)
      fit_pred = fit_pred[order(fit_pred[,1]),]

      colnames(fit_pred) = c(colnames(model_pred)[i],"fit")
      results[[i]] = fit_pred
    }
    names(results) = colnames(model_pred)

    myplot <- function(data){
      if(class(data[,1]) =="factor"){
        g1 = ggplot(data, aes(x = data[,1], y = data[,2])) + geom_point(pch = 13,cex=5,col="red") + theme_classic()
        g1 = g1 + xlab(names(data)[1]) + ylab(str_to_title(paste0(names(data)[1]," centered effect")))
        if(eq==1){g1 <-g1+ ggtitle(paste0(names(model$model)[1]," mean"))}
        if(eq==2){g1 <-g1+ ggtitle(paste0(names(model$model)[1]," mean"))}
        if(eq==3){g1 <-g1+ ggtitle(paste0(names(x$mu1$model)[1]," variance"))}
        if(eq==4){g1 <-g1 + ggtitle(paste0(names(x$mu2$model)[1]," variance"))}
        if(eq==5){g1 <-g1 + ggtitle(paste0(names(x$mu1$model)[1],"-",names(x$mu2$model)[1]," correlation"))}
        g1
        }else{
        g1 = ggplot(data, aes(x = data[,1], y = data[,2]))
        g1 = g1+geom_line() + theme_classic()+ geom_rug(aes(x=data[,1],y = NULL),col="grey")

        if( names(data)[1] %in% do.call(c,lapply(model$smooth, function(x) x$label)) == TRUE){
          aa = which(do.call(c,lapply(model$smooth, function(x) x$label)) %in% names(data[1]))
          first <- model$smooth[[aa]]$first.para
          last <- model$smooth[[aa]]$last.para
          edf <- format(round(sum(model$edf[first:last]),2),nsmall=2)
        }

        g1 = g1 + xlab(model$smooth[[aa]]$term) + ylab(paste0("s(",model$smooth[[aa]]$term,",",edf,")"))
        if(eq==1){g1 <-g1+ ggtitle(paste0(names(model$model)[1]," mean"))}
        if(eq==2){g1 <-g1+ ggtitle(paste0(names(model$model)[1]," mean"))}
        if(eq==3){g1 <-g1+ ggtitle(paste0(names(x$mu1$model)[1]," variance"))}
        if(eq==4){g1 <-g1 + ggtitle(paste0(names(x$mu2$model)[1]," variance"))}
        if(eq==5){g1 <-g1 + ggtitle(paste0(names(x$mu1$model)[1],"-",names(x$mu2$model)[1]," correlation"))}

         g1
      }
    }

    plots = lapply(results, myplot)
    if(length(plots)<4){return(grid.arrange(grobs = plots,ncol=length(plots)))}
    if(length(plots)>=4){return(grid.arrange(grobs = plots,ncol=length(plots)/2,nrow=length(plots)/2))}
}
