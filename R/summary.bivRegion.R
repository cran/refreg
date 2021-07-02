#' bivRegion summary method
#'
#' This function takes an bivRegion object and indicates which observations are
#' located outside the estimated reference region and why.
#'
#' @param object A bivRegion object.
#' @param tau The data coverage proportion previously obtained by the bivRegion function.
#' @param ... Additional summary options.
#' @return This function indicates the region apparent coverage, and which
#' patients are located outside the estimated reference region.
#' @export

summary.bivRegion = function(object,tau = 0.95,...){
if(sum(tau%in%object$tau)!=length(tau)){return("WARNING: this coverage probability should be previously estimated with bivregion")}

  bb = which((object$tau) %in% tau)
  Out_points = list()
  Coverage = as.numeric()

## Sin modelo e sin newdata
if(is.null(object$fit)){
  for(i in 1:length(bb)){
    pts_test = SpatialPoints(object$Y)
    spol = SpatialPolygons(list(Polygons(list(Polygon(cbind(object$region[[bb[i]]][,1],object$region[[bb[i]]][,2]))),ID = " ")))
    out_res = as.data.frame(object$Y[which(is.na(over(pts_test,spol))),])
    Coverage[i] = 1-(nrow(out_res)/nrow(object$Y))

    cuadrant_1 = which(out_res[,1] > apply(object$Y,2,mean)[1] & out_res[,2] > apply(object$Y,2,mean)[2])
    cuadrant_2 = which(out_res[,1] < apply(object$Y,2,mean)[1] & out_res[,2] > apply(object$Y,2,mean)[2])
    cuadrant_3 = which(out_res[,1] < apply(object$Y,2,mean)[1] & out_res[,2] < apply(object$Y,2,mean)[2])
    cuadrant_4 = which(out_res[,1] > apply(object$Y,2,mean)[1] & out_res[,2] < apply(object$Y,2,mean)[2])

    Out_points[[i]] = list(out_res[cuadrant_1,],out_res[cuadrant_2,],out_res[cuadrant_3,],out_res[cuadrant_4,])
    names(Out_points[[i]]) = c("Both high",paste0(names(out_res)[1]," low and ",names(out_res)[2]," high"),
                               "Both low",paste0(names(out_res)[1]," high and ",names(out_res)[2]," low"))
  }
  names(Out_points) = paste0("Tau = ",object$tau[bb])
}

## Con modelo e sin newdata
if(!is.null(object$fit)){
  for(i in 1:length(bb)){
    pts_test = SpatialPoints(object$Y)
    spol = SpatialPolygons(list(Polygons(list(Polygon(cbind(object$region[[bb[i]]][,1],object$region[[bb[i]]][,2]))),ID = " ")))
    data_out = as.data.frame(object$data[which(is.na(over(pts_test,spol))),])
    out_res = as.data.frame(object$Y[which(is.na(over(pts_test,spol))),])
    names(out_res) = c(paste0(names(data_out)[1],"-res"),
                       paste0(names(data_out)[2],"-res"))
    Coverage[i] = 1-(nrow(out_res)/nrow(object$data))

    cuadrant_1 = which(out_res[,1] > apply(object$Y,2,mean)[1] & out_res[,2] > apply(object$Y,2,mean)[2])
    cuadrant_2 = which(out_res[,1] < apply(object$Y,2,mean)[1] & out_res[,2] > apply(object$Y,2,mean)[2])
    cuadrant_3 = which(out_res[,1] < apply(object$Y,2,mean)[1] & out_res[,2] < apply(object$Y,2,mean)[2])
    cuadrant_4 = which(out_res[,1] > apply(object$Y,2,mean)[1] & out_res[,2] < apply(object$Y,2,mean)[2])

    data_out = cbind(data_out,out_res)
    Out_points[[i]] = list(data_out[cuadrant_1,],data_out[cuadrant_2,],
                           data_out[cuadrant_3,],data_out[cuadrant_4,])
    names(Out_points[[i]]) = c("Both high",paste0(names(data_out)[1]," low and ",names(data_out)[2]," high"),
                               "Both low",paste0(names(data_out)[1]," high and ",names(data_out)[2]," low"))
  }
  names(Out_points) = paste0("Tau = ",object$tau[bb])
}
  return(list(Out_points = Out_points,Coverage = Coverage))
}

