#' Bivariate reference region estimation
#'
#' This functions estimate a probabilistic/reference region for bivariate data. It is
#' based on a kernel density estimation. It may be applied to a
#' set of bivariate data points, or to a bivRegr object. In the former case, the
#' function will estimate a bivariate reference region for the model standarized
#' residuals.
#'
#' @param Y A set of bivariate data points, or a bivRegr object.
#' @param tau A number or vector defining the desired coverage(s) of the bivariate
#' reference region.
#' @param H_choice Kernel bandwidth selection method: "plug.in" for plug.in method,
#' "LSCV" for least squate cross valiation, "SCV" for smooth cross validation,
#' and "Hcov" for a bandwidth selection method which optimize the region coverage.
#' @param display_plot A logical indicating if plot must be displayed during "Hcov"
#' bandwidht estimation procedure. The plot depicts region's coverage, evaluated
#' with k fold cross validation, depending on kernel bandwidth value.
#' @param k In case of using "Hcov" the number of k fold cross validations
#' replicates to be performed.
#' @param shape Shape parameter modulating the final shape of the bivariate
#' probabilistic/reference region by hand.
#' @param ... Additional parameters to be modified in KernSmooth::bkde2D()
#' function by the user (e.g. gridsize).
#' @return This function return a region or a set of regions containing a given
#' percentage of bivariate data points.
#' @references Duong, T. (2019) ks: Kernel Smoothing. R package version 1.11.6. https://CRAN.R--project.org/package=ks.
#' @references Matt Wand (2020). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995). R package version 2.23--18. https://CRAN.R--project.org/package=KernSmooth
#' @export
#' @examples
#' Y <- cbind(rnorm(100), rnorm(100))
#' Y <- as.data.frame(Y)
#' names(Y) <- c("y1", "y2")
#' region <- bivRegion(Y, tau = 0.95, shape = 2)
#' plot(region)
#' @importFrom "KernSmooth" "dpik" "bkde2D"
#' @importFrom "ks" "Hlscv.diag" "Hscv.diag"
#' @importFrom "grDevices" "contourLines"
#' @importFrom "graphics" "abline" "lines" "plot"
#' @importFrom "methods" "is"


bivRegion <- function(Y = fit, H_choice = "Hcov", tau = 0.95, k = 20,
                      display_plot = TRUE, shape = NULL, ...) {
  if (is(object=Y,class2 = "bivRegr")) {
    fit <- Y
    data <- Y$data
    Y <- Y$bivres
  } else {
    fit <- NULL
    data <- NULL
  }

  ## Bandwidht selection
  if (is.null(shape)) {
    if (H_choice == "Hcov") {
      HCV <- Hcov(Y = Y, shape = seq(1, 10, 0.25), k = k, tau = max(tau), display_plot = display_plot)
      bandwidth <- c(dpik(Y[, 1]), dpik(Y[, 2])) * HCV$best_shape
    }

    if (H_choice == "plug.in") {
      bandwidth <- c(dpik(Y[, 1]), dpik(Y[, 2]))
    }

    if (H_choice == "LSCV") {
      bandwidth <- diag(Hlscv.diag(Y))
    }

    if (H_choice == "SCV") {
      bandwidth <- diag(Hscv.diag(Y))
    }
  }

  if (!is.null(shape)) {
    bandwidth <- c(dpik(Y[, 1]), dpik(Y[, 2])) * shape
  }

  # Kernel estimation
  surf <- bkde2D(Y, bandwidth = bandwidth, ...)

  # Region definition
  cl <- contourLines(surf[[1]], surf[[2]], surf[[3]], nlevels = 500)
  pts <- SpatialPoints(Y)
  pin <- numeric(length(cl))

  regions <- list()
  for (i in 1:length(cl)) {
    region <- Polygon(cbind(cl[[i]]$x, cl[[i]]$y))
    region <- Polygons(list(region), ID = " ")
    region <- SpatialPolygons(list(region))
    pin[i] <- sum(!is.na(over(pts, region)))
    regions[[i]] <- region
  }
  pin <- pin / nrow(Y)

  tau_region <- list()
  cob <- as.numeric()
  for (i in 1:length(tau)) {
    best <- which.min(abs(pin - tau[i])) # podo adaptalos a cantos queira
    best.region <- regions[[best]]@polygons[[1]]@Polygons[[1]]
    tau_region[[i]] <- cbind(xcord = best.region@coords[, 1], ycord = best.region@coords[, 2])
    cob[i] <- pin[best]
  }
  l <- list(region = tau_region, tau = tau, Y = Y, fit = fit, bandwidth = bandwidth, data = data)
  class(l) <- "bivRegion"
  return(l)
}







## If you are reading this code. In the following I have an alternative way
# of estimating the bivariate region
##' arroba/importFrom "stats" "predict" "quantile" "as.formula" "na.omit"

# bivRegion <- function(Y, tau = 0.95, shape = 2){
# if(class(Y)=="bivRegr"){fit <- Y;data <- Y$data;Y <- Y$bivres}else{fit=NULL;data=NULL}

# bandwidth <- c(dpik(Y[,1]), dpik(Y[,2]))

# f <- kde(Y, H = diag(bandwidth * shape), approx.cont = T, gridsize = c(50,50))
# fhat <- predict(f, x = Y)
# k <- quantile(fhat, probs = 1-tau)

# region <- list()
# for(i in 1:length(tau)){
# cl1 <- contourLines(f$eval.points[[1]], f$eval.points[[2]], f$estimate,levels = k[i])[[1]]
# region[[i]] <- as.data.frame(cbind(cl1$x, cl1$y))
# names(region[[i]]) = names(Y)
# }

# l <- list(region = region, ks = k, f = fhat, tau = tau, Y = Y,fit = fit,data = data)
# class(l) <- "bivRegion"
# return(l)
# }
