#' Prediction for a bivRegion object
#'
#' This function takes a bivRegion object and allow to obtain region limits
#' for a given covariate value if cond=TRUE. If not, it can be applied
#' to a new dataset to evaluate which points will be included inside the
#' standarized region (for instance, if we estimate a reference region with
#' healthy patients results, we can know where non healthy patients will be located).
#'
#' @param object A bivRegion object.
#' @param tau A number or vector defining the desired coverage(s) of the
#' bivariate reference/probabilistic region.
#' @param newdata An optional data frame containing new covariate values, or
#' new observations.
#' @param cond A logical argument, if TRUE, a matrix of values defining the
#' conditional reference region limits is given.
#' @param ... Additional prediction parameters.
#' @return This function returns reference region limits values for
#' a given covariate value (if cond=TRUE), or which observations falls
#' outside the estimated region (if cond=FALSE).
#' @export
#' @importFrom "sp" "SpatialPoints" "SpatialPolygons" "Polygons" "Polygon" "over"

predict.bivRegion <- function(object, tau = 0.95, newdata = NULL, cond = FALSE, ...) {
  if (sum(tau %in% object$tau) != length(tau)) {
    return("WARNING: this coverage probability should be previously estimated with bivregion")
  }
  bb <- which((object$tau) %in% tau)
  Out_points <- list()
  Coverage <- as.numeric()

  #-----------------------
  ## ---- Bivariate points
  #-----------------------
  if (is.null(object$fit) & !is.null(newdata)) {
    for (i in 1:length(bb)) {
      pts_test <- SpatialPoints(newdata[, names(object$Y)])
      spol <- SpatialPolygons(list(Polygons(list(Polygon(cbind(object$region[[bb[i]]][, 1], object$region[[bb[i]]][, 2]))), ID = " ")))
      out_res <- as.data.frame(newdata[which(is.na(over(pts_test, spol))), names(object$Y)])
      Coverage[i] <- 1 - (nrow(out_res) / nrow(newdata))

      cuadrant_1 <- which(out_res[, 1] > apply(object$Y, 2, mean)[1] & out_res[, 2] > apply(object$Y, 2, mean)[2])
      cuadrant_2 <- which(out_res[, 1] < apply(object$Y, 2, mean)[1] & out_res[, 2] > apply(object$Y, 2, mean)[2])
      cuadrant_3 <- which(out_res[, 1] < apply(object$Y, 2, mean)[1] & out_res[, 2] < apply(object$Y, 2, mean)[2])
      cuadrant_4 <- which(out_res[, 1] > apply(object$Y, 2, mean)[1] & out_res[, 2] < apply(object$Y, 2, mean)[2])

      Out_points[[i]] <- list(out_res[cuadrant_1, ], out_res[cuadrant_2, ], out_res[cuadrant_3, ], out_res[cuadrant_4, ])
      names(Out_points[[i]]) <- c(
        "Both high", paste0(names(out_res)[1], " low and ", names(out_res)[2], " high"),
        "Both low", paste0(names(out_res)[1], " high and ", names(out_res)[2], " low")
      )
    }
    names(Out_points) <- paste0("Tau = ", object$tau[bb])
    return(list(Values_out_of_reference_region_tau = Out_points, Coverage = Coverage))
  }

  #-----------------
  ## bivRegr object
  #-----------------
  if (!is.null(object$fit) & !is.null(newdata) & cond == FALSE) {
    newdata <- na.omit(newdata[, unique(c(
      all.vars(object$fit$formula[[1]])[1], all.vars(object$fit$formula[[2]])[1],
      do.call(c, lapply(object$fit$formula, all.vars))
    ))])

    pred_mu1 <- predict(object$fit$mu1, newdata = newdata)
    pred_mu2 <- predict(object$fit$mu2, newdata = newdata)
    pred_sd1 <- sqrt(exp(predict(object$fit$var1, newdata = newdata)))
    pred_sd2 <- sqrt(exp(predict(object$fit$var2, newdata = newdata)))
    pred_rho <- tanh(predict(object$fit$rho, newdata = newdata))

    newdata$res1 <- (eval(parse(text = paste0("newdata$", all.vars(object$fit$mu1$formula)[1]))) - pred_mu1) / pred_sd1
    newdata$res2 <- (eval(parse(text = paste0("newdata$", all.vars(object$fit$mu2$formula)[1]))) - pred_mu2) / pred_sd2

    YC <- cbind(newdata$res1, newdata$res2)

    for (i in 1:nrow(newdata)) {
      Sigma <- matrix(c(1, pred_rho[i], pred_rho[i], 1), ncol = 2, nrow = 2)
      P <- chol(solve(Sigma))
      YC[i, ] <- P %*% YC[i, ]
    }

    YC <- as.data.frame(YC)
    names(YC) <- c(all.vars(object$fit$mu1$formula)[1], all.vars(object$fit$mu2$formula)[1])

    for (i in 1:length(bb)) {
      pts_test <- SpatialPoints(YC)
      spol <- SpatialPolygons(list(Polygons(list(Polygon(cbind(object$region[[bb[i]]][, 1], object$region[[bb[i]]][, 2]))), ID = " ")))
      newdata_out <- as.data.frame(newdata[which(is.na(over(pts_test, spol))), ])
      out_res <- as.data.frame(YC[which(is.na(over(pts_test, spol))), ])
      names(out_res) <- c(paste(names(newdata_out)[1], " residuals"), paste(names(newdata_out)[2], " residuals"))
      Coverage[i] <- 1 - (nrow(out_res) / nrow(newdata))

      cuadrant_1 <- which(out_res[, 1] > apply(object$Y, 2, mean)[1] & out_res[, 2] > apply(object$Y, 2, mean)[2])
      cuadrant_2 <- which(out_res[, 1] < apply(object$Y, 2, mean)[1] & out_res[, 2] > apply(object$Y, 2, mean)[2])
      cuadrant_3 <- which(out_res[, 1] < apply(object$Y, 2, mean)[1] & out_res[, 2] < apply(object$Y, 2, mean)[2])
      cuadrant_4 <- which(out_res[, 1] > apply(object$Y, 2, mean)[1] & out_res[, 2] < apply(object$Y, 2, mean)[2])

      newdata_out <- cbind(newdata_out, out_res)
      Out_points[[i]] <- list(newdata_out[cuadrant_1, ], newdata_out[cuadrant_2, ], newdata_out[cuadrant_3, ], newdata_out[cuadrant_4, ])
      names(Out_points[[i]]) <- c(
        "Both high", paste0(names(newdata_out)[1], " low and ", names(newdata_out)[2], " high"),
        "Both low", paste0(names(newdata_out)[1], " high and ", names(newdata_out)[2], " low")
      )
    }
    names(Out_points) <- paste0("Tau = ", object$tau[bb])
    return(list(Values_out_of_reference_region_tau = Out_points, Coverage = Coverage))
  }

  #---------------------------
  ## bivRegr conditional case
  #---------------------------
  if (!is.null(object$fit) & !is.null(newdata) & cond == TRUE) {
    vars <- unique(c(
      names(object$fit$mu1$model)[-1], names(object$fit$mu2$model)[-1],
      names(object$fit$var1$model)[-1], names(object$fit$var2$model)[-1], names(object$fit$rho$model)[-1]
    ))
    vars <- vars[-which(vars == "(weights)")]

    if (sum(tau %in% object$tau) != length(tau)) {
      return("WARNING: this coverage probability should be previously estimated with bivregion")
    }

    bb <- which((object$tau) %in% tau)

    if (is.null(newdata)) {
      medias <- cbind(predict(object$fit$mu1), predict(object$fit$mu2))
    } else {
      medias <- cbind(predict(object$fit$mu1, newdata), predict(object$fit$mu2, newdata))
    }

    if (is.null(newdata)) {
      sds <- sqrt(cbind(exp(predict(object$fit$var1)), exp(predict(object$fit$var2))))
    } else {
      sds <- sqrt(cbind(exp(predict(object$fit$var1, newdata)), exp(predict(object$fit$var2, newdata))))
    }

    if (is.null(newdata)) {
      R <- predict(object$fit$rho)
      R[R > 1.8] <- 1.8
      R[R < (-1.8)] <- -1.8
      R <- tanh(R)
    } else {
      R <- predict(object$fit$rho, newdata = newdata)
      R[R > 1.8] <- 1.8
      R[R < (-1.8)] <- -1.8
      R <- tanh(R)
    }

    pred_tau <- list()
    for (j in 1:length(bb)) {
      pred_regions <- list()

      for (i in 1:dim(medias)[1]) {
        Sigma <- matrix(c(1, R[i], R[i], 1), ncol = 2, nrow = 2)
        P <- chol(solve(Sigma))
        P <- solve(P)

        pred_regions[[i]] <- cbind(
          medias[i, 1] + (t(apply(object$region[[bb[j]]], 1, function(x) P %*% x)))[, 1] * sds[i, 1],
          medias[i, 2] + (t(apply(object$region[[bb[j]]], 1, function(x) P %*% x)))[, 2] * sds[i, 2]
        )
      }

      if (!is.null(newdata)) {
        for (k in 1:dim(medias)[1]) names(pred_regions)[k] <- paste(names(newdata), " = ", newdata[k, ], collapse = "&")
      } else {
        for (k in 1:dim(medias)[1]) names(pred_regions)[k] <- paste(vars, "=", format(round(object$data[k, vars], 4)), collapse = " & ")
      }

      pred_tau[[j]] <- pred_regions
    }

    names(pred_tau) <- paste("Coverage probability = ", format(object$region$tau[bb], nsmall = 2))
    return(pred_tau)
  }
}
