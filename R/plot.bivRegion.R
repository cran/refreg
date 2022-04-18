#' Plot a bivRegion object
#'
#' This function allow to depict the estimated bivariate reference/probabilistic region, in
#' the estandarized residuals scale (cond=FALSE), or for any covariate value (cond=TRUE).
#'
#' @param x A bivRegion object.
#' @param tau A number, or vector, defining the desired coverage(s) of the bivariate reference region.
#' @param cond A logical argument, if TRUE a conditional reference region is depicted.
#' @param newdata If cond=FALSE, a data.frame with new values to be depicted in
#' the standarized residuals scale. If cond=TRUE, a data frame containing covariate
#' values for which the reference/probabilistic region will be depicted.
#' @param add A logical argument, if TRUE the conditional reference region is depicted over a pre existing plot.
#' @param legend A logical argument, if TRUE a legend is given along with the reference region.
#' @param reg.col Region line colour, in case of more than one tau it can be a vector.
#' @param reg.lwd Region line width, in case of more than one tau it can be a vector.
#' @param reg.lty Region line type, in case of more than one tau it can be a vector.
#' @param axes Logical; if TRUE (and cond=FALSE), vertical and horizontal lines
#' are added indicating four quadrants in the model residuals scale.
#' @param axes.col Axes colour.
#' @param axes.lwd Axes line width.
#' @param ... Further plot parameters.
#' @return This function return a graphical representation for a bivRegion object.
#' @examples
#' Y <- cbind(rnorm(100), rnorm(100))
#' Y <- as.data.frame(Y)
#' names(Y) <- c("y1", "y2")
#' reg <- bivRegion(Y, tau = 0.95, shape = 2)
#' plot(reg)
#' @importFrom "stats" "median"
#' @export

plot.bivRegion <- function(x, tau = 0.95, newdata = NULL, reg.col = NULL, reg.lwd = 1, reg.lty = NULL, axes = TRUE,
                           axes.col = "black", axes.lwd = 2L, cond = FALSE, add = FALSE, legend = TRUE, ...) {
  if (sum(tau %in% x$tau) != length(tau)) {
    return("WARNING: this coverage probability should be previously estimated with bivregion")
  }

  if (length(reg.col) == 1) {
    reg.col <- rep(reg.col, length(tau))
  }
  if (length(reg.lwd) == 1) {
    reg.lwd <- rep(reg.lwd, length(tau))
  }
  if (length(reg.lty) == 1) {
    reg.lty <- rep(reg.lty, length(tau))
  }


  colours_seq <- c(
    brewer.pal(n = 8, name = "Dark2"),
    brewer.pal(n = 8, name = "Accent"),
    brewer.pal(n = 12, name = "Paired"),
    brewer.pal(n = 8, name = "Pastel1"),
    brewer.pal(n = 8, name = "Pastel2"),
    brewer.pal(n = 8, name = "Set1"),
    brewer.pal(n = 8, name = "Set2"),
    brewer.pal(n = 8, name = "Set3")
  )

  if (is.null(reg.lty)) {
    reg.lty <- 1:length(tau)
  }
  if (is.null(reg.col)) {
    reg.col <- colours_seq
  }


  if (is.null(x$fit)) {
    if (is.null(newdata)) {
      bb <- which((x$tau) %in% tau)
      plot(x$Y, ...)
      for (i in 1:length(bb)) {
        lines(x$region[[bb[i]]], col = reg.col, lwd = reg.lwd)
      }
    } else {
      plot(newdata[, names(x$Y)], ...)
      bb <- which((x$tau) %in% tau)
      for (i in 1:length(bb)) {
        lines(x$region[[bb[i]]], col = reg.col[i], lwd = reg.lwd[i], lty = reg.lty[i])
      }
    }
  }


  if (!is.null(x$fit) & cond == FALSE) {
    if (is.null(newdata)) {
      bb <- which((x$tau) %in% tau)
      plot(x$Y, ...)
      for (i in 1:length(bb)) {
        lines(x$region[[bb[i]]], col = reg.col[i], lwd = reg.lwd[i], lty = reg.lty[i])
      }
      if (legend == TRUE) legend("bottomright", col = reg.col[1:length(tau)], legend = tau, lwd = 2.5, bty = "n", lty = reg.lty[1:length(tau)])
    } else {
      pred_mu1 <- predict(x$fit$mu1, newdata = newdata)
      pred_mu2 <- predict(x$fit$mu2, newdata = newdata)
      pred_sd1 <- sqrt(exp(predict(x$fit$var1, newdata = newdata)))
      pred_sd2 <- sqrt(exp(predict(x$fit$var2, newdata = newdata)))
      pred_rho <- tanh(predict(x$fit$rho, newdata = newdata))

      newdata$res1 <- (eval(parse(text = paste0("newdata$", all.vars(x$fit$mu1$formula)[1]))) - pred_mu1) / pred_sd1
      newdata$res2 <- (eval(parse(text = paste0("newdata$", all.vars(x$fit$mu2$formula)[1]))) - pred_mu2) / pred_sd2

      YC <- cbind(newdata$res1, newdata$res2)

      for (i in 1:nrow(newdata)) {
        Sigma <- matrix(c(1, pred_rho[i], pred_rho[i], 1), ncol = 2, nrow = 2)
        P <- chol(solve(Sigma))
        YC[i, ] <- P %*% YC[i, ]
      }
      YC <- as.data.frame(YC)
      names(YC) <- c(
        all.vars(x$fit$mu1$formula)[1],
        all.vars(x$fit$mu2$formula)[1]
      )

      bb <- which((x$tau) %in% tau)
      plot(YC, ...)
      for (i in 1:length(bb)) {
        lines(x$region[[bb[i]]], lty = reg.lty[i], col = reg.col[i], lwd = reg.lwd[i])
      }
      if (legend == TRUE) legend("bottomright", col = reg.col[1:length(tau)], legend = tau, lwd = 2.5, bty = "n", lty = reg.lty[1:length(tau)])
    }
    if (axes == TRUE) abline(h = mean(x$Y[, 1]), v = mean(x$Y[, 2]), col = axes.col, lwd = axes.lwd)
  }




  if (cond == TRUE) {
    if (is.null(newdata)) {
      stop("newdata = NULL. This function was developed to represent the bivariate reference region at an specific covariate values. Please supply those values in the newdata argument.")
    }

    reg_pred <- predict(x, cond = TRUE, newdata = newdata, tau = tau)

    if (add == TRUE) {
      colours_seq2 <- c(
        brewer.pal(n = 8, name = "Accent"),
        brewer.pal(n = 8, name = "Set1"),
        brewer.pal(n = 12, name = "Paired")
      )


      for (i in 1:length(reg_pred)) {
        for (j in 1:length(reg_pred[[i]])) {
          lines(reg_pred[[i]][[j]], col = colours_seq2[j], lwd = reg.lwd, lty = reg.lty[i])
        }
      }

      if (legend == TRUE) {
        legend("topleft", legend = names(reg_pred[[1]]), lwd = 2.5, bty = "n", col = colours_seq2[1:length(reg_pred[[1]])], ncol = 2)
        legend("bottomright", legend = tau, lwd = 2.5, bty = "n", lty = reg.lty[1:length(tau)])
      }
    }

    if (add == FALSE) {
      y1 <- x$fit$mu1$model[, 1]
      y2 <- x$fit$mu2$model[, 1]
      plot(y1, y2, ...)
      for (i in 1:length(reg_pred)) {
        for (j in 1:length(reg_pred[[i]])) lines(reg_pred[[i]][[j]], col = colours_seq[j], lwd = reg.lwd, lty = reg.lty[i])
      }
      if (legend == TRUE) {
        legend("topleft", legend = names(reg_pred[[1]]), lwd = 2.5, bty = "n", col = colours_seq[1:length(reg_pred[[1]])], ncol = 2)
        legend("bottomright", legend = tau, lwd = 2.5, bty = "n", lty = reg.lty[1:length(tau)])
      }
    }
  }
}
