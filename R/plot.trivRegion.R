#' Default trivRegion plotting
#'
#' This function allow to depict in an interactive rgl plot the estimated trivariate reference/probabilistic
#' region. If cond=FALSE it showes trivariate standarized residuals, while with
#' cond=TRUE it represents the region shape for any covariate(s) value.
#'
#' @param x A trivRegion object.
#' @param cond A logical argument, if TRUE a conditional reference region is depicted.
#' @param newdata If cond==TRUE, a data frame containing covariate
#' values for which the reference/probabilistic region will be depicted.
#' @param add A logical argument, if TRUE the conditional reference region is depicted over a pre existing plot.
#' @param legend A logical argument, if TRUE a legend is given along with the reference region.
#' @param reg.col Region contour colour.
#' @param incol Colours for the points included inside the reference region.
#' @param planes Logical; if TRUE, planes are added indicating (x=0,y=0,z=0).
#' @param ... Further rgl plot parameters.
#' @return This function return an interactive rgl plot of a bivRegion object.
#' @importFrom "misc3d" "computeContour3d" "makeTriangles" "drawScene.rgl"
#' @importFrom "rgl" "plot3d" "points3d" "legend3d" "planes3d"
#' @importFrom "grDevices" "adjustcolor"
#' @export

plot.trivRegion <- function(x, cond = FALSE, planes = FALSE, newdata = NULL, add = FALSE, reg.col = NULL, incol = "grey", legend = FALSE, ...) {
  fhat <- x$kernel_fit
  standarized_region <- computeContour3d(fhat$estimate,
    maxvol = max(fhat$estimate),
    level = x$k_limit,
    x = fhat$eval.points[[1]],
    y = fhat$eval.points[[2]],
    z = fhat$eval.points[[3]],
    mask = NULL
  )


  if (cond == FALSE) {
    plot3d(x$trivres[x$which_out, ], ...)
    points3d(x$trivres, col = incol)

    if (is.null(reg.col)) reg.col <- "orange"
    drawScene.rgl(makeTriangles(standarized_region, color = reg.col, alpha = 0.2), add = TRUE)

    if (planes == TRUE) {
      planes3d(10, 0, 0, col = "azure")
      planes3d(0, 10, 0, col = "azure")
      planes3d(0, 0, 10, col = "azure")
    }
  }

  if (cond == TRUE) {
    fhat <- x$kernel_fit
    standarized_x <- computeContour3d(fhat$estimate,
      maxvol = max(fhat$estimate),
      level = x$k_limit,
      x = fhat$eval.points[[1]],
      y = fhat$eval.points[[2]],
      z = fhat$eval.points[[3]],
      mask = NULL
    )
    means <- cbind(
      predict(x$fit$mean1, newdata = newdata),
      predict(x$fit$mean2, newdata = newdata),
      predict(x$fit$mean3, newdata = newdata)
    )

    sds <- cbind(
      sqrt(exp(predict(x$fit$var1$fit, newdata = newdata))),
      sqrt(exp(predict(x$fit$var2$fit, newdata = newdata))),
      sqrt(exp(predict(x$fit$var3$fit, newdata = newdata)))
    )

    rhos <- cbind(
      tanh(predict(x$fit$rho1$fit, newdata = newdata)),
      tanh(predict(x$fit$rho2$fit, newdata = newdata)),
      tanh(predict(x$fit$rho3$fit, newdata = newdata))
    )



    for (i in 1:nrow(newdata)) {
      Sigma <- matrix(c(1, rhos[i, 1], rhos[i, 2], rhos[i, 1], 1, rhos[i, 3], rhos[i, 2], rhos[i, 3], 1), ncol = 3, nrow = 3, byrow = TRUE)
      P <- chol(solve(Sigma))
      P <- solve(P)

      pol_pred <- cbind(
        means[i, 1] + (t(apply(standarized_region, 1, function(x) P %*% x)))[, 1] * sds[i, 1],
        means[i, 2] + (t(apply(standarized_region, 1, function(x) P %*% x)))[, 2] * sds[i, 2],
        means[i, 3] + (t(apply(standarized_region, 1, function(x) P %*% x)))[, 3] * sds[i, 3]
      )

      colours_seq <- c(
        brewer.pal(n = 8, name = "Dark2"),
        brewer.pal(n = 8, name = "Accent"),
        brewer.pal(n = 8, name = "Set2"),
        brewer.pal(n = 8, name = "Set3")
      )

      if (is.null(reg.col)) {
        reg.col <- colours_seq
        reg.col_aux <- colours_seq[i]
      }
      if (length(reg.col) >= 1) {
        reg.col_aux <- reg.col[i]
      }

      if (add == FALSE) {
        if (i == 1) plot3d(pol_pred, col = reg.col_aux, size = 1, ...) else plot3d(pol_pred, col = reg.col_aux, add = TRUE, size = 1)
        drawScene.rgl(makeTriangles(pol_pred, color = reg.col_aux, alpha = 0.2), add = TRUE)
      }
      if (add == TRUE) {
        drawScene.rgl(makeTriangles(pol_pred, color = reg.col_aux, alpha = 0.2), add = TRUE)
      }
    }

    if (legend == TRUE) {
      vars <- names(newdata)
      nomes <- as.numeric()
      for (k in 1:nrow(newdata)) {
        nomes[k] <-
          paste(vars, "=", newdata[k, vars],
            collapse = " & "
          )
      }

      legend3d("bottomright",
        legend = nomes,
        col = adjustcolor(reg.col, alpha.f = 0.3), lwd = 10, cex = 2,
        bty = "n"
      )
    }
  }
}
