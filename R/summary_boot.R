#' bivRegr summary function
#'
#' This function perform a bootstrap procedure in order to obtain
#' confidence intervals for the bivRegr estimated effects. The function allow
#' to paralelize the bootstrap resampling scheme using doParalell, and
#' foreach libraries.
#'
#' @param object A bivRegr fit.
#' @param B A number indicating the bootstrap iterations.
#' @param parallel A logical indicating if bootstrap parallelization must be
#' applied.
#' @param cores If parallel = TRUE, a number indicating computer cores
#' to be used during paralellization. If NULL the function use all cores available
#' but one.
#' @return This function returns the bootstrap replicates of the bivRegr
#' sub models. Results might be checked applying plot.summary_boot().
#' @importFrom "utils" "setTxtProgressBar" "txtProgressBar"
#' @importFrom "parallel" "detectCores" "makeCluster" "stopCluster"
#' @importFrom "doParallel" "registerDoParallel"
#' @importFrom "foreach" "foreach" "%dopar%"
#' @export

summary_boot <- function(object, B = 100, parallel = FALSE, cores = NULL) {

  # Pillo as covariables
  covar <- unique(c(
    names(object$mu1$model)[-1], names(object$mu2$model)[-1],
    names(object$var1$model)[-1], names(object$var2$model)[-1],
    names(object$rho$model)[-1]
  ))
  covar <- covar[which(covar != "(weights)")]


  if (parallel == FALSE) {
    boot_mu1 <- boot_mu2 <- boot_var1 <- boot_var2 <- boot_rho <- list()
    pb <- txtProgressBar(min = 0, max = B, style = 3)

    for (k in 1:B) {
      Sys.sleep(0.1)

      erros <- as.matrix(object$bivres[sample(1:nrow(object$bivres), nrow(object$bivres), replace = TRUE), ])

      # Xero unha base de datos por bootstrap
      dat_star <- matrix(0, ncol = 2, nrow = nrow(object$data))
      for (i in 1:nrow(object$data)) {
        Sigma <- matrix(c(1, object$correlation[i], object$correlation[i], 1), ncol = 2, nrow = 2)
        P <- chol(solve(Sigma))
        P <- solve(P)

        dat_star[i, ] <- P %*% erros[i, ]
        dat_star[i, 1] <- object$means[i, 1] + dat_star[i, 1] * object$std.variations[i, 1]
        dat_star[i, 2] <- object$means[i, 2] + dat_star[i, 2] * object$std.variations[i, 2]
      }

      dat_star <- as.data.frame(cbind(dat_star, object$data[, covar]))
      names(dat_star) <- c(names(object$mu1$model)[1], names(object$mu2$model)[1], covar)

      fit_star <- bivRegr(f = object$formula, data = dat_star)

      boot_mu1[[k]] <- predict(fit_star$mu1, newdata = object$data, type = "terms")
      boot_mu2[[k]] <- predict(fit_star$mu2, newdata = object$data, type = "terms")
      boot_var1[[k]] <- predict(fit_star$var1, newdata = object$data, type = "terms")
      boot_var2[[k]] <- predict(fit_star$var2, newdata = object$data, type = "terms")
      boot_rho[[k]] <- predict(fit_star$rho, newdata = object$data, type = "terms")
      setTxtProgressBar(pb, k)
    }
    close(pb)
  }


  if (parallel == TRUE) {
    if (is.null(cores)) cores <- detectCores() - 1
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    resampling <- foreach(i = 1:B) %dopar% {
      erros <- as.matrix(object$bivres[sample(1:nrow(object$bivres), nrow(object$bivres), replace = TRUE), ])

      # Xero unha base de datos por bootstrap
      dat_star <- matrix(0, ncol = 2, nrow = nrow(object$data))
      for (i in 1:nrow(object$data)) {
        Sigma <- matrix(c(1, object$correlation[i], object$correlation[i], 1), ncol = 2, nrow = 2)
        P <- chol(solve(Sigma))
        P <- solve(P)

        dat_star[i, ] <- P %*% erros[i, ]
        dat_star[i, 1] <- object$means[i, 1] + dat_star[i, 1] * object$std.variations[i, 1]
        dat_star[i, 2] <- object$means[i, 2] + dat_star[i, 2] * object$std.variations[i, 2]
      }

      dat_star <- as.data.frame(cbind(dat_star, object$data[, covar]))
      names(dat_star) <- c(names(object$mu1$model)[1], names(object$mu2$model)[1], covar)

      fit_star <- bivRegr(f = object$formula, data = dat_star)

      boot_mu1 <- predict(fit_star$mu1, newdata = object$data, type = "terms")
      boot_mu2 <- predict(fit_star$mu2, newdata = object$data, type = "terms")
      boot_var1 <- predict(fit_star$var1, newdata = object$data, type = "terms")
      boot_var2 <- predict(fit_star$var2, newdata = object$data, type = "terms")
      boot_rho <- predict(fit_star$rho, newdata = object$data, type = "terms")

      list(boot_mu1, boot_mu2, boot_var1, boot_var2, boot_rho)
    }
    stopCluster(cl)
    boot_mu1 <- lapply(resampling, function(x) x[[1]])
    boot_mu2 <- lapply(resampling, function(x) x[[2]])
    boot_var1 <- lapply(resampling, function(x) x[[3]])
    boot_var2 <- lapply(resampling, function(x) x[[4]])
    boot_rho <- lapply(resampling, function(x) x[[5]])
  }

  l <- list(
    boot_mu1 = boot_mu1, boot_mu2 = boot_mu2,
    boot_var1 = boot_var1, boot_var2 = boot_var2, boot_rho = boot_rho, fit = object
  )
  class(l) <- "summary_boot"
  return(l)
}
