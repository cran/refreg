#' Default summary_boot plotting
#'
#' This function takes the bivRegr bootstrap replicates obtained with summary_boot
#' function, and plots the parametric, and smooth effects for each model.
#'
#' @param x A summary_boot object.
#' @param eq  A number indicating the model effects to be depicted; 1 = first response mean,
#' 2 =  second response mean, 3 = first response variance, 4 = second response variance, and
#' 5 = correlation model.
#' @param select  An optional parameter to represent an specific effect for each equation.
#' @param ... Additional plot parameters, not yet implemented.
#' @return This function returns a ggplot2 plot for the estimated effects along with
#' bootstrap 95\% confidence intervals.
#' @importFrom "stringr" "str_detect"
#' @export


plot.summary_boot <- function(x, eq = 1, select = NULL, ...) {
  results_boot <- list()
  if (eq == 1) model <- x$fit$mu1
  if (eq == 2) model <- x$fit$mu2
  if (eq == 3) model <- x$fit$var1
  if (eq == 4) model <- x$fit$var2
  if (eq == 5) model <- x$fit$rho

  model_pred <- predict(model, type = "terms")


  for (i in 1:ncol(model_pred)) {
    if (eq == 1) boot_res <- x$boot_mu1
    if (eq == 2) boot_res <- x$boot_mu2
    if (eq == 3) boot_res <- x$boot_var1
    if (eq == 4) boot_res <- x$boot_var2
    if (eq == 5) boot_res <- x$boot_rho

    fit_boot <- do.call(cbind, lapply(boot_res, function(x) x[, colnames(model_pred)[i], drop = FALSE]))

    if (str_detect(colnames(model_pred)[i], "s") == TRUE) {
      covar_data <- x$fit$data[, which(str_detect(colnames(model_pred)[i], names(x$fit$data))), drop = FALSE]
    } else {
      covar_data <- x$fit$data[, colnames(model_pred)[i], drop = FALSE]
    }

    fit_pred <- cbind(
      covar_data,
      model_pred[, colnames(model_pred)[i], drop = FALSE],
      t(apply(fit_boot, 1, quantile, probs = c(0.05, 0.95)))
    )
    fit_pred <- unique(fit_pred)
    fit_pred <- as.data.frame(fit_pred)
    fit_pred <- fit_pred[order(fit_pred[, 1]), ]

    colnames(fit_pred) <- c(colnames(model_pred)[i], "fit", "q2.5", "q97.5")
    results_boot[[i]] <- fit_pred
  }
  names(results_boot) <- colnames(model_pred)


  myplot <- function(data) {
    if (is(data[, 1],"factor")) {
      g1 <- ggplot(data, aes(x = data[, 1], y = data[, 2])) +
        geom_point(pch = 13, cex = 5, col = "#E41A1C") +
        theme_classic()
      g1 <- g1 + geom_errorbar(aes(ymin = data[, 3], ymax = data[, 4]), col = "#377EB8", size = 0.75, width = 0.25)
      g1 <- g1 + xlab(names(data)[1]) + ylab(str_to_title(paste0(names(data)[1], " centered effect")))

      if (eq == 1) {
        g1 <- g1 + ggtitle(paste0(names(model$model)[1], " mean"))
      }
      if (eq == 2) {
        g1 <- g1 + ggtitle(paste0(names(model$model)[1], " mean"))
      }
      if (eq == 3) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu1$model)[1], " variance"))
      }
      if (eq == 4) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu2$model)[1], " variance"))
      }
      if (eq == 5) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu1$model)[1], "-", names(x$fit$mu2$model)[1], " correlation"))
      }
      g1
    } else {
      g1 <- ggplot(data, aes(x = data[, 1], y = data[, 2]))
      g1 <- g1 + geom_line() + theme_classic() + geom_rug(aes(x = data[, 1], y = NULL), col = "grey")
      g1 <- g1 + geom_ribbon(data = data, aes(ymin = data[, 3], ymax = data[, 4]), alpha = 0.3)

      if (names(data)[1] %in% do.call(c, lapply(model$smooth, function(x) x$label)) == TRUE) {
        aa <- which(do.call(c, lapply(model$smooth, function(x) x$label)) %in% names(data[1]))
        first <- model$smooth[[aa]]$first.para
        last <- model$smooth[[aa]]$last.para
        edf <- format(round(sum(model$edf[first:last]), 2), nsmall = 2)
      }

      g1 <- g1 + xlab(model$smooth[[aa]]$term) + ylab(paste0("s(", model$smooth[[aa]]$term, ",", edf, ")"))
      if (eq == 1) {
        g1 <- g1 + ggtitle(paste0(names(model$model)[1], " mean"))
      }
      if (eq == 2) {
        g1 <- g1 + ggtitle(paste0(names(model$model)[1], " mean"))
      }
      if (eq == 3) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu1$model)[1], " variance"))
      }
      if (eq == 4) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu2$model)[1], " variance"))
      }
      if (eq == 5) {
        g1 <- g1 + ggtitle(paste0(names(x$fit$mu1$model)[1], "-", names(x$fit$mu2$model)[1], " correlation"))
      }
      g1
    }
  }

  plots <- lapply(results_boot, myplot)

  if (is.null(select)) {
    if (length(plots) < 4) {
      return(grid.arrange(grobs = plots, ncol = length(plots)))
    }
    if (length(plots) >= 4) {
      return(grid.arrange(grobs = plots, ncol = length(plots) / 2, nrow = length(plots) / 2))
    }
  } else {
    return(plots[[select]])
  }
}
