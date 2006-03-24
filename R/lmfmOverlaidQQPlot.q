lmfmOverlaidQQPlot <- function(x, ...)
{
  n.models <- length(x)
  mod.names <- names(x)

	res <- na.omit(sapply(x, residuals))
  n <- length(res)
	px <- py <- matrix(0, n, n.models)

	for(i in 1:n.models) {
		tmp <- qqnorm(res[, i], plot.it = FALSE)
		px[, i] <- tmp$x
		py[, i] <- tmp$y
	}

  matplot(px, py,
    pch = 1:n.models,
    col = 1:n.models,
    xlab = "Quantiles of Standard Normal",
    ylab = "Residuals",
    main = "Normal QQ-Plot of Residuals",
    ...)

  key(min(px), max(py),
    text = list(mod.names),
    lines = list(type = "p", col = 1:n.models, pch = 1:n.models),
    transparent = TRUE)

  invisible(x)
}


