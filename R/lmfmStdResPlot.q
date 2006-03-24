lmfmStdResPlot <- function(x, level = 0.95, id.n = 3, ...)
{
  n.models <- length(x)
  mod.names <- names(x)
  n <- length(residuals(x[[1]]))
  x.sum <- summary(x)
  threshold <- qnorm(level)

  std.resids <- matrix(0.0, n, n.models)
  for(i in 1:n.models) {
    if(is.null(x.sum[[i]]$sigma) || is.na(x.sum[[i]]$sigma))
      std.resids[, i] <- residuals(x[[i]])
    else
      std.resids[, i] <- residuals(x[[i]]) / x.sum[[i]]$sigma
  }

  y.range <- range(std.resids)
  y.range[1] <- 1.05 * min(y.range[1], -threshold)
  y.range[2] <- 1.05 * max(y.range[2], threshold)

  panel.special <- function(x, y, threshold = 1.645, id.n = 3)
  {
    n <- length(y)
    type <- ifelse(n > 40, "l", "b")
    panel.xyplot(x, y, type = type, col = 6, pch = 16)
    outliers <- which(abs(y) > threshold)
    if(length(outliers) > id.n)
      outliers <- order(abs(y))[(n - id.n + 1):n]
    if(id.n > 0 && any(outliers))
      panel.text(x[outliers], y[outliers], paste(" ", outliers, sep = ""), adj = 0)
    panel.abline(h = threshold, lty = 2)
    panel.abline(h = -threshold, lty = 2)
    invisible()
  }

  df <- data.frame(indicies = rep(1:n, n.models),
    std.resid = as.vector(std.resids),
    mod = rep(mod.names, each = n))

  print(xyplot(std.resid ~ indicies | mod,
    data = df,
    xlab = "Index (Time)",
    panel = panel.special,
    id.n = id.n,
    ylim = y.range,
    ylab = "Standardized Residuals",
    main = "Standardized Residuals vs. Index (Time)",
    strip = function(...) strip.default(..., style = 1),
    threshold = threshold,
    layout = c(n.models, 1, 1),
    ...))

  invisible(x)
}


