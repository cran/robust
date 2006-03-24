lmfmSqrtResVsFittedPlot <- function(x, smooths = FALSE, rugplot = FALSE,
  id.n = 3, ...)
{
  n.models <- length(x)
  mod.names <- names(x)
  n <- length(residuals(x[[1]]))

  panel.special <- function(x, y, smooths = FALSE, rugplot = FALSE, id.n = 3)
  {
    panel.xyplot(x, y, pch = 16, col = 6)
    panel.addons(x, y, smooths = smooths,
      rugplot = rugplot, id.n = id.n)
    invisible()
  }
    
  df <- data.frame(f = as.vector(sapply(x, fitted)),
    r = as.vector(sqrt(abs(sapply(x, residuals)))),
    mod = rep(mod.names, each = n))
      
  print(xyplot(r ~ f | mod,
    data = df,
    xlab = "Fitted Values",
    ylab = "Sqrt(abs(Residuals))",
    main = "sqrt(abs(Residuals)) vs Fitted Values",
    panel = panel.special,
    smooths = smooths,
    rugplot = rugplot,
    id.n = id.n,
    strip = function(...) strip.default(..., style = 1),
    layout = c(n.models, 1, 1),
    ...))

  invisible(x)
}


