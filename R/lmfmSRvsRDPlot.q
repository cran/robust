lmfmSRvsRDPlot <- function(x, level = 0.95, id.n = 3, ...)
{
  n.models <- length(x)
  mod.names <- names(x)
  n <- length(residuals(x[[1]]))
  x.sum <- summary(x)

  std.resids <- matrix(0.0, n, n.models)
  for(i in 1:n.models) {
    if(is.null(x.sum[[i]]$sigma) || is.na(x.sum[[i]]$sigma))
      std.resids[, i] <- residuals(x[[i]])
    else
      std.resids[, i] <- residuals(x[[i]]) / x.sum[[i]]$sigma
  }

  model <- sapply(x, function(u) !is.null(u$model))
  if(!any(model))
    stop("none of the fitted models in ", sQuote(deparse(substitute(x))),
          "contain a model frame component")
  model <- x[[(1:n.models)[model][1]]]$model

  model.terms <- attributes(model)$terms
  term.labels <- attr(model.terms, "term.labels")
  dataClasses <- attr(model.terms, "dataClasses")
  numeric.vars <- names(dataClasses == "numeric")
  term.labels <- intersect(term.labels, numeric.vars)

  if(length(term.labels)) {

    model <- model[term.labels]
    p <- dim(model)[2]
    dist <- covRob(model, distance = TRUE)$dist

    res.thresh <- qnorm(level)
    dist.thresh <- qchisq(level, df = p)

    y.range <- range(std.resids)
    y.range[1] <- 1.05 * min(y.range[1], -res.thresh)
    y.range[2] <- 1.05 * max(y.range[2], res.thresh)

    x.range <- c(0.0, max(dist))
    x.range[2] <- 1.05 * max(x.range[2], res.thresh)

    panel.special <- function(x, y, res.thresh = 1.0, dist.thresh = 1.0, id.n = 3)
    {
      panel.xyplot(x, y, col = 6, pch = 16)
      panel.addons(x, y, smooths = FALSE, rugplot = FALSE, id.n = id.n)
      panel.abline(v = dist.thresh, lty = 2)
      panel.abline(h = res.thresh, lty = 2)
      panel.abline(h = -res.thresh, lty = 2)
      invisible()
    }

    df <- data.frame(RD = rep(dist, n.models),
      RR = as.vector(std.resids),
      mod = rep(mod.names, each = n))

    print(xyplot(RR ~ RD | mod,
      data = df,
      xlab = "Robust Distances",
      panel = panel.special,
      xlim = x.range,
      ylim = y.range,
      ylab = "Standardized Residuals",
      main = "Standardized Residuals vs Robust Distances",
      strip = function(...) strip.default(..., style = 1),
      res.thresh = res.thresh,
      dist.thresh = dist.thresh,
      id.n = id.n,
      layout = c(n.models, 1, 1),
      ...))
    }

  else
    warning("robust distances could not be computed because there are no numeric variables in the model frame")

  invisible(x)
}


