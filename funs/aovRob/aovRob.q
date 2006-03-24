##
## Robust ANOVA 
## Author: Jeffrey Wang
## Date  : 08/09/2000
##

aovRob <- function(formula, data=sys.parent(), contrasts=NULL, ...)
{
  neval <- sys.parent()
	m <- model.frame(formula, data)
	m <- m[, -1, drop = F]
	m <- apply(m, 1, paste, collapse = "")
	balanced <- length(unique(table(m))) == 1
  if (missing(data))
    Terms <- terms(formula, "Error") 
  else 
    Terms <- terms(formula, "Error", data = data)
  if (!inherits(formula, "formula"))
    formula <- attr(Terms, "formula")
  if (!is.null(attr(Terms, "specials")$Error))
    warning("Error is not supported in aovRob().")
  lmcall <- call <- match.call()
  lmcall[[1]] <- as.name("lmRob")
  lmcall$formula <- Terms
  result <- eval(lmcall, neval)
  oldClass(result) <- c("aovRob", "lmRob")
  result$call <- call
	result$balanced <- balanced
  result
}

print.aovRob <- function(x, intercept=F, ...)
{
  ss <- as.data.frame(anova.lmRob(x, test="RF"))
  tmp.d <- dim(ss)
  tmp.n <- dimnames(ss)
  ss <- ss[-1,-tmp.d[2]]
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nTerms:\n")
  prmatrix(t(array(	c(format(zapsmall(ss[,2])), format(ss[,1])), 
										c(tmp.d[1]-1, 2), list(tmp.n[[1]][-1],
										rev(tmp.n[[2]][-tmp.d[2]])))), quote = F,
										right = T)
  cat("\nRobust residual scale:", format(x$scale), "\n")
  invisible(x)
}

summary.aovRob <- function(object, ...)
{
  ss <- anova.lmRob(object, test="RF")
  attr(ss, "heading") <- NULL
	if(object$balanced) {
		cf.ss <- dummy.coef.lm(object)

		if(!is.null(cf.ss[["(Intercept)"]])) {
			cf.ss[["(Intercept)"]] <- NULL
			intercept <- 1
		}
		else
			intercept <- 0

		terms <- names(cf.ss)

		term.sets <- list()
		n.levels <- list()
		dtf <- list()

		n <- length(object$residuals)

		for(term in terms) {
			term.sets[[term]] <- unlist(unpaste(term, sep = ":"))
			n.levels[[term]] <- length(cf.ss[[term]])
			dtf[[term]] <- prod(unlist(n.levels[term.sets[[term]]]) - 1)
			cf.ss[[term]] <- sum(cf.ss[[term]]^2)
			cf.ss[[term]] <- (n / n.levels[[term]]) * cf.ss[[term]]
		}

	ss <- data.frame(ss[[1]],c(NA, unlist(cf.ss)),c(NA, unlist(cf.ss))/ss[[1]],ss[[2]],ss[[3]])

	names(ss) <- c("Df", "Sum of Sq", "Mean Sq", "RobustF", "Pr(F)")
	}
  ss[-1, , drop = F]
}

#"plot.aovRob" <- function(x, smooths=F, rugplot=F, id.n=3, 
#                          ask=T, which.plots=NULL, ...)
#{
#  add.ons <- function(x, y, smooths=T, rugplot=T, id.n=3) {
#    if(smooths) {
#      prediction <- loess.smooth(x, y, span = 1, degree = 1)
#      lines(prediction)
#    }
#    if(rugplot) {
#      jx <- jitter(x[!is.na(x)])
#      xlim <- range(jx)
#      rug(jx)
#    }
#    if(id.n) {
#      n <- length(y)
#      oy <- order(abs(y))
#      which <- oy[(n - id.n + 1):n]
#      text(x[which], y[which], as.character(which), adj = 0)
#    }
#  }
#  choices <- c("All", 
#               "Treatment Means",
#               "Response vs Fitted Values", 
#               "Normal QQplot of Residuals")
#  tmenu <- paste("plot:", choices)
#  if (!is.null(which.plots)) {
#    ask <- F
#    which.plots <- c(which.plots+2, 1)
#  }
#  else {
#    if (ask == F) which.plots <- c(3:5, 1)
#  }
#  while (T) {
#    if(ask) {
#      which.plots <- lmRob.checkbox(tmenu, title = 
#                     "\nMake plot selections (or 0 to exit):\n")
#      if (any(which.plots == 1))
#        which.plots <- 2:4
#      which.plots <- 1+which.plots
#    }
#    for (idx in 1:length(which.plots)) {
#      pick <- which.plots[idx]
#      switch(pick,
#             invisible(return(x)),
#### Plot all choices one by one ###
#             {
#               ask.now <- F
#             },
#### Treatment Means ###
#             {
#               plot.design(x)
#             },
#### Response vs Fitted Values ###
#             {
#               Residuals <- x$residuals
#               fits <- x$fitted
#               response <- fits + Residuals
#               form <- formula.lmRob(x)
#               x.lab <- paste("fitted(",deparse(substitute(x)),")",sep="")
#               response.name <- deparse(form[[2]])
#               plot(fits, response, xlab=x.lab, ylab=response.name, ...)
#               abline(0, 1, lty=2)
#             add.ons(fits, response, smooths=smooths, rugplot=rugplot, 
#                     id.n=id.n)
#             },
#### Normal QQplot of Residuals ###
#             {
#               qqxy <- qqnorm(x$residuals)
#               add.ons(qqxy$x, qqxy$y, smooths=F, rugplot=F, id.n=id.n)
#               qqline(x$residuals, lty=2)
#             }
#      )
#    }
#  }
#  invisible(x)
#}

plot.aovRob <- function(x, which.plots = "ask", vertical.outlier = 0.99,
				smooths = F, rugplot = F, id.n = 3, envelope = T, half.normal = F,
				robustQQline = T, mc.samples = 100, level = 0.95, seed = 289, ...)
{
	x.name <- deparse(substitute(x))
	model.list <- list(x$call)
	names(model.list) <- x.name
	x <- list(x = x)
	names(x) <- x.name
	attr(x, "model.list") <- model.list

	plot.aovfm(x, which.plots = which.plots, vertical.outlier =
		vertical.outlier, smooths = smooths, rugplot = rugplot,
		id.n = id.n, envelope = envelope, half.normal = half.normal,
		robustQQline = robustQQline, mc.samples = mc.samples,
		level = level, seed = seed, ...)

	invisible(x[[1]])
}


