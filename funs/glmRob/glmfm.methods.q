print.glmfm <- function(x, ...)
{
	n.models <- length(x)
	mod.names <- names(x)

	cat("\nCalls: \n")
	for(i in 1:n.models) {
		cat(mod.names[i], ": ")
		print(x[[i]]$call)
	}

	mod.coefs <- lapply(x, coef)
	coef.names <- unique(unlist(lapply(mod.coefs, names)))
	n.coefs <- length(coef.names)
	tmp <- matrix(NA, n.coefs, n.models)
	dimnames(tmp) <- list(coef.names, mod.names)
	for(i in 1:n.models)
		tmp[match(names(mod.coefs[[i]]), coef.names), i] <- mod.coefs[[i]]

	cat("\nCoefficients:\n")
	tmp.idx <- is.na(tmp)
	tmp <- format(round(tmp, options()$digits - 3))
	tmp[tmp.idx] <- ""
	print(tmp, quote = F, ...)
	cat("\n")

	cat("Residual Deviance Estimates:\n")
	for(i in 1:n.models)
	  cat(mod.names[i], ":", format(signif(x[[i]]$deviance)),
			"on", x[[i]]$df.resid, "degrees of freedom\n")

	invisible(x)
}


summary.glmfm <- function(object, correlation = T, ...)
{
	append.p.values <- function(obj)
	{
		cbind(obj, 2 * (1 - pnorm(abs(obj[, ncol(obj)]))))
	}

	n.models <- length(object)
	model.list <- attr(object, "model.list")

	ans <- list()
	ans$mod.names <- mod.names <- format(names(model.list))
	ans$calls <- lapply(object, function(x) x$call)

	ans$restats <- t(sapply(object, function(x)
		quantile(residuals(x), na.rm = T)))

	dimnames(ans$restats) <- list(paste(ans$mod.names, ":"),
		c("Min", "1Q", "Median", "3Q", "Max"))

	mod.coefs <- lapply(object, coef)
	coef.names <- unique(unlist(lapply(mod.coefs, names)))
	n.coefs <- length(coef.names)

	ans$coefs <- array(NA, dim = c(n.coefs, 4, n.models))
	ans$devs <- ans$df <- rep(NA, n.models)

	if(correlation)
		ans$correlations <- array(NA, dim = c(n.coefs, n.coefs, n.models))

	for(i in 1:n.models) {
		ans$devs[i] <- object[[i]]$deviance
		tmp <- summary(object[[i]], correlation = correlation)
			coef.idx <- match(names(mod.coefs[[i]]), coef.names)

	##  add a column with p-values for the glm / glmRobEst summary  ##

		if(ncol(tmp$coefficients) == 3) {
			old.dimnames <- dimnames(tmp$coefficients)
			tmp$coefficients <- append.p.values(tmp$coefficients)
			old.dimnames[[2]] <- c(old.dimnames[[2]], "Pr(>|t|)")
			dimnames(tmp$coefficients) <- old.dimnames
		}

		ans$coefs[coef.idx, , i] <- tmp$coefficients

		if(correlation)
			ans$correlations[coef.idx, coef.idx,i] <- tmp$correlation

		ans$sigmas[i] <- tmp$sigma
		ans$r.squared[i] <- tmp$r.squared
		ans$df[i] <- tmp$df[2]
	}

	dimnames(ans$coefs) <- list(coef.names, dimnames(tmp$coefficients)[[2]], NULL)
	if(correlation)
			dimnames(ans$correlations) <- list(coef.names, coef.names, NULL)
	
	oldClass(ans) <- "summary.glmfm"
	ans
}



print.summary.glmfm <- function(x, ...)
{
	n.models <- length(x$mod.names)
	digits <- options()$digits - 3

	cat("\nCalls: \n")
	for(i in 1:n.models) {
		cat(x$mod.names[i], ": ")
		print(x$calls[[i]])
	}

	cat("\nResidual Statistics:\n")
	print(format(round(x$restats, digits = digits)), quote = F, ...)
	
	cat("\nCoefficients:\n")
	np <- dim(x$coefs)[1]
	tmp <- numeric()
	tmp.names <- character()

	if(np == 1) {
		for(j in 1:n.models) {
			tmp <- rbind(tmp, x$coefs[,  , j])
			tmp.names <- c(tmp.names, dimnames(x$coefs[,  , j, drop = F])[[1]])
		}
	}

	else {
		for(i in 1:np) {
			for(j in 1:n.models) {
			  tmp <- rbind(tmp, (x$coefs[,  , j])[i,  ])
			  tmp.names <- c(tmp.names, dimnames(x$coefs[,  , j])[[1]][i])
			}
		}
	}

	tmp.names <- paste(rep(paste(x$mod.names, ":"), np), format(tmp.names))
	dimnames(tmp)[[1]] <- tmp.names
	tmp.idx <- is.na(tmp)
	tmp <- format(round(tmp, digits = digits))
	tmp[tmp.idx] <- ""
	print(tmp, quote = F, ...)

	cat("\nResidual Deviance of model(s):\n")
	for(i in 1:n.models)
		cat(x$mod.names[i], ":", format(signif(x$devs[i], digits)), "\n")

	if(!is.null(x$correlations)) {
		n.coefs <- dim(x$correlations)[1]
		if(n.coefs > 1) {
			cat("\nCorrelations:\n")
			for(i in 1:n.models) {
				cat(paste(x$mod.names[i], ":", sep = ""), "\n")
				tmp <- x$correlations[,  , i]
				tmp.idx <- is.na(tmp)
				tmp <- format(round(tmp, digits = digits))
				tmp[col(tmp) > row(tmp)] <- tmp[tmp.idx] <- ""
				print(tmp, quote = F, ...)
				cat("\n")
			}
		}
	}

	invisible(x)
}


plot.glmfm <- function(x, which.plots = "ask", type = "pearson",
	chisq.percent = 0.99, vertical.outlier = .99, smooths = F,
	rugplot = F, id.n = 3,	...)
{

	model.list <- attr(x, "model.list")
	n.models <- length(x)
	xsum <- summary.glmfm(x, correlation = F)
	if(!length(xsum$mod.names))
		xsum$mod.names <- paste("Model", 1:n.models)
	mod.names <- xsum$mod.names

	##
	##	Define colors, line styles and bandwidth fnct.
	##

	if(n.models < 5)
		colors <- c(1,6,4,8)[1:n.models]
	else colors <- 1:n.models

	if(n.models < 5)
		styles <- c(1,4,6,8)[1:n.models]
	else styles <- 1:n.models

	bandwidth.nrd <- function(x) {
		r <- quantile(x, c(0.25, 0.75), na.rm = T)
		h <- (r[2] - r[1])/1.34
		4 * 1.06 * min(sqrt(var(x, na.method = "omit")), h) * length(x)^{-1/5}
	}

	vertical.outlier <- qnorm(vertical.outlier)
	r <- na.omit(sapply(x, residuals, type = "deviance"))
	p.r <- na.omit(sapply(x, residuals, type = "pearson"))
	f <- na.omit(sapply(x, fitted))
	n <- nrow(f)
	pr <- na.omit(sapply(x, predict))
	idx <- c(1, 3:(n.models + 1))

	##
	##  make denx/deny and px/py on the same plotting scale
	##

	denx <- deny <- matrix(0, 100, n.models)
	ppx <- ppy <- px <- py <- matrix(0, n, n.models)
	for(i in 1:n.models) {
		b <- bandwidth.nrd(r[, i])
		den <- density(r[, i], width = b, n = 100, na.rm = T)
		denx[, i] <- den$x
		deny[, i] <- den$y
		tmp <- qqnorm(r[, i], plot = F)
		px[, i] <- tmp$x
		py[, i] <- tmp$y
		tmp <- qqnorm(p.r[, i], plot = F)
		ppx[, i] <- tmp$x
		ppy[, i] <- tmp$y
	}
	den.range <- c(min(py, denx), max(py, denx))

	##
	##  menu choices
	##

	choices <- c("All",
		"Deviances vs Fitted Values", 
		"Response vs Fitted Values", 
		"QQ-Plot of Pearson Residuals",
		"Deviances QQ-Plot",
		"Standardized Deviances vs Robust Distances", 
		"Standardized Deviances vs Index (Time)",
		"Sqrt of abs(Deviances) vs Fitted Values")
	tmenu <- paste("plot:", choices)

	ask <- T

	if(is.integer(which.plots)) {
		ask <- F
		which.plots <- c(which.plots + 1, 1)
	}

	else if(which.plots == "all") {
		ask <- F
		which.plots <- c(3:9, 1)
	}

	while(T) {
		if(ask) {
			which.plots <- lmRob.checkbox(tmenu, title = 
				"\nMake plot selections (or 0 to exit):\n")
			if(any(which.plots == 1))
				which.plots <- 2:8
			which.plots <- 1 + which.plots
		}

		graph.number <- 1
		if(dev.cur() == 1 && which.plots[1] != 1)
			trellis.device()
		for(iwhich in 1:length(which.plots)) {
			pick <- which.plots[iwhich]
			switch(pick,
				invisible(return(x)),
				{
					ask.now <- F
				}
				,
				{
					panel.special <- function(x, y, smooths = smooths,
															rugplot = rugplot, id.n = id.n)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths,
												rugplot = rugplot, id.n = id.n)
						abline(h = 0, lty = 2)
						invisible()
					}

					df <- data.frame(	f = as.vector(f),
														r = as.vector(r),
														mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(r ~ f | mod,
						data = df,
						xlab = "Fitted Values",
						ylab = "Deviances",
						main = "Deviances vs. Fitted Values",
						panel = panel.special,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						strip = function(...)
							strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Deviances vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					y <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						m <- x[[i]]$call
						m <- call("model.frame", formula = m$formula, data = m$data,
											na.action = m$na.action)
						m <- eval(m, sys.parent())

						y[,i] <- model.extract(m, response)
					}

					panel.special <- function(x, y, smooths = smooths,
															rugplot = rugplot, id.n = id.n)
					{
						panel.xyplot(x, y, col = 6, pch = 16)
						panel.addons(x, y, smooths = smooths,
										rugplot = rugplot, id.n = id.n)
						abline(0, 1, lty = 2)
						invisible()
					}

					df <- data.frame(	y = as.vector(y),
														f = as.vector(f),
														mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(y ~ f | mod,
					data = df,
						xlab = "Fitted Values",
						ylab = "Response",
						main = "Response vs. Fitted Values",
						panel = panel.special,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						strip = function(...)
							strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Response vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					df <- data.frame(	py = as.vector(ppy),
														px = as.vector(ppx),
														mod = rep(mod.names, rep(n, n.models)))

					panel.special <- function(x, y, id.n)
						{
							panel.xyplot(x, y, pch = 16, col = 6)
							panel.addons(x, y, smooths = F, rugplot = F, id.n = id.n)
							invisible()
						}

					print(xyplot(py ~ px | mod,
						data = df,
						ylab = "Pearson Residuals",
						xlab = "Quantiles of Standard Normal",
						main = "QQ-Plot of Pearson Residuals",
						panel = panel.special,
						id.n = id.n,
						strip = function(...)
							strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "QQ-Plot of Pearson Residuals")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					fam <- family(x[[1]])$family["name"]
					if(fam == "Binomial")
						dist <- 0
					else dist <- 1
					qq.b <- qq.a <- y <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						m <- x[[i]]$call
						m <- call("model.frame", formula = m$formula, data = m$data)
						m <- eval(m, sys.parent())
						y[,i] <- model.extract(m, response)
						tmp <- qqplot.glmRob(y[,i], f[,i], dist)
						qq.a[order(r[, i]), i] <- tmp$quantiles
					}

					df <- data.frame(	qq.a = as.vector(qq.a),
														qq.b = as.vector(r),
														mod = rep(mod.names, rep(n, n.models)))

					panel.special <- function(x, y, id.n = id.n) {
						panel.xyplot(x, y, col = 6, pch = 16)
						panel.addons(x, y, smooths = F, rugplot = F, id.n = id.n)
						invisible()
					}

					print(xyplot(qq.a ~ qq.b | mod,
						data = df,
						ylab = "Deviances",
						xlab = "Estimated Quantiles",
						main = "Deviances QQ-Plot",
						panel = panel.special,
						id.n = id.n,
						strip = function(...)
							strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Deviances QQ-Plot")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					RR <- RD <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						RR[, i] <- r[, i]
						tmp.mod <- numeric.model.matrix(x[[i]])

						if(is.null(ncol(tmp.mod)))
							stop("Robust distances can not be computed: all regressors are factors.")
						RD[, i] <- sqrt(covRob(tmp.mod)$dist)
					}

					x.range <- range(RD)
					x.chisq <- sqrt(qchisq(chisq.percent, df = ncol(tmp.mod)))
					if(x.range[2] < x.chisq)
						x.range[2] <- x.chisq
					y.range <- range(RR)
					if(y.range[1] > -vertical.outlier)
						y.range[1] <-  -vertical.outlier
					if(y.range[2] < vertical.outlier)
						y.range[2] <- vertical.outlier

					panel.special <- function(x, y, vo = 2.5, id.n = 3, chi = 0.975) {

						panel.xyplot(x, y, col = 6, pch = 16)
						n <- length(y)
						y.index <- (1:n)[abs(y) > vo]

						if(length(y.index) > id.n)
							y.index <- order(abs(y))[(n - id.n + 1):n]

						if(any(y.index)) {
							labels <- paste("  ", as.character((1:n)[y.index]), sep = "")
							text(x[y.index], y[y.index], labels, adj = 0)
						}

						abline(v = chi, lty = 2)
						abline(h = vo, lty = 2)
						abline(h = -vo, lty = 2)
						invisible()
					}

					df <- data.frame(	RD = as.vector(RD),
														RR = as.vector(RR),
														mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(RR ~ RD | mod,
						data = df,
						xlab = "Robust Distances",
						panel = panel.special,
						id.n = id.n,
						xlim = x.range,
						ylim = y.range,
						ylab = "Standardized Deviances",
						main = "Standardized Deviances vs. Robust Distances",
						strip = function(...)
							strip.default(..., style = 1),
						vo = vertical.outlier,
						chi = x.chisq,
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Deviances vs. Robust Distances")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					SR <- IDX <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						if(is.null(xsum$sigmas[i]) || is.na(xsum$sigmas[i]))
							SR[, i] <- r[, i]
						else
							SR[, i] <- r[, i] / xsum$sigmas[i]
						IDX[, i] <- 1:n
					}
					y.range <- range(SR)
					if(y.range[1] > -vertical.outlier)
						y.range[1] <- -vertical.outlier
					if(y.range[2] < vertical.outlier)
						y.range[2] <- vertical.outlier

					panel.special <- function(x, y, vo = 2.5, chi = 0.975, id.n = 3)
					{
						if(length(y) > 40) type <- "l"
						else type <- "b"
						panel.xyplot(x, y, type = type, col = 6, pch = 16)
						n <- length(y)
						y.index <- (1:n)[abs(y) > vo]
						if(length(y.index) > id.n)
							y.index <- order(abs(y))[(n - id.n + 1):n]
						if(any(y.index)) {
							labels <- paste(" ", as.character((1:n)[y.index]), sep = "")
							text(x[y.index], y[y.index], labels, adj = 0)
						}
						abline(h = vo, lty = 2)
						abline(h = -vo, lty = 2)
						invisible()
					}

					df <- data.frame(IDX = as.vector(IDX),
						SR = as.vector(SR),
						mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(SR ~ IDX | mod,
						data = df,
						xlab = "Index (Time)",
						panel = panel.special,
						id.n = id.n,
						ylim = y.range,
						ylab = "Standardized Deviances",
						main = "Standardized Deviances vs. Index (Time)",
						strip = function(...) strip.default(..., style = 1),
						vo = vertical.outlier,
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Deviances vs. Index")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					df <- data.frame(	r = as.vector(sqrt(abs(r))),
														pr = as.vector(pr),
														mod = rep(mod.names, rep(n, n.models)))
						
					panel.special <- function(x, y, smooths = smooths,
															rugplot = rugplot, id.n = id.n)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths,
											rugplot = rugplot, id.n = id.n)
						invisible()
					}

					print(xyplot(r ~ pr | mod,
						data = df,
						panel = panel.special,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						xlab = "Fitted Values",
						ylab = "Sqrt(abs(Deviances))",
						main = "sqrt(abs(Deviances)) vs. Fitted Values",
						strip = function(...)
							strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "sqrt(abs(Deviances)) vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
			)
		}
	}
	invisible(x)
}


coef.glmfm <- function(object)
{
	coef.lmfm(object)
}

