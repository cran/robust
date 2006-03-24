print.lmfm <- function(x, ...)
{
	n.models <- length(x)
	model.list <- attr(x, "model.list")
	mod.names <- format(names(model.list))

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
	tmp <- format(round(tmp, digits = options()$digits))
	tmp[tmp.idx] <- ""
	print(tmp, quote = F, ...)
	cat("\n")

	sigmas <- devs <- rep(NA, n.models)
	for(i in 1:n.models) {
		if(is.null(x[[i]]$scale))
		  sigmas[i] <- sqrt(sum(x[[i]]$residuals^2)/x[[i]]$df.resid)
		else
			sigmas[i] <- x[[i]]$scale
		devs[i] <- x[[i]]$deviance
	}

	if(any(!is.na(sigmas))) {
		cat("Residual standard errors:\n")
		for(i in 1:n.models)
			cat(mod.names[i], ":", format(signif(sigmas[i], options()$digits)),
				"on", x[[i]]$df.resid, "degrees of freedom\n")
	}

	if(any(!is.na(devs))) {
		cat("Residual Deviance Estimates:\n")
		for(i in 1:n.models)
		  cat(mod.names[i], ":", format(signif(devs[i], options()$digits)),
				"on", x[[i]]$df.resid, "degrees of freedom\n")
	}

	invisible(x)
}


summary.lmfm <- function(object, correlation = T, ...)
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
	ans$sigmas <- ans$r.squared <- ans$df <- rep(NA, n.models)

	if(correlation)
		ans$correlations <- array(NA, dim = c(n.coefs, n.coefs, n.models))

	for(i in 1:n.models) {
		ans$devs[i] <- object[[i]]$deviance
		tmp <- summary(object[[i]], correlation = correlation)
		coef.idx <- match(names(mod.coefs[[i]]), coef.names)

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

	test <- lapply(object, function(x) grep("lmRob", as.character(x$call[[1]])))
	if(length(test)) {
		test.names <- names(test)
		test <- as.logical(unlist(lapply(test, length)))
		test <- test.names[test]
		bias.test <- list()
		if(length(test)) for (i in test)
			bias.test[[i]] <- test.lmRob(object[[i]])
		ans$biasTest <- bias.test
	}

	oldClass(ans) <- "summary.lmfm"
	ans
}


print.summary.lmfm <- function(x, ...)
{
	digits <- options()$digits - 3
	n.models <- length(x$calls)

	cat("\nCalls: \n")
	for(i in 1:n.models) {
		cat(x$mod.names[i], " : ")
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
	tmp.names <- paste(rep(x$mod.names, np), format(tmp.names))
	dimnames(tmp)[[1]] <- tmp.names
	tmp.idx <- is.na(tmp)
	tmp <- format(round(tmp, digits = digits))
	tmp[tmp.idx] <- ""
	print(tmp, quote = F, ...)

	cat("\nResidual Scale Estimates:\n")
	for(i in 1:n.models)
		cat(x$mod.names[i], ":", format(signif(x$sigmas[i],
			digits)), "on", x$df[i], "degrees of freedom\n")


	cat("\nProportion of variation in response(s) explained by model(s):\n")
	for(i in 1:n.models)
		cat(x$mod.names[i], ":", format(signif(x$r.squared[i], digits)), "\n")

	if(!is.null(x$correlations)) {
		n.coefs <- dim(x$correlations)[1]
		if(n.coefs > 1) {
			cat("\nCorrelations:\n")
			for(i in 1:n.models) {
				cat(x$mod.names[i], "\n")
				tmp <- x$correlations[,  , i]
				tmp.idx <- is.na(tmp)
				tmp <- format(round(tmp, digits = digits))
				tmp[col(tmp) > row(tmp)] <- tmp[tmp.idx] <- ""
				print(tmp, quote = F, ...)
				cat("\n")
			}
		}
	}

	if(length(x$biasTest)) {
		cat("\nBias Tests for Robust Models:")
		the.names <- names(x$biasTest)
		for(i in the.names) {
			cat(paste("\n", i, ":", sep = ""))
			print(x$biasTest[[i]])
		}
		cat("\n")
	}

	invisible(x)
}


panel.addons <- function (x, y, smooths = F, rugplot = F, id.n = 3, ...)
{
	if(smooths)
		panel.loess(x, y, col = 1, ...)
	if(rugplot)
		rug(x)
	if(id.n) {
		n <- length(y)
		id.n <- order(abs(y))[(n - id.n + 1):n]
		labels <- paste("  ", as.character((1:n)[id.n]), sep = "")
		text(x[id.n], y[id.n], labels, adj = 0)
	}
	invisible()
}


plot.lmfm <- function(x, which.plots = "ask", chisq.percent = 0.99,
	vertical.outlier = .99, smooths = F, rugplot = F, id.n = 3,
	envelope = T, half.normal = F, robustQQline = T, mc.samples = 100,
	level = .95, seed = 289, cutoff = T, ...)
{
	##
	## RESIDUAL QQPLOT WITH CONFIDENCE ENVELOPE
	##

	qqplot.with.envelope <- function(px, py, envelope = T, half.normal = F,
														 mc.samples = 100, level = .95, n.models = 1,
														 mod.names = c("a", "b"), id.n = 3,
														 robustQQline = T, seed = 289, ...)

	{

	##
	##  calculate alpha / 2
	##

		level <- (1 - level) / 2

	##
	##  which order statistics to use for envelope
	##

		sample.ci <- function(x, mc.samples, level) {
			lower <- round(level * mc.samples) + 1
			upper <- round((1 - level) * mc.samples)
			sort(x)[c(lower, upper)]
		}

		n <- nrow(px)

		if(envelope && half.normal) {
			temp.seed <- .Random.seed
			set.seed(seed)
			py <- abs(py)
				px <- .5 + (0:(n-1)) / (2*n)
			px <- matrix(rep(qnorm(px), n.models), n, n.models)
			for(i in 1:n.models)
				px[order(py[,i]),i] <- px[,i]
			envelope <- matrix(abs(rnorm(mc.samples * n)), n, mc.samples)
			set.seed(temp.seed)
			envelope <- apply(envelope, 2, function(x) sort(x))
			envelope <- t(apply(envelope, 1, sample.ci,
										mc.samples = mc.samples, level = level))
			den.range <- c(min(py, envelope), max(py, envelope))
		}

		else if(envelope) {
			temp.seed <- .Random.seed
			set.seed(seed)
			envelope <- matrix(rnorm(mc.samples * n), n, mc.samples)
			set.seed(temp.seed)
			envelope <- apply(envelope, 2, function(x) sort(x))
			envelope <- t(apply(envelope, 1, sample.ci,
										mc.samples = mc.samples, level = level))
			den.range <- c(min(py, envelope), max(py, envelope))
		}

		else
			den.range <- c(min(py), max(py))

 		if(length(envelope) > 1 && half.normal) {
			ordered.px <- matrix(0, n, n.models)

			for(i in 1:n.models)
				ordered.px[order(py[,i]),i] <- px[,i]

		 	df <- data.frame(
	 			py = c(	as.vector(py),
								rep(envelope[,1], n.models),
								rep(envelope[,2], n.models)),
				px = c(as.vector(px), rep(as.vector(ordered.px), 2)),
				grp = c(rep("data", n * n.models),
								rep("min", n * n.models),
								rep("max", n * n.models)),
				mod = rep(rep(mod.names, rep(n, n.models)), 3))

			 	panel.special <- function(x, y, id.n, robQQln, ...) {
				 	dat.idx <- 1:(length(x)/3)

				 	panel.xyplot(	x[dat.idx],
				 								y[dat.idx],
				 								pch = 16,
												col = 6,
												...)

				 	if(robQQln)
						abline(coef(lmRob(y[dat.idx] ~ x[dat.idx])))

					panel.addons(	x[dat.idx],
					 							y[dat.idx],
							 					smooths = F,
							 					rugplot = F,
						 						id.n = id.n)

				 	dat.idx <- ((length(x)/3)+1):(2*length(x)/3)

				 	panel.xyplot(	sort(x[dat.idx]),
			 									sort(y[dat.idx]),
			 									type = "l",
	 											col = 1,
 												lty = 2,
 												...)

					dat.idx <- (2*(length(x)/3)+1):(length(x))

					panel.xyplot(	sort(x[dat.idx]),
					 							sort(y[dat.idx]),
								 				type = "l",
												col = 1,
												lty = 2,
												...)
					invisible()
				}
      }

	  	else if(length(envelope) > 1) {
				df <- data.frame(
				py = c(	as.vector(py),
								rep(envelope[,1], n.models),
								rep(envelope[,2], n.models)),
				px = rep(as.vector(px), 3),
			  grp = c(rep("data", n * n.models),
								rep("min", n * n.models),
								rep("max", n * n.models)),
			  mod = rep(rep(mod.names, rep(n, n.models)), 3))

				panel.special <- function(x, y, id.n, robQQln, ...) {
					dat.idx <- 1:(length(x)/3)

				  panel.xyplot(	x[dat.idx],
												y[dat.idx],
						 						pch = 16,
							 					col = 6,
												...)

					panel.addons(	x[dat.idx],
						 						y[dat.idx],
							 					smooths = F,
						 						rugplot = F,
					 							id.n = id.n)

					if(robQQln)
						abline(coef(lmRob(y[dat.idx] ~ x[dat.idx])))

					dat.idx <- ((length(x)/3)+1):(2*length(x)/3)

					panel.xyplot(	sort(x[dat.idx]),
						 						sort(y[dat.idx]),
										 		type = "l",
									 			col = 1,
								 				lty = 2,
							 					...)

					dat.idx <- (2*(length(x)/3)+1):(length(x))

					panel.xyplot(	sort(x[dat.idx]),
					 							sort(y[dat.idx]),
					 							type = "l",
						 						col = 1,
					 							lty = 2,
				 								...)
			 		invisible()
	  		}
		 	}

			else {
				df <- data.frame(
					py = as.vector(py),
					px = as.vector(px),
					mod = rep(mod.names, rep(n, n.models)))

				panel.special <- function(x, y, id.n, robQQln) {
					panel.xyplot(x, y, pch = 16, col = 6)
					panel.addons(x, y, smooths = F, rugplot = F, id.n = id.n)

					if(robQQln) abline(coef(lmRob(y ~ x)))

					invisible()
				}
			}

		print(xyplot(	py ~ px | mod,
							data = df,
							ylab = "Residuals",
							xlab = "Quantiles of Standard Normal",
							main = "Normal QQ-Plot of Residuals",
							id.n = id.n,
							robQQln = robustQQline,
							panel = panel.special,
							strip = function(...)
								strip.default(..., style = 1),
							layout = c(n.models, 1, 1),
							...))
	}

	## END OF qqplot.with.envelope

	##
	## One variable regression plot
	##

	one.variable <- function(x, cutoff = T, contrasts = NULL, ...) 

	{

		n.models <- length(x)
		mod.names <- names(x)
		Call <- x[[1]]$call

		xvar <- as.character(Call$formula[[3]])
		yvar <- as.character(Call$formula[[2]])

		main <- paste(yvar, "~", xvar)

		m <- call("model.frame", formula = Call$formula, data = Call$data)
		m <- eval(m)
		X <- attr(m, "terms")
		X <- model.matrix(X, m, contrasts)

		if(dim(X)[2] == 2)
			X <- X[,2]

		else if(dim(X)[2] > 2)
			stop("This method only supports one explanatory variable")

		Y <- model.extract(m, "response")
		n <- length(Y)

		plot(X, Y, type = "n", xlab = xvar, ylab = yvar, main = main)
		points(X, Y, pch = 16, col = 6)

		the.lines <- integer(n.models)
		the.wd <- integer(n.models)

		for(i in 1:n.models) {
			if(length(grep("Rob", x[[i]]$call))) {
				a <- x[[i]]$yc * x[[i]]$scale
				if(cutoff) abline(coef(x[[i]]) + c(-a, 0), lty = 2)
				abline(coef(x[[i]]), lwd = 2)
				if(cutoff) abline(coef(x[[i]]) + c(a, 0), lty = 2)
				the.lines[i] <- 1
				the.wd[i] <- 2
			}

			else {
				abline(coef(x[[i]]), lty = 4)
				the.lines[i] <- 4
				the.wd[i] <- 1
			}
		}

		key(text = mod.names,
			lines = list(lty = the.lines, lwd = the.wd),
			transparent = T,
			corner = c(.5, 0))

		invisible()
	}

	## END OF One variable regression plot


	##
	## always use comparison plots
	##

	model.list <- attr(x, "model.list")
	n.models <- length(x)
	xsum <- summary.lmfm(x, correlation = F)
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
	r <- na.omit(sapply(x, residuals))
	f <- na.omit(sapply(x, fitted))
	n <- nrow(f)
	idx <- c(1, 3:(n.models + 1))

	##
	##  make denx/deny and px/py on the same plotting scale
	##

	denx <- deny <- matrix(0, 100, n.models)
	px <- py <- matrix(0, n, n.models)

	for(i in 1:n.models) {
		b <- bandwidth.nrd(r[, i])
		den <- density(r[, i], width = b, n = 100, na.rm = T)
		denx[, i] <- den$x
		deny[, i] <- den$y
		tmp <- qqnorm(r[, i], plot = F)
		px[, i] <- tmp$x
		py[, i] <- tmp$y
	}

	den.range <- c(min(py, denx), max(py, denx))

	##
	##  menu choices
	##

	choices <- c("All",
		"Normal QQ-Plot of Residuals", 
		"Estimated Kernel Density of Residuals", 
		"Robust Residuals vs Robust Distances", 
		"Residuals vs Fitted Values", 
		"Sqrt of abs(Residuals) vs Fitted Values", 
		"Response vs Fitted Values", 
		"Standardized Residuals vs Index (Time)", 
		"Overlaid Normal QQ-Plot of Residuals", 
		"Overlaid Estimated Density of Residuals")

	##
	## if only one variable, append one variable plot
	##

	if(length(attr(x[[1]]$terms, "term.labels")) == 1)
		choices <- c(choices, "Data with Robust and LS Fits")
		
	tmenu <- paste("plot:", choices)

	if(is.integer(which.plots)) {
		ask <- F
		which.plots <- c(which.plots + 1, 1)
	}

	else if(which.plots == "all") {
		which.plots <- c(3:11, 1)
		ask <- F
	}

	else
		ask <- T

	while(T) {
		if(ask) {
			which.plots <- lmRob.checkbox(tmenu,
				title = "\nMake plot selections (or 0 to exit):\n")
			if(any(which.plots == 1))
				which.plots <- 2:10
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

				##
				##	Normal QQplot of Residuals
				##

					mpy <- py
					for (i in 1:n.models) {
						if(!is.null(x[[i]]$scale)) 
							mpy[,i] <- mpy[,i]/x[[i]]$scale
						else
							mpy[,i] <- mpy[,i]/stdev(x[[i]]$residuals)
					}

					qqplot.with.envelope(	px, mpy,
																envelope = envelope,
																half.normal = half.normal,
																mc.samples = mc.samples,
																level = level,
																mod.names = mod.names,
																id.n = id.n,
																n.models = n.models,
																robustQQline = robustQQline,
																...)

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residuals QQ")

						graph.number <- graph.number + 1
					}
				}

				,
				{
					panel.special <- function(x, y)
					{
						panel.xyplot(x, y, type = "l", col = 6)
						abline(v = 0, lty = 2)
						invisible()
					}
						
					df <- data.frame(denx = as.vector(denx),
							deny = as.vector(deny),
							mod = rep(mod.names, rep(100, n.models)))
								
					print(xyplot(deny ~ denx | mod,
						data = df,
						xlab = "Residuals",
						ylab = "Kernel Density",
						panel = panel.special,
						main = "Kernel Density of Residuals",
						strip = function(...) strip.default(..., style = 1),
						layout = c(n.models, 1, 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residual Density")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					if(!any(as.vector(sapply(x, function(u)
						dim(numeric.model.matrix(u)))) == 0)) {

						RR <- RD <- matrix(0, n, n.models)

						for(i in 1:n.models) {
							if(is.na(xsum$sigmas[i]) || is.null(xsum$sigmas[i]))
								RR[, i] <- r[, i]
							else RR[, i] <- r[, i] / xsum$sigmas[i]
							tmp.mod <- numeric.model.matrix(x[[i]])
							if (ncol(tmp.mod)>1){
								RD[, i] <- covRob(tmp.mod)$dist
							}
						}

						x.range <- range(RD)
						x.chisq <- sqrt(qchisq(chisq.percent, df = ncol(tmp.mod)))
						if(x.range[2] < x.chisq)
						x.range[2] <- x.chisq
						y.range <- range(RR)
						if(y.range[1] > -vertical.outlier)
							y.range[1] <- -vertical.outlier
						if(y.range[2] < vertical.outlier)
							y.range[2] <- vertical.outlier

						panel.special <- function(x, y, vo = 2.5, chi = 0.975, id.n = 3)
						{
							panel.xyplot(x, y, col = 6, pch = 16)
							n <- length(y)
							y.index <- (1:n)[abs(y) > vo]
							if(length(y.index) > id.n)
								y.index <- order(abs(y))[(n - id.n + 1):n]
							if(length(y.index))
							{
								labels <- paste("  ", as.character((1:n)[y.index]), sep = "")
									text(x[y.index], y[y.index], labels, adj = 0)
							}
							abline(v = chi, lty = 2)
							abline(h = vo, lty = 2)
							abline(h = -vo, lty = 2)
							invisible()
						}

						df <- data.frame(RD = as.vector(RD),
							RR = as.vector(RR),
							mod = rep(mod.names, rep(n, n.models)))

						print(xyplot(RR ~ RD | mod,
							data = df,
							xlab = "Robust Distances",
							panel = panel.special,
							id.n = id.n,
							xlim = x.range,
							ylim = y.range,
							ylab = "Standardized Residuals",
							main = "Standardized Residuals vs Robust Distances",
							strip = function(...) strip.default(..., style = 1),
							vo = vertical.outlier,
							chi = x.chisq,
							layout = c(n.models, 1, 1),
							...))
						}

					else {
						cat("\nMessage:")
						cat("\n    All regressors are factors.  Robust Residuals vs. Robust")
						cat("\n    Distances skipped.\n")
					}

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residuals vs. Robust Distances")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					panel.special <- function(x, y, smooths = F, rugplot = F, id.n = 3)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths,
							rugplot = rugplot, id.n = id.n)
						abline(h = 0, lty = 2)
						invisible()
					}
						
					df <- data.frame(f = as.vector(f),
						r = as.vector(r),
						mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(r ~ f | mod,
						data = df,
						xlab = "Fitted Values",
						ylab = "Residuals",
						main = "Residuals vs. Fitted Values",
						panel = panel.special,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						strip = function(...) strip.default(..., style = 1),
						layout = c(n.models, 1, 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residuals vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					panel.special <- function(x, y, smooths = F, rugplot = F, id.n = 3)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths,
							rugplot = rugplot, id.n = id.n)
						invisible()
					}
						
					df <- data.frame(f = as.vector(f),
						r = as.vector(sqrt(abs(r))),
						mod = rep(mod.names, rep(n, n.models)))
							
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

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Sqrt(abs(Residuals)) vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					y <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						m <- call("model.frame",
											formula = xsum$calls[[i]]$formula,
											data = xsum$calls[[i]]$data,
											na.action = xsum$calls[[i]]$na.action)
						if(!is.null(xsum$calls[[i]]$subset))
							m$subset <- xsum$calls[[i]]$subset
						if(!is.null(xsum$calls[[i]]$weights))
							m$weights <- xsum$calls[[i]]$weights
						m <- eval(m, sys.parent())
						y[, i] <- model.extract(m, response)
					}

					panel.special <- function(x, y, smooths = F, rugplot = F)
					{
						panel.xyplot(x, y, col = 6, pch = 16)
						panel.addons(x, y, smooths = smooths, rugplot = rugplot, id.n = 0)
						abline(0, 1, lty = 2)
						invisible()
					}
						
					df <- data.frame(y = as.vector(y),
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
						strip = function(...) strip.default(..., style = 1),
						layout = c(n.models, 1, 1),
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
					SR <- IDX <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						if(is.na(xsum$sigmas[i]) || is.null(xsum$sigmas[i]))
							SR[, i] <- r[, i]
						else SR[, i] <- r[, i] / xsum$sigmas[i]
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
						ylab = "Standardized Residuals",
						main = "Standardized Residuals vs. Index (Time)",
						strip = function(...) strip.default(..., style = 1),
						vo = vertical.outlier,
						layout = c(n.models, 1, 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Standardized Residuals vs. Index")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					x.range <- range(c(py, denx))
					matplot(py, px,
						pch = c(16, 2:n.models)[1:n.models],
						col = c(1, 6, 2:n.models)[1:n.models],
						xlim = x.range,
						xlab = "Residuals",
						ylab = "Quantiles of Standard Normal",
						main = "Normal QQ-Plot of Residuals",
						lty = 1,
						...)
							
					key(x.range[1], max(px),
						text = list(mod.names),
						lines = list(type = "p", col = c(1, 6, 2:n.models)[1:n.models],
							pch = c(16, 2:n.models)[1:n.models]),
						transparent = T)
	
					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Overlaid Residuals QQ")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					x.range <- range(c(py, denx))
					matplot(denx, deny,
						type = "n",
						xlim = x.range,
						xlab = "Residuals",
						ylab = "Kernel Density",
						main = "Kernel Density of Residuals",
						...)

					for(i in 1:n.models)
						lines(denx[,i], deny[,i], lwd = 2,
							lty = (3 * ((1:n.models) - 1) + 1)[i],
							col = c(1, 6, 2:n.models)[1:n.models][i])

					key(x.range[1], max(deny),
						text = list(mod.names),
						lines = list(type = "l",
											col = c(1, 6, 2:n.models)[1:n.models],
											lty = 3 * ((1:n.models) - 1) + 1,
											lwd = 2),
						transparent = T)
	
					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Overlaid Residual Density")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					##
					## Robust and LS Fits + Data
					##

					one.variable(x, cutoff = cutoff, ...)
	
					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Scatter Plot")

						graph.number <- graph.number + 1
					}
				}
			)
		}
	}
}


partial.plot.lmfm <- function(object, fit = F, scale = T, rugplot = F, ...)
{
	models <- names(object)
	n.models <- length(models)
	fits <- list()
	partial.res <- list()
	print.warning <- F

	for(model in models) {
		if(!is.null(object[[model]]$na.action))
			object[[model]]$na.action <- NULL
		terms <- labels.lm(object[[model]])
		Terms <- object[[model]]$terms
		a <- attributes(Terms)
		Call <- object[[model]]$call
		all.terms <- labels(Terms)
		xvars <- as.vector(Terms)
		names(xvars) <- all.terms
		terms <- match.arg(terms, all.terms)
		Interactions <- a$order > 1

		if(any(Interactions)) {
			all.terms <- all.terms[!Interactions]
			TM <- match(terms, all.terms, 0)
			if(!all(TM)) {
				terms <- terms[TM > 0]
				print.warning <- T
			}
		}

		xvars <- xvars[terms]
		xnames <- as.list(terms)
		names(xnames) <- terms
		modes <- sapply(xvars, mode)

		for(term in terms[modes != "name"]) {
			evars <- all.names(xvars[term], functions = F, unique = T)
			if(!length(evars))
				next
			xnames[[term]] <- evars
			evars <- parse(text = evars)
			if(length(evars) == 1)
				evars <- evars[[1]]
			else {
				evars <- c(as.name("list"), evars)
				mode(evars) <- "call"
			}
			xvars[[term]] <- evars
		}

		xvars <- c(as.name("list"), xvars)
		mode(xvars) <- "call"

		if(!is.null(Call$subset) | !is.null(Call$na.action) |
			!is.null(options("na.action")[[1]])) {

			Rownames <- names(object[[model]]$fitted)
			if(!(Rl <- length(Rownames)))
				stop("need to have names for fitted.values when	call has a subset or na.action argument")
			Mcall <- c(as.name("model.frame"), list(formula = 
				terms.inner(parse(text = unlist(xnames))),
				subset = Rownames, na.action = function(x) x))
			mode(Mcall) <- "call"
			Mcall$data <- Call$data
			xvars <- eval(xvars, eval(Mcall))
		}

		else {
			ecall <- substitute(eval(expression(xvars)))
			ecall$local <- Call$data
			xvars <- eval(ecall)

		}
		fits[[model]] <- predict(object[[model]], type = "terms",
			terms = terms, se.fit = T)$fit
		partial.res[[model]] <- fits[[model]] + residuals(object[[model]])
	}

	if(scale) {
		y.limits <- matrix(0, ncol = 2, nrow = n.models)
		for(i in 1:n.models)
			y.limits[i,  ] <- range(partial.res[[i]])
		y.limits <- range(y.limits)
	}

	else
		y.limits <- NULL

	n.obs <- length(xvars[[1]])
	df <- list()

	if(dev.cur() == 1)
		trellis.device()

	for(term in terms) {
		partials <- predicted <- vector()

		for(model in models) {
			partials <- c(partials, partial.res[[model]][, term])
			predicted <- c(predicted, fits[[model]][, term])
		}

		x <- xvars[[term]]

		if(is.factor(x)) {
			factor.levels <- levels(x)
			store <- matrix(nrow = length(factor.levels), ncol = n.models)
			dimnames(store) <- list(factor.levels, models)
			for(model in models) {
				temp <- data.frame(fac = x, x = fits[[model]][, term])
				for(level in factor.levels) {
					store[level, model] <- temp[temp$fac == level, ]$x[1]
				}
			}

			x <- as.numeric(x)
			x.limits <- range(x) + c(-0.5, 0.5)
			temp.seed <- .Random.seed
			set.seed(317)
			x <- x + runif(length(x), min = -0.2, max = 0.2)
			set.seed(temp.seed)
			x.var <- y.var <- vector()

			for(model in models) {
				x.var <- c(x.var, x, 1:dim(store)[1])
				y.var <- c(y.var, partial.res[[model]][, term],
					store[, model])
			}

			df <- data.frame(x = x.var, y = y.var, mod = rep(models,
				each = n.obs + length(factor.levels)))
			scales <- list(x = list(relation = "same",
				at = 1:length(factor.levels), labels = factor.levels, crt = 30))

			panel.special <- function(x, y, rugplot, n.levels, fit, ...) {
				resids <- 1:(length(x) - n.levels)
				panel.xyplot(x[resids], y[resids], pch = 16,
					col = 6)
				if(fit) {
					resids <- (length(resids) + 1):length(
						x)
					segments(x[resids] - 0.2, y[resids],
						x[resids] + 0.2, y[resids],
						lwd = 2)
				}
			}

			print(xyplot(y ~ x | mod,
				data = df,
				xlab = term,
				ylab = paste("Partial for", term),
				xlim = x.limits,
				ylim = y.limits,
				scales = scales,
				panel = panel.special,
				fit = fit,
				n.levels = length(factor.levels),
				rugplot = rugplot,
				strip = function(...)
					strip.default(..., style = 1), ...))
		}

		else {
			if(fit) {
				df <- data.frame(x = rep(x, 2 * n.models),
					y = c(partials, predicted),
					type = c(rep("partials", length(partials)),
						rep("predicted", length(predicted))),
					mod = c(rep(models, each = n.obs),
						rep(models, each = n.obs)))
			}

			else {
				df <- data.frame(x = rep(x, n.models),
					y = partials,
					type = rep("partials", length(partials)),
					mod = c(rep(models, each = n.obs)))
			}

			panel.special <- function(x, y, groups, rugplot, ...) {
				panel.superpose(x, y, groups = groups, pch = 16,
					lwd = c(1, 2), type = c("p", "l"), col = c(6, 1), ...)
				if(rugplot)
					rug(x)
			}

			print(xyplot(y ~ x | mod,
				data = df,
				groups = type,
				xlab = term,
				ylab = paste("Partial for", term),
				panel = panel.special,
				rugplot = rugplot,
				strip = function(...)
					strip.default(..., style = 1), ...))
		}
	}

	if(print.warning)
		warning("No terms saved for \"a:b\" style interaction terms")

	invisible(object)
}


coef.lmfm <- function(object)
{
	n.models <- length(object)
	the.coefs <- lapply(object, coef)
	the.names <- list(names(the.coefs), names(the.coefs[[1]]))
	the.coefs <- matrix(unlist(the.coefs), nrow = n.models, byrow = T)
	dimnames(the.coefs) <- the.names
	the.coefs
}

