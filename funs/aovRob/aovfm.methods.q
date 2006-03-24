print.aovfm <- function(x, intercept = F, ...)
{
	n.models <- length(x)
	digits <- options()$digits - 3
	model.list <- attr(x, "model.list")
	mod.names <- format(names(model.list))

	cat("\nCalls: \n")
	for(i in 1:n.models) {
		cat(mod.names[i], ": ")
		print(x[[i]]$call)
	}

	effects <- lapply(x, function(x) x$effects)
	asgn <- lapply(x, function(x) x$assign)
	rl <- unique(unlist(lapply(asgn, names)))
	nterms <- length(rl)
	df <- ss <- matrix(NA, ncol = nterms, nrow = n.models)
		if(nterms) {
			for(j in 1:n.models) {
			  tmp <- match(names(asgn[[j]]), rl)
			  for(i in 1:length(tmp)) {
			    ai <- asgn[[j]][[i]]
			    df[j, tmp[i]] <- length(ai)
			    ss[j, tmp[i]] <- sum(effects[[j]][ai]^2)
			  }
			}
			keep <- df > 0
			for(j in 1:n.models) {
			  if(!intercept && names(asgn[[j]])[1] == "(Intercept)")
			    keep[j, 1] <- F
			}
			df[as.vector(!keep)] <- NA
			ss[as.vector(!keep)] <- NA
			keep <- apply(keep, 2, any)
			rl <- rl[keep]
			nterms <- length(rl)
			df <- df[, keep]	
			ss <- ss[, keep]
		}

		cat("\nTerms:\n")

		if(nterms == 0) { #empty model
			prmatrix(array(0, c(1, 2), list("<empty>", c("Sum of Squares", 
			  "Deg. of Freedom"))))
			return(invisible(x))
		}

		df.res <- sapply(x, function(x) x$df.resid)

		if(any(df.res > 0)) {
			nterms <- nterms + 1
			rl[nterms] <- "Residuals"
			df <- cbind(df, df.res)
			ss <- cbind(ss, sapply(x, function(x)
			sum((x$resid)^2)))
		}

		ss <- zapsmall(ss)
		tmp <- rbind(format(ss), format(df))
		dimnames(tmp) <- list(c(paste(mod.names, "Sum of Squares"),
			paste(mod.names, "Deg. ofFreedom")), rl)
		prmatrix(tmp, quote = F, right = T)
		cat("\n")

# residual standard errors

	sigmas <- rep(NA, n.models)
	for(i in 1:n.models) {
		if(model.list[[i]][[1]] == "lm" || model.list[[i]][[1]] == "aov")
		  sigmas[i] <- sqrt(sum(x[[i]]$residuals^2)/x[[i]]$df.resid)
		else sigmas[i] <- x[[i]]$scale
	}

	if(any(!is.na(sigmas))) {
		cat("Residual standard errors:\n")
		for(i in 1:n.models)
		  cat(mod.names[i], ":", format(signif(sigmas[i], digits)), "on", 
		    x[[i]]$df.resid, "degrees of freedom\n")
	}

	invisible(x)
}


summary.aovfm <- function(object, ...)
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
	ans$aovtable <- lapply(object, summary)
	n.columns <- sapply(ans$aovtable, function(x) dim(x)[[2]])
	if(length(unique(n.columns)) != 1) {
			
		tmp.fun <- function(tmp) {
			nc <- dim(tmp)[[2]]
			tmp <- tmp[, c(1, nc - 1, nc), drop = F]
			tmp
		}

		ans$aovtable <- lapply(ans$aovtable, tmp.fun)
	}

	tmp.fun <- function(tmp) {
		tmpnames <- dimnames(tmp)[[1]]
		nnames <- length(tmpnames)
		if(any(strip.blanks(tmpnames) == "Residuals"))
			tmp <- tmp[-nnames, , drop = F]
		tmp
	}

	ans$aovtable <- lapply(ans$aovtable, tmp.fun)

	oldClass(ans) <- "summary.aovfm"
	ans
}


print.summary.aovfm <- function(x, ...)
{
	n.models <- length(x$calls)
	digits <- options()$digits - 3

	coln <- unique(strip.blanks(unlist(lapply(x$aovtable,
							function(x) dimnames(x)[[1]]))))
	np <- length(coln)
	colnames <- matrix(rep(coln, n.models), nrow = n.models, byrow = T)
	colnames <- paste(x$mod.names, colnames)
		n.columns <- dim(x$aovtable[[1]])[[2]]
	tmp <- matrix(NA, ncol = n.columns, nrow = length(colnames))
	dimnames(tmp) <- list(colnames, dimnames(x$aovtable[[1]])[[2]])
	for(i in 1:n.models) {
	  idx <- match(strip.blanks(dimnames(x$aovtable[[i]])[[1]]), coln)
	  tmp[n.models*(idx-1)+i, ] <- as.matrix(x$aovtable[[i]])
	}

	cat("\nComparison of ANOVA Tables:\n")
	tmp.idx <- is.na(tmp)
	tmp[, 2:n.columns] <- format(round(tmp[, 2:n.columns], digits = digits))
	tmp[tmp.idx] <- ""
	print(tmp, quote = F, right = T, ...)

	invisible(x)
}


plot.aovfm <- function(x, which.plots = "ask", vertical.outlier = .99,
	smooths = F, rugplot = F, id.n = 3, envelope = T, half.normal = F,
	robustQQline = T, mc.samples = 100, level = .95, seed = 289, ...)
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
							...))
	}

	## END OF qqplot.with.envelope ##

	model.list <- attr(x, "model.list")
	n.models <- length(x)
	xsum <- summary.aovfm(x, correlation = F)
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
		"Residuals vs Fitted Values", 
		"Sqrt of abs(Residuals) vs Fitted Values", 
		"Response vs Fitted Values", 
		"Residual-Fit Spread",
		"Standardized Residuals vs Index (Time)",
		"Overlaid Normal QQ-Plot of Residuals", 
		"Overlaid Estimated Density of Residuals")
	tmenu <- paste("plot:", choices)

	ask <- T

	if(is.integer(which.plots)) {
		ask <- F
		which.plots <- c(which.plots + 1, 1)
	}

	else if(which.plots == "all") {
		ask <- F
		which.plots <- c(3:11, 1)
	}

	while(T) {
		if(ask) {
			which.plots <- lmRob.checkbox(tmenu, title = 
				"\nMake plot selections (or 0 to exit):\n")
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

					qqplot.with.envelope(	px, py,
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
					df <- data.frame(
						denx =	as.vector(denx),
						deny = as.vector(deny),
						mod = rep(mod.names, rep(100, n.models)))
							
					print(xyplot(deny ~ denx | mod,
						data = df,
						xlab = "Residuals",
						ylab = "Kernel Density",
						panel = panel.special,
						main = "Kernel Density of Residuals",
						strip = function(...) strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Kernel Density of Residuals")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					panel.special <- function(x, y, smooths, rugplot, id.n)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths, rugplot = rugplot,
								id.n = id.n)
						invisible()
					}
						
					df <- data.frame(
						f = as.vector(f),
						r = as.vector(r),
						mod = rep(mod.names,rep(n, n.models)))
						print(xyplot(r ~ f | mod,
						data = df,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						xlab = "Fitted Values",
						ylab = "Residuals",
						main = "Residuals vs. Fitted Values",
						panel = panel.special,
						strip = function(...) strip.default(..., style = 1),
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
					panel.special <- function(x, y, smooths, rugplot, id.n)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths, rugplot = rugplot,
							id.n = id.n)
						invisible()
					}
						
					df <- data.frame(
						f = as.vector(f),
						r = as.vector(sqrt(abs(r))),
						mod = rep(mod.names, rep(n, n.models)))
							
					print(xyplot(r ~ f | mod,
						data = df,
						smooths = smooths,
						rugplot = rugplot,
						id.n = id.n,
						panel = panel.special,
						xlab = "Fitted Values",
						ylab = "Sqrt(abs(Residuals))",
						main = "sqrt(abs(Residuals)) vs. Fitted Values",
						strip = function(...) strip.default(..., style = 1),
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "sqrt(abs(Residuals)) vs. Fitted Values")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					panel.special <- function(x, y, smooths, rugplot)
					{
						panel.xyplot(x, y, pch = 16, col = 6)
						panel.addons(x, y, smooths = smooths, rugplot = rugplot,
								id.n = F)
						invisible()
					}

					y <- matrix(0, n, n.models)

					for(i in 1:n.models) {
						m <- x[[i]]$call
						m <- call("model.frame", formula = m$formula, data = m$data)
						m <- eval(m, sys.parent())
						y[, i] <- model.extract(m, response)
					}

					df <- data.frame(	y = as.vector(y),
														f = as.vector(f),
														mod = rep(mod.names, rep(n, n.models)))

					print(xyplot(y ~ f | mod,
						data = df,
						smooths = smooths,
						rugplot = rugplot,
						xlab = "Fitted Values", 
						ylab = "Response",
						main = "Response vs. Fitted Values",
						panel = panel.special,
						strip = function(...) strip.default(..., style = 1),
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
					resfit <- numeric()
					title <- character()
					frac <- numeric()

					for(i in 1:n.models) {
						resfit <- c(resfit, sort(x[[i]]$fitted - mean(x[[i]]$fitted)))
						resfit <- c(resfit, sort(x[[i]]$residuals))
						
						data.len <- length(x[[i]]$fitted)
						
						title <- c(title, c(rep(paste(mod.names[i], " Fitted", sep = ":"), data.len),
							rep(paste(mod.names[i], " Residuals", sep = ":"), data.len)))
							
						frac <- c(frac, rep((1:data.len)/data.len,2))
					}

					print(xyplot(resfit ~ frac | as.factor(title),
						pch = 16,
						col = 6,
						main = "Residual-Fit Spread",
						xlab = "f-value",
						ylab = "",
						between = list(x = 0, y = 1),	
						strip = function(...) strip.default(..., style = 1)))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residual-Fit Spread")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					SR <- IDX <- matrix(0, n, n.models)
					for(i in 1:n.models) {
						if(!is.null(x[[i]]$scale))
							SR[, i] <- r[, i] / x[[i]]$scale
						else {
							x.scale <- sqrt(sum(x[[i]]$residuals^2)/x[[i]]$df.residual)
							SR[, i] <- r[, i] / x.scale
						}
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
						...))

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Residuals vs. Index")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					plot.char <- c(16, 2:n.models)[1:n.models]

					x.range <- range(c(py, denx))
					matplot(py, px,
						pch = plot.char,
						col = colors,
						xlim = x.range,
						xlab = "Residuals",
						ylab = "Quantiles of Standard Normal",
						main = "Normal QQ-Plot of Residuals",
						lty = c(1,1),
						...)

					key(x.range[1], max(px),
						text = list(mod.names),
						lines = list(	type = "p",
													col = colors,
													pch = plot.char),
						transparent = T)

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "QQ-Plot of Residuals")

						graph.number <- graph.number + 1
					}
				}
				,
				{
					wts <- 2

					x.range <- range(c(py, denx))
					matplot(denx, deny,
						type = "l",
						lwd = wts,
						lty = styles,
						col = colors,
						xlim = x.range,
						xlab = "Residuals",
						ylab = "Kernel Density",
						main = "Kernel Density of Residuals",
						...)

					abline(v = 0, lty = 2)
					key(x.range[1], max(deny),
						text = list(mod.names),
						lines = list( type = "l",
													col = colors,
													lwd = wts,
													lty = styles), 
						transparent = T)

					if(names(dev.cur()) == "graphsheet") {
						guiModify("GraphSheetPage",
							Name = paste("$", graph.number, sep = ""),
							NewName = "Overlaid Residual Density")

						graph.number <- graph.number + 1
					}
				}
			)
		}
	}

	invisible(x)
}


coef.aovfm <- function(object)
{
	coef.lmfm(object)
}


dummy.coef.aovfm <- function(object)
{
	n.models <- length(object)
	mod.names <- names(object)
	table.list <- list()
	for(i in 1:n.models) {
		oldClass(object[[i]]) <- "aov"
		table.list[[i]] <- dummy.coef(object[[i]])
	}
	n.coef <- length(table.list[[1]])
	if(n.models > 1) {
		for(i in 1:n.coef) {
			for(j in 2:n.models) {
				table.list[[1]][[i]] <- rbind(table.list[[1]][[i]],
					table.list[[j]][[i]])
			}
			dimnames(table.list[[1]][[i]]) <- list(mod.names, 
				dimnames(table.list[[1]][[i]])[[2]])
		}
	}
	table.list[[1]]
}
