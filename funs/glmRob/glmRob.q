# We now call S_GICSTP
# and S_GYTSTP 

# have to fix glmRob.modified.cubif.null


##
## file: glmrob.q
##


"glmRob.cubif.control" <- function(epsilon = 0.001, maxit = 50, 
         bpar=2, cpar=1.5, trc=F, ...)
{
        list(epsilon = epsilon, maxit = maxit,
             bpar = bpar, cpar=cpar, trc=trc) 
}

##
## glmRob
##

glmRob <- function(formula = formula(data), family = binomial, data = sys.parent(),
	subset, na.action, start = eta, fit.method = 'cubif', model = F, x = F,
  y = T, contrasts = NULL, cubif.control = glmRob.cubif.control(...), 
	misclass.control = glmRob.misclass.control(...),
  mallows.control = glmRob.mallows.control(...), method = "glmRob.fit",
	estim = "mcd", robust.cov.control = covRob.control(estim = estim, quan = .75, ...),
  ...)
{
	Call <- match.call()
	robust.cov.control$estim <- estim

#	m <- match.call(expand = F)
#	m$family <- m$fit.method <- m$model <- m$x <- m$y <- m$control <- NULL 
#	m$cubif.control <- NULL
#	m$misclass.control <- m$mallows.control <- NULL
#	m$contrasts <- m$... <- NULL
#	m$method <- NULL
#	m$robust.cov.control <- NULL
#	m[[1]] <- as.name("model.frame")
  
  m <- call("model.frame", formula = substitute(formula))
  
  if(!missing(data))
    m$data <- substitute(data)

  if(!missing(start))
    m$start <- substitute(start)

  if(!missing(na.action))
    m$na.action <- substitute(na.action)

	m <- eval(m, sys.parent())
	Terms <- attr(m, "terms")

	if(method == "model.frame")
		return(m)

#	xvars <- as.character(attr(Terms, "variables"))
#	if(length(xvars) > 0) {
#		xlevels <- lapply(m[xvars], levels)
#		xlevels <- xlevels[!sapply(xlevels, is.null)]
#		if(length(xlevels) == 0)
#			xlevels <- NULL
#	}
#
#	else
#		xlevels <- NULL

	a <- attributes(m)
	Y <- model.extract(m, response)
	X <- model.matrix(Terms, m, contrasts)
	start <- model.extract(m, start)
	offset <- model.extract(m, offset)
	family <- as.family(family)

	if(fit.method == "cubif")
		fit <- glmRob.cubif(x = X, y = Y, intercept = F, offset = offset, family = family,
      null.dev = T, control = cubif.control, robust.cov.control = robust.cov.control,
      ...)

	else if(family$family["name"] != "Binomial") 
		stop(paste("method ", fit.method, " not implemented for family ",
			family$family["name"], sep=''))

	else if(fit.method == "misclass") 
		fit <- glmRob.misclass(x = X, y = Y, control = misclass.control, offset = offset,
      null.dev = T, family = family, Terms = Terms, ...) 

	else if(fit.method == "mallows") 
		fit <- glmRob.mallows(x = X, y = Y, control = mallows.control, offset = offset,
    null.dev = T, family = family, Terms = Terms, ...)

	else 
		stop(paste('method: ', fit.method, ' not implemented'))

	oldClass(fit) <- "glmRob"
	fit$terms <- Terms
	fit$formula <- as.vector(attr(Terms, "formula"))
	fit$fit.method <- fit.method
	fit$call <- Call
	fit$robust.cov.control <- robust.cov.control
	if(model)
		fit$model <- m
	if(x)
		fit$x <- X
	if(!y)
		fit$y <- NULL 
	fit$R <- qr.R(qr(X))

	fit
}





glmRob.gintac <- function(x, y, ni, oi, icase, maxtt, maxta, tolt, tola, b, c) 
{
	mdx <- nrow(x)
	n <- length(y)
	np <- ncol(x)
	mdt <- n
	ncov <- np*(np+1)/2
	nitt <- integer(1)
	nita <- integer(1)
	sigma <- double(1)
	a <- double(ncov)
	if (length(oi)==1) oi <- rep(0,n)
	theta <- double(mdt)
	ci <- double(n)
	dist <- double(n)
	rw1 <- double(5*ncov+3*n)
	rw2 <- matrix(double(1),mdx,np)
	iw1 <- integer(np)
	dw1 <- double(2*ncov+np+n)
	f.res <- .Fortran( "s_gintac",
		x=as.double(x),
		y=as.double(y),
		ni=as.integer(ni),
		oi=as.double(oi),
		mdx=as.integer(mdx),
		mdt=as.integer(mdt),
		n=as.integer(n),
		np=as.integer(np),
		ncov=as.integer(ncov),
		icase=as.integer(icase),
		maxtt=as.integer(maxtt),
		maxta=as.integer(maxta),
		tau = as.double(tolt),
		tolt=as.double(tolt),
		tola=as.double(tola),
		b=as.double(b),
		c=as.double(c),
		nitt=as.integer(nitt),
		nita=as.integer(nita),
		sigma=as.double(sigma),
		a=as.double(a),
		theta=as.double(theta),
		ci=as.double(ci),
		dist=as.double(dist),
		rw1=as.double(rw1),
		rw2=as.double(rw2),
		iw1=as.integer(iw1),
		dw1=as.double(dw1),
		ips=as.integer(3),
		xk=as.double(1.5))
	list(nitt=f.res$nitt, nita=f.res$nita, sigma=f.res$sigma,
		a=f.res$a, theta=f.res$theta, ci=f.res$ci, dist=f.res$dist)
}


glmRob.gfedca <- function(vtheta,ci,wa,ni,oi=0,icase) 
{
	n <- length(vtheta)
	d <- double(n)
	e <- double(n)
	if (length(oi)==1) oi <- rep(0,n)
	f.res_.Fortran("s_gfedca",
		vtheta=as.double(vtheta),
		ci=as.double(ci),
		wa=as.double(wa),
		ni=as.integer(ni),
		oi=as.double(oi),
		n=as.integer(n),
		icase=as.integer(icase),
		d=as.double(d),
		e=as.double(e))
	list(d=f.res$d,e=f.res$e)
}


glmRob.ktaskw <- function(x,d,e,tau,ia,f,f1,
			iainv,a) {
	n_length(d)
	np_ncol(x)
	mdx_nrow(x)
	mdsc_np
	ncov_np*(np+1)/2
	if (missing(a)) a_double(ncov)
	s1inv_double(ncov)
	s2_double(ncov)
	ainv_double(ncov)
	cov_double(ncov)
	sc_matrix(double(1),mdsc,np)
	f.res_.Fortran("s_ktasbi",
		x=as.double(x),
		d=as.double(d),
		e=as.double(e),
		n=as.integer(n),
		np=as.integer(np),
		mdx=as.integer(mdx),
		mdsc=as.integer(mdsc),
		ncov=as.integer(ncov),
		tau=as.double(tau),
		ia=as.integer(ia),
		f=as.double(f),
		f1=as.double(f1),
		iainv=as.integer(iainv),
		a=as.double(a),
		s1inv=as.double(s1inv),
		s2=as.double(s2),
		ainv=as.double(ainv),
		cov=as.double(cov),	
		sc=as.double(sc))
list(a=f.res$a,s1inv=f.res$s1inv,s2=f.res$s2,ainv=f.res$ainv,cov=f.res$cov)
}

glmRob.gymain <- function(x, y, ni, cov, a, theta, oi=0, b, 
	gam, tau, icase, iugl, iopt, ialg, icnvt, icnva, 
	maxit, maxtt, maxta, maxtc, 
	tol, tolt, tola, tolc, trc=F) 
{
	mdx_nrow(x)
	n_length(y)
	np_ncol(x)
	ncov_length(cov)
	if (length(oi)==1) oi <- rep(0,n)
	if (missing(a)) a_double(ncov)
	if (missing(theta)) theta_double(np)
	nit_integer(1)
	ci_double(n)
	wa_double(n)
	vtheta_double(n)
	delta_double(np)
	grad_double(np)
	hessnv_double(ncov)
	rw1_double(5*ncov+3*n)
	rw2_matrix(double(1),mdx,np)
	iw1_integer(np)
	dw1_double(2*ncov+np+n)

	if(trc) trc <- 1
	else trc <- 0
	f.res <- .Fortran("s_gymain",
		x=as.double(x),
		y=as.double(y),
		ni=as.integer(ni),
		cov=as.double(cov),
		a=as.double(a),
		theta=as.double(theta),
		oi=as.double(oi),
		mdx=as.integer(mdx),
		n=as.integer(n),
		np=as.integer(np),
		ncov=as.integer(ncov),
		b=as.double(b),
		gam=as.double(gam),
		tau=as.double(tau),
		icase=as.integer(icase),
		iugl=as.integer(iugl),
		iopt=as.integer(iopt),
		ialg=as.integer(ialg),
		icnvt=as.integer(icnvt),
		icnva=as.integer(icnva),
		maxit=as.integer(maxit),
		maxtt=as.integer(maxtt),
		maxta=as.integer(maxta),
		maxtc=as.integer(maxtc),
		tol=as.double(tol),
		tolt=as.double(tolt),
		tola=as.double(tola),
		tolc=as.double(tolc),
		nit=as.integer(nit),
		ci=as.double(ci),
		wa=as.double(wa),
		vtheta=as.double(vtheta),
		delta=as.double(delta),
		grad=as.double(grad),
		hessnv=as.double(hessnv),
		rw1=as.double(rw1),
		rw2=as.double(rw2),
		iw1=as.integer(iw1),
		dw1=as.double(dw1),
		trc=as.integer(trc))
	list(a=f.res$a,theta=f.res$theta,nit=f.res$nit,ci=f.res$ci,wa=f.res$wa,
		vtheta=f.res$vtheta,delta=f.res$delta,grad=f.res$grad,
		hessnv=f.res$hessnv, rw1=f.res$rw1, tol=tol)
}

glmRob.glmdev <- function(y, ni, ci, wa, vtheta, offset=0, icase = ics)
{
        n <- length(y)
        dev <- double(1)
        thetas <- double(n)
        li  <- double(n)      
        sc  <- double(n)  
        f.res <- .Fortran("s_glmdev",
                y = as.double(y),
                ni = as.integer(ni),
                ci = as.double(ci),
                wa = as.double(wa),
                vtheta = as.double(vtheta),
                offset = as.double(offset),
                n = as.integer(n),
                icase = as.integer(icase),
                dev = as.double(dev),
                thetas = as.double(thetas),
                li = as.double(li),
                sc = as.double(sc))  
        sc <- f.res$sc 
        list(dev = f.res$dev, thetas = f.res$thetas, li = f.res$li, 
		sc = f.res$sc)
}


"print.glmRob" <- function(x, ai=F, ci=F, A.mat=F, 
		digits=.Options$digits-3,...)
{
        if(!is.null(cl <- x$call)) {
                cat("Call:\n")
                dput(cl)
        }
        coef <- x$coef
        if(any(nas <- is.na(coef))) {
                if(is.null(names(coef))) names(coef) <- paste("b", 1:length(
                                coef), sep = "")        
                cat("\nCoefficients: (", sum(nas),
                        " not defined because of singularities)\n", sep = "")
        }
        else cat("\nCoefficients:\n")
        print(coef, ...)
        rank <- x$rank
        if(is.null(rank))
                rank <- sum(!nas)
        nobs <- length(x$residuals)
        rdf <- x$df.resid
        if(is.null(rdf))
                rdf <- nobs - rank
        cat("\nDegrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
        cat("Residual Deviance:", format(x$deviance), "\n")
        if (A.mat) {
          A <- x$A
          if(!is.null(A)) {
                  p <- dim(A)[2]
                  if(p > 1) {
                          cat("\nLower-triangular A matrix :\n")
                          ll <- row(A) >= col(A)
                          A  <- format(round(A, digits), ...)
                          A[!ll] <- ""
                          print(A[drop = F], quote = F)
                  }
          }
        }
        if (ai) {
             cat("\nContants a_1,..., a_n :\n")
             print(x$ai, quote=F)
        }
        if (ci) {
             cat("\nContants c_1,..., c_n :\n")
             print(x$ci, quote=F)
        }
        invisible(x)
}

"summary.glmRob" <- function(object, dispersion = NULL,
		correlation = T, tl=1.e-10)
{
        coef <- object$coef
        resid <- object$residuals
        dresid <- residuals(object, "deviance")
        n <- length(resid)
        p <- object$rank
        if(is.null(p))
                p <- sum(!is.na(coef))
        if(!p) {
                warning("This model has zero rank --- no summary is provided")
                return(object)
        }
        nsingular <- length(coef) - p
        rdf <- object$df.resid
        if(is.null(rdf))
                rdf <- n - p
        famname <- object$family["name"]
        dispersion <- 1
        names(dispersion) <- famname

        covun <- object$cov
        var <- diag(covun)
        nas <- is.na(coef)
        cnames <- names(coef[!nas])
        coef <- matrix(rep(coef[!nas], 3), ncol = 3)
        dimnames(coef) <- list(cnames, c("Value", "Std. Error", "t value"))
        coef[, 2] <- sqrt(var)
        coef[, 3] <- coef[, 1]/coef[, 2]
        if(correlation) {
                cor <- covun
                for (i in 1:nrow(cor)) {
                  if (var[i]<1.e-10)
                    {str <- paste("Variance number",i,"smaller than 1.e-10",
                                "(set to 1.e-10)")
                     print(str)}
                  cor[i,1:i] <- cor[i,1:i]/sqrt(var[i]*var[1:i])
                  cor[1:i,i] <- cor[i,1:i]
                }
                dimnames(cor) <- list(cnames, cnames)
        }
        else cor <- NULL
        dimnames(covun) <- list(cnames, cnames)
        ocall <- object$call
        if(!is.null(form <- object$formula)) {
                if(is.null(ocall$formula))
                        ocall <- match.call(get("glm"), ocall)
                ocall$formula <- form
        }
        structure(list(call = ocall, terms = object$terms, coefficients = coef,
                dispersion = dispersion, df = c(p, rdf), deviance.resid =
                dresid, cov.unscaled = covun, correlation = cor, deviance =
                deviance(object), null.deviance = object$null.deviance, iter =
                object$iter, nas=nas), class = "summary.glmRob" )
}

"print.summary.glmRob" <-
function(x, digits = .Options$digits-3, quote = T, prefix = "")
{
        nas <- x$nas
        coef <- x$coef
        correl <- x$correl
        if(any(nas)) {
                nc <- length(nas)
                cnames <- names(nas)
                coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
                coef1[!nas,  ] <- coef
                coef <- coef1
                if(!is.null(correl)) {
                        correl1 <- matrix(NA, nc, nc, dimnames = list(cnames,
                                cnames))
                        correl1[!nas, !nas] <- correl
                        correl <- correl1
                }
        }
        cat("\nCall: ")
        dput(x$call)
        dresid <- x$deviance.resid
        n <- length(dresid)
        rdf <- x$df[2]
        if(rdf > 5) {
                cat("Deviance Residuals:\n")
                rq <- quantile(as.vector(dresid), na.rm=T)
                names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
                print(rq, digits = digits)
        }
        else if(rdf > 0) {
                cat("Deviance Residuals:\n")
                print(dresid, digits = digits)
        }
        if(any(nas))
                cat("\nCoefficients: (", sum(nas),
                        " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        print(coef, digits = digits)
        cat(paste("\n(Dispersion Parameter for", names(x$dispersion),
                "family taken to be", format(round(x$dispersion, digits)),
                ")\n"))
        int <- attr(x$terms, "intercept")
        if(is.null(int))
                int <- 1
        cat("\n    Null Deviance:", format(round(x$null.deviance, digits)),
                "on", n - int, "degrees of freedom\n")
        cat("\nResidual Deviance:", format(round(x$deviance, digits)), "on",
                round(rdf, digits), "degrees of freedom\n")
        cat("\nNumber of Iterations:", format(trunc(x$iter)),
                "\n")
        if(!is.null(correl)) {
                p <- dim(correl)[2]
                if(p > 1) {
                        cat("\nCorrelation of Coefficients:\n")
                        ll <- lower.tri(correl)
                        correl[ll] <- format(round(correl[ll], digits))
                        correl[!ll] <- ""
                        print(correl[-1,  - p, drop = F], quote = F, digits =
                                digits)
                 }
        }
        invisible(NULL)
}

"residuals.glmRob" <- function(obj, type = c("deviance", 
	"pearson", "working", "response"))
{
	type <- match.arg(type)
	mu <- obj$fitted
	y <- obj$y
	family <- as.family(obj)
	switch(type,
		working = obj$y - obj$ci - obj$fitted.values,
		pearson = (y - mu)/sqrt(family$variance(mu)),
		deviance = {
			w <- obj$prior.w
			if(is.null(w))
				w <- rep(1, length(mu))
			family$deviance(mu, y, w, residuals = T)
			},
		response = obj$y - obj$fitted
	)
}

plot.glmRob <- function(x, which.plots = "ask", type = "pearson",
	chisq.percent = 0.99, vertical.outlier = 0.99, smooths = F,
	rugplot = F, id.n = 3, ...)
{
	x.name <- deparse(substitute(x))
	model.list <- list(x$call)
	names(model.list) <- x.name
	x <- list(x = x)
	names(x) <- x.name
	attr(x, "model.list") <- model.list

	plot.glmfm(x, which.plots = which.plots, type = type,
		chisq.percent = chisq.percent, vertical.outlier = vertical.outlier,
		smooths = smooths, rugplot = rugplot, id.n = id.n, ...)

	invisible(x[[1]])
}



#plot.glmRob <- 
#function(x, residuals = NULL, smooths = F, rugplot = F, 
#	id.n = F, ask = T, which.plots = NULL, 
#	...)
#{
#	glm.obj <- x
#	Residuals <- resid(glm.obj, type = 'deviance')
#	tmp.x <- numeric.model.matrix(glm.obj)
#	resid.name <- "Deviance Residuals"
#	if(!is.null(residuals)) {
#		if(length(residuals) == 1 && residuals)
#			residuals <- Residuals
#		else {
#			Residuals <- residuals
#			resid.name <- deparse(substitute(residuals))
#		}
#	}
#	fits <- fitted(glm.obj)
#	preds <- predict(glm.obj)
#	response <- glm.obj$y
#	if(!is.null(glm.obj$na.action))
#		response <- nafitted(glm.obj$na.action, response)
#	form <- formula(glm.obj)
#	response.name <- deparse(form[[2]])
#	model <- deparse(form[[3]])
#	fit.lab <- paste("Fitted :", model, sep = " ")
#	pred.lab <- paste("Predicted :", model, sep = " ")
#	add.ons <- function(x, y, smooths = T, rugplot = T, id.n = 3)
#	{
#		if(smooths) {
#			which <- !is.na(x) & !is.na(y)
#			prediction <- loess.smooth(x[which], y[which], span = 1, degree = 1)
#			lines(prediction)
#		}
#		if(rugplot) {
#			jx <- jitter(x[!is.na(x)])
#			xlim <- range(jx)
#			rug(jx)
#		}
#		if(id.n) {
# Identify id.n greatest y-values (in absolute value)
#			n <- length(y)
#			id.n <- id.n + sum(is.na(y))
#			oy <- order(abs(y))
#			which <- oy[(n - id.n + 1):n]
#			text(x[which], y[which], as.character(which), adj = -0.05)
#		}
#	}
#	plotfun1 <- function(fits, Residuals, xlab, ylab, smooths, rugplot, id.n, add.ons, ...)
#	{
#
# Residuals vs Fitted Values
#		plot(fits, Residuals, xlab = xlab, ylab = ylab, ...)
#		abline(h = 0, lty = 2)
#		add.ons(fits, Residuals, smooths = smooths, rugplot = rugplot, id.n = id.n)
#	}
#	plotfun2 <- function(preds, Residuals, xlab, ylab, smooths, rugplot, id.n, add.ons, ...)
#	{
#
# Sqrt of abs(Residuals) vs Fitted Values
#		y <- sqrt(abs(Residuals))
#		plot(preds, y, xlab = xlab, ylab = ylab, ...)
#		add.ons(preds, y, smooths = smooths, rugplot = rugplot, id.n = id.n)
#	}
#	plotfun3 <- function(fits, response, xlab, ylab, smooths, rugplot, id.n, add.ons, ...)
#	{
#
# Response vs Fitted Values
#		xylims <- range(fits[!is.na(fits)], response[!is.na(response)])
#		plot(fits, response, xlab = xlab, ylab = ylab, xlim = xylims, ylim = xylims, ...)
#		abline(0, 1, lty = 2)
#		add.ons(fits, response, smooths = smooths, rugplot = rugplot, id.n = F)
#	}
#	plotfun4 <- function(glm.obj)
#	{
#
# Normal QQplot of Residuals
#		x <- resid(glm.obj, type = "pearson")
#		qqnorm(x, ylab = "Pearson Residuals")
#		qqline(x, lty = 2)
#	}
#	plotfun5 <- function(x, r) plot.glmRob.RvsRD(r,x)
#	plotfun6 <- function(glm.obj)
#	{
#
# NOT-NORMAL QQplot of Deviance Residuals
#
#		x <- resid(glm.obj, type = "deviance")
#		#qqnorm(x, ylab = "Deviance Residuals")
#		#qqline(x, lty = 2)
#
#		fam <- family(glm.obj)$family["name"]
#		if(fam=="Binomial") dist <- 0
#			else dist <- 1
#
#		resp <- glm.obj$y
#		param <- fitted(glm.obj)
#
#		#qq <- qqplot.glmRob(resp, param, dist)$quantiles
#		qq <- qqplot.glmRob(resp, param, dist)
#
#		a <- qq$quantiles
#		b <- qq$deviance.residuals
#
#		#plot(qq, sort(x), ylab = 'Deviance Residuals',
#		plot(a, b, ylab = 'Deviances',
#			xlab='Estimated quantiles') 
#		abline(0,1, lty = 2)
#
#
#	}
#	if(is.null(which.plots)) {
#		choices <- c("All", "Residuals vs Fitted Values", 
#			"Sqrt of abs(Residuals) vs Predictions", "Response vs Fitted Values", 
#			"Normal QQplot of Pearson Residuals", "RR vs RD",
#			"QQplot of Deviances")
#		choices <- substring(choices, 1, 40)	#truncate long names
#		tmenu <- paste("plot:", choices)
#		if (ask == F) which.plots <- c(3:8, 1)
#		while (T) {
#			if(ask) {
#				which.plots <- lmRob.checkbox(tmenu, title = 
#					  "\nMake a plot selection (or 0 to exit):\n")
#				if (any(which.plots == 1))
#					which.plots <- 2:6
#				which.plots <- which.plots+1
#			}
#			for (idx in 1:length(which.plots)) {
#				pick <- which.plots[idx]
#				switch(pick,
#				invisible(return(glm.obj)),
#				{
# Plot all choices one by one
#				ask.now <- F
#				}
#				,
#				{
# Residuals vs Fitted Values
#				  plotfun1(fits, Residuals, fit.lab, ylab = resid.name, smooths, 
#				    rugplot, id.n, add.ons, ...)
#				}
#				,
#				{
# Sqrt of abs(Residuals) vs Predicted Values
#				  plotfun2(preds, Residuals, xlab = pred.lab, ylab = paste(
#				    "sqrt(abs(", resid.name, "))", sep = ""), smooths, rugplot, 
#				    id.n, add.ons, ...)
#				}
#				,
#				{
# Response vs Fitted Values
#				  plotfun3(fits, response, xlab = fit.lab, ylab = response.name, 
#				    smooths, rugplot, id.n, add.ons, ...)
#				}
#				,
#				{
# Normal QQplot of Pearson Residuals
#				  plotfun4(glm.obj)
#				}
#				,
#				{
#				  plotfun5(tmp.x, Residuals)
#				}
#				,
#				{ 
# NON-Normal QQplot of Deviance Residuals
#				  plotfun6(glm.obj)
#				}
#				)
#			}
#		}
#	}
#	else {
#		for(i in which.plots) {
#			switch(i,
#				{
# Residuals vs Fitted Values
#				  plotfun1(fits, Residuals, fit.lab, ylab = if(is.null(residuals)
#				    ) "Deviance Residuals" else deparse(substitute(residuals)), 
#				    smooths, rugplot, id.n, add.ons, ...)
#				}
#				,
#				{
# Sqrt of abs(Residuals) vs Fitted Values
#				  plotfun2(preds, Residuals, xlab = pred.lab, ylab = paste(
#				    "sqrt(abs(", resid.name, "))", sep = ""), smooths, rugplot, 
#				    id.n, add.ons, ...)
#				}
#				,
#				{
# Response vs Fitted Values
#				  plotfun3(fits, response, xlab = fit.lab, ylab = response.name, 
#				    smooths, rugplot, id.n, add.ons, ...)
#				}
#				,
#				{
# Normal QQplot of Residuals
#				  plotfun4(glm.obj)
#				}
#				,
#				{
#				  plotfun5(tmp.x, Residuals)
#				}
#				,
#				{ 
# NON-Normal QQplot of Deviance Residuals
#				  plotfun6(glm.obj)
#				}
#				,
#				warning(paste("There is no plot number", i)))
#		}
#	}
#	invisible(glm.obj)
#}


# predict.glmRob
# adapted from predict.glm

predict.glmRob <-
function(object, newdata, type = c("link", "response", "terms"), 
		se.fit = F, terms = labels(object), dispersion = NULL, ...)
{
	type <- match.arg(type)
	save.na.action <- object$na.action
	object$na.action <- NULL
	out <- if(!se.fit) {
#No standard errors
		if(missing(newdata)) switch(type,
				link = object$linear.predictors,
				response = object$fitted,
				terms = { 
				  #class(object) <- c(class(object), "lm")
				  oldClass(object) <- "lm"
				  NextMethod("predict")
				}
		)		
		else switch(type, 
				response =  {
				  #class(object) <- c(class(object), "lm")
				  oldClass(object) <- "lm"
				  family(object)$inverse(NextMethod("predict"))
				}
				,
				link = {
				  type <- "response"
				  #class(object) <- c(class(object), "lm")
				  oldClass(object) <- "lm"
				  NextMethod("predict")
				}
				,
				{
				  #class(object) <- c(class(object), "lm")
				  oldClass(object) <- "lm"
				  NextMethod("predict")
				} )
	}
	else {
#With standard errors
		if(is.null(dispersion) || dispersion == 0)
			dispersion <- summary(object, dispersion = dispersion)$dispersion
		switch(type,
			response = {

				Terms <- object$terms
				offset <- attr(Terms, "offset")
				intercept <- attr(Terms, "intercept")
				xbar <- NULL
				if(missing(newdata)) 
					x <- model.matrix(object)
				else if(!((is.atomic(newdata) && length(newdata) == 1
				  && length(object$coef) != 1 && newdata > 0 
				  && (newdata - trunc(newdata) < .Machine$single.eps))
				  | is.list(newdata))) {
				  	if(!is.null(offset)) {
						warning("Offset not included")
						offset <- NULL
					}
					TT <- length(object$coef)
					if(is.matrix(newdata) && ncol(newdata) == TT)
						x <- newdata
					else if(length(newdata)==TT) 
						x <- matrix(x, 1, TT)
					else stop("Argument \"newdata\" cannot be coerced")
				} else {
					x <- model.matrix(delete.response(Terms),
						newdata, contrasts = object$contrasts,
						xlevels = attr(object, "xlevels"))
					if(!is.null(offset))
						offset <- eval(attr(Terms, "variables")[offset], 
							newdata)
				}
				coefs <- coef(object)
				asgn <- attr(coefs, "assign")
				nac <- is.na(objects$coef)
				if(any(nac)) {
				 	xbar <- xbar[!nac]
					x <- x[,!nac]
				}
				attr(x, "constant") <- xbar
				fit.summary <- summary(object)
				#pred <- Build.terms(x, coefs, ??, asgn, collapse = T)
				cov <- fit.summary$cov.unscaled
				ft <- drop(x %*% coefs)
				vv <- rowSums((x %*% cov) * x)
				pred <- list(fit = ft, se.fit = drop(sqrt(vv)))
				pred$df <- object$df.resid
				if(!is.null(offset)) {
					if(missing(newdata))
						warning("Offset not included")
					else pred$fit <- pred$fit + offset
				}
				if(missing(newdata) && !is.null(object$na.action)) {
					pred$fit <- napredict(object$na.action, pred$fit)
					pred$se.fit <- napredict(object$na.action, pred$se.fit)
				}
				famob <- family(object)
				pred$fit <- famob$inverse(pred$fit)
				pred$se.fit <- (pred$se.fit * sqrt(dispersion))
				pred$residual.scale <- as.vector(sqrt(dispersion))
				pred$se.fit <- pred$se.fit/abs(famob$deriv(pred$fit))
				pred
			}
			,
			link = {
				type <- "response"
				Terms <- object$terms
				offset <- attr(Terms, "offset")
				intercept <- attr(Terms, "intercept")
				xbar <- NULL
				if(missing(newdata)) 
					x <- model.matrix(object)
				else if(!((is.atomic(newdata) && length(newdata) == 1
				  && length(object$coef) != 1 && newdata > 0 
				  && (newdata - trunc(newdata) < .Machine$single.eps))
				  | is.list(newdata))) {
				  	if(!is.null(offset)) {
						warning("Offset not included")
						offset <- NULL
					}
					TT <- length(object$coef)
					if(is.matrix(newdata) && ncol(newdata) == TT)
						x <- newdata
					else if(length(newdata)==TT) 
						x <- matrix(x, 1, TT)
					else stop("Argument \"newdata\" cannot be coerced")
				} else {
					x <- model.matrix(delete.response(Terms),
						newdata, contrasts = object$contrasts,
						xlevels = attr(object, "xlevels"))
					if(!is.null(offset))
						offset <- eval(attr(Terms, "variables")[offset], 
							newdata)
				}
				coefs <- coef(object)
				asgn <- attr(coefs, "assign")
				nac <- is.na(objects$coef)
				if(any(nac)) {
				 	xbar <- xbar[!nac]
					x <- x[,!nac]
				}
				attr(x, "constant") <- xbar
				fit.summary <- summary(object)
				#pred <- Build.terms(x, coefs, ??, asgn, collapse = T)
				cov <- fit.summary$cov.unscaled
				ft <- drop(x %*% coefs)
				vv <- rowSums((x %*% cov) * x)
				pred <- list(fit = ft, se.fit = drop(sqrt(vv)))
				pred$df <- object$df.resid
				if(!is.null(offset)) {
					if(missing(newdata))
						warning("Offset not included")
					else pred$fit <- pred$fit + offset
				}
				if(missing(newdata) && !is.null(object$na.action)) {
					pred$fit <- napredict(object$na.action, pred$fit)
					pred$se.fit <- napredict(object$na.action, pred$se.fit)
				}
				pred
			}
			, # type = "terms"
			{
				Terms <- object$terms
				offset <- attr(Terms, "offset")
				intercept <- attr(Terms, "intercept")
				xbar <- NULL
				if(missing(newdata)) {
					x <- model.matrix(object)
					if(intercept) {
						xbar <- colMeans(x)
						x <- sweep(x, 2, xbar)
					}
				}
				else if(!((is.atomic(newdata) && length(newdata) == 1
				  && length(object$coef) != 1 && newdata > 0 
				  && (newdata - trunc(newdata) < .Machine$single.eps))
				  | is.list(newdata))) {
				  	if(!is.null(offset)) {
						warning("Offset not included")
						offset <- NULL
					}
					TT <- length(object$coef)
					if(is.matrix(newdata) && ncol(newdata) == TT)
						x <- newdata
					else if(length(newdata)==TT) 
						x <- matrix(x, 1, TT)
					else stop("Argument \"newdata\" cannot be coerced")
				} else {
					x <- model.matrix(delete.response(Terms),
						newdata, contrasts = object$contrasts,
						xlevels = attr(object, "xlevels"))
					if(!is.null(offset))
						offset <- eval(attr(Terms, "variables")[offset], 
							newdata)
				}
				if(!missing(newdata) && intercept) {
					xold <- model.matrix(object)
					xbar <- colMeans(xold)
					x <- sweep(x, 2, xbar)
				}
				coefs <- coef(object)
				asgn <- attr(coefs, "assign")
				terms <- match.arg(terms, labels(object))
				asgn <- asgn[terms]
				nac <- is.na(object$coef)
				if(any(nac)) {
					xbar <- xbar[!nac]
					x <- x[, !nac]
				}
				attr(x, "constant") <- xbar
				fit.summary <- summary(object)
				#pred <- Build.terms(x, coefs, ??, asng, collapse = F)
				cov <- fit.summary$cov.unscaled
				constant <- attr(x, "constant")
				if(!is.null(constant))
					constant <- sum(constant*coefs)
				assign <- asgn
				se <- fit <- array(0, c(nrow(x), length(assign)), 
					list(dimnames(x)[[1]], names(assign)))
				TL <- sapply(assign, length)
				simple <- TL == 1
				complex <- TL > 1
				if(any(simple)) {
						asss <- unlist(assign[simple])
						ones <- rep(1, nrow(x))
						fit[, simple] <- x[, asss]*outer(ones, coefs[asss])
						se[,simple] <- abs(x[,asss])*outer(ones, 
							sqrt(diag(cov))[asss])
				}
				if(any(complex)) {
					assign <- assign[complex]
					for(term in names(assign)) {
						TT <- assign[[term]]
						xt <- x[,TT]
						fit[, term] <- xt %*% coefs[TT]
						se[, term] <- sqrt(rowSums(((xt %*% cov[TT, TT]) * xt)))
					}
				}
				attr(fit, "constant") <- constant
				pred <- list(fit=fit, se.fit = se)
				pred$df <- object$df.resid
				if(missing(newdata) && !is.null(object$na.action)) {
						pred$fit <- napredict(object$na.action, pred$fit)
						pred$se.fit <- napredict(object$na.action, pred$se.fit)
				}
				pred
			}
			)
	}
	if(missing(newdata) && !is.null(save.na.action))
		if(is.list(out)) {
			out$fit <- napredict(save.na.action, out$fit)
			out$se.fit <- napredict(save.na.action, out$se.fit)
		}
		else out <- napredict(save.na.action, out)
	out
}


anova.glmRob <- function(object, ..., test = c("none", "Chisq", "F", "Cp"))
{
	test <- match.arg(test)

	margs <- function(...)
		nargs()

	if(margs(...))
		return(anova.glmRoblist(list(object, ...), test = test))

	Terms <- object$terms
	term.labels <- attr(Terms, "term.labels")
	nt <- length(term.labels)
	m <- model.frame(object)
	x <- model.matrix(Terms, m, contrasts = object$contrasts)
	ass <- attr(x, "assign")
	control <- object$control
	robust.cov.control <- object$robust.cov.control

	if(is.null(control)) {
		fit.method <- object$fit.method
		if(fit.method=="cubif") 
			control <- glmRob.cubif.control()
		else if(fit.method == "mallows")
			control <- glmRob.mallows.control()
		else if(fit.method == "misclass")
			control <- glmRob.misclass.control()
		else 
			stop(paste("method ",fit.method," does not exist"))
	}

	Family <- as.family(object)
	a <- attributes(m)
	y <- model.extract(m, "response")
	w <- model.extract(m, "weights")

	if(!length(w))
		w <- rep(1, nrow(m))

	offset <- attr(Terms, "offset")

	if(is.null(offset))
		offset <- 0

	else
		offset <- m[[offset]]

	dev.res <- double(nt)
	df.res <- dev.res
	nulld <- object$null.deviance

	if(is.null(nulld))
		nulld <- sum(w * (y - weighted.mean(y, w))^2)

	dev.res[1] <- nulld
	df.res[1] <- nrow(x) - attr(Terms, "intercept")

	if(nt > 1) {
		for(iterm in seq(nt, 2)) {
			x <- x[,  - (ass[[(term.labels[iterm])]]), drop = F]

			fit.call <- object$call
			fit.call[[1]] <- as.name(paste('glmRob.',object$fit.method, 
					sep = ''))
			fit.call$x <- x
			fit.call$y <- y
			fit.call$formula <- NULL
			fit.call$control <- control
			fit.call$family <- family(object)
			fit.call$offset <- offset
			fit.call$Terms <- Terms
			fit.call$null.dev <- T
			fit.call$robust.cov.control <- robust.cov.control
			fit <- eval(fit.call, sys.parent())

			dev.res[iterm] <- deviance(fit)
			df.res[iterm] <- fit$df.resid
		}
	}

	if(nt) {
		dev.res <- c(dev.res, deviance(object))
		df.res <- c(df.res, object$df.resid)
		dev <- c(NA,  - diff(dev.res))
		df <- c(NA,  - diff(df.res))
	}

	else dev <- df <- as.numeric(NA)	

	heading <- c("Analysis of Deviance Table\n", 
			paste(Family$family[1], "model\n"), 
			paste("Response: ", 
				as.character(formula(object))[2], 
				"\n", sep = ""), 
			"Terms added sequentially (first to last)")

	aod <- data.frame(Df = df, Deviance = dev, 
			"Resid. Df" = df.res, "Resid. Dev" = dev.res, 
			row.names = c("NULL", term.labels), 
			check.names = F)

	attr(aod, "heading") <- heading

	oldClass(aod) <- c("anova", "data.frame")

	if(test == "none")
		return(aod)

	else
		stat.anova(aod, test, 
			deviance.lm(object)/object$df.resid, 
			object$df.resid, nrow(x))
}



anova.glmRoblist <-
function(object, ..., test = c("none", "Chisq", "F", "Cp"))
{
	diff.term <- function(term.labels, i)
	{
		t1 <- term.labels[[1]]
		t2 <- term.labels[[2]]
		m1 <- match(t1, t2, F)
		m2 <- match(t2, t1, F)
		if(all(m1)) {
			if(all(m2))
				return("=")
			else return(paste(c("", t2[ - m1]), collapse = "+"))
		}
		else {
			if(all(m2))
				return(paste(c("", t1[ - m2]), collapse = "-"))
			else return(paste(i - 1, i, sep = " vs. "))
		}
	}
	test <- match.arg(test)
	rt <- length(object)
	if(rt == 1) {
		object <- object[[1]]
		UseMethod("anova")
	}
	forms <- sapply(object, function(x)
	as.character(formula(x)))
	subs <- as.logical(match(forms[2,  ], forms[2, 1], F))
	if(!all(subs))
		warning("Some fit objects deleted because response differs from the first model")
	if(sum(subs) == 1)
		stop("The first model has a different response from the rest")
	forms <- forms[, subs]
	object <- object[subs]
	dfres <- sapply(object, "[[", "df.resid")
	dev <- sapply(object, "[[", "deviance")
	tl <- lapply(object, labels)
	rt <- length(dev)
	effects <- character(rt)
	for(i in 2:rt)
		effects[i] <- diff.term(tl[c(i - 1, i)], i)
	ddev <-  - diff(dev)
	ddf <-  - diff(dfres)
	heading <- c("Analysis of Deviance Table", paste("\nResponse: ", forms[2, 1], "\n", sep = ""
		))
	aod <- data.frame(Terms = forms[3,  ], "Resid. Df" = dfres, "Resid. Dev" = dev, Test = 
		effects, Df = c(NA, ddf), Deviance = c(NA, ddev), check.names = F)
	aod <- as.anova(aod, heading)
	if(test != "none") {
		n <- length(object[[1]]$residuals)
		o <- order(dfres)
		stat.anova(aod, test, deviance.lm(object[[o[1]]])/dfres[o[1]], dfres[o[1]], n)
	}
	else aod
}


"family.glmRob" <- function(object)
{
  family(object$call$family)
}

"labels.glmRob" <- labels.lm

"coef.glmRob" <- coef.lm



model.frame.glmRob <-
function(formula, data = NULL, na.action = NULL, ...)
{
	m <- formula$model
	if(!is.null(m))
	return(m)
	oc <- formula$call
	oc$method <- "model.frame"
	oc[[1]] <- as.name("glmRob")
	if(length(data)) {
		oc$data <- substitute(data)
		eval(oc, sys.parent())
	}
	else eval(oc, list())
}

model.matrix.glmRob <-
function(object, ...)
{
	if(n <- match("x", names(object), 0))
		object[[n]]
	else {
		data <- model.frame(object, xlevels = attr(object, "xlevels"), ...)	#print(data)
		NextMethod("model.matrix", data = data, contrasts = object$contrasts)
	}
}

plot.glmRob.RvsRD <-
function(RR, x, method="mve", chisq.percent=0.975, 
    vertical.outlier=2.5, seed=1313)
{
#
# RR -> Robust Residuals? 
#
  x <- as.matrix(x)
  np <- ncol(x)
  tmp.n <- as.character(1:length(RR))
  set.seed(seed)
  if (is.null(method)) {
    if (ncol(x) <= 10)
      method <- "mve"
    else
      method <- "covRob"
  }
  if (method == "mcd")
    x.cov <- cov.mcd(x, print=F)
  else if (method == "mve")
    x.cov <- cov.mve(x, print=F)
  else if (method == "covRob")
    x.cov <- covRob(x)
  else
    stop("Other methods are not implemented yet.")
  x.center <- x.cov$center
  x.cov <- x.cov$cov
  RD <- sqrt(mahalanobis(x, x.center, x.cov))
 x.range <- range(RD)
  x.chisq <- sqrt(qchisq(chisq.percent, df=np))
  if (x.range[2] < x.chisq)
    x.range[2] <- x.chisq
  y.range <- range(RR)
  if (y.range[1] > -vertical.outlier) 
    y.range[1] <- -vertical.outlier
  if (y.range[2] < vertical.outlier)
    y.range[2] <- vertical.outlier
  plot(RD, RR, xlim=x.range, ylim=y.range, xlab="Robust Distance",
       ylab="Standardized Residuals")
  yp.idx <- (RR >  vertical.outlier)
  ym.idx <- (RR < -vertical.outlier)
  hjit <- (x.range[2]-x.range[1])/100
  if (any(yp.idx)) 
    text(RD[yp.idx]+hjit, RR[yp.idx], tmp.n[yp.idx], adj=0)
  if (any(ym.idx)) 
    text(RD[ym.idx]+hjit, RR[ym.idx], tmp.n[ym.idx], adj=0)
  abline(v=x.chisq, lty=2)
  abline(h= vertical.outlier, lty=2)
  abline(h=-vertical.outlier, lty=2)
  invisible()
}



glmRob.misclass.control <- function( mc.gamma = 0.01, mc.maxit = 30,
		mc.trc = F, mc.tol = 0.001, mc.initial = NULL, ...)
{
	list(mc.gamma = mc.gamma , mc.maxit = mc.maxit, mc.trc = mc.trc,
		mc.tol = mc.tol, mc.initial = mc.initial)
}

glmRob.mallows.control <- function( wt.fn = wt.carroll, wt.tuning = 8, ...)
{
	list(wt.fn = wt.fn, wt.tuning = wt.tuning)
}



glmRob.mallows <- 
	function(x, y, control, offset, null.dev, family, Terms, ...)  
{

			#fam.name <- if(is.null(family$family["name"]))
			#			family["name"]
			#		else
			#			family$family["name"]

			fam.name <- family$family["name"]

			if(fam.name != "Binomial") 
				stop(paste("Method mallows not implemented for family ", 
						fam.name)) 

			if (is.category(y)) y <- y != levels(y)[1]
              else y <- as.vector(y)

	       if( any(y>1) || is.matrix(y) ) 
		   	stop("Response doesn't look Bernoulli. Method mallows only implemented for Bernoulli responses")


			wt.fn <- control$wt.fn
			b <- control$wt.tuning
			
			x0 <- x
			tmp <- dimnames(x0)[[2]] == "(Intercept)"
			if(any(tmp)) x0 <- x0[,!tmp, drop = F]

			n <- dim(x0)[1]
			p <- dim(x0)[2]

			tmp <- cov.mcd(x0, print.it =F)
			mu <- tmp$center
			V <- tmp$cov
			x1 <- scale(x0, center = mu, scale = rep(1,p))
			
			#V <- var(x0)
			#x1 <- scale(x0, scale = rep(1,p))

			tmp2 <- solve(V) %*% t(x1)
			d1 <- diag(x1 %*% tmp2)

			d1 <- sqrt( d1 / p )

			w <- wt.fn(d1, b)

			# now use the weights


			w.glm.fit <- glm.fit(x=x, y=y, family=binomial, 
					offset = offset, w=w, null.dev = T)
			w.glm.fit$call <- match.call()
			w.glm.fit$control <- control
			w.glm.fit$prior.weights <- NULL

			## we need null.deviance here!!!


			if( any(offset) && attr(Terms, "intercept") ) {
				null.deviance <- if(length(Terms))
						glm.fit(x[, "(Intercept)", drop = F], y, w, 
							offset = offset, family = family,
						 	null.dev = NULL)$deviance
					else
				 		w.glm.fit$deviance
				w.glm.fit$null.deviance <- null.deviance
			}


			# cov matrix

			p <- dim(x)[2]
			n <- dim(x)[1]

			tmp1 <- tmp2 <- matrix(0, p, p)

			tmp3 <- x %*% w.glm.fit$coef

			for(i in 1:n) {
				tmp <- x[i,] %*% t(x[i,]) 
				tmp <- tmp * glmRob.mc.f( tmp3[i] )
				tmp1 <- tmp1 + tmp * w[i]
				tmp2 <- tmp2 + tmp * w[i] * w[i]
			}

			tmp1 <- tmp1 / n
			tmp2 <- tmp2 / n

			cov <- solve(tmp1) %*% tmp2 %*% solve(tmp1)
			cov <- cov / n

			xn <- dimnames(x)[[2]]

			dimnames(cov) <- list(xn, xn)
			w.glm.fit$cov <- cov

			c(w.glm.fit, list(mallows.weights = w))
}

wt.carroll <- function(u, b) 
	( 1 - (u/b)^2)^3 * (abs(u) <= b)


glmRob.misclass <- 
	function(x, y, control, offset,
		family, Terms, ...)  
{

			fam.name <- family$family["name"]

			if(fam.name != "Binomial") 
				stop(paste("Method mallows not implemented for family ", 
						fam.name)) 


			if (is.category(y)) y <- y != levels(y)[1]
              else y <- as.vector(y)

	       if( any(y>1) || is.matrix(y) ) 
		   	stop("Response doesn't look Bernoulli. Method misclass only implemented for Bernoulli responses")

			mc.gamma <- control$mc.gamma
			maxit <- control$mc.maxit
			mc.trc <- control$mc.trc
			mc.tol <- control$mc.tol
			initial <- control$mc.initial
		
			if(is.null(initial)) initial <- coef(glm.fit(x,y, family=binomial))


			
			n <- dim(x)[1]
			p <- dim(x)[2]
			mc.beta <- s.logistic.mc.fit(x=x, y=y, mc.gamma = mc.gamma, 
							maxit=maxit, mc.trc=mc.trc, mc.tol=mc.tol, initial=initial) 

			# now use the weights

			beta <- mc.beta$coefficients
			eta <- drop(x %*% beta)
			mu <- binomial()$inverse(eta)
			w <- glmRob.mc.w(eta, mc.gamma)


			w.glm.fit <- glm.fit(x=x, y=y, family=binomial, 
							w=w, null.dev = T)

			#w.glm.fit <- glm(y~x-1, family=binomial, w=w)
			w.glm.fit$call <- match.call()
			w.glm.fit$control <- control
			w.glm.fit$prior.weights <- NULL

			if( any(offset) && attr(Terms, "intercept") ) {
				null.deviance <- if(length(Terms))
						glm.fit(x[, "(Intercept)", drop = F], y, w, 
							offset = offset, family = family,
						 	null.dev = NULL)$deviance
					else
				 		w.glm.fit$deviance
				w.glm.fit$null.deviance <- null.deviance
			}




			# cov matrix

			p <- dim(x)[2]
			n <- dim(x)[1]

			tmp1 <- tmp2 <- matrix(0, p, p)

			tmp3 <- x %*% w.glm.fit$coef

			for(i in 1:n) {
				tmp <- x[i,] %*% t(x[i,]) 
				tmp <- tmp * glmRob.mc.f( tmp3[i] )
				tmp1 <- tmp1 + tmp * w[i]
				tmp2 <- tmp2 + tmp * w[i] * w[i]
			}

			tmp1 <- tmp1 / n
			tmp2 <- tmp2 / n

			cov <- solve(tmp1) %*% tmp2 %*% solve(tmp1)
			cov <- cov / n


			xn <- dimnames(x)[[2]]

			dimnames(cov) <- list(xn, xn)
			w.glm.fit$cov <- cov


			c(w.glm.fit, 
				list(mc.iter = mc.beta$iter,
				mc.beta = mc.beta$coefficients, mc.weights = w, 
				mc.converged = mc.beta$converged))
}



s.logistic.mc.fit <- function(x,y,mc.gamma,maxit,mc.trc,mc.tol,initial) {
		
	#mc.fit = misclassification fit
	# solves eq 2.4
	
			

	
	n <- dim(x)[1]
	p <- dim(x)[2]

	converged <- T
	
	beta1 <- initial

	if(mc.trc) cat("\n")
		
	v <- 10*mc.tol
	j <- 0
	
	# newton raphson
	
	while( (sum(abs(v)) >  mc.tol) & (j<maxit) ) { 
				
			j <- j + 1			
			beta0 <- beta1
			a <- matrix(0, p, p)
			v <- rep(0, p)
			for(i in 1:n) {
				tmp <- drop(x[i,] %*% beta0)
				tmp2 <- glmRob.mc.w(tmp, mc.gamma) * glmRob.mc.g(tmp, mc.gamma)
				tmp3 <- (tmp2 * x[i,]) %*% t(x[i,])
				a <- a + tmp3
				tmp4 <- glmRob.mc.w(tmp, mc.gamma) * x[i,] * (
						y[i] - glmRob.mc.G(tmp, mc.gamma) )
				v <- v + tmp4
			}
			a <- -a
			beta1 <- as.vector(beta0 - solve(a) %*% v)
			if(mc.trc) cat(j, beta1, "\n", sep = ' - ')
	}
	
	beta <- as.vector(beta0)
	eta <- drop(x %*% beta)
	mu <- binomial()$inverse(eta)
	
	fit <- list(coefficients = as.vector(beta0), 
			iter = j, converged = sum(abs(v)) < mc.tol )

	fit
			
}



glmRob.mc.F <- function(u) 1/(1+exp(-u))

glmRob.mc.G <-  function(u, lambda) 
				glmRob.mc.F(u) + lambda*(1-2*glmRob.mc.F(u))

glmRob.mc.w <- function(u, lambda) 
	( (1-2*lambda)*glmRob.mc.F(u)*(1 - glmRob.mc.F(u)) ) / 
			( glmRob.mc.G(u, lambda) * (1 - glmRob.mc.G(u, lambda)) )

	
glmRob.mc.f <- function(u)  exp(-u) / ((1+exp(-u))^2)

glmRob.mc.g <- function(u, lambda) 
			glmRob.mc.f(u) - 2 * lambda * glmRob.mc.f(u)


###### VY, Aug 2000
	
qqplot.glmRob<-function(y, par, dist)
{
# y is the response variable, 
# par the parameter values
# dist=0 binomial, dist=1 Poisson
# the  output  are the quantiles (out$quantiles) and 
#the  deviances residuals (out$deviance.residuals)
	n <- length(y)
	dev <- vector("numeric", n)
	if(dist == 0) {
                        par<-pmax(par,.000000000001)
                        par<-pmin(par,.999999999999)
		uu <- glmRob.binom.dev(par)
		dev <- y * sqrt(-log(par)) - (1 - y) * sqrt(-log(1 - par))
                        dev<-sort(dev)
	}
	if(dist > 0) {
                        par<-pmax(par,.000000000001)
		uu <- glmRob.poiss.dev(par)
		v1 <- (y == 0)
		v2 <- (y > 0)
		v3 <- 1:n
		v1 <- v3[v1]
		v2 <- v3[v2]
		dev[v1] <- -sqrt(par[v1])
		dev[v2] <- par[v2] - y[v2] * log(par[v2]) -
y[v2] + y[v2] * log(y[v2])
                        dev[v2]<-sqrt(dev[v2])*sign((y[v2]-par[v2]))
                        dev<-sort(dev)
	}
            dev<-sqrt(2)*dev
	uu <- glmRob.group.values(uu[[1]], uu[[2]])
	q <- glmRob.quant(uu[[1]], uu[[2]],n)
	out <- list(q, dev)
            names(out)<-c("quantiles","deviance.residuals")
	out
}
 
glmRob.binom.dev<-function(p)
# this function computes the posibles  desviances  and its probabilities
#  when the responses  have distribution binomial with parameter p
#  The output is out: first element are the  the values  of the deviances 
#  (there may be
# repeated elements) 
# and the probabilities the second element. 
{p<-pmax(p,.000000000001)
p<-pmin(p,.999999999999)
n<-length(p)
z1<-sqrt(-log(p))
u1<-p/n
z2<-sqrt(-log(1-p))
u2<-(1-p)/n
d<-sqrt(2)*c(z1,-z2)
p<-c(u1,u2)
out<-list(d,p)
out}

glmRob.poiss.dev<-function(lan)
# this function computes the posibles  desviances  and its probabilities
#  when the responses  have distribution Poisson with parameter lan
#  The output is out: first element are the  the values  of the deviances 
#  (there may be
# repeated elements) 
# and the probabilities the second element.
# Warning if some elements
# of lan are very large the computing time may be very large
{lan<-pmax(lan,.000000000001)
n<-length(lan)
qp1<-qpois(.0001,lan)
qp1<-pmax(qp1-1,0)
qp11<-pmax(qp1,1)
qp2<-qpois(.9999,lan)+1

# qp[i] gives the number of values of the Poisson distribution that
# will be consider for the i value of the parameter

y1<-qp11[1]:qp2[1]

# In d we acumulate all he values of the deviances
# In pp we acumulate all the probabilities

d<-lan[1]-(log(lan[1])*y1)-y1+(y1*log(y1))
d<-sqrt(d)*sign(y1-lan[1])
if (qp1[1]>0)
{
    pp<-dpois(y1,lan[1])
}
if (qp1[1]==0)
{
    d<-c(-sqrt(lan[1]),d)
    pp<-dpois(c(0,y1),lan[1])
                            } 
for(i in 2:n)
{
  y1<-qp11[i]:qp2[i]
  d1<-lan[i]-(log(lan[i])*y1)-y1+(y1*log(y1))
  d1<-sqrt(d1)*sign(y1-lan[i])

   if(qp1[i]>0)
{
     pp1<-dpois(y1,lan[i])
}
   if(qp1[i]==0)
{      
      d1<-c(-sqrt(lan[i]),d1)
      pp1<-dpois(c(0,y1),lan[i])
}
  d<-c(d,d1)
  pp<-c(pp,pp1)
 
}
d<-sqrt(2)*d
ppp<-sum(pp)
pp<-pp/ppp
out<-list(d,pp)
out}


glmRob.quant<- function(x,p,n)
# Given a discrete distribution x the values (ordered and different)
# with probabilities p, this subroutine compute the quantiles 
# (i-.5)/n, i=1,...,n

{
z<-vector('numeric',n)
pp<-vector('numeric',n)
xx<-sort(x)
u<-order(x)
po<-p[u]
j<-0
aa<-0
for( i in 1:n)
{
    hh<-0
    while(hh==0)
    {
      j<-(j+1)
      pp[j]<-aa+po[j]
      aa<-pp[j]
      hh<-(pp[j]>=((i-.5)/n))
      }
    z[i]<-xx[j]
    }
z}





glmRob.group.values<-function(x,p)
# given a vector x of possible no ordered values (may be repeated values)
#  with probabilities p the function group 
# the  same values adding their probabilities
# the output is out with first element the ordered values (all
# differents) and second element
# the probabilities
{n<-length(x)
xxx<-vector("numeric", n)
pp<-vector("numeric",n)
xx<-sort(x)
u<-order(x)
po<-p[u]
u<-diff(xx)
v<-(abs(u)>.0000000001)
r<-1:(n-1)
vv<-r[v]
h<-length(vv)+1
xxx[1]<-xx[1]
pp[1]<-sum(po[1:vv[1]])
for (i in 2:(h-1))
{
  l<-vv[i-1]+1
  xxx[i]<-xx[l]
  pp[i]<-sum(po[(vv[i-1]+1):vv[i]])
}
pp[h]<-sum(po[(vv[h-1]+1):n])
xxx[h]<-xx[n]
out<-list(xxx[1:h],pp[1:h])
out
}

##### end VY, Aug 2000


## April 2001
##
## glmRob.modified.cubif
##

glmRob.cubif <- function(x, y, intercept=F, offset = 0, family = binomial(), null.dev = T,
  control, robust.cov.control, ...) 
{
	family <- as.family(family)
	nn <- family$family["name"]
  
	if(nn == "Gaussian")
    stop("Use lmRob(formula, ...) for the gaussian case")

  ics <- 0

  if(nn == "Binomial") {
    ics <- 2 
    if(is.matrix(y)) {

      if(dim(y)[2]>2)
        stop("only binomial response matrices (2 columns)")

      ni <- as.vector(y %*% c(1,1))
      y <- y[,1]
      ly <- length(y)
    }

    else {
      ics <- 1

      if(is.category(y))
        y <- y != levels(y)[1]
      else
        y <- as.vector(y)

      if(any(y>1)) stop("Response doesn't look Binomial")
        ly <- length(y)
        ni <- rep(1,ly)
    } 
  }

  if(nn == "Poisson") {

    if(!is.null(dim(y)[2]))
      stop("Poisson response cannot be a matrix")

    eval(family$initialize)

    ics <- 3 
    ly <- length(y)
		ni <- rep(1,ly)
  }

  if(ics == 0)
    stop(paste(nn, "Family not implemented in glmRob yet", sep=" "))

  eta <- ci <- ai <- rsdev <- y 
	yy  <- y 
	dn  <- dimnames(x)
	xn  <- dn[[2]]
	yn  <- dn[[1]]

	if(intercept)
    x <- cbind(1,x)

	if(length(offset) <= 1)
    offset <- rep(0, ly) 

	p <- ncol(x)
	n <- nrow(x)

  zero <- rep(F, n)

#
# Initializations
#

	qrx <- qr(x)[c("qr", "rank", "pivot", "qraux")]
	rank <- qrx$rank
	piv <- 1:rank
	epsilon <- control$epsilon
	tua <- control$epsilon
	maxit <- control$maxit
	mxt <- maxit
	mxf <- mxt
	gma <- control$gma
	iug <- control$iug
	ipo <- control$ipo
	ilg <- control$ilg
	icn <- control$icn
	icv <- control$icv
	ia1 <- control$ia1
	trc <- control$trc

	gma <- 1
	iug <- 1
	ipo <- 1
	ilg <- 1
	icn <- 1
	icv <- 1
	ia1 <- 1

	tmp <- list(...)$singular.ok
	if (!is.null(tmp)) singular.ok <- tmp else singular.ok <- F
	tmp <- list(...)$qr.out

	if(!is.null(tmp))
    qr.out <- tmp
  else
    qr.out <- F

	if(rank < p) {
 
    if(!singular.ok)
      stop(paste("x is singular, rank(x)=", rank))

    else {
      piv <- qrx$pivot[1:rank]
      x <- x[, piv, drop=F]

      if (any(zero))
        x0 <- x0[,piv,drop=F]

      xn <- dimnames(x)[[1]] 
    }
  }

	bpar <- control$bpar
  cpar <- control$cpar

# Deviance for the model reduced to the constant term.

	if(null.dev)
    null.dev <- glmRob.cubif.null(x, y, ni, offset, ics, family, control, robust.cov.control)

  else
    null.dev <- NULL 

# Initial theta, A (A0) and c (c0)

	zi <- glmRob.cubif.Ini(X = x, y = y, ni = ni, icase = ics, offset = offset, b = bpar,
    zmin = 0.001, epsilon = epsilon, ia1 = ia1, control = robust.cov.control)

	ci <- zi$ci
	cov <- zi$cov
	ai <- zi$ai
	theta0 <- theta <- zi$theta

	if(icn != 1) 
    stop("icn !=1, in glmRob.cubif")

	if(trc)
    cat("\nFull model.\n")

	sqcov <- rep(1,p)

	for(i in 1:p) {
    j <- (i*(i+1))/2; 
    sqcov[i] <- sqrt(cov[j])
  }

# Iterations

	nit <- 1
	repeat {

# theta-step

    zt <- glmRob.gytstp(x, y, ci, theta, ai, cov, ni, offset, tol = epsilon,
      icase = ics, maxit = mxt)

    theta  <- zt$theta[1:p]
    vtheta <- zt$vtheta
    nitt <- zt$nit

# Check convergence

			if(nit == maxit)
        break

			delta <- theta - theta0

			if(all(epsilon*sqcov - abs(delta) > 0))
        break

			theta0 <- theta

# c-step

    zc <- glmRob.gicstp(icase = ics, ialg = ilg, ni = ni, vtheta = vtheta,
      ai = ai, oi = offset, tol = epsilon, maxit = mxt)

    ci <- zc$ci
    nit <- nit+1
	}

# End Iterations

# Final covariance matrix of estimated coefficients

	z <- glmRob.gfedca(vtheta, ci, ai, ni, offset, icase = ics)
	zc <- glmRob.ktaskw(x = x, d = z$d, e = z$e, f = 1/n, f1 = 1, iainv = 0,
    ia = ia1, tau = epsilon)

  covf <- zc$cov
	zf <- list(theta = theta, ai = ai, vtheta = vtheta, ci = ci, cov = covf,
    nit = nit, delta = delta)

	coefs <- zf$theta

#  Residual Deviance

	zd <- glmRob.glmdev(y, ni, zf$ci, zf$ai, zf$vtheta, offset, icase = ics)

	cov <- matrix(0, nrow=rank, ncol=rank)
	i2 <- 0

  for(i in 1:rank) {
    i1 <- i2 + 1
    i2 <- i1 + i - 1
    cov[i,1:i] <- zf$cov[i1:i2]
    cov[1:i,i] <- zf$cov[i1:i2]
  }


  xn <- dimnames(x)[[2]]
	names(coefs) <- xn
	dimnames(cov) <- list(xn,xn)
	asgn <- attr(x, "assign")

# Compute eta, mu and residuals.

	ind <- 1:n
	dni <- c(ind[1],diff(ind)) 
	iii <- cumsum(dni[dni!=0])
	jjj <- (1:n)*(dni!=0) 
	eta[iii] <- zf$vtheta[jjj]

  if(any(offset))
    offset[iii] <- offset[jjj]

  ci[iii] <- zf$ci[jjj]
	ai[iii] <- zf$ai[jjj] 
	rsdev[iii] <- zd$li[jjj] - zd$sc[jjj]
	dni <- ni 
	ni <- rep(1,length(eta))
	ni[iii] <- dni[jjj] 
	mu <- family$inverse(eta+offset)
	names(eta) <- yn

	if(rank < p) {
    coefs[piv[ - seq(rank)]] <- NA
    pasgn <- asgn
    newpos <- match(1:p, piv)
    names(newpos) <- xn

    for(j in names(asgn)) {
      aj <- asgn[[j]]
      aj <- aj[ok <- (nj <- newpos[aj]) <= rank]

      if(length(aj)) {
        asgn[[j]] <- aj
        pasgn[[j]] <- nj[ok]
      }
      else
        asgn[[j]] <- pasgn[[j]] <- NULL
    }

    cnames <- xn[piv]
	}
  
	new.dev <- family$deviance(mu, yy/ni, w=rep(1,n))
  resp <- yy/ni - mu
	names(ai) <- yn
	names(mu) <- yn
	df.residual <- ly - rank

	fit <- list(coefficients = coefs, fitted.values = mu, ci = ci, rank = rank,
    assign = asgn, df.residual = df.residual, control = control)

	if(rank < p) {
    if(df.residual > 0)
      fit$assign.residual <- (rank + 1):n
	}
  
	if(qr.out)
    fit$qr <- qrx

	fit$ai  <- ai
	fit$cov <- cov
	fit$ni <- ni
	rsdev <- sign(y-mu)*sqrt(2*abs(rsdev)) 
	fit$weights <- pmin(1, fit$ai / abs(yy/ni-ci-mu) )

  this.call <- match.call()
	names(fit$fitted) <- names(eta)
  
	c(fit, list(family = family$family, ics=ics, linear.predictors = eta, deviance = new.dev,
    null.deviance = null.dev, call = this.call, iter = nit, y = yy/ni,
    contrasts = attr(x,"contrasts"), rsdev = rsdev, gradient = zf$grad, inv.hessian = zf$hessnv,
    residuals = yy/ni - ci - mu, wa = zf$ai, vtheta = zf$vtheta))
}



# glmRob( y ~ x, fit.method='modified.cubif')


glmRob.cubif.Ini <- function(X,y,ni,icase,
		offset,b,zmin,epsilon,ia1, control)
{

	# Compute distances (dist), weights (ai), initial coefficients,
	# intial ci.  If an intercept is present, assume that the first
	# column of X is made of ones

	n <- length(y)
	np <- dim(as.matrix(X))[2]

	# Mahalanobis distances and weights

	if(np == 1 & all(X == 1)) {
		dist <- rep(1, n)
		ai <- rep(b, n)
	}
  
	else {
  
		if(all(X[, 1] == 1)) {
			Z <- as.matrix(X[, -1])
			Zmcd <- covRob(Z, estim = control$estim, control = control, distance = F)
			Minv <- solve(Zmcd$cov)
			Zc <- sweep(Z, 2, Zmcd$center)
		}

		else {
			Zc <- as.matrix(X)
			Zmcd <- covRob(Zc, estim = control$estim, control = control, distance = F)
			mu <- as.matrix(Zmcd$center)
			Mu <- Zmcd$cov + mu %*% t(mu)
			Minv <- solve(Mu)
		}

		ai <- dist <- rep(1, n)
		for(i in 1:n) {
			z <- as.matrix(Zc[i,  ])
			dist[i] <- sqrt((t(z) %*% Minv) %*% z)
			ai[i] <- b/max(zmin, dist[i])
		}
	}

#	# Mahalanobis distances and weights
#
#	dist <- rep(1,n)
#	ai <- rep(b,n)
#
#	if(!(dim(X)[[2]] == 1 || all(X == 1))) {
#
#		# remove interaction columns from X	#
#		interaction.cols <- attr(X, "order") > 1
#		interaction.cols <- (1:length(interaction.cols))[interaction.cols]
#		if(length(interaction.cols))
#			Z <- X[, -interaction.cols, drop = F]
#		else
#			Z <- X
#
#		# remove factor columns from Z #
#		factor.cols <- names(attr(X, "contrasts"))
#		good.cols <- setdiff(dimnames(Z)[[2]], factor.cols)
#		Z <- Z[, good.cols, drop = F]
#
#		if(length(good.cols)) {
#
#			if(dimnames(Z)[[2]][1] == "(Intercept)") {
#				intercept <- 1
#				Zc   <- as.matrix(Z[, -1, drop = F])
#			}
#
#			else {
#				intercept <- 0
#				Zc   <- as.matrix(Z)
#			}
#
#			Zcov <- covRob(Zc, estim = control$estim, control = control)
#			dist <- ((np - intercept) / dim(Z)[2])*Zcov$dist
#			ai <- dist
#			ai[ai < zmin] <- zmin
#			ai <- b / ai
#		}
#	}

	# Initial value of theta and vartheta; set c_i=0 #

	zi <- glmRob.Initial.LMS(X,y,ni,dist,offset,icase) 

	theta0 <- zi$theta0[1:np]
	ci <- rep(0,n)
	vtheta <- as.vector(X %*% theta0)

	# Initial covariance matrix of estimated theta #

	z <- glmRob.gfedca(vtheta,ci,ai,ni,offset,icase)
	zc <- glmRob.ktaskw(x = X, d = z$d, e = z$e, f = 1/n, f1 = 1,
											iainv = 0, ia = ia1, tau = epsilon)

	covi <- zc$cov

	list(theta = theta0, ci = ci, cov = covi, dist = dist, ai = ai)
}


glmRob.Initial.LMS <- function(X,y,ni,dist,offset,icase) {
		# LMS estimate using transformed y-s. 
		n <- length(y); 
		np <- ncol(X)
		if (length(offset)==0) 
			offset <- rep(0,n)
		# Transform the y-s
		if (icase!=3) {
			nip1  <- ni+1
			sy    <- (y+0.5)/nip1 
			ytild <- log(sy/(1-sy))-offset
		} else {
			sy    <- y
			sy[y<=0] <- 0.5
			ytild <- log(sy)-offset
		}
		# Initial estimate of theta (LMS)
		if (all(X[,1]==1)) {
			if( np ==1) {
			theta0 <- median(ytild)
			sigma0 <- mad(ytild)
			return(list(theta0=theta0, sigma0 =sigma0))
		} else {
			XX <- X[,-1]; 
			itc <- T
		} } else {
			XX <- X; 
			itc <- F
		}
		zt0    <- lmsreg(XX,ytild,intercept=itc)
		theta0 <- zt0$coef;
		sigma0 <- zt0$scale
		list(theta0=theta0,sigma0=sigma0)
	}



"glmRob.cubif.null" <- 
function(x, y, ni, offset, ics, family, control, robust.cov.control)
{
        tua <- epsilon <- control$epsilon
        mxf <- mxt <- maxit <- control$maxit
		trc <- control$trc
		gma <- 1
		iug <- 1
		ipo <- 1
		ilg <- 1
		icn <- 1
		icv <- 1
		ia1 <- 1
        rank <- 1

		bpar <- control$bpar
        cpar <- control$cpar

        n <- ly   <- length(y)
        w    <- rep(1,ly)
        ai   <- ci <- rep(0,ly)
        intl <- attr(x, "term.labels")
        int  <- if(is.null(intl)) 
			F 
		else 
			as.logical(match(intl[1], c("(Int.)",
                        "(Intercept)"), F))
        if(!int) {
		  eta <- rep(0,ly)
		  mu <- family$inverse(eta+offset) 
          cval <- 0.5
		  if (ics==3) cval <- 1
          ci <- rep(cval,ly)
		  ai <- rep(9999.,ly) 
		} else 
		{
          X    <- as.matrix(x[,1])

	 zi   <- glmRob.cubif.Ini(X=X, y=y, ni=ni, icase=ics,
	 	offset=offset, b=bpar, zmin = 0.001, epsilon = epsilon,
		ia1 = ia1, control = robust.cov.control)
	 ci <- zi$ci
	 cov <- zi$cov
	 ai <- zi$ai
	 theta <- theta0 <- zi$theta

	 p <- 1

		sqcov  <- rep(1,p)
		for (i in 1:p) {
			j <- (i*(i+1))/2; 
			sqcov[i] <- sqrt(cov[j])
		} 

# Iterations
		nit <- 1
		repeat {
# theta-step
			zt     <- glmRob.gytstp(X,y,ci,theta,ai,cov,
					ni,offset,tol=epsilon,icase=ics,maxit=mxt)
			theta  <- zt$theta[1:p]
			vtheta <- zt$vtheta
			nitt <- zt$nit
# Check convergence
			if (nit==maxit) break
			delta  <- theta-theta0
			if (all(epsilon*sqcov-abs(delta) >0)) break
			theta0 <- theta
# c-step
			zc     <- glmRob.gicstp(icase=ics,ialg=ilg,ni=ni,
						vtheta=vtheta,
						ai=ai,oi=offset,tol=epsilon,maxit=mxt)
			ci     <- zc$ci
			nit    <- nit+1
		}
#
# end of iterations
#
# Final covariance matrix of estimated coefficients
		z      <- glmRob.gfedca(vtheta,ci,ai,ni,offset,icase=ics)
		zc     <- glmRob.ktaskw(x=X,d=z$d,e=z$e,f=1/n, f1=1, iainv=0,
					ia=ia1, tau=epsilon)
		covf   <- zc$cov

		zf <- list(theta=theta,ai=ai,vtheta=vtheta,ci=ci,
			cov=covf,nit=nit,delta=delta)

#          z    <- glmRob.gintac(X, y, ni, offset, icase = ics, 
#					tolt=10*epsilon,
#					tola=10*epsilon, b = upar, c = cpar, maxtt=mxt, 
#					maxta=mxf)
#          t0   <- z$theta[1] 
#	   	  A0 <- z$a
#	   	  c0 <- z$ci
#          wa   <- upar/pmax(1.e-4,z$dist)
#          vtheta <- rep(t0,ly)
#          z    <- glmRob.gfedca(vtheta, c0, wa, ni, offset, ics)
#          zc   <- glmRob.ktaskw(X, z$d, z$e, f=1/ly, f1=1, iainv=0,
#				ia = ia1, tau = epsilon)
#          covi <- zc$cov
#          if (icn != 1) covi <- 1/covi
#          zf   <- glmRob.gymain(X, y, ni, covi, A0, t0, offset, 
#				b=upar, gam=gma, 
#				tau=epsilon, icase=ics, iugl=iug, iopt=ipo, 
#				ialg=ilg, icnvt=icn,
#				icnva=icv, maxit=maxit, maxtt=mxt, maxta=mxf, 
#				maxtc=mxt, 
#				tol=epsilon, tolt=10*epsilon, tola=10*epsilon, tolc=10*epsilon,
#				trc=trc)
#           ai   <- zf$wa
#           ci   <- zf$ci

           eta  <- zf$vtheta
		   #print(eta)
           mu   <- family$inverse(eta+offset)
        } # end if(!int...
		#print(mu)
		family$deviance(mu, y/ni, w)
}	




#			zt     <- glmRob.gytstp(X,y,ci,theta,ai,cov,
#					ni,offset,tol=epsilon,icase=ics,maxit=mxt)

"glmRob.gytstp"<-
function(x, y, ci, theta, wa, cov, ni, offset, 
	tol, 
	icase, 
	maxit)
{

	gam <- 1
	tau <- tol
	iopt <- 1
	icnv <- 1
	nitmon <- 0

	n <- length(y)
	np <- ncol(x)
	mdx <- nrow(x)
	ncov <- length(cov)
	if(length(offset) == 1)
		offset <- rep(0, n)
	nit <- integer(1)
	q0 <- double(1)
	delta <- double(np)
	f0 <- double(n)
	f1 <- double(n)
	f2 <- double(n)
	vtheta <- double(n)
	grad <- double(np)
	hessnv <- double(ncov)
	rw1 <- double(5 * np)
	rw2 <- matrix(double(1), mdx, np)
	iw1 <- integer(np)
	f.res <- .Fortran("s_gytstp",
		x = as.double(x),
		y = as.double(y),
		ci = as.double(ci),
		theta = as.double(theta),
		wa = as.double(wa),
		cov = as.double(cov),
		ni = as.integer(ni),
		oi = as.double(offset),
		n = as.integer(n),
		np = as.integer(np),
		mdx = as.integer(mdx),
		ncov = as.integer(ncov),
		gam = as.double(gam),
		tol = as.double(tol),
		tau = as.double(tau),
		iopt = as.integer(iopt),
		icase = as.integer(icase),
		icnv = as.integer(icnv),
		maxit = as.integer(maxit),
		nitmon = as.integer(nitmon),
		nit = as.integer(nit),
		q0 = as.double(q0),
		delta = as.double(delta),
		f0 = as.double(f0),
		f1 = as.double(f1),
		f2 = as.double(f2),
		vtheta = as.double(vtheta),
		grad = as.double(grad),
		hessnv = as.double(hessnv),
		rw1 = as.double(rw1),
		rw2 = as.double(rw2),
		iw1 = as.integer(iw1))
	list(theta = f.res$theta, nit = f.res$nit, q0 = f.res$q0, delta = f.res$
		delta, f0 = f.res$f0, f1 = f.res$f1, f2 = f.res$f2, vtheta = 
		f.res$vtheta, grad = f.res$grad, hessnv = f.res$hessnv)
}

glmRob.gicstp <- function(icase,ialg,ni,vtheta,ai,oi,tol,maxit){
# This is the interface to an existing auxiliary robeth subroutine

	n <- length(vtheta)
	if (length(oi)==1) oi <- rep(0,n)
	ci <- double(n)
	f.res<-.Fortran("s_gicstp",
		icase=as.integer(icase),
		ialg=as.integer(ialg),
		nn=as.integer(ni),
		vtheta=as.double(vtheta),
		wa=as.double(ai),
		oi=as.double(oi),
		n=as.integer(n),
		tol=as.double(tol),
		maxit=as.integer(maxit),
		c=as.double(ci))
	list(ci=f.res$c)
}










