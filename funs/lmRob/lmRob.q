#----------------------------------------------------------------#
# Robust Linear Regression                                       #
# Author: Jeffrey Wang & Kjell Konis                             #
# Date:   02/06/2002                                             #
# Insightful Corporation                                         #
#----------------------------------------------------------------#

#----------------------------------------------------------------#
# LM-like functions for lmRob                                    #
#----------------------------------------------------------------#

"lmRob" <- function(formula, data, weights, subset,
		na.action, model = F, x = F, y = F, contrasts = NULL,
		nrep = NULL, robust.control = lmRob.robust.control(...),
		genetic.control = NULL, ...)
{
  call <- match.call()
  m <- match.call(expand = F)
  m$model <- m$x <- m$y <- m$contrasts <- m$nrep <- NULL
  m$robust.control <- m$genetic.control <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
  weights <- model.extract(m, weights)
  Y <- model.extract(m, response)

  ##	In this section we make the model.matrix X and the sub
  ##	model.matices X1 and X2 such that X1 contains all of the
  ##	columns of X that correspond to factor variables and X2
  ##	contains all the columns of X that correspond to numeric
  ##	variables.  x1.idx is the (column) indicies of X that
  ##	are in X1.
  
  factor.names <- names(m)[unlist(lapply(m, is.factor))]
  X <- model.matrix(Terms, m, contrasts)
  asgn <- attr(X, "assign")
  x1.idx <- unlist(asgn[factor.names])
  names(x1.idx) <- NULL
  X1 <- X[, x1.idx, drop = F]
  X2 <- X[, setdiff(1:(dim(X)[2]), x1.idx), drop = F]

  ##	If there are both factor and numeric variables
  ##	then put the intercept term in X1

  if(dim(X1)[2] > 0 && dim(X2)[2] > 1 && dimnames(X2)[[2]][1] == "(Intercept)") {
    X1 <- cbind(X2[, 1], X1)
    dimnames(X1)[[2]][1] <- "(Intercept)"
    X2 <- X2[, 2:(dim(X2)[2]), drop = F]
    x1.idx <- c(1, x1.idx)
  }
  
  ##	If the only column in X2 is the intercept
  ##	move it to X1 and set X2 to NULL.
  
  else if(dim(X2)[2] == 1 && dimnames(X2)[[2]] == "(Intercept)") {
    X1 <- X
    x1.idx <- 1:(dim(X)[2])
    X2 <- NULL
  }
  
  ##	If X1 is empty then set it to NULL.
  
  if(length(dim(X1)) && dim(X1)[2] == 0) {
    X1 <- NULL
    x1.idx <- NULL
  }

  ##	If X2 is empty then set it to NULL.
  
  if(length(dim(X2)) && dim(X2)[2] == 0)
    X2 <- NULL

  if (nrow(X) <= ncol(X)) 
    stop("Robust method is inappropriate: not enough observations.")

  if(length(weights))
    fit <- lmRob.wfit(X, Y, weights, x1=X1, x2=X2, x1.idx=x1.idx, 
                         nrep=nrep, robust.control=robust.control,
                         genetic.control=genetic.control, ...) 
  else
    fit <- lmRob.fit(X, Y, x1=X1, x2=X2, x1.idx=x1.idx, nrep=nrep,
                        robust.control=robust.control, 
                        genetic.control=genetic.control, ...)
  if(is.null(fit)) 
    return(NULL)
  fit$terms <- Terms
  fit$call <- call

  x.names <- dimnames(X)[[2]]
  pasgn <- asgn
  qrx <- qr(X)
  Rk <- qrx[["rank"]]
  piv <- qrx[["pivot"]][1:Rk]
  newpos <- match(1:Rk, piv)
  if(length(x.names)) 
    names(newpos) <- x.names
  for(j in names(asgn)) {
    aj <- asgn[[j]]
    aj <- aj[ok <- (nj <- newpos[aj]) <= Rk]
    if(length(aj)) {
      asgn[[j]] <- aj
      pasgn[[j]] <- nj[ok]
    }
    else 
      asgn[[j]] <- pasgn[[j]] <- NULL
  }

  effects <- X * matrix(fit$coefficients, byrow=T, nrow=nrow(X), ncol=ncol(X))
  fit$effects <- sqrt(colSums(effects^2))
  fit$assign <- asgn
  if (model) 
    fit$model <- m
  if (x)
    fit$x <- X
  if (y)
    fit$y <- Y
  attr(fit, "na.message") <- attr(m, "na.message")
  if (!is.null(attr(m, "na.action")))
    fit$na.action <- attr(m, "na.action")
  fit
}


lmRob.fit <- function(x, y, x1 = NULL, x2 = NULL, x1.idx = NULL, nrep = NULL,
							robust.control = NULL, genetic.control = NULL, ...)
{
  if(!is.numeric(x))
    stop("model matrix must be numeric")
  if(!is.numeric(y))
    stop("response must be numeric")
  fit <- lmRob.fit.compute(x2, y, x1=x1, x1.idx=x1.idx, nrep=nrep,
                           robust.control=robust.control, 
                           genetic.control=genetic.control, ...)
  if (is.null(fit)) return(NULL)
  fit$contrasts <- attr(x, "contrasts")
  fit
}


"lmRob.wfit" <- function(x, y, w, x1=NULL, x2=NULL, x1.idx=NULL, nrep=NULL, 
                         robust.control=NULL, genetic.control=NULL, ...)
{
  if(!is.numeric(x))
    stop("model matrix must be numeric")
  if(!is.numeric(y))
    stop("response must be numeric")
  if(any(w < 0))
    stop("negative weights not allowed")
  contr <- attr(x, "contrasts")
  zero <- w == 0

  if(any(zero)) {
    pos <- !zero
    r <- f <- y
    ww <- w
    x0 <- x[zero,  , drop = F]
    y0 <- y[zero]
    x2 <- x2[pos,  , drop = F]
    if (!is.null(x1))
      x1 <- x1[pos,  , drop = F]
    y <- y[pos]
    w <- w[pos]
  }

  w.factor <- sqrt(w)
  if(!is.null(x2))
    x2 <- x2 * w.factor
  if(!is.null(x1))
    x1 <- x1 * w.factor
  y <- y * w.factor

  fit <- lmRob.fit.compute(x2, y, x1=x1, x1.idx=x1.idx, nrep=nrep,
                           robust.control=robust.control,
                           genetic.control=genetic.control, ...)

  if(is.null(fit)) 
    return(NULL)
  fit$residuals <- fit$residuals/w.factor
  fit$fitted.values <- fit$fitted.values/w.factor

  if(any(zero)) {
    nas <- is.na(fit$coef)
    if(any(nas))
      f0 <- x0[, !nas] %*% fit$coef[!nas]
    else 
      f0 <- x0 %*% fit$coef
    r[pos] <- fit$resid
    f[pos] <- fit$fitted
    r[zero] <- y0 - f0
    f[zero] <- f0
    fit$residuals <- r
    fit$fitted.values <- f
    w <- ww
  }

  fit$weights <- w
  fit$contrasts <- contr
  fit
}


lmRob.fit.compute <- function(x2, y, x1=NULL, x1.idx=NULL, nrep=NULL,
                       robust.control=NULL, genetic.control=NULL, ...)
{

##
## Step 0. Set Control Parameters
##

  tua <- robust.control$tua
  tlo <- robust.control$tlo
  mxr <- robust.control$mxr
  mxf <- robust.control$mxf
  mxs <- robust.control$mxs
  tl  <- robust.control$tl
  psi <- robust.control$weight
  sed <- robust.control$seed
  ini <- robust.control$initial.alg
  fnl <- robust.control$final.alg
  est <- robust.control$estim
  eff <- robust.control$efficiency
  trc <- robust.control$trace

  if(is.null(x1)) {
    x <- as.matrix(x2)
    xnames <- dimnames(x)[[2]]
    n <- nrow(x)
    p <- p2 <- ncol(x)
  }

  else {
    x1 <- as.matrix(x1)
    p1 <- ncol(x1)

    if(is.null(x1.idx))
      x1.idx <- 1:p1

    if (!is.null(x2)) {
      x2 <- as.matrix(x2)
      p2 <- ncol(x2)
    }

    else
      p2 <- 0

    p <- p1+p2
    n <- nrow(x1)
    x <- matrix(NA, nrow=n, ncol=p)
    xnames <- rep(NA, p)
    x[, x1.idx] <- x1
    xnames[x1.idx] <- dimnames(x1)[[2]]

    if (p2) {
      x[,-x1.idx] <- x2
      xnames[-x1.idx] <- dimnames(x2)[[2]]
    }

  }

  if(casefold(ini) == "auto") {
    if((n <= 250 && p == 2) || (n <= 80 && p == 3)) 
      ini <- "exhaustive"
    else if (p2 > 15)
      ini <- "fast"
    else 
      ini <- "random"
  }

  if(casefold(ini) == "random") {
    iop <- 2

    if(is.null(nrep)) 
      nrep <- round(4.6*2^p2)

    if (p2 > 15 && trc) {
      if (exists("guiDisplayDialog")) {
        guiDisplayDialog("Function", "confirmLmRob", bModal=T)
        if (yes.p == F)
          return(NULL)
      }

      else {
        yesorno <- c("Yes", "No")
        title <- paste("Random resampling may take a long time.",
                       "Do you want to continue?")
        wmenu <- menu(yesorno, title=title)
        if (wmenu != 1) 
          return(NULL)
      }
    }
  }

  else if(casefold(ini) == "exhaustive") {

    if (n > 300 || p2 > 10)
      stop("The data set is large. ",
           "Try using random resampling or the fast procedure.")

    iop <- 3
    nrep <- choose(n,p2)
  }

  else if(casefold(ini) == "genetic")
    iop <- 4

  else if(casefold(ini) == "fast")
    iop <- 6

  else
    stop("Invalid choice of resampling method.")

  if (p2 == 0)
    iop <- 5

  msg.UCV <- paste("Sum(psi.weight(wi)) less than",tl,"in lmRob.ucovcoef.")

  if(match("(Intercept)", xnames, nomatch = 0) == 0)
    intercept <- F

  else 
    intercept <- T

  if(length(y) != n)
    stop(paste("Length of response vector must be the same", 
         " as the number of rows in the independent variables",sep=""))

  qrx <- qr(x)[c("rank", "pivot")]
  rank <- qrx$rank
  if(rank < p) stop(paste("x is singular, rank(x)=", rank))

##
## Step 1. Compute Initial Estimates
##

  if(casefold(psi[1]) == "bisquare") {
    ipsi <- 2
		xk <- 1.5477
		beta <- 0.5
  }

  else if(casefold(psi[1]) == "optimal") {
    ipsi <- 1
		xk <- 0.4047
		beta <- 0.2661
  }

  else
    stop("Invalid choice of weight function.")

  bet0 <- 0.773372647623 #bet0 = pnorm(0.75)
	
  if(iop == 5) {
    tmp <- lmRob.lar(x, y, tol=tl)
    coeff0 <- tmp$coef
    resid0 <- tmp$resid
    scale0 <- mad(resid0)
    tmpn <- double(n)
    z1 <- .Fortran("s_rsigm2",
                   as.double(resid0),
                   as.double(resid0),
                   SIGMAI=as.double(scale0),
                   as.integer(n),
                   as.integer(p),
                   as.double(tlo),
                   ITYPE=as.integer(1),
                   ISIGMA=as.integer(1),
                   MAXIS=as.integer(mxs),
                   NIT=integer(1),
                   SIGMAF=double(1),
                   SW=tmpn,
                   SC=tmpn,
                   as.integer(ipsi),
                   as.double(xk),
                   as.double(beta),
                   as.double(bet0))
    scale0 <- z1$SIGMAF
  }

  else if(iop == 4) { 
    set.seed(sed)

    if(is.null(genetic.control))
      genetic.control <- lmRob.genetic.control()

    z1 <- lmRob.ga(x, y, rank, genetic.control, ipsi, xk, beta, tolr=tlo,
                   tau=tua, maxs1=mxs)

    coeff0 <- z1$theta[1:rank]
    scale0 <- z1$smin
    resid0 <- z1$rs
  }

  else if(!is.null(x1)) {
    tmp <- lmRob.lar(x1, y, tol=tl)
    y.tilde <- tmp$resid
    t1 <- tmp$coef
    x2.tilde <- x2
    T2 <- matrix(0, nrow=p1, ncol=p2)
		
##	linux differs from everything else at least here

    for(i in 1:p2) {
      tmp <- lmRob.lar(x1, x2[,i], tol=tl)
      x2.tilde[,i] <- tmp$resid
      T2[,i] <- tmp$coef
    }

    tmpn <- double(n)
    tmp1 <- double(p1)
    tmpp <- double(p2)
    SFGH1 <- matrix(double(1), nrow=p1, ncol=3)
    SFGH2 <- matrix(double(1), nrow=p2, ncol=3)
    xx <- matrix(double(1),p2,p2)
    storage.mode(x1) <- storage.mode(x2.tilde) <- "double"
    z1 <- .Fortran("s_frstml",
                   X1=x1,
                   X2=x2.tilde,
                   Y=as.double(y.tilde),
                   as.integer(n),
                   as.integer(p1),
                   as.integer(p2),
                   as.integer(p2),
                   as.integer(n),
                   T1=tmp1,
                   T2=tmpp,
                   x1c=x1,
                   XTHETA1=tmpn,
                   XTHETA2=tmpp,
                   XX=xx,
                   YY=tmpp,
                   iopt=as.integer(iop),
                   INTCH=as.integer(1),
                   as.integer(nrep),
                   TOLS=as.double(tlo),
                   TOLR=as.double(tlo),
                   TAU=as.double(tua),
                   MAXS1=as.integer(mxs),
                   as.integer(sed),
                   IERR=integer(1),
                   SMIN=double(1),
                   tmpn,
                   tmpn,
                   SFGH1,
                   SFGH2,
                   integer(p2),
                   tmpn,
                   integer(p2),
                   as.integer(ipsi),
                   as.double(xk),
                   as.double(beta),
                   as.double(bet0),
                   as.integer(trc))
    b2 <- z1$T2
    b1 <- t1 + z1$T1 -T2 %*% b2
    coeff0 <- rep(NA, p)
    coeff0[x1.idx]  <- b1
    coeff0[-x1.idx] <- b2
    scale0 <- z1$SMIN
    resid0 <- y - x1 %*% b1 - x2 %*% b2
  }

  else if(iop == 6) {
    storage.mode(x) = "double"
    storage.mode(y) = "double"
    tmpn = double(n)
    tmpp = double(p)
    tmpi = integer(p)
    tmpm = matrix(double(1), p, p)
    n = as.integer(n)
    p = as.integer(p)
    z1 <- .Fortran("s_fastse",
                   x,
                   x,
                   n,
                   p,
                   y,
                   y,
                   n,
                   RES=tmpn,
                   TAU=as.double(tua),
                   K=integer(1),
                   SF=tmpp,
                   SG=tmpp,
                   SH=tmpp,
                   IP=tmpi,
                   XPXH=tmpm,
                   XPXI=tmpm,
                   HDIAG=tmpn,
                   Q=tmpm,
                   U=tmpm,
                   Z=x,
                   ITER=integer(1),
                   IERR=integer(1),
                   x,
                   y,
                   SMIN=as.double(0),
                   as.integer(ipsi),
                   as.double(xk),
                   as.double(beta),
                   as.double(bet0),
                   as.integer(mxs),
                   as.double(tlo),
                   as.double(tlo),
                   as.integer(mxs),
                   THETA=tmpn,
                   XTHETA1=tmpp,
                   XTHETA2=tmpp,
                   IPERM=integer(n),
                   C1=as.double(2))
    if(z1$IERR == 1) 
      stop("Singular matrix encountered.")
    coeff0 = z1$XTHETA2
    scale0 = z1$SMIN
    resid0 = y - x %*% coeff0
  }

  else {
    mdi <- p + rank
    mdw <- (p+2)*rank + 3*p + n
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    tmpn <- double(n)		
    z1 <- .Fortran("s_hsesm2",
                   x,
                   y,
                   as.integer(n),
                   as.integer(p),
                   nq=as.integer(rank),
                   mdx=as.integer(n),
                   as.integer(mdw),
                   as.integer(mdi),
                   iopt=as.integer(iop),
                   intch=as.integer(1),
                   as.integer(nrep),
                   tols=as.double(tlo),
                   tolr=as.double(tlo),
                   tau=as.double(tua),
                   maxs1=as.integer(mxs),
                   iseed=as.integer(sed),
                   ierr=integer(1),
                   smin=double(1),
                   theta=tmpn,
                   rs=tmpn,
                   it1=integer(rank),
                   work=double(mdw),
                   iwork=integer(mdi),
                   as.integer(ipsi),
                   as.double(xk),
                   as.double(beta),
                   as.double(bet0),
                   as.integer(trc))
    coeff0 <- z1$theta[1:rank]
    scale0 <- z1$smin
    resid0 <- z1$rs
  }
	
  if(scale0 < tl) 
    warning(paste("Initial scale less than tl = ", tl, "."))

  z1 <- lmRob.ucovcoef(x, y, resid0, scale0, rank, ipsi=ipsi, xk=xk,
                       tau=tua, tl=tl) 
  if(z1$ierr==1)
    warning(paste(msg.UCV," during initial estimation.")) 

  ucov0 <- z1$ucov
  M.weights0 <- z1$wi

##
## Step 2. Refine Initial Estimates
##

  if(iop == 5) {
    tmpp <- double(p)
    storage.mode(x) <- "double"
    z1 <- .Fortran("s_ywagm2",
                   x,
                   as.double(y),
                   theta=as.double(coeff0),
                   sigmaf=as.double(scale0),
                   as.integer(n),
                   as.integer(p),
                   as.integer(n),
                   as.double(tlo),
                   GAM=as.double(1),
                   as.double(tua),
                   as.integer(mxr),
                   nit=integer(1),
                   rs=as.double(resid0),
                   tmpp,
                   tmpn,
                   tmpp,
                   tmpp,
                   tmpp,
                   integer(p),
                   x)
  }

  else if(is.null(x1)) {
    coeffi <- vector("numeric", max(n, rank))
    coeffi[1:rank] <- coeff0
    z1 <- lmRob.wm(x, y, coeffi, ucov0, scale0, ipsi=ipsi, xk=xk,
                   beta=beta, tlo=tlo, tua=tua, mxr=mxr)
  }

  else { 
    A <- matrix(double(1), p1, p2)
    B <- matrix(double(1), p1, p1)
    storage.mode(x1) <- storage.mode(x2) <- "double"
    z1 <- .Fortran("s_dscnm2",
                  x1,
                  x2,
                  as.double(y),
                  as.integer(n),
                  as.integer(p1),
                  as.integer(p2),
                  as.integer(n),
                  as.double(scale0),
                  sigmaf=double(1),
                  as.double(b1),
                  as.double(b2),
                  b1=tmp1,
                  b2=tmpp,
                  rs=as.double(resid0),
                  RSTMP=tmpn,
                  as.double(tlo),
                  as.double(tua),
                  as.integer(mxr),
                  as.integer(mxs),
                  SFGH1,
                  as.integer(ipsi),
                  as.double(xk),
                  as.double(beta),
                  double(1),
                  IFAIL=integer(1),
                  tmpn,
                  A,
                  B,
                  CC=xx,
                  C2=xx,
                  D=x2,
                  BD=A,
                  H=tmpp,
                  tc=as.double(1.55),
                  x1,
                  IP=integer(p1),
                  IDX=integer(p2),
                  WP1=tmp1,
                  WP2=tmpp,
                  nit=integer(1),
                  MAXK=as.integer(mxr))
    z1$theta <- rep(NA, p)
    z1$theta[ x1.idx] <- z1$b1
    z1$theta[-x1.idx] <- z1$b2
  }

  iter.refinement <- z1$nit
  tmp <- z1$sigmaf

  if(tmp < tl) 
    warning(paste("Refined scale less than tl = ", tl, "."))

  if(iter.refinement == mxr) 
    warning("Max iteration for refinement reached.")

  if(tmp < scale0) {
    coeff0 <- z1$theta[1:rank]
    scale0 <- tmp
    resid0 <- z1$rs

    z1 <- lmRob.ucovcoef(x, y, resid0, scale0, rank, ipsi=ipsi, xk=xk,
						tau=tua,tl=tl)

    if(z1$ierr==1)
      warning(paste(msg.UCV," when refining initial estimates.")) 

    ucov0 <- z1$ucov
    M.weights0 <- z1$wi
  }

  names(coeff0) <- xnames
  attributes(resid0) <- attributes(y)

##
## Step 3. Robust Cov, Dev, and R2 of Initial Estimates
##

  dev0 <- scale0^2
  tmp <- ucov0 * dev0 / n
  scov0 <- matrix(NA, nrow = rank, ncol = rank)
  dimnames(scov0) <- list(xnames, xnames)
  i2 <- 0

  for(i in 1:rank) { 
    i1 <- i2 + 1
    i2 <- i1 + i - 1
    scov0[i, 1:i] <- tmp[i1:i2]
    scov0[1:i, i] <- tmp[i1:i2]
  }

  t0 <- if (intercept) median(y) else 0
  ymt <- abs(y - t0)
  s0 <- median(ymt)/0.6745

  if(s0 == 0) {
    s0 <- min(ymt[ymt != 0])
    while(sum(chi.weight(ymt/s0,ipsi,xk)) > (n - 1) * 0.5) 
      s0 <- 1.5 * s0
  }

  tmp <- matrix(1, ncol = 1, nrow = n)
  ymt[1] <- t0

  z1 <- lmRob.wm(tmp, y, ymt, ucov0[1], s0, isigma=1, ipsi=ipsi,
                 xk=xk, beta=beta, tlo=tlo, tua=tua, mxr=mxr)

  tstar <- if (intercept) z1$theta[1] else 0
  sy <- (z1$sigmaf)^2
  r.squared0 <- ((n-intercept)*sy - (n-rank)*scale0^2)/((n-intercept)*sy)
  fval0 <- y-resid0
  dfres <- length(resid0) - rank

  z <- list(coefficients=coeff0, scale=scale0, residuals=resid0,
            fitted.values=fval0, cov=scov0, rank=rank,
            iter.refinement=iter.refinement,  df.residual=dfres, 
            est="initial", robust.control=robust.control)

##
## Step 4. Compute Final Estimates
##

  if (casefold(est) == "final") {
    coeffi <- vector("numeric", max(n, rank))
    coeffi[1:rank] <- coeff0

    if (casefold(fnl) == "mm" || casefold(fnl) == "m") {
      if (casefold(psi[2]) == "optimal") {
        ipsi2 <- 1
        if (eff == 0.95)      yc <- 1.060158
        else if (eff == 0.9)  yc <- 0.9440982
        else if (eff == 0.85) yc <- 0.8684
        else if (eff == 0.8)  yc <- 0.8097795
        else                  yc <- lmRob.effvy(eff)
      }

      else if (casefold(psi[2]) == "bisquare") {
        ipsi2 <- 2
        if (eff == 0.95)      yc <- 4.685061
        else if (eff == 0.9)  yc <- 3.882646
        else if (eff == 0.85) yc <- 3.443689
        else if (eff == 0.8)  yc <- 3.136909
        else                  yc <- chb(eff)$cb
      }

      else 
        stop("Invalid choice of weight function.")

      z1 <- lmRob.wm(x, y, coeffi, ucov0, scale0, isigma=0, ipsi=ipsi2, 
                    xk=yc, beta=beta, tlo=tlo, tua=tua, mxr=mxf)

      iter.final.coef <- z1$nit

      if(iter.final.coef == mxf) 
        warning("Max iteration for final coefficients reached.")

      ymt <- vector("numeric", n)
      ymt[1] <- tstar
      tmp <- matrix(1, nrow = n, ncol = 1)

      if(intercept) {
        that <- lmRob.wm(tmp, y, ymt, ucov0[1], scale0, isigma=0, 
                         ipsi=ipsi2, xk=yc, beta=beta, tlo=tlo, 
                         tua=tua, mxr=mxf)$theta[1]
      }

      else 
        that <- 0

      ymt <- (y - that)/scale0
      sy <- sum(chi.weight(ymt,ipsi2,yc))
      ymt <- z1$rs/scale0
      s0 <- sum(chi.weight(ymt,ipsi2,yc))
      r.squared1 <- (sy - s0)/sy
      dev1 <- (scale0^2) * s0

      z2 <- lmRob.ucovcoef(x, y, resid0, scale0, rank, ipsi=ipsi2,
                           xk=yc, tau=tua, tl=tl)

      if(z2$ierr == 1) 
        warning(paste(msg.UCV,"during final scale estimation.")) 

      ucov1 <- z2$ucov
      M.weights1 <- z2$wi
      tmp <- ucov1 * scale0^2 / n
      scov1  <- matrix(NA, nrow = rank, ncol = rank)
      dimnames(scov1) <- list(xnames,xnames)
      i2 <- 0

      for(i in 1:rank) {
        i1 <- i2 + 1
        i2 <- i1 + i - 1
        scov1[i, 1:i] <- tmp[i1:i2]
        scov1[1:i, i] <- tmp[i1:i2]
      }

      z2 <- .Fortran("s_rsigm2",
                     z1$rs,
                     wgt=y,
                     sigmai=as.double(2*scale0),
                     as.integer(n),
                     np=as.integer(rank),
                     tol=as.double(tlo),
                     itype=as.integer(1),
                     isigma=as.integer(1),
                     maxis=as.integer(mxs),
                     nit=integer(1),
                     sigmaf=double(1),
                     double(n),
                     double(n),
                     as.integer(ipsi),
                     as.double(xk),
                     as.double(beta),
                     as.double(1.0))
      scale1 <- z2$sigmaf
    }

    else if (casefold(fnl) == "adaptive") {
      tmp <- double(p)
      residi <- resid0
      storage.mode(x) <- storage.mode(y) <- "double"
      eta <- lmRob.effad(eff)
      nit <- 1

      for (i in 1:nit) {
        z1 <- .Fortran("s_finlml",
                       X=x,
                       Y=y,
                       WGT=double(n),
                       rs=as.double(residi),
                       as.integer(n),
                       as.integer(p),
                       as.integer(n),
                       theta=as.double(coeffi),
                       as.double(scale0),
                       tmp,
                       tmp,
                       tmp,
                       integer(p),
                       x,
                       y,
                       as.double(tua),
                       as.double(eta),
                       ierr=integer(1),
                       as.integer(ipsi),
                       as.double(xk),
                       f=double(1),
                       double(n))
        coeffi <- z1$theta
        residi <- z1$rs
      }

      if (z1$ierr == 1) {
        tmp <- paste("The weighted model matrix is singular.",
                     "Only initial estimates are returned.")
        warning(tmp)
        est <- "initial"
      }

      scov1 <- scale0*scale0*solve(t(x) %*% x)*z1$f
      dimnames(scov1) <- list(xnames, xnames)
    }

    else
      stop("Invalid choice of final estimation.")

    coeff1 <- z1$theta[1:rank]
    names(coeff1) <- xnames
    resid1 <- z1$rs
    attributes(resid1) <- attributes(y)
    fval1 <- y - resid1
  }

  if (casefold(est) == "final") {

    z <- list(coefficients=coeff1, scale=scale0, residuals=resid1,
              fitted.values=fval1, cov=scov1, rank=rank,
              iter.refinement=iter.refinement, df.residual=dfres,
              est="final", robust.control=robust.control)

    z$T.coefficients  <- coeff0
    z$T.residuals     <- resid0
    z$T.fitted.values <- fval0
    z$T.cov          <- scov0
    if (fnl == "MM" || fnl == "m" || fnl == "M") {
      z$dev <- dev1;              z$T.dev <- dev0
      z$M.weights <- M.weights1;  z$T.M.weights <- M.weights0
      z$r.squared <- r.squared1;  z$T.r.squared <- r.squared0
      z$T.scale         <- scale1
      z$iter.final.coef <- iter.final.coef
    }
  }

  else {
    if(is.null(x1)) {
      z$dev <- dev0
      z$M.weights <- M.weights0
      z$r.squared <- r.squared0
    }
  } 

  if(!is.null(genetic.control))
    z$genetic.control <- genetic.control 
  z$qr <- qrx
	z$yc <- lmRob.effvy(eff)
  oldClass(z) <- c("lmRob")
	z
}


#--------------------------------------------------------------#
# Control Lists                                                #
#--------------------------------------------------------------#


lmRob.robust.control <- function(tlo=0.0001, tua=1.5e-06, mxr=50,
    mxf=50, mxs=50, tl=1e-6, estim="Final", initial.alg="Auto",
    final.alg="MM", seed=1313, level=0.10, efficiency=0.90,
    weight=c("Optimal","Optimal"), trace=T)
{
  list(tlo=tlo, tua=tua, mxr=mxr, mxf=mxf, mxs=mxs, tl=tl,
       estim=estim, initial.alg=initial.alg, final.alg=final.alg, 
       seed=seed, level=level, efficiency=efficiency, weight=weight,
       trace=trace)
}


lmRob.genetic.control <- function(popsize=NULL, mutate.prob=NULL,
    random.n=NULL, births.n=NULL, stock=list(), maxslen=NULL,
    stockprob=NULL, nkeep=1)
{
  list(popsize=popsize, mutate.prob=mutate.prob, random.n=random.n,
       births.n=births.n, stock=stock, maxslen=maxslen,
       stockprob=stockprob, nkeep=nkeep)
}


#-------------------------------------------------------------#
# Access Functions                                            #
#-------------------------------------------------------------#

coef.lmRob <- function(object,...) {
  oldClass(object) <- "lm"
  UseMethod("coef")
}

 
formula.lmRob <- function(object) {
  attr(object$terms, "formula")
}


predict.lmRob <- function(object, newdata, type=c("response","terms"), 
   se.fit = F, terms = labels.lm(object))
{
  type <- match.arg(type)
  if(missing(newdata) && type != "terms" && !se.fit)
    return(fitted(object))

  Terms <- object$terms

  if(!inherits(Terms, "terms"))
    stop("invalid terms component of  object")

  offset <- attr(Terms, "offset")
  intercept <- attr(Terms, "intercept")
  xbar <- NULL
  if(missing(newdata)) {
    x <- model.matrix(object)
		#center x if terms are to be computed
    if(type == "terms" && intercept) {
      xbar <- colMeans(x)
      x <- sweep(x, 2, xbar)
    }
  }

  else if(!((is.atomic(newdata) && length(newdata) == 1 && 
             length(object$coef) != 1 && newdata > 0 && 
             (newdata - trunc(newdata) < .Machine$single.eps)) | 
             is.list(newdata))) {

		#try and coerce newdata to look like the x matrix
    if(!is.null(offset)) {
      warning("Offset not included")
      offset <- NULL
    }

    TT <- length(object$coef)

    if(is.matrix(newdata) && ncol(newdata) == TT)
      x <- newdata

    else if(length(newdata) == TT)
      x <- matrix(newdata, 1, TT)

    else 
      stop("Argument \"newdata\" is not a data frame, and cannot be coerced to an appropriate model matrix.")
  }

  else {
		#newdata is a list, data frame or frame number
    x <- model.matrix(delete.response(Terms), newdata, contrasts = 
                      object$contrasts, xlevels=attr(object, "xlevels"))

    if(!is.null(offset))
      offset <- eval(attr(Terms, "variables")[offset], newdata)
  }

  if(!missing(newdata) && type == "terms" && intercept) {
		#need to center x 
    xold <- model.matrix(object)
    xbar <- colMeans(xold)
    x <- sweep(x, 2, xbar)
  }

  coefs <- coef(object)
  asgn <- attr(coefs, "assign")

  if(type == "terms") {
    terms <- match.arg(terms, labels.lm(object))
    asgn <- asgn[terms]
  }

  nac <- is.na(object$coef)

  if(any(nac)) {
    xbar <- xbar[!nac]
    x <- x[, !nac]
  }

  attr(x, "constant") <- xbar

  if(se.fit) {
    fit.summary <- summary.lmRob(object)
    pred <- Build.terms(x, coefs, fit.summary$cov * fit.summary$sigma^2, 
                        asgn, collapse = type != "terms")
    pred$residual.scale <- fit.summary$sigma
    pred$df <- object$df.resid
  }

  else 
    pred <- Build.terms(x, coefs, NULL, assign = asgn, 
                        collapse = type != "terms")

  if(!is.null(offset) && type != "terms") {
    if(missing(newdata))
      warning("Offset not included")
    else {
      if(se.fit)
        pred$fit <- pred$fit + offset
      else 
        pred <- pred + offset
    }
  }

  if(missing(newdata) && !is.null(object$na.action)) {
    if(!se.fit)
      pred <- napredict(object$na.action, pred)
    else {
      pred$fit <- napredict(object$na.action, pred$fit)
      pred$se.fit <- napredict(object$na.action, pred$se.fit)
    }
	}
  pred
}


residuals.lmRob <- function(object, ...) {
  oldClass(object) <- "lm"
  UseMethod("residuals")
}

 
model.frame.lmRob <- function(formula, ...) {
  oldClass(formula) <- "lm"
  UseMethod("model.frame")
}
 

model.matrix.lmRob <- function(object, ...) {
  oldClass(object) <- "lm"
  UseMethod("model.matrix")
}


scale.lmRob <- function(x) {
  x$scale
} 


rsquared.lmRob <- function(x) {
  z <- x$r.squared
  if (is.null(z))
    warning("R^2 is only defined for S-estimates and MM-estimates.")
  z
}


deviance.lmRob <- function(object, ...) {
  z <- object$dev
  if (is.null(z))
    warning("Deviance is only defined for S-estimates and MM-estimates.")
  z
}


weights.lmRob <- function(x) {
  z <- x$M.weights
  if (is.null(z))
    warning("M.weights is only defined for S-estimates and MM-estimates.")
  z
}


cov.lmRob <- function(x) {
  x$cov
}


cor.lmRob <- function(x, tl = 1.e-10) {
  z <- cov.lmRob(x)
  zvar <- diag(z)
  for(i in 1:nrow(z)) {
    if(zvar[i]<tl) {
      str <- paste("Variance number",i,"smaller than",tl,"(set to",tl,")")
      warning(str)
    }
    z[i,1:i] <- z[i,1:i]/sqrt(zvar[i]*zvar[1:i])
    z[1:i,i] <- z[i,1:i]
  }
  z
} 


#--------------------------------------------------------------#
# Methods for Generic Functions                                #
#--------------------------------------------------------------#


test.lmRob <- function(object, type="bias", level=NULL, n.permute=99)
{

##
## Step 1. Construct Response Vector and Model Frame
##

  call <- object$call
  contrasts <- object$contrasts

  m <- list()
  m[[1]] <- as.name("model.frame")
  m$formula <- call$formula
  m$data <- call$data
  m$weights <- call$weights
  m$subset <- call$subset
  m$na.action <- call$na.action
  mode(m) <- "call"

  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
  weights <- model.extract(m, weights)
  y <- model.extract(m, response)
  x <- model.matrix(Terms, m, contrasts)
  n <- length(y)
  p <- length(coef(object))
  control <- object$robust.control

  if(casefold(type) == "permutation") {
    x.nam <- dimnames(x)[[2]]

    if(x.nam[1] == "(Intercept)")
      x <- x[,-1]

    if(!is.null(dim(x)))
      stop("Permutation test is only available for a straight line fit.")

    control$initial.alg <- "Random"
    gen.list <- object$genetic.control
    perm.x <- samp.permute(1:n, n.permute)
    beta.all <- rep(0, n.permute)

    for(i in 1:n.permute) {
      tmp <- lmRob(y~x[perm.x[,i]],robust.control=control,
                   genetic.control=gen.list)
      beta.all[i] <- coef(tmp)[2]
    }

    beta.all <- abs(c(0, beta.all))
    k <- (1:(n.permute+1))[order(beta.all) == 1] - 1
    return((n.permute+1-k)/(n.permute+1))
  }

  if (casefold(type) != "bias")
    stop("Wrong type for test.lmRob")


  tmp <- control$final.alg
  if (casefold(tmp) == "adaptive")
    stop("Tests for bias are only available for MM-estimates.")

##
## Step 2. Extract Components from Model Object
##

  rank <- object$rank
  resid0 <- object$T.residuals
  scale0 <- object$scale
  scale1 <- object$T.scale
  tl <- control$tl
  psi <- control$weight

  if(is.null(level))
    level <- 1-control$level

  eff <- control$efficiency

  if(casefold(psi[1]) == "bisquare") {
    ipsi <- 2
    xk <- 1.5477
  }

  else if(casefold(psi[1]) == "optimal") {
    ipsi <- 1
    xk <- 0.4047
  }

  else 
    stop("Invalid choice of weight function.")

  if(casefold(psi[2]) == "optimal") {
    ipsi2 <- 1
    if (eff == 0.95)      yc <- 1.060158
    else if (eff == 0.9)  yc <- 0.9440982
    else if (eff == 0.85) yc <- 0.8684
    else if (eff == 0.8)  yc <- 0.8097795
    else                  yc <- lmRob.effvy(eff)
  }

  else if(casefold(psi[2]) == "bisquare") {
    ipsi2 <- 2
    if (eff == 0.95) yc <- 4.685061
    else if (eff == 0.9)  yc <- 3.882646
    else if (eff == 0.85) yc <- 3.443689
    else if (eff == 0.8)  yc <- 3.136909
    else                  yc <- chb(eff)$cb
  }

  else 
    stop("Invalid choice of weight function.")

##
## Step 3. Compute Test Statistics
##

  b.OLS <- solve(x, y)
  resid.OLS <- y - x %*% b.OLS
  scale.OLS <- sqrt(sum(resid.OLS * resid.OLS)/(n-p))
  tmp <- resid0/scale0
  sc0.OLS <- tmp
  sc0 <- psi.weight(tmp, ipsi2, yc)
  s1p <- sum(psp.weight(tmp, ipsi2, yc))/n
  s0p <- sum(psp.weight(tmp, ipsi, xk))/n
  sc1 <- psi.weight(tmp, ipsi, xk)
  tmp <- tmp * psi.weight(tmp, ipsi, xk)
  s0r <- (sum(tmp) * scale0)/n

  if(s0r < tl || s0p < tl || s1p < tl) 
    warning(paste("Denominator smaller than tl=",tl," in test for bias."))

  tmp <- sc0/s1p - sc1/s0p
  d2 <- sum(tmp * tmp)/n
  tmp <- sc0.OLS - sc1/s0p
  d2.OLS <- sum(tmp * tmp)/n

  if(d2 < tl) 
    warning(paste("Denominator smaller than tl=",tl," in test for bias."))

  tbias <- (2*n*(scale1-scale0)*s0r)/(s0p*d2*scale0*scale0)
  tbias.OLS <- (2*n*(scale.OLS-scale0)*s0r)/(s0p*d2.OLS*scale0*scale0)
  qchi <- qchisq(level, rank)
  pchi <- 1-pchisq(tbias, rank)
  pchi.OLS <- 1-pchisq(tbias.OLS,rank)
  mm.bias <- list(stat=tbias, pchi=pchi, qchi=qchi)
  ls.bias <- list(stat=tbias.OLS, pchi=pchi.OLS)

  ans <- list(mm=mm.bias, ls=ls.bias, level=level)
  oldClass(ans) <- "biasMM" 
  ans
}


print.biasMM <- function(object) {
  mmbias <- object$mm
  lsbias <- object$ls
  level  <- object$level
  cat("\nTest for Bias:\n")
  tmp <- rbind(c(mmbias$stat, mmbias$pchi),
               c(lsbias$stat, lsbias$pchi))
  dimnames(tmp) <- list(c("M-estimate","LS-estimate"),
                        c("Statistics","P-value"))
  print(tmp,digits=3)

  if (mmbias$stat > mmbias$qchi) {
    warning(paste("Significant test at level ", 100*(1-level), "%. ",
      "The bias is high, and inference based on final estimates is not ",
      "recommended. Use initial estimates as exploratory tools.", sep=""))
  }

  invisible(object)
}


drop1.lmRob <- function(object, scope, scale, keep, fast = F, ...)
{
  rfpe.compute <- function(res, scale, ipsi, yc, p) {
		res <- res / scale
		a <- sum(rho.weight(res, ipsi, yc))
		b <- (p * sum(psi.weight(res, ipsi, yc)^2))
		d <- sum(psp.weight(res, ipsi, yc))

		if(d <= 0)
			return(NA)

		a + b/d
  }

  tmp <- object$robust.control$final.alg

  if(casefold(tmp) == "adaptive")
    stop("drop1 is only available for MM-estimates.")

  b <- coef(object)
  singular <- attr(b, "singular")

  if(singular)
    stop("The model matrix is singular.")

  psif <- object$robust.control$weight
  x <- model.matrix(object)
  asgn <- attr(x, "assign")
  tl <- attr(object$terms, "term.labels")

  if(missing(scope))
    scope <- drop.scope(object)

  else {
    if(!is.character(scope))
      scope <- attr(terms(update.formula(object,scope)),"term.labels")

    if(!all(match(scope, tl, F)))
      stop("scope is not a subset of term labels")
  }

  p <- length(asgn)
  asgn <- asgn[scope]
  k <- length(scope)

  if(missing(scale)) 
    scale <- object$scale

  if(object$est == "initial")
		warning("Inference based on initial estimates is not recommended.")

  rfpe.none <- lmRob.RFPE(object, scale)

  if(!missing(keep)) {
    max.keep <- c("coefficients", "fitted", "residuals")

    if(is.logical(keep) && keep) 
      keep <- max.keep

    else {
      if(!all(match(keep, max.keep, F)))
        stop(paste("Can only keep one or more of: \"",
                   paste(max.keep, collapse = "\", \""), "\"", sep = ""))
    }

    value <- array(vector("list", 3 * k), c(k, 3),
                   list(scope, c("coefficients", "fitted", "residuals")))
  }

  else
    keep <- character(0)

  dfs <- double(k)
  rfpe <- double(k)

  if(fast) {
    Weights <- object$M.weights
    ipsi <- 1
    xk <- .9440982
    itype <- 1
    isigma <- -1
    yc <- object$yc
    mxr <- object$robust.control$mxr
		mxs <- object$robust.control$mxs
    tlo <- object$robust.control$tlo
    tua <- object$robust.control$tua
    tl <- object$robust.control$tl
    y <- object$fitted.values + object$residuals
    n <- length(y)
		tmpn <- double(n)
 		beta <- sum(Weights)/(2*n)
		bet0 <- 1

		for(i in 1:k) {
      ii <- asgn[[i]]
      pii <- length(ii)
      dfs[i] <- pii
      curfrm <- update.formula(object, formula = paste(".~.-", scope[[i]]), evaluate = F)
      Terms <- terms(curfrm)
      x.columns <- unlist(asgn[attr(Terms, "term.labels")])
      curobj <- lsfit(x = x[, x.columns, drop = F], y = y, wt = Weights, intercept = F)
      coeff0 <- curobj$coef

			scale0 <- .Fortran("s_rsigm2",
													as.double(curobj$residuals),
													as.double(curobj$residuals),
													as.double(mad(curobj$residuals)),
													as.integer(n),
													as.integer(dim(x[, x.columns, drop = F])[2]),
													as.double(tlo),
													as.integer(1),
													as.integer(1),
													as.integer(mxs),
													as.integer(1),
													SIGMAF = double(1),
													as.double(tmpn),
													as.double(tmpn),
													as.integer(ipsi),
													as.double(xk),
													as.double(beta),
													as.double(bet0))$SIGMAF

			ucov0 <- length(y) * solve(t(sqrt(Weights) * x[, x.columns, drop = F]) %*%
				(sqrt(Weights) * x[, x.columns, drop = F])) %*% t(Weights * x[, x.columns, drop = F]) %*%
				(Weights * x[, x.columns, drop = F]) %*% solve(t(sqrt(Weights) * x[, x.columns, drop = F]) %*%
				(sqrt(Weights) * x[, x.columns, drop = F]))

			ucov0 <- ucov0[row(ucov0) <= col(ucov0)]

      curobj <- lmRob.wm(x = x[, x.columns, drop = F], y = y, coeff0 = coeff0,
								ucov0 = ucov0, scale0 = scale0, itype = itype, isigma = isigma,
								ipsi = ipsi, xk = xk, beta = beta, wgt = y, tlo = tlo, tua = tua,
								mxr = mxr)

      rfpe[i] <- rfpe.compute(curobj$rs, scale, ipsi, yc, p - 1)
    }
  }

  else {
    for(i in 1:k) {
      ii <- asgn[[i]]
      pii <- length(ii)
      dfs[i] <- pii
      curfrm <- as.formula(paste(".~.-", scope[[i]]))
      curobj <- update(object, curfrm)
      rfpe[i] <- lmRob.RFPE(curobj, scale)
      if(length(keep)) {
        value[i,1] <- list(curobj$coefficients)
        value[i,2] <- list(curobj$fitted)
        value[i,3] <- list(curobj$residuals)
      }
    }
  }
  
  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  rfpe <- c(rfpe.none, rfpe)
  dfs[1] <- NA
  aod <- data.frame(Df=dfs, RFPE=rfpe, row.names=scope, check.names=F)

  head <- c("\nSingle term deletions", "\nModel:",
            deparse(as.vector(formula(object))))

  if(!missing(scale))
    head <- c(head, paste("\nscale: ", format(scale), "\n"))

  oldClass(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head

  if(length(keep))
    list(anova=aod, keep=structure(value[, keep, drop=F], class="matrix"))
  else
    aod
}


step.lmRob <- function(object, scope, scale, direction =
	c("both", "backward", "forward"), trace = T, keep = NULL,
	steps = 1000, control = NULL, fast = F, ...)
{
	if(missing(direction))
		direction <- "backward"

	else direction <- match.arg(direction)

	if(direction != "backward")
		stop("Presently step.lmRob only supports backward model selection.")

	sub.assign <- function(terms, assign)
	{
		a <- attributes(terms)
		tl <- a$term.labels
		if(a$intercept)
			tl <- c(names(assign)[1], tl)
		asgn <- assign[tl]
		poi <- 0
		for(i in tl) {
			la <- length(asgn[[i]])
			asgn[[i]] <- seq(poi + 1, poi + la)
			poi <- poi + la
		}
		asgn
	}

	re.arrange <- function(keep)
	{
		namr <- names(k1 <- keep[[1]])
		namc <- names(keep)
		nc <- length(keep)
		nr <- length(k1)
		array(unlist(keep, recursive = F), c(nr, nc), list(namr, namc))
	}

	make.step <- function(models, fit, scale, object)
	{
		change <- sapply(models, "[[", "change")
		rdf <- sapply(models, "[[", "df.resid")
		ddf <- c(NA, diff(rdf))
		RFPE <- sapply(models, "[[", "RFPE")
		heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
			"\nInitial Model:", deparse(as.vector(formula(object))),
			"\nFinal Model:", deparse(as.vector(formula(fit))),
			"\n")
		aod <- data.frame(Step = change, Df = ddf, "Resid. Df" = rdf,
			RFPE = RFPE, check.names = F)
		attr(aod, "heading") <- heading
		oldClass(aod) <- c("anova", "data.frame")
		fit$anova <- aod
		fit
	}

	backward <- direction == "both" || direction == "backward"
	forward <- direction == "both" || direction == "forward"

	if(missing(scope)) {
		fdrop <- numeric(0)
		fadd <- NULL
	}


	else {
		if(is.list(scope)) {
			fdrop <- if(!is.null(fdrop <- scope$lower)) attr(terms(
					update.formula(object, fdrop)), 
					"factor") else numeric(0)
			fadd <- if(!is.null(fadd <- scope$upper)) attr(terms(
					update.formula(object, fadd)), "factor")
		}
		else {
			fadd <- if(!is.null(fadd <- scope))
				attr(terms(update.formula(object, scope)), "factor")
			fdrop <- numeric(0)
		}
	}

	if(is.null(fadd)) {
		backward <- T
		forward <- F
	}

	m <- model.frame(object)
	obconts <- object$contrasts
	objectcall <- object$call
	robust.control <- object$robust.control

	#build the big model matrix
	if(forward) {
		add.rhs <- paste(dimnames(fadd)[[2]], collapse = "+")
		add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
		new.form <- update.formula(object, add.rhs, evaluate = F)
		fc <- objectcall
		Terms <- terms(new.form)
		fc$formula <- Terms
		fobject <- list(call = fc)
		oldClass(fobject) <- oldClass(object)
		m <- model.frame(fobject)
		x <- model.matrix(Terms, m, contrasts = obconts)
	}

	else {
		Terms <- object$terms
		x <- model.matrix(Terms, m, contrasts = obconts)
	}

	Asgn <- attr(x, "assign")
	a <- attributes(m)
	y <- model.extract(m, "response")
	w <- model.extract(m, "weights")

	if(is.null(w))
		w <- rep(1, nrow(m))

	models <- vector("list", steps)

	if(!is.null(keep)) {
		keep.list <- vector("list", steps)
		nv <- 1
	}

	n <- length(object$fitted)
	scale <- object$scale
	fit <- object
	cf <- attributes(coef(object))
	#check if any terms have zero df
	if(cf$singular) {
		TT <- !match(TL <- attr(object$terms, "term.labels"), names(cf$assign), F)
		if(any(TT)) {
			upd <- eval(parse(text = paste(c(".~.", TL[TT]), collapse = "-")))
			fit <- update(fit, upd)
		}
	}

	bRFPE <- lmRob.RFPE(fit)
	nm <- 1
	Terms <- fit$terms
	if(trace)
		cat("Start:  RFPE=", format(round(bRFPE, 4)), "\n",
      deparse(as.vector(formula(fit))), "\n\n")

	models[[nm]] <- list(df.resid = fit$df.resid, change = "", RFPE = bRFPE)

	if(!is.null(keep))
		keep.list[[nm]] <- keep(fit, bRFPE)

	RFPE <- bRFPE + 1

	while(bRFPE < RFPE & steps > 0) {
		steps <- steps - 1
		RFPE <- bRFPE
		bfit <- fit
		ffac <- attr(Terms, "factor")
		scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
		aod <- NULL
		change <- NULL

		if(backward && (ndrop <- length(scope$drop))) {
			aod <- drop1.lmRob(fit, scope$drop, scale, fast = fast)
			if(trace)
				print(aod)
			change <- rep("-", ndrop + 1)
		}

		if(forward && (nadd <- length(scope$add))) {
			aodf <- add1.lmRob(fit, scope$add, scale, x = x)
			if(trace)
				print(aodf)
			change <- c(change, rep("+", nadd + 1))
			if(is.null(aod))
				aod <- aodf
			else {
				ncaod <- dim(aod)[1]
				aod[seq(ncaod + 1, ncaod + nadd + 1),  ] <- aodf
			}
		}

		if(is.null(aod))
			break

		o <- order(aod[, "RFPE"])[1]

		## If the original model has minimum RFPE then break
		if(o[1] == 1) break

		change <- paste(change[o], dimnames(aod)[[1]][o])
		Terms <- terms(update(formula(fit), eval(parse(text = paste("~ .", change)))))
		attr(Terms, "formula") <- new.formula <- rebld.formula(Terms)
		asgn <- sub.assign(Terms, Asgn)
		tx <- x[, unlist(Asgn[names(asgn)]), drop = F]
		newfit <- lmRob(new.formula, data = m, robust.control = robust.control)
		bRFPE <- aod[, "RFPE"][o]

		if(trace)
			cat("\nStep:  RFPE =", format(round(bRFPE, 4)), "\n",
				deparse(as.vector(formula(Terms))), "\n\n")

		if(bRFPE >= RFPE)
			break

		nm <- nm + 1
		models[[nm]] <- list(df.resid = newfit$df.resid, change = change, RFPE = bRFPE)
		fit <- c(newfit, list(formula = new.formula))
		oc <- objectcall
		oc$formula <- as.vector(fit$formula)
		fit$call <- oc
		oldClass(fit) <- oldClass(object)
		if(!is.null(keep))
			keep.list[[nm]] <- keep(fit, bRFPE)
	}

	if(!is.null(keep))
		fit$keep <- re.arrange(keep.list[seq(nm)])

	make.step(models = models[seq(nm)], fit, scale, object)
}


plot.lmRob <- function(x, which.plots = "ask", chisq.percent = 0.99,
	vertical.outlier = 0.99, smooths = F, rugplot = F, id.n = 3,
	envelope = T, half.normal = F, robustQQline = T, mc.samples = 100,
	level = 0.95, seed = 289, cutoff = T, ...)

{
	x.name <- deparse(substitute(x))
	model.list <- list(x$call)
	names(model.list) <- x.name
	x <- list(x = x)
	names(x) <- x.name
	attr(x, "model.list") <- model.list

	plot.lmfm(x, which.plots = which.plots, chisq.percent = chisq.percent,
		vertical.outlier = vertical.outlier, smooths = smooths, rugplot = rugplot,
		id.n = id.n, envelope = envelope, half.normal = half.normal,
		robustQQline = robustQQline, mc.samples = mc.samples, level = level,
		seed = seed, cutoff = cutoff, ...)

	invisible(x[[1]])
}


print.lmRob <- function(x, ...) 
{
  if (x$est == "initial")
    cat("Initial Estimates.\n\n")

  if(!is.null(cl <- x$call))
    cat("Call:\n"); dput(cl)

  coef <- coefficients(x)

  if(ns <- attr(coef, "singular"))
    cat("\nCoefficients: (", ns, 
        " not defined because of singularities)\n", sep = "")

  else 
    cat("\nCoefficients:\n")

  print(coef, ...)
  rank <- x$rank
  nobs <- length(x$residuals)
  rdf <- x$df.resid

  if(is.null(rdf))
    rdf <- nobs - rank

  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")

  if (!is.null(x$na.action))
    cat(naprint(x$na.action),"\n")

  if(rdf > 0) {
    if(is.null(w <- x$weights))
      cat("Residual standard error:", format(x$scale), "\n")

    else 
      cat("Residual standard error (on weighted scale):", 
          format(x$scale), "\n")
  }
  invisible(x)
}


print.summary.lmRob <- function(x, digits=max(3,.Options$digits-3), ...) 
{
  if (x$est == "initial") 
    cat("Initial Estimates.\n")

  cat("\nCall: ")
	dput(x$call)
  resid <- x$residuals
  attr(resid, ".guiColInfo") <- NULL
  df <- x$df
  rdf <- df[2]

  if(rdf > 5) {
    cat("\nResiduals:\n")
    if(length(dim(resid)) == 2) {
      rq <- apply(t(resid), 1, quantile)
      dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"), 
                           dimnames(resid)[[2]])
    }

    else {
      rq <- quantile(resid)
      names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print(rq, digits = digits, ...)
  }

  else if(rdf > 0) {
    cat("\nResiduals:\n")
    print(resid, digits = digits, ...)
  }

  if(nsingular <- df[3] - df[1])
    cat("\nCoefficients: (", nsingular, 
        " not defined because of singularities)\n", sep = "")

  else 
    cat("\nCoefficients:\n")

	if(!is.null(x$bootstrap)) {
		coef.names <- dimnames(x$coef)
		coef.names[[2]] <- c(coef.names[[2]][1:2],
			"Bootstrap SE", coef.names[[2]][3:4])
		the.coef <- cbind(x$coef[,1:2], x$bootstrap.se, x$coef[,3:4])
		dimnames(the.coef) <- coef.names
	}
	else
		the.coef <- x$coef

	print(format(round(the.coef, digits = digits)), quote = F, ...)

  cat("\nResidual standard error:", format(signif(x$sigma, digits)), 
      "on",rdf, "degrees of freedom\n")

  if (!is.null(x$na.action))
    cat(naprint(x$na.action),"\n")

  if (!is.null(x$r.squared))
    cat("Multiple R-Squared:", 
        format(signif(x$r.squared, digits)),"\n")

  correl <- x$correlation

  if(!is.null(correl)) {
    p <- dim(correl)[2]
    if(p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      correl <- format(round(correl, digits), ...)
      correl[col(correl) > row(correl)] <- ""
      print(correl, quote=F, ...)
    }
  }

  if (!is.null(x$biasTest)) {
    print(x$biasTest)
  }

  invisible(x)
}


summary.lmRob <- function(object, correlation = T, bootstrap.se = F) {
  wt <- object$M.weights
  wt1 <- object$weights

  if(!is.null(wt1) && !is.null(wt))
    wt <- wt * wt1

  coef <- coefficients(object)
  cnames <- labels(coef)
  ctotal <- object$coef
  ptotal <- length(ctotal)
  resid <- object$residuals
  fv <- object$fitted
  n <- length(resid)
  p <- object$rank

  if(is.null(p)) 
    p <- sum(!is.na(ctotal))

  if(any(na <- is.na(coef))) {
    coef <- coef[!na]
    p <- length(coef)
  }

  rdf <- object$df.resid

  if(is.null(rdf)) 
    rdf <- n - p

  if(!is.null(wt1)) {
    wt1 <- wt1^0.5
    resid <- resid * wt1
    fv <- fv * wt1
    excl <- wt1 == 0

    if(any(excl)) {
      warning(paste(sum(excl),"rows with zero weights not counted"))
      resid <- resid[!excl]
      fv <- fv[!excl]
      wt1 <- wt1[!excl]
      if(is.null(object$df.residual))
        rdf <- rdf - sum(excl)
      wt <- wt * wt1
    }
  }

  stddev <- object$scale
  cov <- object$cov
  var <- diag(cov)

  if(p < ptotal)
    R <- R[1:p, 1:p, drop = F]

  if(correlation) {
    correl <- cov/sqrt(var)
    correl <- t(correl)/sqrt(var)
  }

  else 
    correl <- NULL

  coef <- array(coef, c(p, 4))
  dimnames(coef) <- list(cnames, 
                         c("Value", "Std. Error", "t value", "Pr(>|t|)"))
  coef[, 2] <- sqrt(var)
  coef[, 3] <- coef[, 1]/coef[, 2]
  coef[, 4] <- if(rdf > 0) 2 * (1 - pt(abs(coef[, 3]), rdf)) else NA
  yy <- fv + resid

  if (is.null(object$robust.control))
    testbias <- T

  else {
    est <- object$est
    fnl <- object$robust.control$final.alg
    if ((casefold(est) == "final") && 
        (casefold(fnl) == "mm" || casefold(fnl) == "m"))
      testbias <- T
    else
      testbias <- F
  }

  if(testbias)
    biasTest <- test.lmRob(object)

	if(bootstrap.se) {
		bootstrap.se <- rb(object)$se
	}

  int <- attr(object$terms, "intercept")
  object <- object[c("call", "terms", "iter.final.coef", "iter.refinement", 
                     "M.weights", "r.squared", "est","robust.control",
                     "genetic.control","na.action")]

  object$residuals <- resid
  object$coefficients <- coef
  object$sigma <- stddev
	if(bootstrap.se[1])
		object$bootstrap.se <- bootstrap.se
  object$df <- c(p, rdf, ptotal)
  object$cov.unscaled <- cov/stddev^2
  object$correlation <- correl

  if (testbias)
    object$biasTest <- biasTest
  oldClass(object) <- c("summary.lmRob", "summary.lm")

  object
}


add1.lmRob <- function(object, scope=.~., scale, keep, x = NULL, ...)
{
  tmp <- object$robust.control$final.alg

  if (tmp == "Adaptive" || tmp == "adaptive")
    stop("add1 is only available for final MM-estimates.")

  p <- length(object$coef)

  if(!is.character(scope))
    scope <- add.scope(object, 
             update.formula(object, scope, evaluate = F))

  if(!length(scope))
    stop("no terms in scope for adding to object")

  if(is.null(x)) {
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(object, add.rhs, evaluate = F)
    fc <- object$call
    Terms <- terms(new.form)
    fc$formula <- Terms
    fob <- list(call = fc)
    oldClass(fob) <- oldClass(object)
    m <- model.frame(fob)
    x <- model.matrix(Terms, m, contrasts = object$contrasts)
  }

  y <- object$fitted+object$residuals
  cnames <- dimnames(x)[[2]]
  asgn <- attr(x, "assign")
  tl <- names(asgn)

  if(!all(match(scope, tl, F)))
    stop("scope is not a subset of term labels of the supplied x")

  xasgn <- unlist(asgn[names(object$assign)])
  asgn <- asgn[scope]
  k <- length(scope)

  if(missing(scale)) 
    scale <- object$scale

  if(!missing(keep)) {
    max.keep <- c("coefficients", "fitted", "residuals")

    if(is.logical(keep) && keep) 
      keep <- max.keep

    else {
      if(!all(match(keep, max.keep, F)))
        stop(paste("Can only keep one or more of: \"",
             paste(max.keep, collapse="\", \""), "\"", sep=""))
    }
  }

  else
    keep <- character(0)

  value <- array(vector("list", 3*k), c(k, 3),
                 list(scope, c("coefficients", "fitted", "residuals")))
  psif <- object$robust.control$weight
  efficiency <- object$robust.control$efficiency

  if(casefold(psif[2]) == "optimal") {
    ipsi <- 1
    if (efficiency == 0.95)      yc <- 1.060158
    else if (efficiency == 0.9)  yc <- 0.9440982
    else if (efficiency == 0.85) yc <- 0.8683765
    else if (efficiency == 0.8)  yc <- 0.8097795
    else                         yc <- lmRob.effvy(efficiency)
  }

  else {
    ipsi <- 2
    if (efficiency == 0.95)      yc <- 4.685061
    else if (efficiency == 0.9)  yc <- 3.882646
    else if (efficiency == 0.85) yc <- 3.443689
    else if (efficiency == 0.8)  yc <- 3.136909
    else                         yc <- chb(efficiency)$cb
  }

  if (object$est == "initial")
    warning("Inference based on initial estimates is not recommended.")
	rfpe.none <- lmRob.RFPE(object, scale)
  dfs <- double(k)
  rfpe <- double(k)

  if(length(xasgn)) {
    oldx <- x[, xasgn, drop = F]
    old.names <- dimnames(oldx)[[2]]
  }

  else 
    stop("need a term or an intercept in initial model")

	newnames <- scope

  for(i in 1:k) {
    ii <- asgn[[i]]
    pii <- length(ii)
    dfs[i] <- pii
    cur.name <- newnames[i]
		curobj <- update(object, as.formula(paste(". ~ . + ", cur.name, sep = "")))
		rfpe[i] <- lmRob.RFPE(curobj, scale)

    if(length(keep)) {
      cur.coef <- curobj$coefficients
      names(cur.coef) <- c(old.names,newnames[i])
      cur.fitted <- curobj$fitted
      cur.resid <- curobj$residuals
      value[i,1] <- list(cur.coef)
      value[i,2] <- list(cur.fitted)
      value[i,3] <- list(cur.resid)
    }
  }

  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  rfpe <- c(rfpe.none, rfpe)
  dfs[1] <- NA
  aod <- data.frame(Df=dfs, RFPE=rfpe, row.names=scope, check.names=F)
  head <- c("Double term additions", "\nModel:",
            deparse(as.vector(formula(object))))

  if(!missing(scale))
    head <- c(head, paste("\nscale: ", format(scale), "\n"))

  oldClass(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head

  if(length(keep))
    list(anova=aod, keep=structure(value[,keep,drop=F],class="matrix"))
  else 
    aod
}


anova.lmRob <- function(object, ..., test=c("RF","RWald")) 
{
  margs <- function(...) {nargs()}
  fun <- function(assign, coeff) {sum(!is.na(coeff[assign]))}
  
  test  <- match.arg(test) 
  psif <- object$robust.control$weight
  efficiency <- object$robust.control$efficiency
  tmp <- object$robust.control$final.alg

  if(casefold(tmp) == "adaptive") {
    if (test == "RF")
      stop("Robust F-test is only available for final MM-estimates.")
  }

  else {
    if (casefold(psif[2]) == "optimal") {
      ipsi <- 1
      if (efficiency == 0.95) {
        cst <- 0.976; yc <- 1.060158 
      }
      else if (efficiency == 0.9) {
        cst <- 0.963; yc <- 0.9440982
      }
      else if (efficiency == 0.85) {
        cst <- 0.953; yc <- 0.8684 
      }
      else if (efficiency == 0.8) {
        cst <- 0.944; yc <- 0.8097795 
      }
      else {
        cst <- lmRob.const(efficiency, ipsi)
        yc <- lmRob.effvy(efficiency) 
      }
    }

    else {
      ipsi <- 2
      if (efficiency == 0.95) {
        cst <- 0.218; yc <- 4.685061 
      }
      else if (efficiency == 0.9) {
        cst <- 0.295; yc <- 3.882646 
      }
      else if (efficiency == 0.85) {
        cst <- 0.357; yc <- 3.443689 
      }
      else if (efficiency == 0.8) {
        cst <- 0.413; yc <- 3.136909 
      }
      else {
        cst <- lmRob.const(efficiency, ipsi)
        yc <- chb(efficiency)$cb 
      }
    }
  }

  if(margs(...))
    return(anova.lmRob.list(list(object, ...),cst,ipsi,yc,test=test))

  if (object$est == "initial")
    warning("Inference based on initial estimates is not recommended.")

  cov <- object$cov
  coef <- object$coef
  np <- length(coef)
  assg <- object$assign

  if(is.null(assg))
    assg <- attributes(object$terms)$assign

  df <- sapply(assg, fun, object$coef)
  nassg <- names(assg)
  names(df) <- nassg

  if (test == "RWald") {
    aod <- as.matrix(round(df, 1))
    rnames <- c("Chisq Df", "Wald", "P(>Wald)")
    aod <- cbind(aod, NA, NA)
    j <- length(df)
    if (j>1) {
      for (i in 2:j) {
        ch.val   <- coef[i]*coef[i]/cov[i,i]
        aod[i,2] <- ch.val
        aod[i,3] <- 1 - pchisq(ch.val,1) 
      } 
    }
  } 

  else {
    aod <- as.matrix(round(df, 1))
    rnames <- c("Chisq Df", "RobustF","Pr(F)")
    aod <- cbind(aod, NA, NA)
    j <- length(df)

    if (nassg[1] == "(Intercept)") 
      frmcar <- ".~1"

    else 
      frmcar <- paste(".~ -1 +",nassg[1])

    curfrm <- as.formula(frmcar)
    curobj <- update(object, curfrm)

    if (curobj$est == "final") 
      res <- curobj$residuals

    else  
      res <- curobj$T.residuals  

    if (j>2) {
      for (i in 2:(j-1)) {
        if (frmcar==".~1") 
          frmcar <- paste(".~",nassg[i])
        else 
          frmcar <- paste(frmcar,nassg[i],sep="+")
        curfrm <- as.formula(frmcar)
        curobj <- update(object,curfrm)
        if (curobj$est == "final") 
          Res <- curobj$residuals
        else  
          Res <- curobj$T.residuals  
        Scale <- curobj$scale 
        FTau <- 2*sum(chi.weight(res/Scale,ipsi,yc)-
                      chi.weight(Res/Scale,ipsi,yc)) 
        aod[i,2] <- FTau
        aod[i,3] <- 1 - pchisq(FTau/cst,1) 
        res <- Res
      }
    }

    if (object$est == "final") 
      Res <- object$residuals

    else  
      Res <- object$T.residuals
    Scale <- object$scale 
    FTau <- 2*sum(chi.weight(res/Scale,ipsi,yc)-
                  chi.weight(Res/Scale,ipsi,yc)) 
    aod[j,2] <- FTau  
    aod[j,3] <- 1 - pchisq(FTau/cst,1) 
  } 

  dimnames(aod) <- list(names(df), rnames)
  heading <- "\nTerms added sequentially (first to last)\n"
  aod <- as.anova(data.frame(aod, check.names = F), heading)
  aod
}


anova.lmRob.list <- function(object, const, ipsi, yc, 
    test=c("RWald","RF")) 
{
  diff.Rn2 <- function(term.coefs, term.cov) {
    t1 <- attr(term.coefs[[1]],"names")
    t2 <- attr(term.coefs[[2]],"names")
    m1 <- match(t1, t2, F)
    m2 <- match(t2, t1, F)
    if(all(m1)) {
      if(all(m2)) 
        return(list(effects="No sub-model",
                    df=NA, Chi.val=NA, Prob.Chi=NA))
      else {
        lab0 <- t2[-m1]
        cov0 <- term.cov[[2]][lab0,lab0]
        coef0 <- term.coefs[[2]][lab0]
        inv0 <- solve(cov0)
        coef0 <- as.matrix(coef0)
        val <- as.numeric( t(coef0) %*% inv0 %*% coef0 )
        df <- length(lab0)
        P <- 1 - pchisq(val, df)

        return(list(effects=paste(c("", t2[ - m1]), collapse = "+"), 
               df=df, Chi.val=val, Prob.Chi=P))
      }
    }

    else {
      if(all(m2)) {
        lab0 <- t1[-m2]
        cov0 <- term.cov[[1]][lab0,lab0]
        coef0 <- term.coefs[[1]][lab0]
        inv0 <- solve(cov0)
        coef0 <- as.matrix(coef0)
        val <- as.numeric( t(coef0) %*% inv0 %*% coef0 )
        df <- length(lab0)
        P <- 1 - pchisq(val, df)

        return(list(effects=paste(c("", t1[ - m2]), collapse = "-"), 
               df=df, Chi.val=val, Prob.Chi=P))
      }

      else return(list(effects="No sub-model",
                  df=NA, Chi.val=NA, Prob.Chi=NA))
    }
  }

  diff.Tau <- function(res, scale, cname, ipsi, yc, const) {     
    t1 <- cname[[1]]
    t2 <- cname[[2]]
    m1 <- match(t1, t2, F)
    m2 <- match(t2, t1, F)
    if(all(m1)) {
      if(all(m2)) 
        return(list(effects="No sub-model",
                    df=NA, Tau.val=NA, Prob.Tau=NA))
      else {
        lab0 <- t2[-m1]
        Scale <- scale[[2]]
        df <- length(lab0)
        FTau <- (2/df)*sum(chi.weight(res[[1]]/Scale,ipsi,yc)-
                            chi.weight(res[[2]]/Scale,ipsi,yc))
        P <- 1 - pchisq(FTau/const, df)

        return(list(effects=paste(c("", t2[ - m1]), collapse = "+"), 
               df=df, Tau.val=FTau, Prob.Tau=P))
      }
    }

    else {
      if(all(m2)) {
        lab0 <- t1[-m2]
        Scale <- scale[[1]]
        df <- length(lab0)
        FTau <- (2/df)*sum(chi.weight(res[[2]]/Scale,ipsi,yc)-
                            chi.weight(res[[1]]/Scale,ipsi,yc))
        P <- 1 - pchisq(FTau/const, df)

        return(list(effects=paste(c("", t1[ - m2]), collapse = "-"), 
               df=df, Tau.val=FTau, Prob.Tau=P))
      }
      else return(list(effects="No sub-model",
                  df=NA, Tau.val=NA, Prob.Tau=NA))
    }
  }

  forms <- sapply(object, function(x) as.character(formula(x)))
  subs <- as.logical(match(forms[2,  ], forms[2, 1], F))

  if(!all(subs)) 
    warning(
    "Some fit objects deleted because response differs from the first model")

  if(sum(subs) == 1)
    stop("The first model has a different response from the rest")

  forms <- forms[, subs]
  object <- object[subs]
  estype <- sapply(object, "[[", "est")
  subs <- estype == "final"  

  if(!all(subs)) {
    warning(paste("Inference based on initial estimates is followed by",
                 " a (*) sign, and not recommended.", sep=""))
  }

  rt <- length(object)
  if(rt == 1) {
    object <- object[[1]]
    UseMethod("anova")
  }

  twrn <- rep("   ", rt)
  effects <- character(rt)
  Probchi <- rep(NA,rt)
  Df <- rep(NA,rt)
  Statval <- rep(NA,rt)
  heading <- paste("\nResponse: ", forms[2, 1], sep = "")

  if(test == "RWald") {
    tc <- list()
    tk <- list()
    for (i in 1:rt) {
      tc[[i]] <- object[[i]]$coef
      tk[[i]] <- object[[i]]$cov
      if (object[[i]]$est != "final") 
        twrn[i] <- "(*)"
    }

    for(i in 2:rt) {
      j <- c(i-1, i)
      diff.val <- diff.Rn2(tc[j],tk[j]) 
      effects[i] <- diff.val$effects
      Df[i] <- diff.val$df
      Statval[i] <- diff.val$Chi.val
      Probchi[i] <- diff.val$Prob.Chi 
    }
    aod <- data.frame(Terms=forms[3,  ], "   "=twrn, Df=Df,
           "Wald"=Statval, "P(>Wald)"=Probchi, check.names=F)
  }

  else {#Tau - test, Robust F
    rs <- list()
    sc <- list()
    cn <- list()
    for (i in 1:rt) {
      sc[[i]] <- object[[i]]$scale 
      cn[[i]] <- names(object[[i]]$coef)
      rs[[i]] <- object[[i]]$residuals
      if (object[[i]]$est != "final") 
        twrn[i] <- "(*)"
    }    

    for(i in 2:rt) {
      j <- c(i-1, i)
      diff.val <- diff.Tau(rs[j],sc[j],cn[j],ipsi,yc,const) 
      effects[i] <- diff.val$effects
      Df[i] <- diff.val$df
      Statval[i] <- diff.val$Tau.val
      Probchi[i] <- diff.val$Prob.Tau
    }
    aod <- data.frame(Terms=forms[3,  ], "   "=twrn, Df=Df,
           "RobustF"=Statval, "Pr(F)"=Probchi, check.names=F)
  }
  aod <- as.anova(aod, heading)
  aod
}


update.lmRob <- function(object, formula, evaluate = T, class, ...)
{
  if(missing(formula)) {

    if(is.null(object$T.coefficients))
      stop("There is nothing to update.")

    ans <- object

    if(object$est == "initial")
      ans$est <- "final"

    else 
      ans$est <- "initial"

    ans$coefficients <- object$T.coefficients
    ans$T.coefficients <- object$coefficients
    ans$cov <- object$T.cov
    ans$T.cov <- object$cov
    ans$residuals <- object$T.residuals
    ans$T.residuals <- object$residuals
    ans$fitted.values <- object$T.fitted.values
    ans$T.fitted.values <- object$fitted.values

    if (casefold(object$robust.control$final) == "mm" ||
        casefold(object$robust.control$final) == "m") {
      ans$dev   <- object$T.dev
      ans$T.dev <- object$dev
      ans$r.squared   <- object$T.r.squared
      ans$T.r.squared <- object$r.squared
      ans$M.weights   <- object$T.M.weights
      ans$T.M.weights <- object$M.weights
    }
  }

  else
    ans <- NextMethod()

  ans
} 


#--------------------------------------------------------------#
# Utility Functions                                            #
#--------------------------------------------------------------#

numeric.model.matrix <- function(obj)
{
	## obj is similar to lm and lmRob.

	m <- model.frame(obj)
	regressors <- attr(obj$terms, "term.labels")
	regressors <- intersect(regressors, names(m))
	m <- m[, regressors, drop = F]
	numeric.idx <- sapply(m, is.numeric)
	m <- m[, numeric.idx, drop = F]
	as.matrix(m)
}


lmRob.lar <- function(x, y, tol=1e-6)
{
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  bet0 <- 0.773372647623  ## bet0 = pnorm(0.75)
  tmpn <- double(n)
  tmpp <- double(p)
  z1 <- .Fortran("s_larsbi",
                 x,
                 y,
                 as.integer(n),
                 as.integer(p),
                 as.integer(n),
                 as.integer(n),
                 as.double(tol),
                 NIT=integer(1),
                 K=integer(1),
                 KODE=integer(1),
                 SIGMA=double(1),
                 THETA=tmpn,
                 RS=tmpn,
                 SC1=tmpn,
                 SC2=tmpp,
                 SC3=tmpp,
                 SC4=tmpp,
                 BET0=as.double(bet0))
  list(coef=z1$THETA[1:p], scale=z1$SIGMA, resid=z1$RS)
}


lmRob.ucovcoef <- function(x, y, resid, sigma, p, ipsi, xk, tau, tl) 
{

##
## Unscaled covariance matrix of cofficients
##

  wi <- resid/sigma
  wcnd <- wi != 0 & !is.na(wi)
  wi[wcnd] <- psi.weight(wi[wcnd],ipsi,xk)/wi[wcnd]
  swi <- sum(wi[wcnd])

  if (abs(swi) <= tl) {
    ans <- list(wi=wi, cov=NA, fact=NA, ierr=1)
    return(ans)
  }

  z  <- as.vector(sqrt(wi))
  sx <- x * z
  np <- ncol(sx)
  n <- length(resid)
 
  storage.mode(sigma) <- "double"
  fact <- .Fortran("s_kffam2",
                   as.double(resid),
                   as.integer(n),
                   as.integer(p),
                   sigma,
                   fh=double(1),
                   as.integer(ipsi),
                   as.double(xk))$fh * swi

  storage.mode(sx) <- "double"
  zz <- .Fortran("s_rmtrm2",
                 x=sx,
                 as.integer(n),
                 np=as.integer(np),
                 as.integer(n),
                 intch=as.integer(1),
                 as.double(tau),
                 k=integer(1),
                 sf=double(np),
                 sg=double(np),
                 sh=double(np),
                 ip=integer(np))

  k  <- zz$k  
  xt <- zz$x
  ncov <- np*(np+1)/2
  zc <- .Fortran("s_kiasm2",
                 xt,
                 k=as.integer(k),
                 as.integer(np),
                 as.integer(n),
                 as.integer(ncov),
                 fu=as.double(1.0),
                 fb=as.double(1.0),
                 cov=double(ncov))$cov

  storage.mode(fact) <- "double"
  cov <- .Fortran("s_kfasm2",
                  xt,
                  cov=zc,
                  k=as.integer(k),
                  as.integer(np),
                  as.integer(n),
                  as.integer(ncov),
                  f=fact,
                  double(np),
                  zz$sg,
                  zz$ip)$cov

  attributes(wi) <- attributes(y)
  list(wi = wi, ucov = cov, fact = fact, ierr=0)
}


lmRob.wm <- function(x, y, coeff0, ucov0, scale0, itype=1, isigma=-1,
    ipsi=1, xk=0.9440982, beta=0.2661, wgt=y, tlo=0.0001, tua=1.5e-06, 
    mxr=50)
{

##
## W-algorithm for robust M-estimation
##

  n <- length(y)
  p <- ncol(x)
  sx <- matrix(double(1), n, p)
  ncov <- length(ucov0)
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"

  rList <- .Fortran("s_rwagm2",
                    x,
                    y,
                    theta=coeff0,
                    wgt=y,
                    cov=as.double(ucov0),
                    psp0=as.double(psp.weight(0,ipsi,xk)),
                    sigmai=as.double(scale0),
                    as.integer(n),
                    as.integer(p),
                    mdx=as.integer(n),
                    as.integer(ncov),
                    tol=as.double(tlo),
                    gam=as.double(1.0),
                    tau=as.double(tua),
                    itype=as.integer(1),
                    as.integer(isigma),
                    icnv=as.integer(1),
                    maxit=as.integer(mxr),
                    maxis=as.integer(1),
                    nit=integer(1),
                    sigmaf=double(1),
                    rs=double(n),
                    delta=double(p),
                    double(n),
                    double(p),
                    double(p),
                    double(p),
                    ip=integer(p),
                    double(n),
                    sx,
                    as.integer(ipsi),
                    as.double(xk),
                    as.double(beta),
                    as.double(1.0))
  rList
}


lmRob.RRvsRD <- function(RR, x, method=NULL, chisq.percent=0.975, 
    vertical.outlier=2.5, seed=1313, ...)
{
  x <- as.matrix(x)
  np <- ncol(x)
  tmp.n <- as.character(1:length(RR))
  set.seed(seed)

  if(is.null(method)) {
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

  plot(RD, RR, xlim=x.range, ylim=y.range, xlab="Robust Distances",
       ylab="Standardized Residuals")

  yp.idx <- (RR >  vertical.outlier)
  ym.idx <- (RR < -vertical.outlier)
  hjit <- (x.range[2]-x.range[1])/100

  if (any(yp.idx)) 
    text(RD[yp.idx]+hjit, RR[yp.idx], tmp.n[yp.idx], adj=0)

  if (any(ym.idx)) 
    text(RD[ym.idx]+hjit, RR[ym.idx], tmp.n[ym.idx], adj=0)

  abline(v = x.chisq, lty = 2)
  abline(h = vertical.outlier, lty = 2)
  abline(h = -vertical.outlier, lty = 2)

  invisible()
}


psi.weight <- function(svals, ips = 1, xk = 1.06) {
  n <- length(svals)
  fvals <- double(n)
  storage.mode(svals) <- "double"
  f.res <- .Fortran("s_psiam2",
                    n=as.integer(n),
                    svals=svals,
                    fvals=fvals,
                    as.integer(ips),
                    as.double(xk))
  f.res$fvals
}


rho.weight <- function(svals, ips = 1, xk = 1.06) {
  n <- length(svals)
  fvals <- double(n)
  storage.mode(svals) <- "double"
  f.res <- .Fortran("s_rhoam2",
                    n=as.integer(n),
                    svals=svals,
                    fvals=fvals,
                    as.integer(ips),
                    as.double(xk))
  f.res$fvals
}


psp.weight <- function(svals, ips = 1, xk = 1.06) {
  n <- length(svals)
  fvals <- double(n)
  storage.mode(svals) <- "double"
  f.res <- .Fortran("s_pspam2",
                    n=as.integer(n),
                    svals=svals,
                    fvals=fvals,
                    as.integer(ips),
                    as.double(xk))
  f.res$fvals
}


chi.weight <- function(svals, ips = 1, xk = 1.06) {
  n <- length(svals)
  fvals <- double(n)
  storage.mode(svals) <- "double"
  f.res <- .Fortran("s_chiam2",
                    n=as.integer(n),
                    svals=svals,
                    fvals=fvals,
                    as.integer(ips),
                    as.double(xk))
  f.res$fvals
}


lmRob.eff0 <- function(mu=1, itype=1, iucv=1, iwww=1, ialfa=1,
    sigmx=1.0, upper=10, til=1.e-4, maxit=150, tol=1.e-5,
    epmach=.Machine$double.eps, uflow=.Machine$double.xmin,
    ta=0, tb=0, tc=0, ipsi=3) 
{

##--------------------------------------------------
## ipsi=3: by default, use Huber's function
##--------------------------------------------------

  index <- c(mu,iwww,iucv,ipsi,itype,1,0)
  tc <- c(0,0,ta,tb,tc,epmach,uflow,upper,til)
  z <- .Fortran("s_ref0bi",
                as.integer(index),
                as.double(tc),
                xlcnst=as.double(-1),
                as.integer(ialfa),
                as.double(sigmx),
                as.integer(maxit),
                as.double(tol),
                nit=as.integer(0),
                alfa=as.double(0),
                beta=as.double(0),
                reff=as.double(0))
  list(nit=z$nit,alfa=z$alfa,beta=z$beta,reff=z$reff)
}


lmRob.effad <- function(eff)
{

##
## Computes the cutoff value for adaptive estimator given eff
##

  eff.func <- function(cc, eff)
    -0.531700164+0.785133187*cc-0.103931349*cc^2+0.000637741*cc^3-eff

  uniroot(eff.func, interval=c(2.0,3.9), eff=eff)$root
}

lmRob.effvy <- function(eff, ipsi = 1)
{

##
## Computes the tuning constant for optimal weight function given eff
##

  eff.func <- function(cc, eff, ipsi=1)
    lmRob.eff0(itype=1,ta=cc,tc=cc,ipsi=ipsi)$reff-eff

  if(ipsi == 1)      ## optimal function is used
    c.inv <- c(0.2, 2.5)

  else if(ipsi == 2) ## bisquare function is used
    c.inv <- c(0.1, 30)

  else if(ipsi == 3) ## huber function is used
    c.inv <- c(0.1, 3.5)

  uniroot(eff.func, interval = c.inv, eff = eff, ipsi = ipsi)$root
}

lmRob.const <- function(eff, ipsi = 1)
{

##
## Computes the factor used for robust tau (RF) test
##

  if(ipsi == 1) {
    if(eff == 0.95) cc <- 1.060158
    else if(eff == 0.9) cc <- 0.9440982
    else if(eff == 0.85) cc <- 0.8684
    else if(eff == 0.8) cc <- 0.8097795
    else cc <- lmRob.effvy(eff)
  }
  else {
    if(eff == 0.95) cc <- 4.685061
    else if(eff == 0.9) cc <- 3.882646
    else if(eff == 0.85) cc <- 3.443689
    else if(eff == 0.8) cc <- 3.136909
    else cc <- chb(eff)$cb
  }

  tmp <- lmRob.eff0(itype = 1, ta = cc, tc = cc, ipsi = ipsi)
  tmp$alfa/tmp$beta
}


lmRob.RFPE <- function(object, scale=NULL)
{
  if (object$est == "initial")
    warning("Inference based on initial estimates is not recommended.")

  tmp <- object$robust.control$final.alg

  if (tmp == "Adaptive" || tmp == "adaptive")
    stop("RFPE is only available for final MM-estimates.")

  p <- length(object$coef)

  if (is.null(scale))
		scale <- object$scale

  res <- object$residuals/scale
  psif <- object$robust.control$weight
  efficiency <- object$robust.control$efficiency

  if (casefold(psif[2]) == "optimal")
		ipsi <- 1
	else
		ipsi <- 2

	yc <- object$yc
  a <- sum(rho.weight(res,ipsi,yc))
  b <- p*sum(psi.weight(res,ipsi,yc)^2)
	d <- sum(psp.weight(res,ipsi,yc))
	if(d <= 0) return(NA)

	a + b/d
}

lmRob.ga <- function(x, y, rk, control, ips, xk, beta, intch=1,
    tolr=1.e-3, tau=1.e-6, maxs1=50) 
{

##
## Step 0. Extract Parameters
##

  n <- length(y)
  p <- ncol(x)
  smin <- double(1)
  theta <- double(p)
  rs <- double(n)
  popsize <- control$popsize
  mutate.prob <- control$mutate.prob
  random.n <- control$random.n
  births.n <- control$births.n
  stock <- control$stock
  maxslen <- control$maxslen
  stockprob <- control$stockprob
  nkeep <- control$nkeep

  if (is.null(popsize))
		popsize <- 10*rk

  if (is.null(mutate.prob))
    mutate.prob <- c(0.15,0.2,0.2,0.2)

  else {
    if (length(mutate.prob) != 4)
      stop("mutate.prob must have length 4.")
    if (any(mutate.prob < 0))
      stop("negative value in mutate.prob.")
    if (sum(mutate.prob[2:4]) > 1)
      stop("sum of last 3 mutation probabilities greater than 1.")
  }

  if (is.null(random.n)) 
    random.n <- 50*rk

  if (is.null(births.n)) 
    births.n <- 50*rk+15*rk*rk

  if (is.null(maxslen)) {
    np2 <- trunc((n-rk)/2)
    if (np2 <= rk) 
      maxslen <- rk
    else 
      maxslen <- min(5*rk, np2)
  }

  else if (maxslen < p)
    stop("maxslen is too small.")

  if (is.null(stockprob))
    stockprob <- cumsum((2*(popsize:1))/popsize/(popsize+1))

  else {
    if (length(stockprob) != popsize)
      stop("length of stockprob must be equal to popsize.")
    if (any(stockprob < 0) || any(diff(stockprob) < 0) || 
        is.character(all.equal.numeric(stockprob[popsize],1)))
      stop("stockprob must be cumulative probabilities.")
  }

##
## Step 1. Obtain the Stock and Others
##

  noldstock <- length(stock)
  stockmat <- matrix(0,maxslen,popsize)
  stock.len <- integer(popsize)

  if (noldstock) {
    noldstock <- min(noldstock, popsize)
    for (i in 1:noldstock) {
      si <- stock[[i]]
      ll <- length(si)
      if (ll > maxslen || ll < rk)
        stop(paste("length of component", i, "of stock is not between",
            rk, "and", maxslen))
      if (any(si > n) || any(si <= 0))
        stop(paste("bad observation number(s) in component",i,"of stock."))
      stockmat[1:ll,i] <- stock[[i]]
      stock.len[i] <- ll
    }
  }

  storage.mode(x) <- "double"
  storage.mode(stockmat) <- "integer"
  xx <- matrix(0,maxslen,p)
  storage.mode(xx) <- "double"
  yy <- double(maxslen)

  z <- .Fortran("s_genem2",
                x,
                as.double(y),
                as.integer(n),
                as.integer(rk),
                as.integer(popsize),
                as.double(mutate.prob),
                as.integer(random.n),
                as.integer(births.n),
                stock=stockmat,
                as.integer(maxslen),
                objective=double(popsize),
                integer(2*maxslen),
                stock.len=as.integer(stock.len),
                as.integer(noldstock),
                as.double(stockprob),
                as.integer(intch),
                as.double(tolr),
                as.double(tau),
                as.integer(maxs1),
                smin=smin,
                theta=theta,
                rs=rs,
                sz=double(n),
                integer(rk),
                double(rk),
                double(rk),
                xtheta=double(maxslen),
                yy=yy,
                double(rk),
                xx=xx,
                integer(maxslen),
                ips=as.integer(ips),
                xk=as.double(xk),
                beta=as.double(beta),
                as.double(1.0))

  ord <- order(z$objective)
  nkeep <- max(1,nkeep)

  if (nkeep > popsize) 
    nkeep <- popsize

  ord <- ord[1:nkeep]
  lengths <- z$stock.len[ord]
  stockmat <- z$stock
  stock <- vector("list",nkeep)

  for (i in 1:nkeep)
    stock[[i]] <- stockmat[1:lengths[i],ord[i]]

  list(theta=z$theta, smin=z$smin, rs=z$rs, objective=z$objective[ord],
       stock=stock, births.n=births.n)
}


gen.data <- function(coeff, n = 100, eps = 0.1, sig = 3,
							snr = 1/20, seed = 837)
{

# coeff : 3 x 1 vector of coefficients
# eps   : the contamination ratio, between 0 and 0.5
# sig   : standard deviation of most observations
# snr   : signal-to-noise ratio, well, not really
# Note  : the regressors are generated as: rnorm(n,1),
#         rnorm(n,1)^3, exp(rnorm(n,1)). It also
#         generates an unused vector x4.

  set.seed(seed)
  x <- cbind(rnorm(n, 1), rnorm(n, 1)^3, exp(rnorm(n, 1)))
  ru <- runif(n)
  n1 <- sum(ru < eps)
  u <- numeric(n)
  u[ru < eps] <- rnorm(n1, sd = sig/snr)
  u[ru > eps] <- rnorm(n - n1, sd = sig)

  data.frame(y = x %*% matrix(coeff, ncol = 1) + u,
    x1 = x[,1], x2 = x[,2], x3 = x[,3], x4 = rnorm(n, 1))
}






























