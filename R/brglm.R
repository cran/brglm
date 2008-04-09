## * 'brglm' and 'brglm.fit' were written using as basis the code
##   or 'glm' and 'glm.fit', respectively.
## * 'print.brglm' is a modification of 'print.glm'
## * 'summary.brglm' is a modification of 'summary.brglm'
## * 'print.summary.brglm' is a modification of 'print.summary.glm'
## Ioannis Kosmidis <I.Kosmidis@warwick.ac.uk> [15/02/2008]
`brglm` <-
function (formula, family = binomial, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control.glm = glm.control1(...), 
    model = TRUE, method = "brglm.fit", pl = FALSE, x = FALSE, 
    y = TRUE, contrasts = NULL, control.brglm = brglm.control(...), 
    ...) 
{
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    br <- method == "brglm.fit"
    ####################
    ## More families to be implemented
    if (br & family$family != "binomial") 
        stop("families other than 'binomial' are not currently implemented")
    ####################
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = fit.proc <- glm.fit, 
        brglm.fit = fit.proc <- brglm.fit, stop("invalid 'method': ", 
            method))
    ####################
    ## Arg control of fit.proc
    if (br) {
        formals(fit.proc)$control.brglm <- control.brglm
    }
    if (pl) 
        formals(fit.proc)$pl <- TRUE
    ####################
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    Xor <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    Xmax <- apply(abs(Xor), 2, max)
    Xmax[Xmax==0] <- 1
    X <- sweep(Xor, 2, Xmax, "/")
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(Y))
        else if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    par.names <- colnames(X)
    fit <- fit.proc(x = X, y = Y, weights = weights, start = start, 
        etastart = etastart, mustart = mustart, offset = offset, 
        family = family, control = control.glm, intercept = attr(mt, 
            "intercept") > 0)
    if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
            y = Y, weights = weights, offset = offset, family = family, 
            control = control.glm, intercept = TRUE)$deviance
    }
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    ## Move back to the original scale
    if (nPars <- ncol(X)) {
        redundant <- if (br) 
            fit$redundant
        else rep.int(0, nPars)
        fit$coefficients <- fit$coefficients/Xmax[!redundant]
        #fit$qr <- qr(sqrt(fit$weights) * Xor[, !redundant])  
        fit$qr <- qr(sqrt(fit$weights) * Xor)
        if (br) {
            fit$FisherInfo <- fit$FisherInfo * tcrossprod(Xmax[!redundant])
            fit$control.brglm <- control.brglm
        }
        ####################
        ## Aliasing
        coefs <- rep(NA, ncol(X))
        names(coefs) <- par.names
        coefs[!redundant] <- fit$coefficients
        fit$coefficients <- coefs
        ####################        
    }
    fit$control.glm <- control.glm
    if (x) 
        fit$x <- Xor
    if (!y) 
        fit$y <- NULL
    if (br) 
        fit$penalized.deviance <- if (all(family$link == "logit") | 
            pl) 
            fit$deviance - log(det(fit$FisherInfo))
        else NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, method = method, pl = pl, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c("brglm", "glm", "lm")
    fit
}

`brglm.fit` <-
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = binomial(), 
    control = glm.control(), control.brglm = brglm.control(), 
    intercept = TRUE, pl = FALSE) 
{
    x <- as.matrix(x)
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkfun <- family$linkfun
    dmu.deta <- family$mu.eta
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object")
    if (EMPTY) {
        return(glm.fit(x = x, y = y, weights = weights, start = start, 
            etastart = etastart, mustart = mustart, offset = offset, 
            family = family, control = control, intercept = intercept))
    }
    valideta <- family$valideta
    if (is.null(valideta)) 
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu)) 
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    options(warn = -1)
    cur.repr <- modifications(family, pl = pl)
    y.count <- y * weights
    wt <- weights + 1
    y.adj <- (y.count + 0.5)/wt
    temp.fit <- glm.fit(x = x, y = y.adj, weights = wt, start = start, 
        etastart = etastart, mustart = mustart, offset = offset, 
        family = family, control = control, intercept = intercept)
    redundant <- is.na(temp.fit$coefficients)
    if (any(redundant)) {
        x <- x[, -which(redundant), drop = FALSE]
        nvars <- nvars - sum(redundant)
    }
    nIter <- 0
    test <- TRUE
    x.t <- t(x)
    while (test & (nIter < control.brglm$br.maxit)) {
        nIter <- nIter + 1
        ps <- temp.fit$fitted.values
        etas <- linkfun(ps)
        ww <- temp.fit$weights/wt * weights
        W.X <- sqrt(ww) * x
        XWXinv <- chol2inv(chol(crossprod(W.X)))
        hats <- gethats(nobs, nvars, x.t, XWXinv, ww)
        #hats <- diag(x%*%XWXinv%*%t(ww * x))
        cur.model <- cur.repr(ps)
        wt <- weights + hats * cur.model$at
        y.adj <- (y.count + hats * cur.model$ar)/wt
        temp.fit <- glm.fit(x = x, y = y.adj, weights = wt, etastart = etas, 
            offset = offset, family = family, control = control, 
            intercept = intercept)
        modscore <- t(dmu.deta(etas)/variance(ps) * x) %*% ((y.adj - 
            ps) * wt)
        if (control.brglm$br.trace) {
            cat("Iteration:", nIter, "\n")
            cat("Modified scores:", modscore, "\n")
        }
        test <- sum(abs(modscore)) > control.brglm$br.epsilon
    }
    options(warn = 0)
    temp.fit$converged <- nIter < control.brglm$br.maxit
    if (!temp.fit$converged) 
        warning("Iteration limit reached")
    temp.fit$ModifiedScores <- c(modscore)
    ww <- temp.fit$weights/wt * weights
    temp.fit$weights <- ww
    W.X <- sqrt(ww) * x
    temp.fit$FisherInfo <- crossprod(W.X)
    XWXinv <- chol2inv(chol(temp.fit$FisherInfo))
    temp.fit$hats <- gethats(nobs, nvars, x.t, XWXinv, ww)
    temp.fit$qr <- qr(W.X)
    temp.fit$nIter <- nIter
    temp.fit$prior.weights <- weights
    temp.fit$y <- y
    temp.fit$deviance <- sum(dev.resids(temp.fit$y, temp.fit$fitted.values, 
        temp.fit$prior.weights))
    temp.fit$cur.model <- cur.model
    temp.fit$redundant <- redundant
    temp.fit
}


`print.brglm` <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (x$method == "glm.fit" | !(nPars <- length(coef(x)))) 
        return(print.glm(x, digits, ...))
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (nPars) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts)) 
            cat("  [contrasts: ", apply(cbind(names(co), co), 
                1, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    if (!is.null(x$penalized.deviance)) 
        cat("Deviance:\t   ", format(round(x$deviance, digits)), 
            "\nPenalized Deviance:", format(round(x$penalized.deviance, 
                digits)), "\tAIC:", format(round(x$aic, digits)), 
            "\n")
    else cat("Deviance:\t   ", format(round(x$deviance, digits)), 
        "\tAIC:", format(round(x$aic, digits)), "\n")
    invisible(x)
}
`print.summary.brglm` <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    }
    else {
        df <- if ("df" %in% names(x)) 
            x[["df"]]
        else NULL
        if (!is.null(df) && (nsingular <- df[3] - df[1])) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null", 
            "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", 
            "deviance")]), digits = max(5, digits + 1)), " on", 
            format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 
            1, paste, collapse = " "), sep = "")
    if (!is.null(x$penalized.deviance)) 
        cat("Penalized deviance:", format(round(x$penalized.deviance, 
            digits = max(5, digits + 1))), "\n")
    if (nchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), 
        "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}
`summary.brglm` <-
function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, 
    ...) 
{
    if (object$method == "glm.fit") 
        return(summary.glm(object, dispersion = NULL, correlation = FALSE, 
            symbolic.cor = FALSE, ...))
    df.r <- object$df.residual
    if (is.null(dispersion)) 
        dispersion <- 1
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(chol(object$FisherInfo))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", 
            "Pr(>|z|)"))
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        df.f <- length(aliased)
    }
    keep <- match(c("call", "terms", "family", "deviance", "aic", 
        "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter", "na.action", "penalized.deviance"), names(object), 
        0)
    ans <- c(object[keep], list(deviance.resid = residuals(object, 
        type = "deviance"), coefficients = coef.table, aliased = aliased, 
        dispersion = dispersion, df = c(object$rank, df.r, df.f), 
        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.brglm"
    return(ans)
}
