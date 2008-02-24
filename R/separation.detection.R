`separation.detection` <-
function (fit, nsteps = 30) 
{
    fit.class <- class(fit)[1]
    if (fit.class != "glm") 
        stop("Only objects of class 'glm' are accepted.")
    eps <- .Machine$double.eps
    stdErrors <- matrix(0, nsteps, length(fit$coef))
    for (i in 1:nsteps) {
        suppressWarnings(temp.fit <- update(fit, control = glm.control(maxit = i, 
            epsilon = eps)))
        stdErrors[i, ] <- summary(temp.fit)$coef[, "Std. Error"]
    }
    res <- sweep(stdErrors, 2, stdErrors[1, ], "/")
    colnames(res) <- names(coef(fit))
    res
}
