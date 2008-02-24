`print.profile.brglm` <-
function (x, ...) 
{
    cat("'level' was set to", attr(x, "level"), "\n")
    cat("Methods that apply:\n")
    cat("'confint'  'plot' 'pairs'\n")
}
`profile.brglm` <-
function (fitted, gridsize = 10, stdn = 5, stepsize = 0.5, level = 0.95, 
    which = 1:length(coef(fitted)), verbose = TRUE, zero.bound = 1e-08, 
    scale = FALSE, ...) 
{
    if (level <= 0 | level >= 1) 
        stop("invalid 'level'.")
    if (fitted$method == "glm.fit") {
        if (verbose) 
            cat("Profiling the ordinary deviance for the supplied fit...\n")
        res1 <- profileModel(fitted, gridsize = gridsize, stdn = stdn, 
            stepsize = stepsize, grid.bounds = NULL, quantile = qchisq(level, 
                1), objective = "ordinaryDeviance", agreement = TRUE, 
            verbose = FALSE, trace.prelim = FALSE, which = which, 
            profTraces = TRUE, zero.bound = zero.bound, scale = scale, 
            stdErrors = summary(fitted)$coefficients[, 2])
        res2 <- NULL
    }
    else {
        fitted1 <- update(fitted, method = "glm.fit")
        if (verbose) 
            cat("Profiling the ordinary deviance for the corresponding ML fit...\n")
        res1 <- profileModel(fitted1, gridsize = gridsize, stdn = stdn, 
            stepsize = stepsize, grid.bounds = NULL, quantile = qchisq(level, 
                1), objective = "ordinaryDeviance", agreement = TRUE, 
            verbose = FALSE, trace.prelim = FALSE, which = which, 
            profTraces = TRUE, zero.bound = zero.bound, scale = scale, 
            stdErrors = summary(fitted1)$coefficients[, 2])
        if (fitted$pl | all(fitted$family$link == "logit")) {
            if (verbose) 
                cat("Profiling the penalized deviance for the supplied fit...\n")
            res2 <- profileModel(fitted, gridsize = gridsize, 
                stdn = stdn, stepsize = stepsize, grid.bounds = NULL, 
                quantile = qchisq(level, 1), objective = "penalizedDeviance", 
                agreement = TRUE, verbose = FALSE, trace.prelim = FALSE, 
                which = which, profTraces = TRUE, zero.bound = zero.bound, 
                scale = scale, stdErrors = summary(fitted)$coefficients[, 
                  2], X = model.matrix(fitted))
        }
        else {
            if (verbose) 
                cat("Profiling the modified score statistic for the supplied fit...\n")
            res2 <- profileModel(fitted, gridsize = gridsize, 
                stdn = stdn, stepsize = stepsize, grid.bounds = NULL, 
                quantile = qchisq(level, 1), objective = "modifiedScoreStatistic", 
                agreement = TRUE, verbose = FALSE, trace.prelim = FALSE, 
                which = which, profTraces = TRUE, zero.bound = zero.bound, 
                scale = scale, stdErrors = summary(fitted)$coefficients[, 
                  2], X = model.matrix(fitted))
        }
    }
    res <- list(profilesML = res1, profilesBR = res2)
    attr(res, "level") <- level
    class(res) <- "profile.brglm"
    res
}
