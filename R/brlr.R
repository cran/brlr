brlr <-
    function (formula, data = NULL, offset, weights, start, ...,
              subset, dispersion = 1, na.action = na.omit,
              contrasts = NULL, x = FALSE, br = TRUE,
              control = list(maxit = 200)) 
{
    glimlog <- function(y) {  
    # log function as in GLIM
    # to cope with zero counts and zero denominators
        ifelse(y == 0 | is.na(y), 0, log(y))
    }
    fmin <- function(beta) {
        eta <- offset. + drop(x %*% beta)
        pr <- plogis(eta)
        w <- wt * denom * pr * (1 - pr)
        detinfo <- det(t(x) %*% sweep(x, 1, w, "*"))
        if (all(pr > 0) && all(pr < 1) && detinfo > 0) {
            sum(wt * (y * glimlog(y/(denom * pr)) +
                     (denom - y) *
                        glimlog((denom - y)/(denom * (1 - pr))))
                ) - 
                br * 0.5 * log(detinfo)
        } else Inf
    }
    gmin <- function(beta) {
        eta <- offset. + drop(x %*% beta)
        pr <- plogis(eta)
        h <- hat(sweep(x, 1, sqrt(wt * denom * pr * (1 - pr)), 
                       "*"), intercept = FALSE)
        -drop((wt * (y + br * h/2 - pr * (denom + br * h))) %*% 
              x)
    }
    if (inherits(formula, "glm")){
        call <- formula$call
        if (is.null(call$family) || call$family != as.name("binomial")) stop(
            "model is not a binomial glm")
        call[[1]] <- as.name("brlr")
        return(eval(call, parent.frame()))}
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$start <- m$br <- m$control <- m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    keep.xmat <- x
    x <- model.matrix(Terms, m, contrasts)
    if (ncol(x) == 0) {  ## no estimation to do
        thecall <- match.call(expand.dots = TRUE)
        thecall[[1]] <- as.name("glm")
        return(eval(thecall, parent.frame()))} 
    xmax <- apply(abs(x), 2, max)
    x.unscaled <- x
    x <- sweep(x, 2, xmax, "/")
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    wt <- model.extract(m, weights)
    n <- nrow(x)
    if (!length(wt)) 
        wt <- rep(1, n)
    offset <- model.extract(m, offset)
    offset. <- if (is.null(offset)) rep(0, n) else offset 
    if (length(offset.) != n) stop("offset has wrong length")
    y <- model.extract(m, response)
    denom.adj <- denom <- rep(1, n)
    if (is.factor(y) && nlevels(y) == 2) y <- as.numeric(y) - 1
    if (is.matrix(y) && ncol(y) == 2 && is.numeric(y)) {
        denom <- as.vector(apply(y, 1, sum))
        denom.adj <- denom + (denom < 0.01)
        y <- as.vector(y[, 1])
    }
    ow <- options("warn")[[1]]
    options(warn = -1)
    fit <- glm.fit(x, y/denom.adj, wt * denom, family = binomial(), 
                   offset = offset.,
                   control = glm.control(maxit = 1))
    pr <- fit$fitted
    eta <- qlogis(pr) - offset.
    w <- wt * denom * pr * (1 - pr)
    leverage <- hat(sweep(x, 1, sqrt(w), "*"), intercept = FALSE)
    z <- eta + (y + leverage/2 - (denom + leverage) * pr)/w
    fit <- lm.wfit(x, z, w)
    options(warn = ow)
    est.start <- fit$coefficients
    resdf <- fit$df.residual
    nulldf <- fit$df.null
    if (missing(start)) start <- est.start
    redundant <- is.na(est.start)
    xstored <- x
    xmax.stored <- xmax
    if (any(redundant)){
        x <- x[, -which(redundant), drop = FALSE]
        xmax <- xmax[-which(redundant)]
        start <- start[-which(redundant)]}
    fstart <- fmin(start)
    parscale <- fstart/(1e-8 + abs(fstart -
                     sapply(seq(along = start), function(r) {
                           fmin(start + (seq(along = start) == r))}
                            )))
    control <- c(control, list(parscale = parscale, fnscale = abs(fstart)))
    res <- optim(start, fmin, gmin,
                 method = "BFGS",
                 control = control)
    beta <- res$par/xmax
    x <- sweep(x, 2, xmax, "*")
    penalized.deviance <- NULL
    if (br) {
        penalized.deviance <- 2 * fmin(beta)
        br <- FALSE
        deviance <- 2 * fmin(beta)
        br <- TRUE
    }
    else deviance <- 2 * fmin(beta)
    niter <- c(f.evals = res$counts[1],
               g.evals = res$counts[2])
    names(beta) <- colnames(x)
    eta <- as.vector(x %*% beta)
    lp <- offset. + eta
    pr <- plogis(lp)
    convergence <- if (res$convergence == 0) 
        TRUE
    else res$convergence
    h <- hat(sweep(x, 1, sqrt(wt * denom * pr * (1 - pr)), "*"), 
             intercept = FALSE)
    coefs <- est.start
    coefs[!is.na(coefs)] <- beta
    fit <- list(coefficients = coefs,
                deviance = deviance,
                penalized.deviance = penalized.deviance, 
                fitted.values = pr,
                linear.predictors = lp,
                call = match.call(), 
                formula = formula,
                convergence = convergence, 
                niter = niter,
                df.residual = resdf,
                df.null = nulldf, 
                model = m,
                y = y/denom,
                family = binomial(),
                offset = offset, 
                prior.weights = wt * denom,
                weights = wt * denom * pr * (1 - pr),
                terms = Terms,
                dispersion = dispersion, 
                bias.reduction = br,
                leverages = h,
                control = control)
    class(fit) <- c("brlr", "glm", "lm")
    W <- fit$weights
    fit$qr <- qr(model.matrix(fit) * sqrt(W))
    fit$rank <- fit$qr$rank
    fit$FisherInfo <- t(x) %*% sweep(x, 1, W, "*")
    attr(fit, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        fit$na.action <- attr(m, "na.action")
    fit$contrasts <- attr(x.unscaled, "contrasts")
    fit$xlevels <- xlev
    if (missing(data)) 
        data <- environment(formula)
    fit$data <- data
    fit$boundary <- FALSE
    fit$residuals <- (y - pr*denom)/(denom*pr*(1-pr))
    if (keep.xmat) fit$x <- x.unscaled
    fit
}

vcov.brlr <- function (object, ...) 
{
    structure(
        object$dispersion * chol2inv(chol(object$FisherInfo)),
        dimnames = dimnames(object$FisherInfo))
}

print.brlr <-
    function (x, digits = max(3, getOption("digits") - 3),
              na.print = "", 
              ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:  ")
        dput(cl)
    }
    if (length(coef(x))) {
        cat("\nCoefficients:\n")
        print(coef(x), digits = digits, ...)
    }
    else {
        cat("\nNo coefficients\n")
    }
    cat("\nDeviance:",
        format(round(x$deviance, digits), nsmall = 2), 
        "\n")
    if (x$bias.reduction) {
        cat("Penalized deviance:",
            format(round(x$penalized.deviance, digits), nsmall = 2),
            "\n")}
    cat("Residual df:", x$df.residual, "\n")
    invisible(x)
}

summary.brlr <-
function (object, dispersion = NULL,
          digits = max(3, .Options$digits - 3), ...) 
{
    cc <- coef(object)
    object$pc <- pc <- length(coef(object))
    coef <- matrix(0, pc, 3, dimnames = list(names(cc),
                     c("Value", "Std. Error", "t value")))
    coef[, 1] <- cc
    vc <- vcov.brlr(object)
    if (is.null(dispersion)) 
        dispersion <- object$dispersion
    else vc <- vc * dispersion/(object$dispersion)
    sd <- cc
    sd[!is.na(cc)] <- sqrt(diag(vc))
    coef[, 2] <- sd
    coef[, 3] <- coef[, 1]/coef[, 2]
    object$coefficients <- coef
    object$digits <- digits
    class(object) <- "summary.brlr"
    object
}

print.summary.brlr <-
    function (x, digits = x$digits, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    coef <- format(round(x$coefficients, digits = digits))
    pc <- x$pc
    if (pc > 0) {
        cat("\nCoefficients:\n")
        print(coef[seq(len = pc), ], quote = FALSE, ...)
    }
    else {
        cat("\nNo coefficients\n")
    }
    cat("\nDeviance:",
        format(round(x$deviance, digits), nsmall = 2), 
        "\n")
    if (x$bias.reduction) {
        cat("Penalized deviance:",
            format(round(x$penalized.deviance, digits), nsmall = 2),
            "\n")}
    cat("Residual df:", x$df.residual, "\n")
    invisible(x)
}

predict.brlr <-
    function (object, newdata = NULL, type = c("link", "response"), 
              dispersion = NULL, terms = NULL, ...) 
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (missing(newdata)) {
        pred <- switch(type, link = object$linear.predictors, 
                       response = object$fitted)
        if (!is.null(na.act)) 
            pred <- napredict(na.act, pred)
    }
    else {
        newdata <- as.data.frame(newdata)
        Terms <- delete.response(object$terms)
        m <- model.frame(Terms, newdata, na.action = function(x) x, 
                         xlev = object$xlevels)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        offset <- if (!is.null(off.num <- attr(Terms, "offset"))) 
            eval(attr(Terms, "variables")[[off.num + 1]], newdata)
        else if (!is.null(object$offset)) 
            eval(object$call$offset, newdata)
        pred <- offset + drop(X %*% object$coef)
        switch(type, response = {
            pred <- plogis(pred)
        }, link = )
    }
    pred
}
