# Asymmetric Champernowne / Losev distribution
#
# Density:
#   f(x) = d / {exp(a * (x - location)) +
#               exp(-b * (x - location)) + c}
#
# with a > 0, b > 0 and
#   c > -1 / {q^q p^p},
# where p = a/(a+b), q = b/(a+b).
#
# Defaults b = a and c = 0 give the symmetric hyperbolic-secant member.
# c = 0 with a != b gives the Losev / asymmetric hyperbolic-secant member.
# a = b and c = 2 gives the ordinary logistic distribution.

.asymcham_softplus <- function(x) {
    pmax(x, 0) + log1p(exp(-abs(x)))
}

.asymcham_check <- function(a, b, c, location) {
    if (length(a) != 1L || length(b) != 1L ||
        length(c) != 1L || length(location) != 1L) {
        stop("a, b, c, and location must be scalar")
    }
    if (!is.finite(a) || a <= 0) stop("a must be finite and > 0")
    if (!is.finite(b) || b <= 0) stop("b must be finite and > 0")
    if (!is.finite(c)) stop("c must be finite")
    if (!is.finite(location)) stop("location must be finite")

    s <- a + b
    p <- a / s
    q <- b / s

    # w(u) = u^q (1-u)^p has its maximum at u = q.
    log_wmax <- q * log(q) + p * log(p)
    wmax <- exp(log_wmax)
    cmin <- -1 / wmax

    if (c <= cmin) {
        stop(sprintf("c must be greater than %.17g for these a and b", cmin))
    }

    list(a = a, b = b, c = c, location = location,
         s = s, p = p, q = q, wmax = wmax, cmin = cmin)
}

.asymcham_logw <- function(v, par) {
    -par$q * .asymcham_softplus(-v) -
        par$p * .asymcham_softplus(v)
}

.asymcham_g <- function(v, par, r = 0) {
    logw <- .asymcham_logw(v, par)
    w <- exp(logw)
    logg <- r * v + logw - log1p(par$c * w)
    exp(logg)
}

.asymcham_z0 <- function(par) {
    if (par$c == 0) {
        return(exp(lbeta(par$q, par$p)))
    }

    integrate(function(v) .asymcham_g(v, par),
              lower = -Inf, upper = Inf,
              subdivisions = 1000L,
              rel.tol = 1e-10, abs.tol = 1e-12,
              stop.on.error = TRUE)$value
}

# Normalizing constant d.
asymcham_const <- function(a, b = a, c = 0) {
    par <- .asymcham_check(a, b, c, 0)
    par$s / .asymcham_z0(par)
}

# Density.
dasymcham <- function(x, a, b = a, c = 0, location = 0,
                      log = FALSE) {
    par <- .asymcham_check(a, b, c, location)
    z0 <- .asymcham_z0(par)
    v <- par$s * (x - par$location)
    logw <- .asymcham_logw(v, par)
    w <- exp(logw)

    ans <- log(par$s) - log(z0) + logw - log1p(par$c * w)
    if (log) ans else exp(ans)
}

.asymcham_tail_one <- function(v, par, z0, lower.tail) {
    if (is.na(v)) return(NA_real_)

    if (lower.tail) {
        if (v == -Inf) return(0)
        if (v == Inf) return(1)
        val <- integrate(function(z) .asymcham_g(z, par),
                         lower = -Inf, upper = v,
                         subdivisions = 1000L,
                         rel.tol = 1e-10, abs.tol = 1e-12,
                         stop.on.error = TRUE)$value / z0
    } else {
        if (v == -Inf) return(1)
        if (v == Inf) return(0)
        val <- integrate(function(z) .asymcham_g(z, par),
                         lower = v, upper = Inf,
                         subdivisions = 1000L,
                         rel.tol = 1e-10, abs.tol = 1e-12,
                         stop.on.error = TRUE)$value / z0
    }

    min(1, max(0, val))
}

# Distribution function.
pasymcham <- function(q, a, b = a, c = 0, location = 0,
                      lower.tail = TRUE, log.p = FALSE) {
    par <- .asymcham_check(a, b, c, location)
    v <- par$s * (q - par$location)

    # Exact beta/logit representation for c = 0.
    if (par$c == 0) {
        u <- plogis(v)
        return(pbeta(u, shape1 = par$q, shape2 = par$p,
                     lower.tail = lower.tail, log.p = log.p))
    }

    z0 <- .asymcham_z0(par)
    ans <- vapply(v, .asymcham_tail_one, numeric(1),
                  par = par, z0 = z0, lower.tail = lower.tail)

    if (log.p) log(ans) else ans
}

.asymcham_to_lower_prob <- function(p, lower.tail, log.p) {
    bad <- rep(FALSE, length(p))

    if (log.p) {
        bad <- !is.na(p) & p > 0
        ans <- if (lower.tail) exp(p) else -expm1(p)
    } else {
        bad <- !is.na(p) & (p < 0 | p > 1)
        ans <- if (lower.tail) p else 1 - p
    }

    if (any(bad)) {
        warning("NaNs produced")
        ans[bad] <- NaN
    }
    ans
}

.asymcham_quantile_one <- function(prob, par, z0) {
    if (is.na(prob)) return(NA_real_)
    if (prob == 0) return(-Inf)
    if (prob == 1) return(Inf)

    # Use the c = 0 quantile as a good initial center.
    u0 <- qbeta(prob, shape1 = par$q, shape2 = par$p)
    center <- qlogis(u0)
    if (!is.finite(center)) center <- log(par$q / par$p)

    use_lower <- prob <= 0.5
    target <- if (use_lower) prob else 1 - prob

    fun <- function(v) {
        .asymcham_tail_one(v, par, z0, lower.tail = use_lower) - target
    }

    step <- 2
    lo <- center - step
    hi <- center + step

    for (i in seq_len(60L)) {
        flo <- fun(lo)
        fhi <- fun(hi)
        if (is.finite(flo) && is.finite(fhi) && flo * fhi <= 0) break
        step <- step * 2
        lo <- center - step
        hi <- center + step
        if (i == 60L) stop("failed to bracket quantile")
    }

    uniroot(fun, interval = c(lo, hi), tol = 1e-10)$root
}

# Quantile function.
qasymcham <- function(p, a, b = a, c = 0, location = 0,
                      lower.tail = TRUE, log.p = FALSE) {
    par <- .asymcham_check(a, b, c, location)
    prob <- .asymcham_to_lower_prob(p, lower.tail, log.p)

    if (par$c == 0) {
        u <- qbeta(prob, shape1 = par$q, shape2 = par$p)
        return(par$location + qlogis(u) / par$s)
    }

    z0 <- .asymcham_z0(par)
    v <- vapply(prob, .asymcham_quantile_one, numeric(1),
                par = par, z0 = z0)
    par$location + v / par$s
}

# Random generation.
#
# Uses the c = 0 beta/logit distribution as a rejection proposal.  If
# U ~ Beta(q,p), the required correction factor is
#   h(U) = 1 / {1 + c U^q (1-U)^p}.
rasymcham <- function(n, a, b = a, c = 0, location = 0) {
    par <- .asymcham_check(a, b, c, location)

    if (length(n) > 1L) {
        n <- length(n)
    } else {
        if (!is.finite(n) || n < 0) stop("invalid n")
        n <- floor(n)
    }
    if (n == 0) return(numeric(0))

    if (par$c == 0) {
        u <- rbeta(n, shape1 = par$q, shape2 = par$p)
        return(par$location + qlogis(u) / par$s)
    }

    z0 <- .asymcham_z0(par)
    beta0 <- exp(lbeta(par$q, par$p))

    # Envelope constant for the correction h(u).
    if (par$c >= 0) {
        M <- 1
    } else {
        M <- 1 / (1 + par$c * par$wmax)
    }

    acc.rate <- min(1, z0 / (beta0 * M))
    out <- numeric(n)
    filled <- 0L

    while (filled < n) {
        need <- n - filled
        m <- ceiling(1.1 * need / max(acc.rate, 1e-6))
        m <- max(100L, min(m, 1000000L))

        u <- rbeta(m, shape1 = par$q, shape2 = par$p)
        logw <- par$q * log(u) + par$p * log1p(-u)
        w <- exp(logw)
        h <- 1 / (1 + par$c * w)
        accept <- runif(m) < h / M
        u <- u[accept]

        take <- min(length(u), need)
        if (take > 0) {
            idx <- filled + seq_len(take)
            out[idx] <- par$location + qlogis(u[seq_len(take)]) / par$s
            filled <- filled + take
        }
    }

    out
}

# Moment-generating function E[exp(t X)].
# It is finite exactly for -b < t < a.
masymcham <- function(t, a, b = a, c = 0, location = 0,
                      log = FALSE) {
    par <- .asymcham_check(a, b, c, location)
    z0 <- .asymcham_z0(par)

    one <- function(tt) {
        if (is.na(tt)) return(NA_real_)
        if (tt <= -par$b || tt >= par$a) return(Inf)
        if (tt == 0) return(if (log) 0 else 1)

        r <- tt / par$s

        if (par$c == 0) {
            logm <- tt * par$location +
                lgamma(par$q + r) + lgamma(par$p - r) -
                lgamma(par$q) - lgamma(par$p)
        } else {
            zr <- integrate(function(v) .asymcham_g(v, par, r = r),
                            lower = -Inf, upper = Inf,
                            subdivisions = 1000L,
                            rel.tol = 1e-10, abs.tol = 1e-12,
                            stop.on.error = TRUE)$value
            logm <- tt * par$location + log(zr) - log(z0)
        }

        if (log) logm else exp(logm)
    }

    vapply(t, one, numeric(1))
}

# Mode of the distribution.
modeasymcham <- function(a, b = a, c = 0, location = 0) {
    par <- .asymcham_check(a, b, c, location)
    par$location + log(b / a) / par$s
}

# Maximum-likelihood fitting.
#
# By default a, b, c, and location are estimated.  Set fixed.c to a numeric
# value to hold c fixed, and symmetric = TRUE to impose b = a.
fitasymcham <- function(x, fixed.c = NULL, symmetric = FALSE,
                        start = NULL, hessian = TRUE,
                        control = list(maxit = 1000, factr = 1e7,
                                       pgtol = 1e-8)) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 5L) stop("need at least five finite observations")
    sx <- sd(x)
    if (!is.finite(sx) || sx <= 0) stop("x must have positive sample standard deviation")

    if (!is.null(fixed.c)) {
        if (length(fixed.c) != 1L || !is.finite(fixed.c)) {
            stop("fixed.c must be NULL or one finite numeric value")
        }
        if (fixed.c <= -2) stop("fixed.c must be greater than -2")
    }

    loc0 <- median(x)
    s0 <- pi / sx
    q0 <- 0.5

    z0 <- (x - mean(x)) / sx
    sk0 <- mean(z0^3)
    if (is.finite(sk0)) q0 <- plogis(0.35 * sk0)
    q0 <- min(0.85, max(0.15, q0))

    cmin_from_q <- function(q) {
        p <- 1 - q
        -exp(-q * log(q) - p * log(p))
    }

    if (is.null(start)) {
        ss0 <- s0
        qq0 <- if (symmetric) 0.5 else q0
        l0 <- loc0
        c0 <- if (is.null(fixed.c)) 0 else fixed.c
    } else {
        if (!is.list(start)) stop("start must be NULL or a named list")
        a0 <- if (!is.null(start$a)) as.numeric(start$a) else s0 * (1 - q0)
        b0 <- if (!is.null(start$b)) as.numeric(start$b) else if (symmetric) a0 else s0 * q0
        l0 <- if (!is.null(start$location)) as.numeric(start$location) else loc0
        c0 <- if (!is.null(start$c)) as.numeric(start$c) else if (is.null(fixed.c)) 0 else fixed.c
        if (length(a0) != 1L || length(b0) != 1L ||
            !is.finite(a0) || a0 <= 0 || !is.finite(b0) || b0 <= 0) {
            stop("start$a and start$b must be finite scalars > 0")
        }
        if (length(l0) != 1L || !is.finite(l0)) stop("start$location must be finite")
        if (symmetric) b0 <- a0
        ss0 <- a0 + b0
        qq0 <- if (symmetric) 0.5 else b0 / ss0
    }

    # Keep the transformed skew parameter away from machine 0 and 1.  These
    # are numerical, not substantive, bounds: q in [~6e-6, 0.999994].
    qlo <- plogis(-12)
    qhi <- plogis(12)
    qq0 <- min(qhi, max(qlo, qq0))

    theta0 <- c(log_s = log(ss0))
    if (!symmetric) theta0 <- c(theta0, logit_q = qlogis(qq0))
    theta0 <- c(theta0, location = l0)

    if (is.null(fixed.c)) {
        cmin0 <- cmin_from_q(qq0)
        if (length(c0) != 1L || !is.finite(c0) || c0 <= cmin0) {
            stop("start$c is outside the admissible range")
        }
        theta0 <- c(theta0, log_c_gap = log(c0 - cmin0))
    } else if (fixed.c < 0) {
        # A negative fixed c is not admissible for every q.  Check the start;
        # later invalid optimizer proposals receive a large objective value.
        if (fixed.c <= cmin_from_q(qq0)) {
            stop("fixed.c is outside the admissible range at the starting asymmetry")
        }
    }

    # optim() is allowed to pass an unnamed parameter vector to fn().  The
    # original implementation indexed theta by name inside fn(), which made
    # all fits fail on R builds where those names were dropped.  Decode by
    # position instead.
    unpack <- function(theta) {
        k <- 1L
        log_s <- theta[[k]]; k <- k + 1L
        logit_q <- if (symmetric) 0 else { z <- theta[[k]]; k <- k + 1L; z }
        location <- theta[[k]]; k <- k + 1L
        log_c_gap <- if (is.null(fixed.c)) theta[[k]] else NA_real_

        ss <- exp(log_s)
        qq <- if (symmetric) 0.5 else plogis(logit_q)
        aa <- ss * (1 - qq)
        bb <- ss * qq
        cc <- if (is.null(fixed.c)) {
            cmin_from_q(qq) + exp(log_c_gap)
        } else {
            fixed.c
        }
        list(a = aa, b = bb, c = cc, location = location, s = ss, q = qq)
    }

    nll <- function(theta) {
        par <- unpack(theta)
        if (!all(is.finite(c(par$a, par$b, par$c, par$location))) ||
            par$a <= 0 || par$b <= 0) return(1e100)

        # For a fixed negative c, reject asymmetries for which the denominator
        # can vanish.  Free c is parameterized to be valid automatically.
        if (!is.null(fixed.c)) {
            cm <- cmin_from_q(par$q)
            if (!is.finite(cm) || fixed.c <= cm) return(1e100)
        }

        val <- tryCatch(
            suppressWarnings(dasymcham(x, a = par$a, b = par$b, c = par$c,
                                       location = par$location, log = TRUE)),
            error = function(e) rep(-Inf, length(x))
        )
        if (length(val) != length(x) || any(!is.finite(val))) return(1e100)
        ans <- -sum(val)
        if (!is.finite(ans)) 1e100 else ans
    }

    # Broad numerical bounds prevent exp()/plogis() overflow while leaving a
    # vastly wider parameter range than is relevant for standardized returns.
    lower <- c(-12)
    upper <- c( 12)
    if (!symmetric) {
        lower <- c(lower, -12)
        upper <- c(upper,  12)
    }
    lower <- c(lower, -Inf)
    upper <- c(upper,  Inf)
    if (is.null(fixed.c)) {
        lower <- c(lower, -20)
        upper <- c(upper,  20)
    }

    # Make control compatible with L-BFGS-B even if a caller supplied reltol,
    # which is a BFGS/Nelder-Mead control parameter.
    ctl <- control
    ctl$reltol <- NULL
    if (is.null(ctl$maxit)) ctl$maxit <- 1000
    if (is.null(ctl$factr)) ctl$factr <- 1e7
    if (is.null(ctl$pgtol)) ctl$pgtol <- 1e-8

    opt <- optim(theta0, nll, method = "L-BFGS-B",
                 lower = lower, upper = upper,
                 hessian = hessian, control = ctl)
    par <- unpack(opt$par)
    k <- length(opt$par)
    ll <- -opt$value

    natural_coef <- c(a = par$a, b = par$b, c = par$c,
                      location = par$location)

    vc <- NULL
    if (hessian && !is.null(opt$hessian) && all(is.finite(opt$hessian))) {
        vc <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    }

    ans <- list(
        coefficients = natural_coef,
        transformed.coefficients = opt$par,
        logLik = ll,
        nobs = length(x),
        df = k,
        aic = -2 * ll + 2 * k,
        bic = -2 * ll + log(length(x)) * k,
        convergence = opt$convergence,
        message = opt$message,
        value = opt$value,
        hessian = if (hessian) opt$hessian else NULL,
        vcov.transformed = vc,
        fixed.c = fixed.c,
        symmetric = symmetric,
        data = x,
        call = match.call()
    )
    class(ans) <- "asymcham_fit"
    ans
}

coef.asymcham_fit <- function(object, ...) object$coefficients

logLik.asymcham_fit <- function(object, ...) {
    ans <- object$logLik
    attr(ans, "df") <- object$df
    attr(ans, "nobs") <- object$nobs
    class(ans) <- "logLik"
    ans
}

nobs.asymcham_fit <- function(object, ...) object$nobs

vcov.asymcham_fit <- function(object, ...) object$vcov.transformed

print.asymcham_fit <- function(x, ...) {
    cat("Asymmetric Champernowne maximum-likelihood fit\n")
    cat("nobs:", x$nobs, "  logLik:", format(x$logLik, digits = 8),
        "  AIC:", format(x$aic, digits = 8),
        "  BIC:", format(x$bic, digits = 8), "\n")
    cat("convergence:", x$convergence, "\n")
    print(x$coefficients)
    invisible(x)
}

summary.asymcham_fit <- function(object, ...) {
    out <- list(
        call = object$call,
        coefficients = object$coefficients,
        logLik = object$logLik,
        nobs = object$nobs,
        df = object$df,
        AIC = object$aic,
        BIC = object$bic,
        convergence = object$convergence,
        message = object$message,
        fixed.c = object$fixed.c,
        symmetric = object$symmetric
    )
    class(out) <- "summary.asymcham_fit"
    out
}

print.summary.asymcham_fit <- function(x, ...) {
    cat("Asymmetric Champernowne fit summary\n\n")
    if (!is.null(x$call)) print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nlogLik:", format(x$logLik, digits = 8),
        "  AIC:", format(x$AIC, digits = 8),
        "  BIC:", format(x$BIC, digits = 8), "\n")
    cat("convergence:", x$convergence, "\n")
    invisible(x)
}
