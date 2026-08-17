# Compare distributions fitted to EWMA-volatility-normalized asset log returns.
#
# Required packages:
#   gamlss.dist
#   GeneralizedHyperbolic
#   ghyp
#
# The script expects asymcham.r and r_utils.r in the same directory.  By
# default it reads asset_class_etf_prices.csv from that directory.
#
# Models:
#   normal
#   Student t
#   Jones-Faddy skew t (gamlss.dist::ST5)
#   Laplace
#   asymmetric Laplace
#   power exponential / GED (gamlss.dist::PE)
#   hyperbolic (GeneralizedHyperbolic::hyperbFit)
#   symmetric hyperbolic (beta fixed at 0)
#   normal inverse Gaussian (ghyp::fit.NIGuv), general and symmetric
#   variance gamma (ghyp::fit.VGuv), general and symmetric
#   generalized hyperbolic, free lambda (ghyp::fit.ghypuv), general and symmetric
#   AsymCham with c = 0
#   AsymCham with c = 2
#   AsymCham with c estimated
#
# Main outputs are CSV files in output_dir.  Set RUN_OOS below to TRUE for a
# chronological train/test log-score comparison in addition to the full-sample
# fits.

# ------------------------------- settings ---------------------------------

script_start_time <- Sys.time()

options(width=10000)
EWMA_LAMBDA <- 0.94
EWMA_WARMUP <- 100L
EWMA_CENTER <- FALSE
TAIL_THRESHOLDS <- c(2, 3, 4)
RUN_OOS <- FALSE
OOS_TRAIN_FRACTION <- 0.70
MIN_FIT_N <- 200L

# Empty means all price columns.  Example: c("SPY", "EFA", "EEM").
SYMBOLS <- character(0)

# Which models to fit, in the order they will be fit and reported.  Remove
# entries from this vector to skip those models.
MODELS <- c("normal", "student_t", "jones_faddy_st5", "laplace",
            "asymmetric_laplace", "power_exponential",
            "hyperbolic", "hyperbolic_sym",
            "nig", "nig_sym", "vg", "vg_sym", "gh", "gh_sym",
            "asymcham_c0", "asymcham_c2") # , "asymcham_free")

# --------------------------- paths and packages ----------------------------

script_path <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_path)) {
    dirname(normalizePath(sub("^--file=", "", script_path[[1L]])))
} else {
    getwd()
}

args <- commandArgs(trailingOnly = TRUE)
price_file <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "asset_class_etf_prices.csv")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(script_dir, "results")

if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("Package 'gamlss.dist' is required. Install it with install.packages('gamlss.dist').")
}
if (!requireNamespace("GeneralizedHyperbolic", quietly = TRUE)) {
    stop("Package 'GeneralizedHyperbolic' is required. Install it with install.packages('GeneralizedHyperbolic').")
}
if (!requireNamespace("ghyp", quietly = TRUE)) {
    stop("Package 'ghyp' is required. Install it with install.packages('ghyp').")
}

source(file.path(script_dir, "r_utils.r"))
source(file.path(script_dir, "asymcham.r"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (exists("print_script_header", mode = "function")) print_script_header()
cat("price file:", price_file, "\n")
cat("output dir:", output_dir, "\n")
cat("EWMA lambda:", EWMA_LAMBDA, " warmup:", EWMA_WARMUP,
    " center:", EWMA_CENTER, "\n\n")

# ------------------------------ fit helpers --------------------------------

.safe_logsum <- function(logd) {
    if (!length(logd) || any(!is.finite(logd))) return(-Inf)
    sum(logd)
}

.safe_optim <- function(starts, objective, method = "BFGS",
                        control = list(maxit = 1000, reltol = 1e-9)) {
    best <- NULL
    for (st in starts) {
        z <- tryCatch(
            optim(st, objective, method = method, control = control),
            error = function(e) NULL,
            warning = function(w) suppressWarnings(
                tryCatch(optim(st, objective, method = method, control = control),
                         error = function(e) NULL))
        )
        if (is.null(z) || !is.finite(z$value)) next
        if (is.null(best) || z$value < best$value) best <- z
    }
    best
}

.make_fit <- function(model, x, coef, loglik, convergence, log_density, cdf,
                      message = NULL) {
    k <- length(coef)
    n <- length(x)
    list(
        model = model,
        coef = coef,
        logLik = loglik,
        nobs = n,
        k = k,
        AIC = -2 * loglik + 2 * k,
        BIC = -2 * loglik + log(n) * k,
        convergence = convergence,
        message = message,
        log_density = log_density,
        cdf = cdf
    )
}

.failed_fit <- function(model, x, message) {
    list(model = model, coef = numeric(0), logLik = NA_real_, nobs = length(x),
         k = NA_integer_, AIC = NA_real_, BIC = NA_real_, convergence = 99L,
         message = message,
         log_density = function(y) rep(NA_real_, length(y)),
         cdf = function(y) rep(NA_real_, length(y)))
}

# ---------------------------- model fitters --------------------------------

.fit_normal <- function(x) {
    mu <- mean(x)
    sigma <- sqrt(mean((x - mu)^2))
    logd <- dnorm(x, mean = mu, sd = sigma, log = TRUE)
    .make_fit(
        "normal", x,
        c(mu = mu, sigma = sigma), .safe_logsum(logd), 0L,
        function(y) dnorm(y, mean = mu, sd = sigma, log = TRUE),
        function(y) pnorm(y, mean = mu, sd = sigma)
    )
}

.fit_student_t <- function(x) {
    mu0 <- median(x)
    sg0 <- sd(x)
    starts <- lapply(c(4, 7, 12, 30), function(df0) c(mu0, log(sg0), log(df0)))
    objective <- function(th) {
        mu <- th[1L]
        sigma <- exp(th[2L])
        df <- exp(th[3L])
        ll <- .safe_logsum(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma))
        if (is.finite(ll)) -ll else 1e100
    }
    opt <- .safe_optim(starts, objective)
    if (is.null(opt)) return(.failed_fit("student_t", x, "optimization failed"))
    mu <- opt$par[1L]
    sigma <- exp(opt$par[2L])
    df <- exp(opt$par[3L])
    .make_fit(
        "student_t", x, c(mu = mu, sigma = sigma, df = df), -opt$value,
        opt$convergence,
        function(y) dt((y - mu) / sigma, df = df, log = TRUE) - log(sigma),
        function(y) pt((y - mu) / sigma, df = df),
        opt$message
    )
}

.dlaplace <- function(x, mu, scale, log = FALSE) {
    ans <- -log(2 * scale) - abs(x - mu) / scale
    if (log) ans else exp(ans)
}

.plaplace <- function(q, mu, scale) {
    z <- (q - mu) / scale
    ifelse(z < 0, 0.5 * exp(z), 1 - 0.5 * exp(-z))
}

.fit_laplace <- function(x) {
    mu <- median(x)
    scale <- mean(abs(x - mu))
    ll <- .safe_logsum(.dlaplace(x, mu, scale, log = TRUE))
    .make_fit(
        "laplace", x, c(mu = mu, scale = scale), ll, 0L,
        function(y) .dlaplace(y, mu, scale, log = TRUE),
        function(y) .plaplace(y, mu, scale)
    )
}

# Asymmetric Laplace parameterized by separate exponential rates:
#   right tail rate = a, left tail rate = b.
.dasym_laplace <- function(x, mu, a, b, log = FALSE) {
    logk <- log(a) + log(b) - log(a + b)
    ans <- ifelse(x < mu, logk + b * (x - mu), logk - a * (x - mu))
    if (log) ans else exp(ans)
}

.pasym_laplace <- function(q, mu, a, b) {
    ifelse(q < mu,
           (a / (a + b)) * exp(b * (q - mu)),
           1 - (b / (a + b)) * exp(-a * (q - mu)))
}

.fit_asym_laplace <- function(x) {
    mu0 <- median(x)
    rate0 <- 1 / max(mean(abs(x - mu0)), 1e-8)
    starts <- list(
        c(mu0, log(rate0), log(rate0)),
        c(mu0, log(1.5 * rate0), log(0.75 * rate0)),
        c(mu0, log(0.75 * rate0), log(1.5 * rate0))
    )
    objective <- function(th) {
        mu <- th[1L]
        a <- exp(th[2L])
        b <- exp(th[3L])
        ll <- .safe_logsum(.dasym_laplace(x, mu, a, b, log = TRUE))
        if (is.finite(ll)) -ll else 1e100
    }
    opt <- .safe_optim(starts, objective)
    if (is.null(opt)) return(.failed_fit("asymmetric_laplace", x, "optimization failed"))
    mu <- opt$par[1L]
    a <- exp(opt$par[2L])
    b <- exp(opt$par[3L])
    .make_fit(
        "asymmetric_laplace", x,
        c(mu = mu, right_rate_a = a, left_rate_b = b), -opt$value,
        opt$convergence,
        function(y) .dasym_laplace(y, mu, a, b, log = TRUE),
        function(y) .pasym_laplace(y, mu, a, b),
        opt$message
    )
}

.fit_power_exponential <- function(x) {
    mu0 <- mean(x)
    sg0 <- sd(x)
    starts <- lapply(c(0.7, 1.0, 1.5, 2.0, 3.0),
                     function(nu0) c(mu0, log(sg0), log(nu0)))
    objective <- function(th) {
        mu <- th[1L]
        sigma <- exp(th[2L])
        nu <- exp(th[3L])
        ld <- suppressWarnings(gamlss.dist::dPE(x, mu = mu, sigma = sigma,
                                                nu = nu, log = TRUE))
        ll <- .safe_logsum(ld)
        if (is.finite(ll)) -ll else 1e100
    }
    opt <- .safe_optim(starts, objective)
    if (is.null(opt)) return(.failed_fit("power_exponential", x, "optimization failed"))
    mu <- opt$par[1L]
    sigma <- exp(opt$par[2L])
    nu <- exp(opt$par[3L])
    .make_fit(
        "power_exponential", x, c(mu = mu, sigma = sigma, nu = nu), -opt$value,
        opt$convergence,
        function(y) gamlss.dist::dPE(y, mu = mu, sigma = sigma, nu = nu, log = TRUE),
        function(y) gamlss.dist::pPE(y, mu = mu, sigma = sigma, nu = nu),
        opt$message
    )
}

# Hyperbolic distribution (Barndorff-Nielsen).  Tails decay like
# exp(-(alpha - beta) x) as x -> +Inf and exp(-(alpha + beta) x) as
# x -> -Inf, so, like AsymCham, it allows different exponential decay
# rates in each tail, but through a different (Bessel-function-based)
# functional form.  Several starting-value heuristics are tried and the
# one with the highest likelihood is kept; "SL" (skew-Laplace starting
# values) is dropped because it is prone to raising an error before the
# optimizer even starts.
.fit_hyperbolic <- function(x) {
    start_values <- c("BN", "MoM", "FN")
    best <- NULL
    errors <- character(0)

    for (sv in start_values) {
        z <- tryCatch(
            suppressWarnings(GeneralizedHyperbolic::hyperbFit(
                x, startValues = sv, plots = FALSE, printOut = FALSE)),
            error = function(e) {
                errors <<- c(errors, conditionMessage(e))
                NULL
            }
        )
        if (is.null(z) || !is.finite(z$maxLik)) next
        if (is.null(best) || z$maxLik > best$maxLik) best <- z
    }

    if (is.null(best)) {
        msg <- if (length(errors)) paste(unique(errors), collapse = " | ") else "optimization failed"
        return(.failed_fit("hyperbolic", x, msg))
    }

    cf <- coef(best)
    mu <- unname(cf[["mu"]])
    delta <- unname(cf[["delta"]])
    alpha <- unname(cf[["alpha"]])
    beta <- unname(cf[["beta"]])
    .make_fit(
        "hyperbolic", x, c(mu = mu, delta = delta, alpha = alpha, beta = beta),
        best$maxLik, best$conv,
        function(y) log(GeneralizedHyperbolic::dhyperb(y, mu = mu, delta = delta,
                                                        alpha = alpha, beta = beta)),
        function(y) GeneralizedHyperbolic::phyperb(y, mu = mu, delta = delta,
                                                    alpha = alpha, beta = beta)
    )
}

# Symmetric hyperbolic distribution: beta fixed at 0, so only (mu, delta,
# alpha) are estimated.  hyperbFit() has no facility for holding a
# parameter fixed, so this is a direct 3-parameter MLE via optim(), mirroring
# .fit_student_t / .fit_power_exponential above.  Having one fewer parameter
# than "hyperbolic" costs nothing on AIC's 2-point-per-parameter penalty but
# can win on BIC's log(n)-point penalty when the data show little skew.
.fit_hyperbolic_sym <- function(x) {
    mu0 <- median(x)
    delta0 <- mad(x, constant = 1)
    if (!is.finite(delta0) || delta0 <= 0) delta0 <- sd(x)
    starts <- lapply(c(0.5, 1, 2, 4),
                     function(alpha0) c(mu0, log(delta0), log(alpha0)))
    objective <- function(th) {
        mu <- th[1L]
        delta <- exp(th[2L])
        alpha <- exp(th[3L])
        ld <- suppressWarnings(log(GeneralizedHyperbolic::dhyperb(
            x, mu = mu, delta = delta, alpha = alpha, beta = 0)))
        ll <- .safe_logsum(ld)
        if (is.finite(ll)) -ll else 1e100
    }
    opt <- .safe_optim(starts, objective)
    if (is.null(opt)) return(.failed_fit("hyperbolic_sym", x, "optimization failed"))
    mu <- opt$par[1L]
    delta <- exp(opt$par[2L])
    alpha <- exp(opt$par[3L])
    .make_fit(
        "hyperbolic_sym", x, c(mu = mu, delta = delta, alpha = alpha), -opt$value,
        opt$convergence,
        function(y) log(GeneralizedHyperbolic::dhyperb(y, mu = mu, delta = delta,
                                                        alpha = alpha, beta = 0)),
        function(y) GeneralizedHyperbolic::phyperb(y, mu = mu, delta = delta,
                                                    alpha = alpha, beta = 0),
        opt$message
    )
}

# NIG, variance gamma, and the full generalized hyperbolic (free lambda),
# each in general and symmetric (gamma fixed at 0) form.  All three are
# special/limiting cases of the generalized hyperbolic family:
#   NIG: lambda = -1/2, semi-heavy tails via the inverse-Gaussian mixture.
#   VG:  delta = 0, lambda > 0, tails set by a Gamma time-change instead.
#   GH:  lambda free, nests hyperbolic (lambda=1), NIG (lambda=-1/2) and,
#        as delta -> 0, VG.
# These are fitted with the ghyp package's own EM-based routine rather than
# a hand-rolled optim() call: generic optimizers are notoriously prone to
# flat ridges and boundary problems for this family (this is also why the
# NIG fit in asymcham.r needed explicit bounds), and ghyp's fitters are
# built specifically to be robust to that.
.fit_ghyp_family <- function(x, model, fitfun, symmetric = FALSE, ...) {
    z <- tryCatch(
        suppressWarnings(fitfun(x, symmetric = symmetric, silent = TRUE, ...)),
        error = function(e) NULL
    )
    if (is.null(z) || !isTRUE(z@converged) || !is.finite(z@llh)) {
        msg <- if (is.null(z)) "optimization failed" else "did not converge"
        return(.failed_fit(model, x, msg))
    }

    all_pars <- c(lambda = z@lambda, alpha.bar = z@alpha.bar,
                  mu = z@mu, sigma = z@sigma, gamma = z@gamma)
    keep <- names(z@fitted.params)[z@fitted.params]

    .make_fit(
        model, x, all_pars[keep], z@llh, 0L,
        function(y) ghyp::dghyp(y, object = z, logvalue = TRUE),
        function(y) ghyp::pghyp(y, object = z)
    )
}

.fit_nig     <- function(x) .fit_ghyp_family(x, "nig",     ghyp::fit.NIGuv,  symmetric = FALSE)
.fit_nig_sym <- function(x) .fit_ghyp_family(x, "nig_sym", ghyp::fit.NIGuv,  symmetric = TRUE)
.fit_vg      <- function(x) .fit_ghyp_family(x, "vg",      ghyp::fit.VGuv,   symmetric = FALSE)
.fit_vg_sym  <- function(x) .fit_ghyp_family(x, "vg_sym",  ghyp::fit.VGuv,   symmetric = TRUE)

# The free-lambda fit has one more parameter to place than NIG/VG/hyperbolic
# and needs more than optim()'s default 500 Nelder-Mead iterations to reach
# its convergence tolerance on data like this.
.fit_gh      <- function(x) .fit_ghyp_family(x, "gh",      ghyp::fit.ghypuv, symmetric = FALSE, control = list(maxit = 2000))
.fit_gh_sym  <- function(x) .fit_ghyp_family(x, "gh_sym",  ghyp::fit.ghypuv, symmetric = TRUE,  control = list(maxit = 2000))

.djones_faddy <- function(x, mu, sigma, alpha, beta, log = FALSE) {
    if (!is.finite(sigma) || sigma <= 0 ||
        !is.finite(alpha) || alpha <= 0 ||
        !is.finite(beta) || beta <= 0) {
        ans <- rep(NaN, length(x))
        return(if (log) ans else exp(ans))
    }

    z <- (x - mu) / sigma
    ab <- alpha + beta
    u <- z / sqrt(ab + z^2)

    # Roundoff alone can put u exactly at +/-1 for enormous |z|.  The true
    # value is strictly inside (-1,1).
    eps <- .Machine$double.eps
    u <- pmin(1 - eps, pmax(-1 + eps, u))

    ans <- (alpha + 0.5) * log1p(u) +
           (beta  + 0.5) * log1p(-u) -
           (alpha + beta - 1) * log(2) -
           0.5 * log(alpha + beta) -
           lbeta(alpha, beta) - log(sigma)
    if (log) ans else exp(ans)
}

.pjones_faddy <- function(q, mu, sigma, alpha, beta) {
    z <- (q - mu) / sigma
    u <- 0.5 * (1 + z / sqrt(alpha + beta + z^2))
    u <- pmin(1, pmax(0, u))
    pbeta(u, shape1 = alpha, shape2 = beta)
}

.jf_ab_to_st5 <- function(alpha, beta) {
    tau <- 2 / (alpha + beta)
    nu <- (alpha - beta) / sqrt(alpha * beta * (alpha + beta))
    c(nu = nu, tau = tau)
}

# Some older gamlss.dist releases contained an ST5 nu,tau reparameterization
# bug (a missing square root).  Detect the installed implementation rather
# than relying on a version number.
.gamlss_st5_ok <- local({
    alpha <- 2.3
    beta <- 3.1
    mu <- 0.2
    sigma <- 1.15
    y <- c(-2.0, -0.3, 0.4, 1.7)
    nt <- .jf_ab_to_st5(alpha, beta)
    ref <- .djones_faddy(y, mu, sigma, alpha, beta, log = TRUE)
    got <- tryCatch(
        suppressWarnings(gamlss.dist::dST5(y, mu = mu, sigma = sigma,
                                           nu = nt[["nu"]], tau = nt[["tau"]],
                                           log = TRUE)),
        error = function(e) rep(NA_real_, length(y))
    )
    length(got) == length(ref) && all(is.finite(got)) &&
        max(abs(got - ref)) < 1e-8
})

.dst5_checked <- function(x, mu, sigma, alpha, beta, log = FALSE) {
    if (.gamlss_st5_ok) {
        nt <- .jf_ab_to_st5(alpha, beta)
        gamlss.dist::dST5(x, mu = mu, sigma = sigma,
                          nu = nt[["nu"]], tau = nt[["tau"]], log = log)
    } else {
        .djones_faddy(x, mu, sigma, alpha, beta, log = log)
    }
}

.pst5_checked <- function(q, mu, sigma, alpha, beta) {
    if (.gamlss_st5_ok) {
        nt <- .jf_ab_to_st5(alpha, beta)
        gamlss.dist::pST5(q, mu = mu, sigma = sigma,
                          nu = nt[["nu"]], tau = nt[["tau"]])
    } else {
        .pjones_faddy(q, mu, sigma, alpha, beta)
    }
}

cat("gamlss.dist version:", as.character(utils::packageVersion("gamlss.dist")), "\n")
if (.gamlss_st5_ok) {
    cat("Jones-Faddy ST5 check: package implementation OK\n\n")
} else {
    cat("Jones-Faddy ST5 check: installed package implementation is inconsistent;\n")
    cat("  using the corrected Jones-Faddy/ST5 formula in this script.\n\n")
}

.fit_jones_faddy <- function(x) {
    # Jones-Faddy shape parameters alpha,beta are fitted directly and kept
    # positive.  They are converted to ST5 nu,tau when the installed
    # gamlss.dist implementation passes the self-test above; otherwise the
    # equivalent corrected density is evaluated directly.
    mu0 <- median(x)
    sg0 <- sd(x)
    starts <- list()
    for (df0 in c(3, 5, 8, 12, 30)) {
        a0 <- df0 / 2
        starts[[length(starts) + 1L]] <- c(mu0, log(sg0), log(a0), log(a0))
        starts[[length(starts) + 1L]] <- c(mu0, log(sg0), log(0.8 * a0), log(1.2 * a0))
        starts[[length(starts) + 1L]] <- c(mu0, log(sg0), log(1.2 * a0), log(0.8 * a0))
    }

    objective <- function(th) {
        mu <- th[1L]
        sigma <- exp(th[2L])
        alpha <- exp(th[3L])
        beta <- exp(th[4L])
        ld <- .dst5_checked(x, mu, sigma, alpha, beta, log = TRUE)
        ll <- .safe_logsum(ld)
        if (is.finite(ll)) -ll else 1e100
    }

    best <- NULL
    lower <- c(-Inf, -12, -6, -6)
    upper <- c( Inf,  12,  7,  7)
    for (st in starts) {
        z <- tryCatch(
            optim(st, objective, method = "L-BFGS-B",
                  lower = lower, upper = upper,
                  control = list(maxit = 1200, factr = 1e7, pgtol = 1e-8)),
            error = function(e) NULL
        )
        if (is.null(z) || !is.finite(z$value)) next
        if (is.null(best) || z$value < best$value) best <- z
    }

    if (is.null(best)) {
        return(.failed_fit("jones_faddy_st5", x, "optimization failed"))
    }

    mu <- best$par[1L]
    sigma <- exp(best$par[2L])
    alpha <- exp(best$par[3L])
    beta <- exp(best$par[4L])

    .make_fit(
        "jones_faddy_st5", x,
        c(mu = mu, sigma = sigma, alpha = alpha, beta = beta),
        -best$value, best$convergence,
        function(y) .dst5_checked(y, mu, sigma, alpha, beta, log = TRUE),
        function(y) .pst5_checked(y, mu, sigma, alpha, beta),
        best$message
    )
}

.fit_asymcham_one <- function(x, fixed.c = NULL) {
    model <- if (is.null(fixed.c)) "asymcham_free" else if (fixed.c == 0) "asymcham_c0" else "asymcham_c2"
    sg <- sd(x)
    loc <- median(x)
    s0 <- pi / sg
    q_starts <- c(0.25, 0.40, 0.50, 0.60, 0.75)
    c_starts <- if (is.null(fixed.c)) c(-0.75, 0, 2) else fixed.c
    best <- NULL
    errors <- character(0)

    for (q0 in q_starts) {
        for (c0 in c_starts) {
            st <- list(a = s0 * (1 - q0), b = s0 * q0,
                       c = c0, location = loc)
            z <- tryCatch(
                suppressWarnings(fitasymcham(x, fixed.c = fixed.c,
                                             symmetric = FALSE, start = st,
                                             hessian = FALSE,
                                             control = list(maxit = 1000,
                                                            factr = 1e7,
                                                            pgtol = 1e-8))),
                error = function(e) {
                    errors <<- c(errors, conditionMessage(e))
                    NULL
                }
            )
            if (is.null(z) || !is.finite(z$logLik)) next
            if (is.null(best) || z$logLik > best$logLik) best <- z
        }
    }

    if (is.null(best)) {
        msg <- if (length(errors)) paste(unique(errors), collapse = " | ") else "optimization failed"
        return(.failed_fit(model, x, msg))
    }
    cf <- coef(best)
    a <- unname(cf[["a"]])
    b <- unname(cf[["b"]])
    cc <- unname(cf[["c"]])
    location <- unname(cf[["location"]])
    .make_fit(
        model, x, c(a = a, b = b, c = cc, location = location),
        best$logLik, best$convergence,
        function(y) dasymcham(y, a = a, b = b, c = cc,
                              location = location, log = TRUE),
        function(y) pasymcham(y, a = a, b = b, c = cc,
                              location = location),
        best$message
    )
}

.fit_all_models <- function(x) {
    all_fitters <- list(
        normal = .fit_normal,
        student_t = .fit_student_t,
        jones_faddy_st5 = .fit_jones_faddy,
        laplace = .fit_laplace,
        asymmetric_laplace = .fit_asym_laplace,
        power_exponential = .fit_power_exponential,
        hyperbolic = .fit_hyperbolic,
        hyperbolic_sym = .fit_hyperbolic_sym,
        nig = .fit_nig,
        nig_sym = .fit_nig_sym,
        vg = .fit_vg,
        vg_sym = .fit_vg_sym,
        gh = .fit_gh,
        gh_sym = .fit_gh_sym,
        asymcham_c0 = function(z) .fit_asymcham_one(z, fixed.c = 0),
        asymcham_c2 = function(z) .fit_asymcham_one(z, fixed.c = 2),
        asymcham_free = function(z) .fit_asymcham_one(z, fixed.c = NULL)
    )

    model_names <- if (length(MODELS)) MODELS else names(all_fitters)
    unknown <- setdiff(model_names, names(all_fitters))
    if (length(unknown)) {
        stop("Unknown model name(s) in MODELS: ", paste(unknown, collapse = ", "),
             "\nAvailable models: ", paste(names(all_fitters), collapse = ", "))
    }
    fitters <- all_fitters[model_names]

    out <- vector("list", length(fitters))
    names(out) <- names(fitters)
    for (nm in names(fitters)) {
        cat("    fitting", nm, "...\n")
        t0 <- Sys.time()
        out[[nm]] <- tryCatch(fitters[[nm]](x),
                              error = function(e) .failed_fit(nm, x, conditionMessage(e)))
        out[[nm]]$elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    out
}

# ------------------------------ diagnostics --------------------------------

.fit_summary_rows <- function(asset, fits) {
    rows <- lapply(fits, function(f) {
        data.frame(
            asset = asset,
            model = f$model,
            nobs = f$nobs,
            npar = f$k,
            logLik = f$logLik,
            AIC = f$AIC,
            BIC = f$BIC,
            convergence = f$convergence,
            message = if (is.null(f$message)) "" else as.character(f$message),
            fit_seconds = if (is.null(f$elapsed)) NA_real_ else f$elapsed,
            stringsAsFactors = FALSE
        )
    })
    z <- do.call(rbind, rows)
    ok <- is.finite(z$AIC)
    z$delta_AIC <- NA_real_
    z$delta_BIC <- NA_real_
    z$AIC_rank <- NA_integer_
    z$BIC_rank <- NA_integer_
    if (any(ok)) {
        z$delta_AIC[ok] <- z$AIC[ok] - min(z$AIC[ok])
        z$delta_BIC[ok] <- z$BIC[ok] - min(z$BIC[ok])
        z$AIC_rank[ok] <- rank(z$AIC[ok], ties.method = "min")
        z$BIC_rank[ok] <- rank(z$BIC[ok], ties.method = "min")
    }
    z
}

# Cross-asset goodness-of-fit summary from the full-sample AIC/BIC rankings.
# One row per model: how it does on average and how often it wins, across
# all assets for which it produced a finite AIC/BIC.
.model_goodness_summary <- function(fit_summary) {
    z <- fit_summary[is.finite(fit_summary$AIC), ]
    if (!nrow(z)) return(data.frame())

    models <- unique(z$model)
    rows <- lapply(models, function(m) {
        zm <- z[z$model == m, ]
        data.frame(
            model = m,
            n_assets = nrow(zm),
            mean_AIC_rank = mean(zm$AIC_rank),
            mean_BIC_rank = mean(zm$BIC_rank),
            aic_wins = sum(zm$AIC_rank == 1L),
            bic_wins = sum(zm$BIC_rank == 1L),
            mean_delta_AIC = mean(zm$delta_AIC),
            mean_delta_BIC = mean(zm$delta_BIC),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$mean_AIC_rank), ]
}

.parameter_rows <- function(asset, fits) {
    rows <- list()
    k <- 0L
    for (f in fits) {
        if (!length(f$coef)) next
        for (nm in names(f$coef)) {
            k <- k + 1L
            rows[[k]] <- data.frame(asset = asset, model = f$model,
                                    parameter = nm, estimate = unname(f$coef[[nm]]),
                                    stringsAsFactors = FALSE)
        }
    }
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
}

.tail_rows <- function(asset, x, fits, thresholds) {
    rows <- list()
    krow <- 0L
    for (thr in thresholds) {
        emp_left <- mean(x <= -thr)
        emp_right <- mean(x >= thr)
        for (f in fits) {
            if (!is.finite(f$logLik)) next
            fit_left <- tryCatch(as.numeric(f$cdf(-thr)), error = function(e) NA_real_)
            fit_right <- tryCatch(1 - as.numeric(f$cdf(thr)), error = function(e) NA_real_)
            krow <- krow + 1L
            rows[[krow]] <- data.frame(
                asset = asset,
                model = f$model,
                threshold = thr,
                empirical_left = emp_left,
                fitted_left = fit_left,
                left_error = fit_left - emp_left,
                empirical_right = emp_right,
                fitted_right = fit_right,
                right_error = fit_right - emp_right,
                stringsAsFactors = FALSE
            )
        }
    }
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
}

.standardized_stats_row <- function(asset, x) {
    mm <- safe_moments(x)
    data.frame(
        asset = asset,
        nobs = length(x),
        mean = mean(x),
        sd = sd(x),
        skew = unname(mm[["skew"]]),
        excess_kurtosis = unname(mm[["excess_kurtosis"]]),
        min = min(x),
        max = max(x),
        stringsAsFactors = FALSE
    )
}

.oos_rows <- function(asset, x, train_fraction) {
    n <- length(x)
    ntrain <- floor(train_fraction * n)
    if (ntrain < MIN_FIT_N || n - ntrain < 50L) return(data.frame())
    train <- x[seq_len(ntrain)]
    test <- x[seq.int(ntrain + 1L, n)]
    cat("  OOS split:", length(train), "train,", length(test), "test\n")
    fits <- .fit_all_models(train)
    rows <- lapply(fits, function(f) {
        ld <- tryCatch(f$log_density(test), error = function(e) rep(NA_real_, length(test)))
        good <- is.finite(ld)
        total <- if (all(good)) sum(ld) else NA_real_
        data.frame(
            asset = asset,
            model = f$model,
            train_n = length(train),
            test_n = length(test),
            test_log_score = total,
            mean_test_log_score = if (is.finite(total)) total / length(test) else NA_real_,
            fit_convergence = f$convergence,
            stringsAsFactors = FALSE
        )
    })
    z <- do.call(rbind, rows)
    ok <- is.finite(z$test_log_score)
    z$log_score_rank <- NA_integer_
    if (any(ok)) z$log_score_rank[ok] <- rank(-z$test_log_score[ok], ties.method = "min")
    z
}

# ------------------------------- main run ----------------------------------

pr <- read_price_file(price_file)
prices <- keep_selected_symbols(pr$prices, SYMBOLS)
dates <- pr$dates
ret <- compute_log_returns(prices)
ret_dates <- dates[-1L]
zret <- standardize_by_ewma_vol(ret, lambda = EWMA_LAMBDA,
                                warmup = EWMA_WARMUP,
                                center = EWMA_CENTER)

# Save the standardized returns with dates for reproducibility.
zout <- data.frame(Date = ret_dates, zret, check.names = FALSE)
write.csv(zout, file.path(output_dir, "ewma_standardized_log_returns.csv"), row.names = FALSE)

summary_all <- list()
param_all <- list()
tail_all <- list()
stats_all <- list()
oos_all <- list()

for (j in seq_along(zret)) {
    asset <- names(zret)[j]
    x <- as.numeric(zret[[j]])
    x <- x[is.finite(x)]
    cat("\n", asset, ": ", length(x), " standardized returns\n", sep = "")

    if (length(x) < MIN_FIT_N || sd(x) <= 0) {
        warning(asset, ": insufficient usable observations; skipping")
        next
    }

    stats_all[[asset]] <- .standardized_stats_row(asset, x)
    fits <- .fit_all_models(x)
    summary_all[[asset]] <- .fit_summary_rows(asset, fits)
    param_all[[asset]] <- .parameter_rows(asset, fits)
    tail_all[[asset]] <- .tail_rows(asset, x, fits, TAIL_THRESHOLDS)

    if (RUN_OOS) oos_all[[asset]] <- .oos_rows(asset, x, OOS_TRAIN_FRACTION)
}

bind_nonempty <- function(x) {
    x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
    if (!length(x)) data.frame() else do.call(rbind, x)
}

fit_summary <- bind_nonempty(summary_all)
parameters <- bind_nonempty(param_all)
tails <- bind_nonempty(tail_all)
stats <- bind_nonempty(stats_all)
oos <- bind_nonempty(oos_all)

model_goodness <- .model_goodness_summary(fit_summary)

write.csv(fit_summary, file.path(output_dir, "distribution_fit_summary.csv"), row.names = FALSE)
write.csv(model_goodness, file.path(output_dir, "model_goodness_of_fit_summary.csv"), row.names = FALSE)
write.csv(parameters, file.path(output_dir, "parameter_estimates.csv"), row.names = FALSE)
write.csv(tails, file.path(output_dir, "tail_diagnostics.csv"), row.names = FALSE)
write.csv(stats, file.path(output_dir, "standardized_return_stats.csv"), row.names = FALSE)
if (RUN_OOS) write.csv(oos, file.path(output_dir, "oos_log_scores.csv"), row.names = FALSE)

cat("\nFull-sample ranking by AIC\n")
if (nrow(fit_summary)) {
    tmp <- fit_summary[order(fit_summary$asset, fit_summary$AIC),
                       c("asset", "model", "logLik", "AIC", "BIC", "delta_AIC",
                         "AIC_rank", "BIC_rank", "fit_seconds")]
    tmp$fit_seconds <- sprintf("%.4f", tmp$fit_seconds)
    print(tmp, row.names = FALSE)
}

cat("\nGoodness-of-fit summary across assets (full-sample AIC/BIC)\n")
if (nrow(model_goodness)) {
    tmp <- model_goodness
    tmp$mean_AIC_rank <- sprintf("%.2f", tmp$mean_AIC_rank)
    tmp$mean_BIC_rank <- sprintf("%.2f", tmp$mean_BIC_rank)
    tmp$mean_delta_AIC <- sprintf("%.3f", tmp$mean_delta_AIC)
    tmp$mean_delta_BIC <- sprintf("%.3f", tmp$mean_delta_BIC)
    print(tmp, row.names = FALSE)
}

if (RUN_OOS && nrow(oos)) {
    cat("\nOut-of-sample ranking by log score\n")
    tmp <- oos[order(oos$asset, oos$log_score_rank), ]
    print(tmp, row.names = FALSE)
}

cat("\nWrote results to:", output_dir, "\n")
cat("Total elapsed time:",
    format(round(as.numeric(difftime(Sys.time(), script_start_time, units = "secs")), 2)),
    "seconds\n")
