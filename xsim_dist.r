#!/usr/bin/env Rscript
# xsim_dist.r
#
# Simulate return series from selected distributions and fit models directly,
# reporting how often the generating model wins AIC/BIC.

options(width = 200)

# ----------------------------
# user options
# ----------------------------

n_sims <- 3 # 50
n_obs <- 100 # 2500
seed <- 123
ret_scale <- 100
fixed_nu <- NA_real_

sim_models <- c(
  "Normal",
  "Logistic",
  "Laplace",
  "ALaplace",
  "T",
  "NCT",
  "SkewNormal",
  "SkewT",
  "Cauchy"
)

fit_models <- c(
  "Normal",
  "Logistic",
  "Laplace",
  "ALaplace",
  "T",
  "NCT",
  "SkewNormal",
  "SkewT",
  "Cauchy"
)

sim_params <- list(
  Normal = list(mean = 0, sd = 1),
  Logistic = list(mu = 0, s = 1),
  Laplace = list(mu = 0, b = 1),
  ALaplace = list(mu = 0, b = 1, kappa = 1.5),
  T = list(xi = 0, omega = 1, nu = 5),
  NCT = list(mu = 0, sigma = 1, nu = 5, ncp = 0.75),
  SkewNormal = list(xi = 0, omega = 1, alpha = 3),
  SkewT = list(xi = 0, omega = 1, alpha = 3, nu = 6),
  Cauchy = list(location = 0, scale = 1)
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_sims <- as.integer(args[1])
if (length(args) >= 2) n_obs <- as.integer(args[2])

#' Stop if a required package is not installed.
require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("package '%s' is required for %s", pkg, purpose))
  }
}

if (any(sim_models %in% c("SkewNormal", "SkewT", "T")) || any(fit_models %in% c("SkewNormal", "SkewT", "T"))) {
  require_pkg("sn", "skew-normal/skew-t simulation or fitting")
  suppressPackageStartupMessages(library("sn"))
}

source("dist_utils.r")
source("return_utils.r")

set.seed(seed)
sim_start <- proc.time()

#' Run sn::selm quietly, suppressing output and warnings.
quiet_selm <- function(expr) {
  out <- NULL
  capture.output({
    out <<- suppressWarnings(try(expr, silent = TRUE))
  })
  out
}

#' Simulate Laplace random variates.
r_laplace <- function(n, mu, b) {
  s <- sample(c(-1, 1), size = n, replace = TRUE)
  e <- rexp(n, rate = 1 / b)
  mu + s * e
}

#' Simulate asymmetric Laplace random variates.
r_alaplace <- function(n, mu, b, kappa) {
  b_pos <- b / kappa
  b_neg <- b * kappa
  p_pos <- b_pos / (b_pos + b_neg)
  u <- runif(n)
  x <- numeric(n)
  pos <- u < p_pos
  x[pos] <- mu + rexp(sum(pos), rate = 1 / b_pos)
  x[!pos] <- mu - rexp(sum(!pos), rate = 1 / b_neg)
  x
}

#' Simulate Azzalini skew-normal random variates.
r_skewnormal <- function(n, xi, omega, alpha) {
  sn::rsn(n, xi = xi, omega = omega, alpha = alpha)
}

#' Simulate Azzalini skew-t random variates.
r_skewt <- function(n, xi, omega, alpha, nu) {
  sn::rst(n, xi = xi, omega = omega, alpha = alpha, nu = nu)
}

#' Simulate symmetric t random variates via skew-t with alpha = 0.
r_symt <- function(n, xi, omega, nu) {
  sn::rst(n, xi = xi, omega = omega, alpha = 0, nu = nu)
}

#' Simulate noncentral t random variates.
r_nct <- function(n, mu, sigma, nu, ncp) {
  mu + sigma * rt(n, df = nu, ncp = ncp)
}

#' Simulate returns for a given model.
sim_returns <- function(model, n, params) {
  if (model == "Normal") {
    return(rnorm(n, mean = params$mean, sd = params$sd))
  }
  if (model == "Logistic") {
    return(rlogis(n, location = params$mu, scale = params$s))
  }
  if (model == "Laplace") {
    return(r_laplace(n, params$mu, params$b))
  }
  if (model == "ALaplace") {
    return(r_alaplace(n, params$mu, params$b, params$kappa))
  }
  if (model == "T") {
    return(r_symt(n, params$xi, params$omega, params$nu))
  }
  if (model == "NCT") {
    return(r_nct(n, params$mu, params$sigma, params$nu, params$ncp))
  }
  if (model == "SkewNormal") {
    return(r_skewnormal(n, params$xi, params$omega, params$alpha))
  }
  if (model == "SkewT") {
    return(r_skewt(n, params$xi, params$omega, params$alpha, params$nu))
  }
  if (model == "Cauchy") {
    return(rcauchy(n, location = params$location, scale = params$scale))
  }
  stop(sprintf("no simulator for model %s", model))
}

#' Negative log-likelihood for logistic distribution.
logistic_nll <- function(par, x) {
  mu <- par[1]
  s <- exp(par[2])
  ll <- logistic_logpdf(x, mu, s)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit logistic distribution by MLE.
fit_logistic <- function(x) {
  mu0 <- median(x)
  s0 <- sd(x) * sqrt(3) / pi
  par0 <- c(mu0, log(s0))
  opt <- try(optim(par0, logistic_nll, x = x, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  s <- exp(par[2])
  list(mu = mu, s = s, logLik = -opt$value)
}

#' Fit asymmetric Laplace by MLE.
fit_alaplace <- function(x) {
  if (length(x) < 3) return(NULL)
  mu0 <- median(x)
  b0 <- mean(abs(x - mu0))
  kappa0 <- 1
  par0 <- c(mu0, log(b0), log(kappa0))
  nll <- function(par, x) {
    mu <- par[1]
    b <- exp(par[2])
    kappa <- exp(par[3])
    ll <- ald_logpdf(x, mu, b, kappa)
    if (any(!is.finite(ll))) return(Inf)
    -sum(ll)
  }
  opt <- try(optim(par0, nll, x = x, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  b <- exp(par[2])
  kappa <- exp(par[3])
  list(mu = mu, b = b, kappa = kappa, logLik = -opt$value)
}

#' Negative log-likelihood for non-central t.
nct_nll <- function(par, x, fixed_nu = NA_real_) {
  mu <- par[1]
  sigma <- exp(par[2])
  if (!is.finite(sigma) || sigma <= 0) return(Inf)
  if (is.finite(fixed_nu)) {
    nu <- fixed_nu
    ncp <- par[3]
  } else {
    nu <- 1 + exp(par[3])
    ncp <- par[4]
  }
  ll <- dt((x - mu) / sigma, df = nu, ncp = ncp, log = TRUE) - log(sigma)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit non-central t by MLE.
fit_nct <- function(x, fixed_nu = NA_real_) {
  mu0 <- mean(x)
  sigma0 <- sd(x)
  ncp0 <- 0
  nu0 <- if (is.finite(fixed_nu)) fixed_nu else 6
  if (is.finite(fixed_nu)) {
    par0 <- c(mu0, log(sigma0), ncp0)
    opt <- try(optim(par0, nct_nll, x = x, fixed_nu = fixed_nu, method = "BFGS"), silent = TRUE)
  } else {
    par0 <- c(mu0, log(sigma0), log(nu0 - 1), ncp0)
    opt <- try(optim(par0, nct_nll, x = x, fixed_nu = fixed_nu, method = "BFGS"), silent = TRUE)
  }
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- if (is.finite(fixed_nu)) fixed_nu else 1 + exp(par[3])
  ncp <- if (is.finite(fixed_nu)) par[3] else par[4]
  list(mu = mu, sigma = sigma, nu = nu, ncp = ncp, logLik = -opt$value)
}

#' Fit a model and return logLik/AIC/BIC.
fit_model <- function(x, model, fixed_nu = NA_real_) {
  n <- length(x)
  if (model == "Normal") {
    mu <- mean(x)
    sigma <- sd(x)
    ll <- sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
    metrics <- fit_metrics(ll, 2, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "Logistic") {
    fit <- fit_logistic(x)
    if (is.null(fit)) return(NULL)
    metrics <- fit_metrics(fit$logLik, 2, n)
    return(list(logLik = fit$logLik, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "Laplace") {
    mu <- median(x)
    b <- mean(abs(x - mu))
    ll <- sum(laplace_logpdf(x, mu, b))
    metrics <- fit_metrics(ll, 2, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "ALaplace") {
    fit <- fit_alaplace(x)
    if (is.null(fit)) return(NULL)
    metrics <- fit_metrics(fit$logLik, 3, n)
    return(list(logLik = fit$logLik, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "T") {
    fixed_list <- list(alpha = 0)
    if (is.finite(fixed_nu)) fixed_list$nu <- fixed_nu
    fit <- quiet_selm(sn::selm(x ~ 1, family = "ST", fixed.param = fixed_list, control = list(trace = 0)))
    if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
    ll <- as.numeric(logLik(fit))
    k <- if (is.finite(fixed_nu)) 2 else 3
    metrics <- fit_metrics(ll, k, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "SkewNormal") {
    fit <- quiet_selm(sn::selm(x ~ 1, family = "SN", control = list(trace = 0)))
    if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
    ll <- as.numeric(logLik(fit))
    metrics <- fit_metrics(ll, 3, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "SkewT") {
    fixed_list <- list()
    if (is.finite(fixed_nu)) fixed_list$nu <- fixed_nu
    fit <- quiet_selm(sn::selm(x ~ 1, family = "ST", fixed.param = fixed_list, control = list(trace = 0)))
    if (is.null(fit) || inherits(fit, "try-error")) return(NULL)
    ll <- as.numeric(logLik(fit))
    k <- if (is.finite(fixed_nu)) 3 else 4
    metrics <- fit_metrics(ll, k, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "NCT") {
    fit <- fit_nct(x, fixed_nu = fixed_nu)
    if (is.null(fit)) return(NULL)
    k <- if (is.finite(fixed_nu)) 3 else 4
    metrics <- fit_metrics(fit$logLik, k, n)
    return(list(logLik = fit$logLik, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  if (model == "Cauchy") {
    nll <- function(par, x) {
      loc <- par[1]
      scale <- exp(par[2])
      ll <- dcauchy(x, location = loc, scale = scale, log = TRUE)
      if (any(!is.finite(ll))) return(Inf)
      -sum(ll)
    }
    par0 <- c(median(x), log(mad(x)))
    opt <- try(optim(par0, nll, x = x, method = "BFGS"), silent = TRUE)
    if (inherits(opt, "try-error")) return(NULL)
    ll <- -opt$value
    metrics <- fit_metrics(ll, 2, n)
    return(list(logLik = ll, AIC = metrics$AIC, BIC = metrics$BIC))
  }
  warning(sprintf("no fit implementation for model %s", model))
  NULL
}

cat(sprintf("simulations per model: %d\n", n_sims))
cat(sprintf("observations per sim: %d\n", n_obs))
cat(sprintf("ret_scale: %g\n", ret_scale))
cat(sprintf("fit models: %s\n", paste(fit_models, collapse = ", ")))

results <- list()

for (model in sim_models) {
  if (!model %in% names(sim_params)) {
    warning(sprintf("no sim params for %s; skipping", model))
    next
  }
  for (i in seq_len(n_sims)) {
    r <- sim_returns(model, n_obs, sim_params[[model]])
    r <- r * ret_scale
    fit_rows <- list()
    for (fm in fit_models) {
      fit <- fit_model(r, fm, fixed_nu = fixed_nu)
      if (is.null(fit)) next
      fit_rows[[fm]] <- data.frame(
        model = fm,
        logLik = fit$logLik,
        AIC = fit$AIC,
        BIC = fit$BIC,
        stringsAsFactors = FALSE
      )
    }
    fits <- do.call(rbind, fit_rows)
    if (is.null(fits) || nrow(fits) == 0) {
      best_aic <- NA_character_
      best_bic <- NA_character_
      best_ll <- NA_character_
    } else {
      best_aic <- fits$model[which.min(fits$AIC)]
      best_bic <- fits$model[which.min(fits$BIC)]
      best_ll <- fits$model[which.max(fits$logLik)]
    }
    results[[length(results) + 1]] <- data.frame(
      sim_model = model,
      sim_id = i,
      best_logLik = best_ll,
      best_aic = best_aic,
      best_bic = best_bic,
      stringsAsFactors = FALSE
    )
  }
}

res <- do.call(rbind, results)
if (is.null(res)) {
  stop("no simulations were run")
}

cat("\nLogLik selection counts\n")
print(table(res$sim_model, res$best_logLik), quote = FALSE)

cat("\nAIC selection counts\n")
print(table(res$sim_model, res$best_aic), quote = FALSE)

cat("\nBIC selection counts\n")
print(table(res$sim_model, res$best_bic), quote = FALSE)

acc <- aggregate(
  data.frame(
    logLik_correct = res$sim_model == res$best_logLik,
    aic_correct = res$sim_model == res$best_aic,
    bic_correct = res$sim_model == res$best_bic
  ),
  by = list(sim_model = res$sim_model),
  FUN = mean
)
acc$logLik_correct <- round(acc$logLik_correct, 3)
acc$aic_correct <- round(acc$aic_correct, 3)
acc$bic_correct <- round(acc$bic_correct, 3)

cat("\nSelection accuracy (share correct)\n")
print(acc, row.names = FALSE)

elapsed <- (proc.time() - sim_start)[["elapsed"]]
cat(sprintf("\nTotal time elapsed: %.2f seconds\n", elapsed))
