#!/usr/bin/env Rscript
# xreturns_dist.r
#
# Read a CSV of prices with first column Date and remaining columns assets,
# compute simple or log returns per asset, fit Azzalini and Fernandez-Steel
# skewed distributions, and optionally plot fitted densities.

# ----------------------------
# user options
# ----------------------------
options(width = 260)

infile <- ""          # if empty, will use first command-line arg
return_type <- "log"  # "log" or "simple"
models_to_fit <- c(
  "Normal",
  "Logistic",
  "EGB2",
  "Cauchy",
  "Laplace",
  "ALaplace",
  "NIG",
  "Hyperbolic",
  "Champernowne",
  "NormalLaplace",
  "GT",
  "SGT",
  "GH",
  "VG",
  "GED",
  "SGED",
  "Sech",
  "GSH",
  "SGSH",
  "NEFGHS",
  "T",
  "NCT",
  "SkewNormal",
  "SkewT",
  "FSSkewNormal",
  "FSSkewT",
  "JFSkewT"
)
max_models <- 2 # NA_integer_ # set to a positive integer to cap models
do_density_plots <- TRUE
do_kde_plot <- FALSE
do_logcondens_plot <- FALSE
do_logspline_plot <- TRUE
fixed_nu <- 5.0 # NA_real_
ret_scale <- 100
start_time <- proc.time()

#' Split a comma-separated environment variable into a trimmed character vector.
split_env_list <- function(x) {
  trimws(unlist(strsplit(x, ",")))
}

models_env <- Sys.getenv("XRETURNS_MODELS", "")
if (nzchar(models_env)) {
  models_to_fit <- split_env_list(models_env)
}
max_models_env <- Sys.getenv("XRETURNS_MAX_MODELS", "")
if (nzchar(max_models_env)) {
  max_models_val <- suppressWarnings(as.integer(max_models_env))
  max_models <- if (is.na(max_models_val)) NA_integer_ else max_models_val
}
plots_env <- Sys.getenv("XRETURNS_DENSITY_PLOTS", "")
if (nzchar(plots_env)) {
  do_density_plots <- tolower(plots_env) %in% c("1", "true", "t", "yes", "y")
}

# ----------------------------
# model fit summary storage
# ----------------------------
aic_rows <- list()
model_times <- list()
fit_moments <- list()
fit_params <- list(
  Normal = list(),
  Logistic = list(),
  EGB2 = list(),
  NIG = list(),
  Hyperbolic = list(),
  GH = list(),
  VG = list(),
  Champernowne = list(),
  NormalLaplace = list(),
  GT = list(),
  SGT = list(),
  Cauchy = list(),
  T = list(),
  NCT = list(),
  Laplace = list(),
  ALaplace = list(),
  GED = list(),
  SGED = list(),
  Sech = list(),
  GSH = list(),
  SGSH = list(),
  NEFGHS = list(),
  SkewNormal = list(),
  SkewT = list(),
  FSSkewNormal = list(),
  FSSkewT = list(),
  JFSkewT = list()
)

#' Track per-asset logLik/AIC/BIC for summary tables.
add_fit_metrics <- function(model, asset, logLik, aic, bic) {
  aic_rows[[length(aic_rows) + 1]] <<- data.frame(
    model = model,
    asset = asset,
    logLik = logLik,
    AIC = aic,
    BIC = bic,
    stringsAsFactors = FALSE
  )
}

#' Track per-asset implied moments.
add_fit_moments <- function(model, asset, mean, sd, skew, kurtosis) {
  fit_moments[[length(fit_moments) + 1]] <<- data.frame(
    model = model,
    asset = asset,
    mean = mean,
    sd = sd,
    skew = skew,
    kurtosis = kurtosis,
    stringsAsFactors = FALSE
  )
}

#' Model metadata for summary tables.
model_meta <- function() {
  data.frame(
    model = c(
      "Normal",
      "Logistic",
      "EGB2",
      "NIG",
      "Hyperbolic",
      "Champernowne",
      "NormalLaplace",
      "GT",
      "SGT",
      "GH",
      "VG",
      "Laplace",
      "ALaplace",
      "GED",
      "SGED",
      "Sech",
      "GSH",
      "SGSH",
      "NEFGHS",
      "T",
      "NCT",
      "SkewNormal",
      "SkewT",
      "FSSkewNormal",
      "FSSkewT",
      "JFSkewT"
    ),
    k = c(2, 2, 4, 4, 4, 3, 3, 4, 5, 5, 4, 2, 3, 3, 4, 2, 3, 4, 4, 3, 4, 3, 4, 3, 4, 4),
    stringsAsFactors = FALSE
  )
}

#' Resolve models enabled by the user list and optional cap.
enabled_models <- function() {
  base <- unique(models_to_fit)
  if (is.finite(max_models) && max_models > 0) {
    base <- base[seq_len(min(length(base), max_models))]
  }
  base
}

#' Check if a model is enabled.
is_enabled <- function(model) {
  model %in% enabled_models()
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
  list(mu = mu, s = s, logLik = -opt$value, convergence = opt$convergence)
}

#' Negative log-likelihood for EGB2 distribution.
egb2_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0) return(Inf)
  ll <- egb2_logpdf(x, mu, sigma, p, q)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit EGB2 distribution by MLE.
fit_egb2 <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  starts <- list(
    c(mu0, log(sigma0), log(2), log(2)),
    c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), log(1.5), log(1.5)),
    c(mu0, log(sigma0), log(1), log(1)),
    c(mu0, log(sigma0), log(5), log(5))
  )
  best <- NULL
  log_sigma_min <- log(sigma0 * 0.01)
  log_sigma_max <- log(sigma0 * 50)
  log_p_min <- log(0.05)
  log_p_max <- log(200)
  log_q_min <- log(0.05)
  log_q_max <- log(200)

  for (par0 in starts) {
    opt <- try(optim(
      par0,
      egb2_nll,
      x = x,
      method = "L-BFGS-B",
      lower = c(-Inf, log_sigma_min, log_p_min, log_q_min),
      upper = c(Inf, log_sigma_max, log_p_max, log_q_max)
    ), silent = TRUE)
    if (inherits(opt, "try-error")) next
    if (!is.finite(opt$value)) next
    if (is.null(best) || opt$value < best$value) best <- opt
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  list(mu = mu, sigma = sigma, p = p, q = q, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for FS skew-normal.
fs_skew_normal_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  gamma <- exp(par[3])
  -sum(fs_skew_normal_logpdf(x, mu, sigma, gamma))
}

#' Fit FS skew-normal distribution by MLE.
fit_fs_skew_normal <- function(x) {
  mu0 <- mean(x)
  sigma0 <- sd(x)
  gamma0 <- 1
  par0 <- c(mu0, log(sigma0), log(gamma0))
  opt <- try(optim(par0, fs_skew_normal_nll, x = x, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  sigma <- exp(par[2])
  gamma <- exp(par[3])
  list(mu = mu, sigma = sigma, gamma = gamma, logLik = -opt$value, convergence = opt$convergence)
}

#' Negative log-likelihood for FS skew-t.
fs_skew_t_nll <- function(par, x, fixed_nu) {
  mu <- par[1]
  sigma <- exp(par[2])
  gamma <- exp(par[3])
  nu <- if (is.finite(fixed_nu)) fixed_nu else 2 + exp(par[4])
  -sum(fs_skew_t_logpdf(x, mu, sigma, gamma, nu))
}

#' Fit FS skew-t distribution by MLE.
fit_fs_skew_t <- function(x, fixed_nu) {
  mu0 <- mean(x)
  sigma0 <- sd(x)
  gamma0 <- 1
  nu0 <- 8
  if (is.finite(fixed_nu)) {
    par0 <- c(mu0, log(sigma0), log(gamma0))
    opt <- try(optim(par0, fs_skew_t_nll, x = x, fixed_nu = fixed_nu, method = "BFGS"), silent = TRUE)
  } else {
    par0 <- c(mu0, log(sigma0), log(gamma0), log(nu0 - 2))
    opt <- try(optim(par0, fs_skew_t_nll, x = x, fixed_nu = fixed_nu, method = "BFGS"), silent = TRUE)
  }
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  sigma <- exp(par[2])
  gamma <- exp(par[3])
  nu <- if (is.finite(fixed_nu)) fixed_nu else 2 + exp(par[4])
  list(mu = mu, sigma = sigma, gamma = gamma, nu = nu, logLik = -opt$value, convergence = opt$convergence)
}

#' Negative log-likelihood for Jones-Faddy skew-t.
jf_skew_t_nll <- function(par, x, fixed_nu) {
  mu <- par[1]
  sigma <- exp(par[2])
  if (is.finite(fixed_nu)) {
    delta <- par[3]
    a <- fixed_nu * exp(delta) / (1 + exp(delta))
    b <- fixed_nu - a
  } else {
    a <- exp(par[3])
    b <- exp(par[4])
  }
  ll <- jf_skew_t_logpdf(x, mu, sigma, a, b)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Negative log-likelihood for asymmetric Laplace.
ald_nll <- function(par, x) {
  mu <- par[1]
  b <- exp(par[2])
  kappa <- exp(par[3])
  ll <- ald_logpdf(x, mu, b, kappa)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit asymmetric Laplace distribution by MLE.
fit_alaplace <- function(x) {
  mu0 <- median(x)
  b0 <- mean(abs(x - mu0))
  k0 <- 1
  starts <- list(
    c(mu0, log(b0), log(k0)),
    c(mean(x), log(sd(x)), 0)
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, ald_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  b <- exp(par[2])
  kappa <- exp(par[3])
  list(mu = mu, b = b, kappa = kappa, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for NIG.
nig_nll <- function(par, x) {
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  eta <- par[4]
  beta <- alpha * tanh(eta)
  ll <- nig_logpdf(x, alpha, beta, delta, mu)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit NIG distribution by MLE.
fit_nig <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  alpha0 <- if (is.finite(sigma0) && sigma0 > 0) 1 / sigma0 else 1
  delta0 <- if (is.finite(sigma0) && sigma0 > 0) sigma0 else 1
  starts <- list(
    c(mu0, log(delta0), log(alpha0), 0),
    c(med0, log(delta0), log(alpha0 * 2), 0)
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, nig_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  beta <- alpha * tanh(par[4])
  list(mu = mu, delta = delta, alpha = alpha, beta = beta, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Hyperbolic.
hyperbolic_nll <- function(par, x) {
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  eta <- par[4]
  beta <- alpha * tanh(eta)
  ll <- hyperbolic_logpdf(x, mu, delta, alpha, beta)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Hyperbolic distribution by MLE.
fit_hyperbolic <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  alpha0 <- if (is.finite(sigma0) && sigma0 > 0) 1 / sigma0 else 1
  delta0 <- if (is.finite(sigma0) && sigma0 > 0) sigma0 else 1
  starts <- list(
    c(mu0, log(delta0), log(alpha0), 0),
    c(med0, log(delta0), log(alpha0 * 2), 0)
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, hyperbolic_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  beta <- alpha * tanh(par[4])
  list(mu = mu, delta = delta, alpha = alpha, beta = beta, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Champernowne.
champernowne_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  lambda <- exp(par[3])
  ll <- champernowne_logpdf(x, mu, sigma, lambda)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Champernowne distribution by MLE.
fit_champernowne <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  lambda0 <- 1
  starts <- list(
    c(mu0, log(sigma0), log(lambda0)),
    c(med0, log(sigma0), log(2))
  )
  methods <- c("L-BFGS-B")
  best <- NULL
  mu_min <- mu0 - 10 * sigma0
  mu_max <- mu0 + 10 * sigma0
  log_sigma_min <- log(sigma0 * 0.05)
  log_sigma_max <- log(sigma0 * 20)
  log_lambda_min <- log(1e-6)
  log_lambda_max <- log(50)

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(
        par0,
        champernowne_nll,
        x = x,
        method = method,
        lower = c(mu_min, log_sigma_min, log_lambda_min),
        upper = c(mu_max, log_sigma_max, log_lambda_max),
        control = list(maxit = 200)
      ), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  lambda <- exp(par[3])
  list(mu = mu, sigma = sigma, lambda = lambda, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Normal-Laplace.
normal_laplace_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  b <- exp(par[3])
  ll <- normal_laplace_logpdf(x, mu, sigma, b)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Normal-Laplace distribution by MLE.
fit_normal_laplace <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  s0 <- sd(x)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  sigma0 <- s0 / sqrt(3)
  b0 <- s0 / sqrt(3)
  starts <- list(
    c(mu0, log(sigma0), log(b0)),
    c(med0, log(sigma0), log(b0))
  )
  methods <- c("L-BFGS-B")
  best <- NULL
  mu_min <- mu0 - 10 * s0
  mu_max <- mu0 + 10 * s0
  log_sigma_min <- log(s0 * 0.05)
  log_sigma_max <- log(s0 * 20)
  log_b_min <- log(s0 * 0.01)
  log_b_max <- log(s0 * 10)

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(
        par0,
        normal_laplace_nll,
        x = x,
        method = method,
        lower = c(mu_min, log_sigma_min, log_b_min),
        upper = c(mu_max, log_sigma_max, log_b_max),
        control = list(maxit = 200)
      ), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  b <- exp(par[3])
  list(mu = mu, sigma = sigma, b = b, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Generalized t (McDonald-Newey).
gt_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  ll <- gt_logpdf(x, mu, sigma, p, q)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Generalized t distribution by MLE.
fit_gt <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  s0 <- sd(x)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  p0 <- 2
  q0 <- 2
  starts <- list(
    c(mu0, log(s0), log(p0), log(q0)),
    c(med0, log(s0), log(1.5), log(1.5))
  )
  methods <- c("L-BFGS-B")
  best <- NULL
  mu_min <- mu0 - 10 * s0
  mu_max <- mu0 + 10 * s0
  log_sigma_min <- log(s0 * 0.05)
  log_sigma_max <- log(s0 * 20)
  log_p_min <- log(0.2)
  log_p_max <- log(20)
  log_q_min <- log(0.2)
  log_q_max <- log(20)

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(
        par0,
        gt_nll,
        x = x,
        method = method,
        lower = c(mu_min, log_sigma_min, log_p_min, log_q_min),
        upper = c(mu_max, log_sigma_max, log_p_max, log_q_max),
        control = list(maxit = 200)
      ), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  list(mu = mu, sigma = sigma, p = p, q = q, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Cauchy.
cauchy_nll <- function(par, x) {
  x0 <- par[1]
  gamma <- exp(par[2])
  ll <- cauchy_logpdf(x, x0, gamma)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Cauchy distribution by MLE.
fit_cauchy <- function(x) {
  x0 <- median(x)
  gamma0 <- mad(x)
  if (!is.finite(gamma0) || gamma0 <= 0) gamma0 <- sd(x)
  if (!is.finite(gamma0) || gamma0 <= 0) gamma0 <- 1
  par0 <- c(x0, log(gamma0))
  opt <- try(optim(par0, cauchy_nll, x = x, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  list(x0 = par[1], gamma = exp(par[2]), logLik = -opt$value, convergence = opt$convergence)
}

#' Negative log-likelihood for Skewed Generalized t (Theodossiou).
sgt_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  lambda <- tanh(par[5])
  ll <- sgt_logpdf(x, mu, sigma, p, q, lambda)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Skewed Generalized t distribution by MLE.
fit_sgt <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  s0 <- sd(x)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  p0 <- 2
  q0 <- 2
  lambda0 <- 0
  starts <- list(
    c(mu0, log(s0), log(p0), log(q0), atanh(lambda0)),
    c(med0, log(s0), log(1.5), log(1.5), 0)
  )
  methods <- c("L-BFGS-B")
  best <- NULL
  mu_min <- mu0 - 10 * s0
  mu_max <- mu0 + 10 * s0
  log_sigma_min <- log(s0 * 0.05)
  log_sigma_max <- log(s0 * 20)
  log_p_min <- log(0.2)
  log_p_max <- log(20)
  log_q_min <- log(0.2)
  log_q_max <- log(20)
  lambda_min <- -3
  lambda_max <- 3

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(
        par0,
        sgt_nll,
        x = x,
        method = method,
        lower = c(mu_min, log_sigma_min, log_p_min, log_q_min, lambda_min),
        upper = c(mu_max, log_sigma_max, log_p_max, log_q_max, lambda_max),
        control = list(maxit = 200)
      ), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  p <- exp(par[3])
  q <- exp(par[4])
  lambda <- tanh(par[5])
  list(mu = mu, sigma = sigma, p = p, q = q, lambda = lambda, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Generalized Hyperbolic.
gh_nll <- function(par, x) {
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  eta <- par[4]
  beta <- alpha * tanh(eta)
  lambda <- par[5]
  ll <- gh_logpdf(x, mu, delta, alpha, beta, lambda)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Generalized Hyperbolic distribution by MLE.
fit_gh <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  alpha0 <- 1 / sigma0
  delta0 <- sigma0
  mu_min <- mu0 - 10 * sigma0
  mu_max <- mu0 + 10 * sigma0
  log_delta_min <- log(sigma0 * 0.05)
  log_delta_max <- log(sigma0 * 20)
  log_alpha_min <- log(1 / (sigma0 * 20))
  log_alpha_max <- log(1 / (sigma0 * 0.05))
  eta_min <- -3
  eta_max <- 3
  lambda_min <- -5
  lambda_max <- 5
  starts <- list(
    c(mu0, log(delta0), log(alpha0), 0, 1),
    c(med0, log(delta0), log(alpha0 * 2), 0, 1)
  )
  best <- NULL

  for (par0 in starts) {
    opt <- try(optim(
      par0,
      gh_nll,
      x = x,
      method = "L-BFGS-B",
      lower = c(mu_min, log_delta_min, log_alpha_min, eta_min, lambda_min),
      upper = c(mu_max, log_delta_max, log_alpha_max, eta_max, lambda_max),
      control = list(maxit = 200)
    ), silent = TRUE)
    if (inherits(opt, "try-error")) next
    if (!is.finite(opt$value)) next
    if (is.null(best) || opt$value < best$value) best <- opt
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  delta <- exp(par[2])
  alpha <- exp(par[3])
  beta <- alpha * tanh(par[4])
  lambda <- par[5]
  list(mu = mu, delta = delta, alpha = alpha, beta = beta, lambda = lambda, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for Variance-Gamma.
vg_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  theta <- par[3]
  nu <- exp(par[4])
  ll <- vg_logpdf(x, mu, sigma, theta, nu)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit Variance-Gamma distribution by MLE.
fit_vg <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  theta0 <- 0
  nu0 <- 0.5
  mu_min <- mu0 - 10 * sigma0
  mu_max <- mu0 + 10 * sigma0
  log_sigma_min <- log(sigma0 * 0.05)
  log_sigma_max <- log(sigma0 * 20)
  theta_min <- -5 * sigma0
  theta_max <- 5 * sigma0
  log_nu_min <- log(0.05)
  log_nu_max <- log(10)
  starts <- list(
    c(mu0, log(sigma0), theta0, log(nu0)),
    c(med0, log(ifelse(is.finite(sigma0) && sigma0 > 0, sigma0, 1)), theta0, log(0.2))
  )
  methods <- c("L-BFGS-B")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(
        par0,
        vg_nll,
        x = x,
        method = method,
        lower = c(mu_min, log_sigma_min, theta_min, log_nu_min),
        upper = c(mu_max, log_sigma_max, theta_max, log_nu_max),
        control = list(maxit = 200)
      ), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  theta <- par[3]
  nu <- exp(par[4])
  list(mu = mu, sigma = sigma, theta = theta, nu = nu, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for GED.
ged_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- exp(par[3])
  ll <- ged_logpdf(x, mu, sigma, nu)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit GED distribution by MLE.
fit_ged <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  nu0 <- 1.5
  starts <- list(
    c(mu0, log(sigma0), log(nu0)),
    c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), log(2))
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, ged_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- exp(par[3])
  list(mu = mu, sigma = sigma, nu = nu, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for skewed GED.
sged_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- exp(par[3])
  kappa <- exp(par[4])
  ll <- sged_logpdf(x, mu, sigma, nu, kappa)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit skewed GED distribution by MLE.
fit_sged <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  nu0 <- 1.5
  k0 <- 1
  starts <- list(
    c(mu0, log(sigma0), log(nu0), log(k0)),
    c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), log(2), 0)
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, sged_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- exp(par[3])
  kappa <- exp(par[4])
  list(mu = mu, sigma = sigma, nu = nu, kappa = kappa, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for GSH.
gsh_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  t <- par[3]
  ll <- gsh_logpdf(x, mu, sigma, t)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit GSH distribution by MLE.
fit_gsh <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  t_starts <- c(-pi / 2, -0.5, 0.5)
  starts <- lapply(t_starts, function(t0) c(mu0, log(sigma0), t0))
  starts <- c(
    starts,
    list(c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), -0.5))
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, gsh_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  t <- par[3]
  list(mu = mu, sigma = sigma, t = t, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for SGSH.
sgsh_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  t <- par[3]
  kappa <- exp(par[4])
  ll <- sgsh_logpdf(x, mu, sigma, t, kappa)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit SGSH distribution by MLE.
fit_sgsh <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  t_starts <- c(-pi / 2, -0.5, 0.5)
  starts <- lapply(t_starts, function(t0) c(mu0, log(sigma0), t0, 0))
  starts <- c(
    starts,
    list(c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), -0.5, 0))
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, sgsh_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  t <- par[3]
  kappa <- exp(par[4])
  list(mu = mu, sigma = sigma, t = t, kappa = kappa, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for NEF-GHS.
nef_ghs_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  lambda <- exp(par[3])
  beta <- par[4]
  ll <- nef_ghs_logpdf(x, mu, sigma, lambda, beta)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit NEF-GHS distribution by MLE.
fit_nef_ghs <- function(x) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  lambda0 <- 0.5
  starts <- list(
    c(mu0, log(sigma0), log(lambda0), 0),
    c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), log(1), 0),
    c(mu0, log(sigma0), log(2), 0.5),
    c(mu0, log(sigma0), log(2), -0.5)
  )
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, nef_ghs_nll, x = x, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  lambda <- exp(par[3])
  beta <- par[4]
  list(mu = mu, sigma = sigma, lambda = lambda, beta = beta, logLik = -best$value, convergence = best$convergence)
}
#' Negative log-likelihood for hyperbolic secant.
sech_nll <- function(par, x) {
  mu <- par[1]
  sigma <- exp(par[2])
  ll <- sech_logpdf(x, mu, sigma)
  if (any(!is.finite(ll))) return(Inf)
  -sum(ll)
}

#' Fit hyperbolic secant distribution by MLE.
fit_sech <- function(x) {
  mu0 <- mean(x)
  sigma0 <- sd(x)
  par0 <- c(mu0, log(sigma0))
  opt <- try(optim(par0, sech_nll, x = x, method = "BFGS"), silent = TRUE)
  if (inherits(opt, "try-error")) return(NULL)
  par <- opt$par
  mu <- par[1]
  sigma <- exp(par[2])
  list(mu = mu, sigma = sigma, logLik = -opt$value, convergence = opt$convergence)
}
#' Fit Jones-Faddy skew-t distribution by MLE.
fit_jf_skew_t <- function(x, fixed_nu) {
  mu0 <- mean(x)
  med0 <- median(x)
  sigma0 <- sd(x)
  mad0 <- mad(x)
  a0 <- 5
  b0 <- 5

  if (is.finite(fixed_nu)) {
    starts <- list(
      c(mu0, log(sigma0), 0),
      c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), 0)
    )
  } else {
    starts <- list(
      c(mu0, log(sigma0), log(a0), log(b0)),
      c(med0, log(ifelse(is.finite(mad0) && mad0 > 0, mad0 * 1.4826, sigma0)), log(8), log(8)),
      c(mu0, log(sigma0), log(2), log(2))
    )
  }
  methods <- c("BFGS", "Nelder-Mead")
  best <- NULL

  for (par0 in starts) {
    for (method in methods) {
      opt <- try(optim(par0, jf_skew_t_nll, x = x, fixed_nu = fixed_nu, method = method), silent = TRUE)
      if (inherits(opt, "try-error")) next
      if (!is.finite(opt$value)) next
      if (is.null(best) || opt$value < best$value) best <- opt
    }
  }

  if (is.null(best)) return(NULL)
  par <- best$par
  mu <- par[1]
  sigma <- exp(par[2])
  if (is.finite(fixed_nu)) {
    delta <- par[3]
    a <- fixed_nu * exp(delta) / (1 + exp(delta))
    b <- fixed_nu - a
  } else {
    a <- exp(par[3])
    b <- exp(par[4])
  }
  list(mu = mu, sigma = sigma, a = a, b = b, logLik = -best$value, convergence = best$convergence)
}

#' Negative log-likelihood for non-central t with location and scale.
nct_nll <- function(par, x, fixed_nu) {
  mu <- par[1]
  sigma <- exp(par[2])
  nu <- if (is.finite(fixed_nu)) fixed_nu else 1 + exp(par[3])
  ncp <- if (is.finite(fixed_nu)) par[3] else par[4]
  -sum(nct_logpdf(x, mu, sigma, nu, ncp))
}

#' Fit non-central t distribution by MLE.
fit_nct <- function(x, fixed_nu) {
  mu0 <- mean(x)
  sigma0 <- sd(x)
  nu0 <- 8
  ncp0 <- 0
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
  list(mu = mu, sigma = sigma, nu = nu, ncp = ncp, logLik = -opt$value, convergence = opt$convergence)
}

# ----------------------------
# helpers
# ----------------------------
source("return_utils.r")
source("dist_utils.r")
source("plot_utils.r")

#' Run sn::selm quietly, suppressing output and warnings.
quiet_selm <- function(expr) {
  out <- NULL
  capture.output({
    out <<- suppressWarnings(try(expr, silent = TRUE))
  })
  out
}

# ----------------------------
# package checks
# ----------------------------
require_pkg <- function(pkg, purpose) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("package '%s' is required for %s", pkg, purpose))
  }
}

require_pkg("sn", "skew-normal/skew-t models")
if (isTRUE(do_density_plots) && isTRUE(do_logcondens_plot)) {
  require_pkg("logcondens", "logcondens density")
}
if (isTRUE(do_density_plots) && isTRUE(do_logspline_plot)) {
  require_pkg("logspline", "logspline density")
}
require_pkg("tseries", "Jarque-Bera p-values")
require_pkg("fBasics", "D'Agostino p-values")
require_pkg("nortest", "Anderson-Darling p-values")

library("sn")

# ----------------------------
# input
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (nchar(infile) == 0) {
  if (length(args) >= 1) infile <- args[1]
}

if (nchar(infile) == 0) {
  stop("provide infile either by setting infile <- 'path.csv' or passing it as the first argument")
}

df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

if (!("Date" %in% names(df))) {
  stop("first column must be named Date")
}

df$Date <- as.Date(df$Date)
df <- df[order(df$Date), , drop = FALSE]

asset_names <- unique(setdiff(names(df), "Date"))
if (length(asset_names) != length(setdiff(names(df), "Date"))) {
  warning("duplicate asset column names detected; using unique names only")
}
if (length(asset_names) == 0) stop("no asset columns found (need at least one column besides Date)")

# coerce price columns to numeric
for (nm in asset_names) {
  df[[nm]] <- as.numeric(df[[nm]])
}

cat(sprintf("file: %s\n", infile))
cat(sprintf("rows: %d  date range: %s to %s\n",
            nrow(df), as.character(min(df$Date, na.rm = TRUE)), as.character(max(df$Date, na.rm = TRUE))))
cat(sprintf("return_type: %s\n", return_type))
cat(sprintf("ret_scale: %g\n", ret_scale))
if (is.finite(fixed_nu)) {
  cat(sprintf("fixed_nu: %g\n", fixed_nu))
}

# ----------------------------
# compute returns per asset
# ----------------------------
ret_list <- list()
for (nm in asset_names) {
  p <- df[[nm]]
  if (return_type == "log" && any(p <= 0, na.rm = TRUE)) {
    warning(sprintf("asset %s has non-positive prices; log returns will be NA for those segments", nm))
  }
  r <- compute_returns(p, type = return_type)
  r <- r * ret_scale
  ret_list[[nm]] <- r
}

# ----------------------------
# summary stats for returns
# ----------------------------
cat("return summary stats by asset (sample; kurtosis is excess)\n")
sumtab <- do.call(rbind, lapply(asset_names, function(nm) {
  st <- summ_stats(ret_list[[nm]])
  cbind(asset = nm, st)
}))
sumtab_print <- sumtab
if ("kurt" %in% names(sumtab_print)) {
  sumtab_print[["kurt"]] <- sumtab[["kurt"]] - 3
}
sumtab_print <- format_table(sumtab_print, digits = 4, integer_cols = "n")
print(sumtab_print, row.names = FALSE)

if (is_enabled("Normal")) {
  t_model <- proc.time()
  # ----------------------------
  # normal fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 2) {
      warning(sprintf("asset %s has too few finite returns for normal fit", nm))
      next
    }
    mu <- mean(r)
    sigma <- sd(r)
    ll <- sum(dnorm(r, mean = mu, sd = sigma, log = TRUE))
    metrics <- fit_metrics(ll, 2, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mean = mu,
      sd = sigma,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Normal[[nm]] <- list(mean = mu, sd = sigma)
    add_fit_metrics("Normal", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("Normal", nm, mu, sigma, 0, 0)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nnormal fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Normal"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Logistic")) {
  t_model <- proc.time()
  # ----------------------------
  # logistic fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 2) {
      warning(sprintf("asset %s has too few finite returns for logistic fit", nm))
      next
    }
    fit <- fit_logistic(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s logistic fit failed", nm))
      next
    }
    moments <- logistic_moments(fit$mu, fit$s)
    metrics <- fit_metrics(fit$logLik, 2, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      s = fit$s,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Logistic[[nm]] <- list(mu = fit$mu, s = fit$s)
    add_fit_metrics("Logistic", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("Logistic", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nlogistic fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Logistic"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("EGB2")) {
  t_model <- proc.time()
  # ----------------------------
  # EGB2 fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for EGB2 fit", nm))
      next
    }
    fit <- fit_egb2(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s EGB2 fit failed", nm))
      next
    }
    moments <- egb2_moments(fit$mu, fit$sigma, fit$p, fit$q)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      p = fit$p,
      q = fit$q,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$EGB2[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, p = fit$p, q = fit$q)
    add_fit_metrics("EGB2", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("EGB2", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nEGB2 fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["EGB2"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Laplace")) {
  t_model <- proc.time()
  # ----------------------------
  # Laplace fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 2) {
      warning(sprintf("asset %s has too few finite returns for Laplace fit", nm))
      next
    }
    mu <- median(r)
    b <- mean(abs(r - mu))
    ll <- sum(laplace_logpdf(r, mu, b))
    metrics <- fit_metrics(ll, 2, length(r))
    moments <- laplace_moments(mu, b)
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = mu,
      b = b,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Laplace[[nm]] <- list(mu = mu, b = b)
    add_fit_metrics("Laplace", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("Laplace", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nLaplace fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Laplace"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("ALaplace")) {
  t_model <- proc.time()
  # ----------------------------
  # asymmetric Laplace fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for asymmetric Laplace fit", nm))
      next
    }
    fit <- fit_alaplace(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s asymmetric Laplace fit failed", nm))
      next
    }
    moments <- ald_moments(fit$mu, fit$b, fit$kappa)
    metrics <- fit_metrics(fit$logLik, 3, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      b = fit$b,
      kappa = fit$kappa,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$ALaplace[[nm]] <- list(mu = fit$mu, b = fit$b, kappa = fit$kappa)
    add_fit_metrics("ALaplace", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("ALaplace", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nasymmetric Laplace fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["ALaplace"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("NIG")) {
  t_model <- proc.time()
  # ----------------------------
  # NIG fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for NIG fit", nm))
      next
    }
    fit <- fit_nig(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s NIG fit failed", nm))
      next
    }
    moments <- nig_moments(fit$alpha, fit$beta, fit$delta, fit$mu)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      delta = fit$delta,
      alpha = fit$alpha,
      beta = fit$beta,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$NIG[[nm]] <- list(mu = fit$mu, delta = fit$delta, alpha = fit$alpha, beta = fit$beta)
    add_fit_metrics("NIG", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("NIG", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nNIG fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["NIG"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Hyperbolic")) {
  t_model <- proc.time()
  # ----------------------------
  # Hyperbolic fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for Hyperbolic fit", nm))
      next
    }
    fit <- fit_hyperbolic(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s Hyperbolic fit failed", nm))
      next
    }
    moments <- hyperbolic_moments(fit$mu, fit$delta, fit$alpha, fit$beta)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      delta = fit$delta,
      alpha = fit$alpha,
      beta = fit$beta,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Hyperbolic[[nm]] <- list(mu = fit$mu, delta = fit$delta, alpha = fit$alpha, beta = fit$beta)
    add_fit_metrics("Hyperbolic", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("Hyperbolic", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nHyperbolic fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Hyperbolic"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Champernowne")) {
  t_model <- proc.time()
  # ----------------------------
  # Champernowne fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for Champernowne fit", nm))
      next
    }
    fit <- fit_champernowne(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s Champernowne fit failed", nm))
      next
    }
    moments <- champernowne_moments(fit$mu, fit$sigma, fit$lambda)
    metrics <- fit_metrics(fit$logLik, 3, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      lambda = fit$lambda,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Champernowne[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
    add_fit_metrics("Champernowne", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("Champernowne", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nChampernowne fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Champernowne"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("NormalLaplace")) {
  t_model <- proc.time()
  # ----------------------------
  # Normal-Laplace fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for Normal-Laplace fit", nm))
      next
    }
    fit <- fit_normal_laplace(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s Normal-Laplace fit failed", nm))
      next
    }
    moments <- normal_laplace_moments(fit$mu, fit$sigma, fit$b)
    metrics <- fit_metrics(fit$logLik, 3, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      b = fit$b,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$NormalLaplace[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, b = fit$b)
    add_fit_metrics("NormalLaplace", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("NormalLaplace", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nNormal-Laplace fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["NormalLaplace"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("GT")) {
  t_model <- proc.time()
  # ----------------------------
  # Generalized t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for GT fit", nm))
      next
    }
    fit <- try(fit_gt(r), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) {
      warning(sprintf("asset %s GT fit failed", nm))
      next
    }
    moments <- gt_moments(fit$mu, fit$sigma, fit$p, fit$q)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      p = fit$p,
      q = fit$q,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$GT[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, p = fit$p, q = fit$q)
    add_fit_metrics("GT", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("GT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nGT fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["GT"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Cauchy")) {
  t_model <- proc.time()
  # ----------------------------
  # Cauchy fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for Cauchy fit", nm))
      next
    }
    fit <- fit_cauchy(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s Cauchy fit failed", nm))
      next
    }
    moments <- cauchy_moments()
    ll <- sum(cauchy_logpdf(r, fit$x0, fit$gamma))
    metrics <- fit_metrics(ll, 2, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      x0 = fit$x0,
      gamma = fit$gamma,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Cauchy[[nm]] <- list(x0 = fit$x0, gamma = fit$gamma)
    add_fit_metrics("Cauchy", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("Cauchy", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nCauchy fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Cauchy"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("SGT")) {
  t_model <- proc.time()
  # ----------------------------
  # Skewed Generalized t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for SGT fit", nm))
      next
    }
    fit <- try(fit_sgt(r), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) {
      warning(sprintf("asset %s SGT fit failed", nm))
      next
    }
    moments <- sgt_moments(fit$mu, fit$sigma, fit$p, fit$q, fit$lambda)
    metrics <- fit_metrics(fit$logLik, 5, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      p = fit$p,
      q = fit$q,
      lambda = fit$lambda,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$SGT[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, p = fit$p, q = fit$q, lambda = fit$lambda)
    add_fit_metrics("SGT", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("SGT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nSGT fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["SGT"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("GH")) {
  t_model <- proc.time()
  # ----------------------------
  # Generalized Hyperbolic fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for GH fit", nm))
      next
    }
    fit <- fit_gh(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s GH fit failed", nm))
      next
    }
    moments <- gh_moments(fit$mu, fit$delta, fit$alpha, fit$beta, fit$lambda)
    metrics <- fit_metrics(fit$logLik, 5, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      delta = fit$delta,
      alpha = fit$alpha,
      beta = fit$beta,
      lambda = fit$lambda,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$GH[[nm]] <- list(mu = fit$mu, delta = fit$delta, alpha = fit$alpha, beta = fit$beta, lambda = fit$lambda)
    add_fit_metrics("GH", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("GH", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nGH fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["GH"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("VG")) {
  t_model <- proc.time()
  # ----------------------------
  # VG fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 4) {
      warning(sprintf("asset %s has too few finite returns for VG fit", nm))
      next
    }
    fit <- try(fit_vg(r), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) {
      warning(sprintf("asset %s VG fit failed", nm))
      next
    }
    moments <- try(vg_moments(fit$mu, fit$sigma, fit$theta, fit$nu), silent = TRUE)
    if (inherits(moments, "try-error")) {
      moments <- list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_)
    }
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      theta = fit$theta,
      nu = fit$nu,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$VG[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, theta = fit$theta, nu = fit$nu)
    add_fit_metrics("VG", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("VG", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nVG fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["VG"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("GED")) {
  t_model <- proc.time()
  # ----------------------------
  # GED fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for GED fit", nm))
      next
    }
    fit <- fit_ged(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s GED fit failed", nm))
      next
    }
    moments <- ged_moments(fit$mu, fit$sigma, fit$nu)
    metrics <- fit_metrics(fit$logLik, 3, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      shape = fit$nu,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$GED[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, nu = fit$nu)
    add_fit_metrics("GED", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("GED", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nGED fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["GED"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("SGED")) {
  t_model <- proc.time()
  # ----------------------------
  # skewed GED fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for skewed GED fit", nm))
      next
    }
    fit <- fit_sged(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s skewed GED fit failed", nm))
      next
    }
    moments <- sged_moments(fit$mu, fit$sigma, fit$nu, fit$kappa)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      shape = fit$nu,
      kappa = fit$kappa,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$SGED[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, nu = fit$nu, kappa = fit$kappa)
    add_fit_metrics("SGED", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("SGED", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nskewed GED fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["SGED"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("Sech")) {
  t_model <- proc.time()
  # ----------------------------
  # hyperbolic secant fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 2) {
      warning(sprintf("asset %s has too few finite returns for sech fit", nm))
      next
    }
    fit <- fit_sech(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s sech fit failed", nm))
      next
    }
    moments <- sech_moments(fit$mu, fit$sigma)
    metrics <- fit_metrics(fit$logLik, 2, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$Sech[[nm]] <- list(mu = fit$mu, sigma = fit$sigma)
    add_fit_metrics("Sech", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("Sech", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nhyperbolic secant fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["Sech"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("GSH")) {
  t_model <- proc.time()
  # ----------------------------
  # GSH fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for GSH fit", nm))
      next
    }
    fit <- fit_gsh(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s GSH fit failed", nm))
      next
    }
    moments <- gsh_moments(fit$mu, fit$sigma, fit$t)
    metrics <- fit_metrics(fit$logLik, 3, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      t = fit$t,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$GSH[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, t = fit$t)
    add_fit_metrics("GSH", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("GSH", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nGSH fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["GSH"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("SGSH")) {
  t_model <- proc.time()
  # ----------------------------
  # SGSH fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for SGSH fit", nm))
      next
    }
    fit <- fit_sgsh(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s SGSH fit failed", nm))
      next
    }
    moments <- sgsh_moments(fit$mu, fit$sigma, fit$t, fit$kappa)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      t = fit$t,
      kappa = fit$kappa,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$SGSH[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, t = fit$t, kappa = fit$kappa)
    add_fit_metrics("SGSH", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("SGSH", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nSGSH fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["SGSH"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("NEFGHS")) {
  t_model <- proc.time()
  # ----------------------------
  # NEF-GHS fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for NEF-GHS fit", nm))
      next
    }
    fit <- fit_nef_ghs(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s NEF-GHS fit failed", nm))
      next
    }
    moments <- nef_ghs_moments(fit$mu, fit$sigma, fit$lambda, fit$beta)
    metrics <- fit_metrics(fit$logLik, 4, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      lambda = fit$lambda,
      beta = fit$beta,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$NEFGHS[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, beta = fit$beta)
    add_fit_metrics("NEFGHS", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("NEFGHS", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nNEF-GHS fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["NEFGHS"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("T")) {
  t_model <- proc.time()
  # ----------------------------
  # symmetric t fits per asset (Azzalini, alpha=0)
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for t fit", nm))
      next
    }
    fixed_list <- list(alpha = 0)
    if (is.finite(fixed_nu)) fixed_list$nu <- fixed_nu
    fit <- quiet_selm(sn::selm(r ~ 1, family = "ST", fixed.param = fixed_list, control = list(trace = 0)))
    if (inherits(fit, "try-error")) {
      warning(sprintf("asset %s t fit failed", nm))
      next
    }
    dp <- fit@param$dp
    nu <- dp[4]
    if (!is.finite(nu) && is.finite(fixed_nu)) nu <- fixed_nu
    dp_mom <- c(dp[1], dp[2], 0, nu)
    moments <- st_moments(dp_mom, nu)
    mean_val <- moments$mean
    sd_val <- moments$sd
    skew_val <- 0
    kurt_val <- moments$kurtosis
    n_obs <- length(r)
    k <- length(dp)
    ll <- as.numeric(logLik(fit))
    metrics <- fit_metrics(ll, k, n_obs)
    nu_out <- if (is.finite(dp[4])) dp[4] else if (is.finite(fixed_nu)) fixed_nu else NA_real_
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      xi = dp[1],
      omega = dp[2],
      nu = nu_out,
      mean = mean_val,
      sd = sd_val,
      skew = skew_val,
      kurtosis = kurt_val,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    dp_full <- if (is.finite(nu_out)) {
      if (length(dp) == 2) {
        c(dp[1], dp[2], 0, nu_out)
      } else if (length(dp) == 3) {
        c(dp[1], dp[2], 0, nu_out)
      } else if (length(dp) >= 4) {
        dp[1:4]
      } else {
        NA_real_
      }
    } else {
      NA_real_
    }
    fit_params$T[[nm]] <- list(dp = dp_full)
    add_fit_metrics("T", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("T", nm, mean_val, sd_val, skew_val, kurt_val)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nsymmetric t fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["T"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("NCT")) {
  t_model <- proc.time()
  # ----------------------------
  # non-central t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for non-central t fit", nm))
      next
    }
    fit <- fit_nct(r, fixed_nu)
    if (is.null(fit)) {
      warning(sprintf("asset %s non-central t fit failed", nm))
      next
    }
    moments <- nct_moments(fit$mu, fit$sigma, fit$nu, fit$ncp)
    k <- if (is.finite(fixed_nu)) 3 else 4
    metrics <- fit_metrics(fit$logLik, k, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      nu = fit$nu,
      ncp = fit$ncp,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$NCT[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, nu = fit$nu, ncp = fit$ncp)
    add_fit_metrics("NCT", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("NCT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nnon-central t fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["NCT"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("SkewNormal")) {
  t_model <- proc.time()
  # ----------------------------
  # Azzalini skew-normal fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for skew-normal fit", nm))
      next
    }
    fit <- quiet_selm(sn::selm(r ~ 1, family = "SN", control = list(trace = 0)))
    dp <- fit@param$dp
    cp <- sn::dp2cp(dp, family = "SN")
    kurt <- cp[4]
    if (is.na(kurt)) {
      kurt <- sn_excess_kurtosis(dp[3])
    }
    n_obs <- length(r)
    k <- length(dp)
    ll <- as.numeric(logLik(fit))
    metrics <- fit_metrics(ll, k, n_obs)
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      xi = dp[1],
      omega = dp[2],
      alpha = dp[3],
      mean = cp[1],
      sd = cp[2],
      skew = cp[3],
      kurtosis = kurt,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$SkewNormal[[nm]] <- list(dp = dp)
    add_fit_metrics("SkewNormal", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("SkewNormal", nm, cp[1], cp[2], cp[3], kurt)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nskew-normal fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["SkewNormal"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("SkewT")) {
  t_model <- proc.time()
  # ----------------------------
  # Azzalini skew-t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for skew-t fit", nm))
      next
    }
    fixed_list <- NULL
    if (is.finite(fixed_nu)) fixed_list <- list(nu = fixed_nu)
    fit <- quiet_selm(sn::selm(r ~ 1, family = "ST", fixed.param = fixed_list, control = list(trace = 0)))
    if (inherits(fit, "try-error")) {
      warning(sprintf("asset %s skew-t fit failed", nm))
      next
    }
    dp <- fit@param$dp
    nu <- dp[4]
    if (!is.finite(nu) && is.finite(fixed_nu)) nu <- fixed_nu
    moments <- st_moments(dp, nu)
    n_obs <- length(r)
    k <- length(dp)
    ll <- as.numeric(logLik(fit))
    metrics <- fit_metrics(ll, k, n_obs)
    nu_out <- if (is.finite(dp[4])) dp[4] else if (is.finite(fixed_nu)) fixed_nu else NA_real_
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      xi = dp[1],
      omega = dp[2],
      alpha = dp[3],
      nu = nu_out,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = ll,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    dp_full <- if (length(dp) < 4 && is.finite(nu_out)) c(dp, nu_out) else dp
    fit_params$SkewT[[nm]] <- list(dp = dp_full)
    add_fit_metrics("SkewT", nm, ll, metrics$AIC, metrics$BIC)
    add_fit_moments("SkewT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nskew-t fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["SkewT"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("FSSkewNormal")) {
  t_model <- proc.time()
  # ----------------------------
  # FS skew-normal fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 3) {
      warning(sprintf("asset %s has too few finite returns for skew-normal fit", nm))
      next
    }
    fit <- fit_fs_skew_normal(r)
    if (is.null(fit)) {
      warning(sprintf("asset %s skew-normal fit failed", nm))
      next
    }
    moments <- fs_skew_normal_moments(fit$mu, fit$sigma, fit$gamma)
    k <- 3
    metrics <- fit_metrics(fit$logLik, k, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      gamma = fit$gamma,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$FSSkewNormal[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, gamma = fit$gamma)
    add_fit_metrics("FSSkewNormal", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("FSSkewNormal", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nFS skew-normal fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["FSSkewNormal"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("FSSkewT")) {
  t_model <- proc.time()
  # ----------------------------
  # FS skew-t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for skew-t fit", nm))
      next
    }
    fit <- fit_fs_skew_t(r, fixed_nu)
    if (is.null(fit)) {
      warning(sprintf("asset %s skew-t fit failed", nm))
      next
    }
    moments <- fs_skew_t_moments(fit$mu, fit$sigma, fit$gamma, fit$nu)
    k <- if (is.finite(fixed_nu)) 3 else 4
    metrics <- fit_metrics(fit$logLik, k, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      gamma = fit$gamma,
      nu = fit$nu,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$FSSkewT[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, gamma = fit$gamma, nu = fit$nu)
    add_fit_metrics("FSSkewT", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("FSSkewT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nFS skew-t fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["FSSkewT"]] <- (proc.time() - t_model)[["elapsed"]]
}

if (is_enabled("JFSkewT")) {
  t_model <- proc.time()
  # ----------------------------
  # Jones-Faddy skew-t fits per asset
  # ----------------------------
  fit_rows <- list()
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 5) {
      warning(sprintf("asset %s has too few finite returns for JF skew-t fit", nm))
      next
    }
    fit <- fit_jf_skew_t(r, fixed_nu)
    if (is.null(fit)) {
      warning(sprintf("asset %s JF skew-t fit failed", nm))
      next
    }
    moments <- jf_skew_t_moments(fit$mu, fit$sigma, fit$a, fit$b)
    k <- 4
    metrics <- fit_metrics(fit$logLik, k, length(r))
    fit_rows[[nm]] <- data.frame(
      asset = nm,
      mu = fit$mu,
      sigma = fit$sigma,
      a = fit$a,
      b = fit$b,
      mean = moments$mean,
      sd = moments$sd,
      skew = moments$skew,
      kurtosis = moments$kurtosis,
      logLik = fit$logLik,
      AIC = metrics$AIC,
      BIC = metrics$BIC,
      stringsAsFactors = FALSE
    )
    fit_params$JFSkewT[[nm]] <- list(mu = fit$mu, sigma = fit$sigma, a = fit$a, b = fit$b)
    add_fit_metrics("JFSkewT", nm, fit$logLik, metrics$AIC, metrics$BIC)
    add_fit_moments("JFSkewT", nm, moments$mean, moments$sd, moments$skew, moments$kurtosis)
  }
  fit_tab <- do.call(rbind, fit_rows)
  rownames(fit_tab) <- NULL
  cat("\nJones-Faddy skew-t fit table\n")
  fit_tab_print <- format_table(fit_tab, digits = 4)
  print(fit_tab_print, row.names = FALSE)
  model_times[["JFSkewT"]] <- (proc.time() - t_model)[["elapsed"]]
}

# ----------------------------
# density plots by asset
# ----------------------------
if (isTRUE(do_density_plots)) {
  logcondens_ok <- requireNamespace("logcondens", quietly = TRUE)
  if (isTRUE(do_logcondens_plot) && !logcondens_ok) {
    warning("package 'logcondens' not available; skipping logcondens density")
  }
  logspline_ok <- requireNamespace("logspline", quietly = TRUE)
  if (isTRUE(do_logspline_plot) && !logspline_ok) {
    warning("package 'logspline' not available; skipping logspline density")
  }
  for (nm in asset_names) {
    r <- ret_list[[nm]]
    r <- r[is.finite(r)]
    if (length(r) < 2) next

    x_grid <- density_x_grid(r, lower = 0.01, upper = 0.99, pad_frac = 0.05, n = 400)

    dens_list <- list()
    if (!is.null(fit_params$Normal[[nm]])) {
      p <- fit_params$Normal[[nm]]
      dens_list$Normal <- dnorm(x_grid, mean = p$mean, sd = p$sd)
    }
    if (!is.null(fit_params$Logistic[[nm]])) {
      p <- fit_params$Logistic[[nm]]
      dens_list$Logistic <- exp(logistic_logpdf(x_grid, p$mu, p$s))
    }
    if (!is.null(fit_params$EGB2[[nm]])) {
      p <- fit_params$EGB2[[nm]]
      dens_list$EGB2 <- exp(egb2_logpdf(x_grid, p$mu, p$sigma, p$p, p$q))
    }
    if (!is.null(fit_params$NIG[[nm]])) {
      p <- fit_params$NIG[[nm]]
      dens_list$NIG <- exp(nig_logpdf(x_grid, p$alpha, p$beta, p$delta, p$mu))
    }
    if (!is.null(fit_params$Hyperbolic[[nm]])) {
      p <- fit_params$Hyperbolic[[nm]]
      dens_list$Hyperbolic <- exp(hyperbolic_logpdf(x_grid, p$mu, p$delta, p$alpha, p$beta))
    }
    if (!is.null(fit_params$Champernowne[[nm]])) {
      p <- fit_params$Champernowne[[nm]]
      dens_list$Champernowne <- exp(champernowne_logpdf(x_grid, p$mu, p$sigma, p$lambda))
    }
    if (!is.null(fit_params$NormalLaplace[[nm]])) {
      p <- fit_params$NormalLaplace[[nm]]
      dens_list$NormalLaplace <- exp(normal_laplace_logpdf(x_grid, p$mu, p$sigma, p$b))
    }
    if (!is.null(fit_params$GT[[nm]])) {
      p <- fit_params$GT[[nm]]
      dens_list$GT <- exp(gt_logpdf(x_grid, p$mu, p$sigma, p$p, p$q))
    }
    if (!is.null(fit_params$Cauchy[[nm]])) {
      p <- fit_params$Cauchy[[nm]]
      dens_list$Cauchy <- exp(cauchy_logpdf(x_grid, p$x0, p$gamma))
    }
    if (!is.null(fit_params$SGT[[nm]])) {
      p <- fit_params$SGT[[nm]]
      dens_list$SGT <- exp(sgt_logpdf(x_grid, p$mu, p$sigma, p$p, p$q, p$lambda))
    }
    if (!is.null(fit_params$GH[[nm]])) {
      p <- fit_params$GH[[nm]]
      dens_list$GH <- exp(gh_logpdf(x_grid, p$mu, p$delta, p$alpha, p$beta, p$lambda))
    }
    if (!is.null(fit_params$VG[[nm]])) {
      p <- fit_params$VG[[nm]]
      dens_list$VG <- exp(vg_logpdf(x_grid, p$mu, p$sigma, p$theta, p$nu))
    }
    if (!is.null(fit_params$T[[nm]])) {
      p <- fit_params$T[[nm]]
      if (length(p$dp) >= 4 && is.finite(p$dp[4])) {
        dens_list$T <- sn::dst(x_grid, dp = p$dp)
      }
    }
    if (!is.null(fit_params$NCT[[nm]])) {
      p <- fit_params$NCT[[nm]]
      dens_list$NCT <- exp(nct_logpdf(x_grid, p$mu, p$sigma, p$nu, p$ncp))
    }
    if (!is.null(fit_params$Laplace[[nm]])) {
      p <- fit_params$Laplace[[nm]]
      dens_list$Laplace <- exp(laplace_logpdf(x_grid, p$mu, p$b))
    }
    if (!is.null(fit_params$ALaplace[[nm]])) {
      p <- fit_params$ALaplace[[nm]]
      dens_list$ALaplace <- exp(ald_logpdf(x_grid, p$mu, p$b, p$kappa))
    }
    if (!is.null(fit_params$GED[[nm]])) {
      p <- fit_params$GED[[nm]]
      dens_list$GED <- exp(ged_logpdf(x_grid, p$mu, p$sigma, p$nu))
    }
    if (!is.null(fit_params$SGED[[nm]])) {
      p <- fit_params$SGED[[nm]]
      dens_list$SGED <- exp(sged_logpdf(x_grid, p$mu, p$sigma, p$nu, p$kappa))
    }
    if (!is.null(fit_params$Sech[[nm]])) {
      p <- fit_params$Sech[[nm]]
      dens_list$Sech <- exp(sech_logpdf(x_grid, p$mu, p$sigma))
    }
    if (!is.null(fit_params$GSH[[nm]])) {
      p <- fit_params$GSH[[nm]]
      dens_list$GSH <- exp(gsh_logpdf(x_grid, p$mu, p$sigma, p$t))
    }
    if (!is.null(fit_params$SGSH[[nm]])) {
      p <- fit_params$SGSH[[nm]]
      dens_list$SGSH <- exp(sgsh_logpdf(x_grid, p$mu, p$sigma, p$t, p$kappa))
    }
    if (!is.null(fit_params$NEFGHS[[nm]])) {
      p <- fit_params$NEFGHS[[nm]]
      dens_list$NEFGHS <- exp(nef_ghs_logpdf(x_grid, p$mu, p$sigma, p$lambda, p$beta))
    }
    if (!is.null(fit_params$SkewNormal[[nm]])) {
      p <- fit_params$SkewNormal[[nm]]
      dens_list$SkewNormal <- sn::dsn(x_grid, dp = p$dp)
    }
    if (!is.null(fit_params$SkewT[[nm]])) {
      p <- fit_params$SkewT[[nm]]
      if (length(p$dp) >= 4 && is.finite(p$dp[4])) {
        dens_list$SkewT <- sn::dst(x_grid, dp = p$dp)
      }
    }
    if (!is.null(fit_params$FSSkewNormal[[nm]])) {
      p <- fit_params$FSSkewNormal[[nm]]
      dens_list$FSSkewNormal <- exp(fs_skew_normal_logpdf(x_grid, p$mu, p$sigma, p$gamma))
    }
    if (!is.null(fit_params$FSSkewT[[nm]])) {
      p <- fit_params$FSSkewT[[nm]]
      dens_list$FSSkewT <- exp(fs_skew_t_logpdf(x_grid, p$mu, p$sigma, p$gamma, p$nu))
    }
    if (!is.null(fit_params$JFSkewT[[nm]])) {
      p <- fit_params$JFSkewT[[nm]]
      dens_list$JFSkewT <- exp(jf_skew_t_logpdf(x_grid, p$mu, p$sigma, p$a, p$b))
    }
    if (length(dens_list) == 0) next

    if (isTRUE(do_kde_plot)) {
      kde <- density(r, n = 512)
      dens_list$KDE <- approx(kde$x, kde$y, xout = x_grid, rule = 2)$y
    }
    if (isTRUE(do_logcondens_plot) && logcondens_ok) {
      lc <- try(logcondens::logConDens(r), silent = TRUE)
      if (!inherits(lc, "try-error") && is.list(lc) && all(c("x", "phi") %in% names(lc))) {
        dens_list$LogCondens <- approx(lc$x, exp(lc$phi), xout = x_grid, rule = 2)$y
      }
    }
    if (isTRUE(do_logspline_plot) && logspline_ok) {
      ls <- try(logspline::logspline(r), silent = TRUE)
      if (inherits(ls, "try-error")) {
        warning(sprintf("logspline fit failed for asset %s", nm))
      } else {
        dens_list$LogSpline <- logspline::dlogspline(x_grid, ls)
      }
    }

    outfile <- sprintf("%s.png", nm)
    plot_density_fits(x_grid, dens_list, outfile, sprintf("%s density fits", nm))
  }
}

# ----------------------------
# AIC/BIC summary tables
# ----------------------------
if (length(aic_rows) > 0) {
  aic_df <- do.call(rbind, aic_rows)
  model_levels <- c(
    "Normal",
    "Logistic",
    "EGB2",
    "NIG",
    "Hyperbolic",
    "Champernowne",
    "NormalLaplace",
    "GT",
    "SGT",
    "Cauchy",
    "GH",
    "VG",
    "Laplace",
    "ALaplace",
    "GED",
    "SGED",
    "Sech",
    "GSH",
    "SGSH",
    "NEFGHS",
    "T",
    "NCT",
    "SkewNormal",
    "SkewT",
    "FSSkewNormal",
    "FSSkewT",
    "JFSkewT"
  )
  aic_df$model <- factor(aic_df$model, levels = model_levels)
  aic_df$model_chr <- as.character(aic_df$model)
  aic_tab <- reshape(aic_df[, c("asset", "model", "AIC")],
                     idvar = "asset", timevar = "model", direction = "wide")
  bic_tab <- reshape(aic_df[, c("asset", "model", "BIC")],
                     idvar = "asset", timevar = "model", direction = "wide")
  best_aic <- by(aic_df, aic_df$asset, function(d) d$model_chr[which.min(d$AIC)])
  best_bic <- by(aic_df, aic_df$asset, function(d) d$model_chr[which.min(d$BIC)])
  aic_tab$best_model <- as.character(best_aic[aic_tab$asset])
  bic_tab$best_model <- as.character(best_bic[bic_tab$asset])
  aic_tab <- aic_tab[, c("asset", "best_model", setdiff(names(aic_tab), c("asset", "best_model")))]
  bic_tab <- bic_tab[, c("asset", "best_model", setdiff(names(bic_tab), c("asset", "best_model")))]
  aic_tab[is.na(aic_tab)] <- NA_real_
  bic_tab[is.na(bic_tab)] <- NA_real_

  cat("\nAIC summary table\n")
  aic_tab_print <- format_table(aic_tab, digits = 4)
  print(aic_tab_print, row.names = FALSE)

  cat("\nBIC summary table\n")
  bic_tab_print <- format_table(bic_tab, digits = 4)
  print(bic_tab_print, row.names = FALSE)

  for (nm in unique(aic_df$asset)) {
    d <- aic_df[aic_df$asset == nm & is.finite(aic_df$logLik), ]
    if (nrow(d) == 0) next

    ord_ll <- d[order(-d$logLik), ]
    ord_aic <- d[order(d$AIC), ]
    ord_bic <- d[order(d$BIC), ]

    rank_tab <- data.frame(
      logLik_model = ord_ll$model_chr,
      AIC_model = ord_aic$model_chr,
      BIC_model = ord_bic$model_chr,
      dLogLik = ord_ll$logLik[1] - ord_ll$logLik,
      dAIC = ord_aic$AIC - ord_aic$AIC[1],
      dBIC = ord_bic$BIC - ord_bic$BIC[1],
      stringsAsFactors = FALSE
    )
    cat(sprintf("\nmodel ranking table: %s\n", nm))
    rank_tab_print <- format_table(rank_tab, digits = 4)
    print(rank_tab_print, row.names = FALSE)
  }

  delta_rows <- lapply(split(aic_df, aic_df$asset), function(d) {
    d <- d[is.finite(d$logLik), ]
    if (nrow(d) == 0) return(NULL)
    best_ll <- max(d$logLik, na.rm = TRUE)
    best_aic <- min(d$AIC, na.rm = TRUE)
    best_bic <- min(d$BIC, na.rm = TRUE)
    data.frame(
      model = d$model_chr,
      dLogLik = best_ll - d$logLik,
      dAIC = d$AIC - best_aic,
      dBIC = d$BIC - best_bic,
      stringsAsFactors = FALSE
    )
  })
  delta_df <- do.call(rbind, delta_rows)
  if (!is.null(delta_df) && nrow(delta_df) > 0) {
    avg_tab <- aggregate(
      delta_df[, c("dLogLik", "dAIC", "dBIC")],
      by = list(model = delta_df$model),
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    moments_df <- if (length(fit_moments) > 0) {
      do.call(rbind, fit_moments)
    } else {
      data.frame(model = character(), mean = numeric(), sd = numeric(), skew = numeric(), kurtosis = numeric(), stringsAsFactors = FALSE)
    }
    avg_moments <- if (nrow(moments_df) > 0) {
      aggregate(
        moments_df[, c("mean", "sd", "skew", "kurtosis")],
        by = list(model = moments_df$model),
        FUN = function(x) mean(x, na.rm = TRUE)
      )
    } else {
      data.frame(model = character(), mean = numeric(), sd = numeric(), skew = numeric(), kurtosis = numeric(), stringsAsFactors = FALSE)
    }
    meta <- model_meta()
    times_df <- if (length(model_times) > 0) {
      data.frame(model = names(model_times), time_sec = unlist(model_times), stringsAsFactors = FALSE)
    } else {
      data.frame(model = character(), time_sec = numeric(), stringsAsFactors = FALSE)
    }
    avg_tab <- merge(avg_tab, meta, by = "model", all.x = TRUE)
    avg_tab <- merge(avg_tab, avg_moments, by = "model", all.x = TRUE)
    avg_tab <- merge(avg_tab, times_df, by = "model", all.x = TRUE)
    avg_tab <- avg_tab[order(avg_tab$dAIC), ]
    keep_cols <- c("model", "dLogLik", "dAIC", "dBIC", "k", "mean", "sd", "skew", "kurtosis", "time_sec")
    for (col in keep_cols) {
      if (!col %in% names(avg_tab)) avg_tab[[col]] <- NA_real_
    }
    avg_tab <- avg_tab[, keep_cols]
    cat("\naverage model ranking table\n")
    avg_tab_print <- format_table(avg_tab, digits = 4, integer_cols = "k")
    print(avg_tab_print, row.names = FALSE)
  }
}

cat(sprintf("\nTotal time elapsed: %.2f seconds\n", (proc.time() - start_time)[["elapsed"]]))
