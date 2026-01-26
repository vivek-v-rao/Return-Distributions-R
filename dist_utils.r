# dist_utils.r
#
# Small utilities for distribution moments and model metrics.

#' Compute excess kurtosis for skew-normal from alpha.
sn_excess_kurtosis <- function(alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  a <- delta * sqrt(2 / pi)
  2 * (pi - 3) * (a^4) / (1 - 2 * delta^2 / pi)^2
}

#' Compute skew-t moments (mean/sd/skew/kurtosis) when they exist.
st_moments <- function(dp, nu) {
  mean_val <- NA_real_
  sd_val <- NA_real_
  skew_val <- NA_real_
  kurt_val <- NA_real_

  if (!is.finite(nu)) {
    return(list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val))
  }

  max_upto <- if (nu > 4) {
    4L
  } else if (nu > 3) {
    3L
  } else if (nu > 2) {
    2L
  } else if (nu > 1) {
    1L
  } else {
    0L
  }

  if (max_upto == 0) {
    return(list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val))
  }

  dp_full <- if (is.finite(nu)) {
    if (length(dp) == 2) {
      c(dp[1], dp[2], 0, nu)
    } else if (length(dp) == 3) {
      c(dp, nu)
    } else {
      dp
    }
  } else {
    dp
  }
  cp <- try(sn:::st.dp2cp(dp_full, cp.type = "proper", upto = max_upto), silent = TRUE)
  if (!inherits(cp, "try-error") && !is.null(cp)) {
    if ("mean" %in% names(cp)) mean_val <- cp[["mean"]]
    if ("s.d." %in% names(cp)) sd_val <- cp[["s.d."]]
    if ("gamma1" %in% names(cp)) skew_val <- cp[["gamma1"]]
    if ("gamma2" %in% names(cp)) kurt_val <- cp[["gamma2"]]
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Compute AIC and BIC from log-likelihood, parameter count, and sample size.
fit_metrics <- function(loglik, k, n) {
  list(
    AIC = 2 * k - 2 * loglik,
    BIC = -2 * loglik + log(n) * k
  )
}

#' Format numeric columns with fixed digits, preserving integer columns.
format_table <- function(df, digits = 4, integer_cols = NULL) {
  out <- df
  if (is.null(integer_cols)) integer_cols <- character()

  for (col in integer_cols) {
    if (col %in% names(df)) {
      out[[col]] <- sprintf("%d", as.integer(round(df[[col]])))
    }
  }

  num_cols <- names(df)[sapply(df, is.numeric)]
  num_cols <- setdiff(num_cols, integer_cols)
  fmt <- paste0("%.", digits, "f")
  out[num_cols] <- lapply(out[num_cols], function(x) sprintf(fmt, x))
  out
}

#' Stable softplus helper.
softplus <- function(z) {
  ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
}

#' Log-density for logistic distribution.
logistic_logpdf <- function(x, mu, s) {
  if (!is.finite(s) || s <= 0) return(rep(-Inf, length(x)))
  z <- (x - mu) / s
  -log(s) - z - 2 * log1p(exp(-z))
}

#' Moments for logistic distribution (excess kurtosis).
logistic_moments <- function(mu, s) {
  list(
    mean = mu,
    sd = (pi / sqrt(3)) * s,
    skew = 0,
    kurtosis = 1.2
  )
}

#' Log-density for EGB2 distribution (McDonald).
egb2_logpdf <- function(x, mu, sigma, p, q) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0) {
    return(rep(-Inf, length(x)))
  }
  z <- (x - mu) / sigma
  logc <- -log(sigma) - lbeta(p, q)
  logc + p * z - (p + q) * softplus(z)
}

#' Closed-form moments for EGB2 distribution (excess kurtosis).
egb2_moments <- function(mu, sigma, p, q) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  k1 <- digamma(p) - digamma(q)
  k2 <- psigamma(p, deriv = 1) + psigamma(q, deriv = 1)
  k3 <- psigamma(p, deriv = 2) - psigamma(q, deriv = 2)
  k4 <- psigamma(p, deriv = 3) + psigamma(q, deriv = 3)

  mean_val <- mu + sigma * k1
  sd_val <- if (is.finite(k2) && k2 > 0) sigma * sqrt(k2) else NA_real_
  skew_val <- if (is.finite(k2) && k2 > 0) k3 / k2^(3/2) else NA_real_
  kurt_val <- if (is.finite(k2) && k2 > 0) k4 / k2^2 else NA_real_

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Normal-Inverse-Gaussian distribution.
nig_logpdf <- function(x, alpha, beta, delta, mu) {
  if (!is.finite(alpha) || alpha <= 0 || !is.finite(delta) || delta <= 0 || !is.finite(beta) || abs(beta) >= alpha) {
    return(rep(-Inf, length(x)))
  }
  y <- x - mu
  gamma <- sqrt(alpha^2 - beta^2)
  z <- sqrt(delta^2 + y^2)
  logk <- log(besselK(alpha * z, 1, expon.scaled = TRUE)) - alpha * z
  logc <- log(alpha) + log(delta) - log(pi) + delta * gamma + beta * y - log(z)
  logc + logk
}

#' Numeric moments for Normal-Inverse-Gaussian distribution (excess kurtosis).
nig_moments <- function(alpha, beta, delta, mu) {
  dens <- function(x) exp(nig_logpdf(x, alpha, beta, delta, mu))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Hyperbolic distribution (GH with lambda = 1).
hyperbolic_logpdf <- function(x, mu, delta, alpha, beta) {
  if (!is.finite(delta) || delta <= 0 || !is.finite(alpha) || alpha <= 0 || !is.finite(beta) || abs(beta) >= alpha) {
    return(rep(-Inf, length(x)))
  }
  y <- x - mu
  gamma <- sqrt(alpha^2 - beta^2)
  z <- sqrt(delta^2 + y^2)
  k1 <- log(besselK(delta * gamma, 1, expon.scaled = TRUE)) - delta * gamma
  k05 <- log(besselK(alpha * z, 0.5, expon.scaled = TRUE)) - alpha * z
  logc <- 0.5 * log(alpha^2 - beta^2) -
    0.5 * log(alpha) -
    0.5 * log(2 * pi) -
    log(delta) -
    k1
  logc + 0.5 * log(z) + k05 + beta * y
}

#' Numeric moments for Hyperbolic distribution (excess kurtosis).
hyperbolic_moments <- function(mu, delta, alpha, beta) {
  dens <- function(x) exp(hyperbolic_logpdf(x, mu, delta, alpha, beta))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Generalized Hyperbolic distribution.
gh_logpdf <- function(x, mu, delta, alpha, beta, lambda) {
  if (!is.finite(delta) || delta <= 0 || !is.finite(alpha) || alpha <= 0 || !is.finite(beta) || abs(beta) >= alpha || !is.finite(lambda)) {
    return(rep(-Inf, length(x)))
  }
  y <- x - mu
  gamma <- sqrt(alpha^2 - beta^2)
  z <- sqrt(delta^2 + y^2)
  k_lambda <- log(besselK(delta * gamma, lambda, expon.scaled = TRUE)) - delta * gamma
  k_order <- log(besselK(alpha * z, lambda - 0.5, expon.scaled = TRUE)) - alpha * z
  logc <- 0.5 * lambda * log(alpha^2 - beta^2) -
    (lambda - 0.5) * log(alpha) -
    0.5 * log(2 * pi) -
    lambda * log(delta) -
    k_lambda
  logc + (lambda - 0.5) * log(z) + k_order + beta * y
}

#' Numeric moments for Generalized Hyperbolic distribution (excess kurtosis).
gh_moments <- function(mu, delta, alpha, beta, lambda) {
  dens <- function(x) exp(gh_logpdf(x, mu, delta, alpha, beta, lambda))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Variance-Gamma distribution.
vg_logpdf <- function(x, mu, sigma, theta, nu) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(nu) || nu <= 0) {
    return(rep(-Inf, length(x)))
  }
  y <- x - mu
  lambda <- 1 / nu
  alpha <- sqrt(theta^2 + 2 * sigma^2 / nu) / sigma^2
  beta <- theta / sigma^2
  if (!is.finite(alpha) || !is.finite(beta) || alpha <= 0 || abs(beta) >= alpha) {
    return(rep(-Inf, length(x)))
  }
  z <- pmax(abs(y), .Machine$double.eps)
  karg <- alpha * z
  order <- lambda - 0.5
  logk <- log(besselK(karg, order, expon.scaled = TRUE)) - karg
  logc <- lambda * log(alpha^2 - beta^2) -
    (lambda - 0.5) * log(2 * alpha) -
    0.5 * log(pi) -
    lgamma(lambda)
  logc + order * log(z) + beta * y + logk
}

#' Closed-form moments for Variance-Gamma distribution (excess kurtosis).
vg_moments <- function(mu, sigma, theta, nu) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(nu) || nu <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  mean_val <- mu + theta
  var_val <- sigma^2 + nu * theta^2
  sd_val <- if (is.finite(var_val) && var_val > 0) sqrt(var_val) else NA_real_
  skew_val <- NA_real_
  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- 2 * theta^3 * nu^2 + 3 * sigma^2 * theta * nu
    m4 <- 3 * (theta^4 * nu^2 * (1 + 2 * nu) + 2 * theta^2 * sigma^2 * (nu + 2 * nu^2) + sigma^4 * (1 + nu))
    skew_val <- m3 / sd_val^3
    kurt_val <- m4 / sd_val^4 - 3
  }
  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Generalized t distribution (McDonald-Newey).
gt_logpdf <- function(x, mu, sigma, p, q) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0) {
    return(rep(-Inf, length(x)))
  }
  z <- abs((x - mu) / sigma)
  logc <- log(p) - log(2 * sigma) - lbeta(1 / p, q)
  logc - (q + 1 / p) * log1p(z^p)
}

#' Moments for Generalized t distribution (excess kurtosis).
gt_moments <- function(mu, sigma, p, q) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  base <- beta(1 / p, q)
  m2 <- NA_real_
  m4 <- NA_real_
  if (q > 2 / p) {
    m2 <- beta(3 / p, q - 2 / p) / base
  }
  if (q > 4 / p) {
    m4 <- beta(5 / p, q - 4 / p) / base
  }
  var_val <- if (is.finite(m2)) sigma^2 * m2 else NA_real_
  sd_val <- if (is.finite(var_val) && var_val > 0) sqrt(var_val) else NA_real_
  kurt_val <- if (is.finite(m4) && is.finite(m2) && m2 > 0) m4 / (m2^2) - 3 else NA_real_
  list(mean = mu, sd = sd_val, skew = 0, kurtosis = kurt_val)
}

#' Log-density for Cauchy distribution.
cauchy_logpdf <- function(x, x0, gamma) {
  if (!is.finite(gamma) || gamma <= 0) return(rep(-Inf, length(x)))
  dcauchy(x, location = x0, scale = gamma, log = TRUE)
}

#' Moments for Cauchy distribution (undefined).
cauchy_moments <- function() {
  list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_)
}

#' Log-density for Skewed Generalized t distribution (Theodossiou).
sgt_logpdf <- function(x, mu, sigma, p, q, lambda) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(p) || p <= 0 || !is.finite(q) || q <= 0 || !is.finite(lambda) || abs(lambda) >= 1) {
    return(rep(-Inf, length(x)))
  }
  z <- (x - mu) / sigma
  c <- p / (2 * beta(1 / p, q))
  s <- 1 + lambda * sign(z)
  if (any(!is.finite(s)) || any(s <= 0)) return(rep(-Inf, length(x)))
  logc <- log(c) - log(sigma)
  logc - (q + 1 / p) * log1p((abs(z) / s)^p)
}

#' Numeric moments for Skewed Generalized t (excess kurtosis).
sgt_moments <- function(mu, sigma, p, q, lambda) {
  dens <- function(x) exp(sgt_logpdf(x, mu, sigma, p, q, lambda))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Fernandez-Steel skew-t distribution.
fs_skew_t_logpdf <- function(x, mu, sigma, gamma, nu) {
  z <- (x - mu) / sigma
  adj <- ifelse(z >= 0, z / gamma, z * gamma)
  logc <- log(2) - log(gamma + 1 / gamma) - log(sigma)
  logc + dt(adj, df = nu, log = TRUE)
}

#' Numeric moments for Fernandez-Steel skew-t (when defined).
fs_skew_t_moments <- function(mu, sigma, gamma, nu) {
  mean_val <- NA_real_
  sd_val <- NA_real_
  skew_val <- NA_real_
  kurt_val <- NA_real_

  if (!is.finite(nu) || nu <= 1) {
    return(list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val))
  }

  dens <- function(x) exp(fs_skew_t_logpdf(x, mu, sigma, gamma, nu))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  if (is.finite(mean_val) && nu > 2) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  if (is.finite(sd_val) && nu > 3) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3) && sd_val > 0) {
      skew_val <- m3 / sd_val^3
    }
  }

  if (is.finite(sd_val) && nu > 4) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4) && sd_val > 0) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Fernandez-Steel skew-normal distribution.
fs_skew_normal_logpdf <- function(x, mu, sigma, gamma) {
  z <- (x - mu) / sigma
  adj <- ifelse(z >= 0, z / gamma, z * gamma)
  logc <- log(2) - log(gamma + 1 / gamma) - log(sigma)
  logc + dnorm(adj, log = TRUE)
}

#' Numeric moments for Fernandez-Steel skew-normal.
fs_skew_normal_moments <- function(mu, sigma, gamma) {
  dens <- function(x) exp(fs_skew_normal_logpdf(x, mu, sigma, gamma))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Laplace (double exponential).
laplace_logpdf <- function(x, mu, b) {
  if (!is.finite(b) || b <= 0) return(rep(-Inf, length(x)))
  -log(2 * b) - abs(x - mu) / b
}

#' Moments for Laplace (excess kurtosis).
laplace_moments <- function(mu, b) {
  list(
    mean = mu,
    sd = sqrt(2) * b,
    skew = 0,
    kurtosis = 3
  )
}

#' Log-density for asymmetric Laplace.
ald_logpdf <- function(x, mu, b, kappa) {
  if (!is.finite(b) || b <= 0 || !is.finite(kappa) || kappa <= 0) {
    return(rep(-Inf, length(x)))
  }
  b_pos <- b / kappa
  b_neg <- b * kappa
  logc <- -log(b_pos + b_neg)
  ifelse(
    x >= mu,
    logc - (x - mu) / b_pos,
    logc + (x - mu) / b_neg
  )
}

#' Moments for asymmetric Laplace (excess kurtosis).
ald_moments <- function(mu, b, kappa) {
  if (!is.finite(b) || b <= 0 || !is.finite(kappa) || kappa <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  b_pos <- b / kappa
  b_neg <- b * kappa
  mean_val <- mu + b_pos - b_neg
  var_val <- b_pos^2 + b_neg^2
  sd_val <- sqrt(var_val)
  m3 <- 2 * (b_pos^3 - b_neg^3)
  m4 <- 9 * (b_pos^4 + b_neg^4)
  skew_val <- if (sd_val > 0) m3 / sd_val^3 else NA_real_
  kurt_val <- if (sd_val > 0) m4 / sd_val^4 - 3 else NA_real_
  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for generalized error distribution (GED).
ged_logpdf <- function(x, mu, sigma, nu) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(nu) || nu <= 0) {
    return(rep(-Inf, length(x)))
  }
  z <- (x - mu) / sigma
  lambda <- sqrt(gamma(1 / nu) / gamma(3 / nu))
  logc <- log(nu) - log(2) - log(lambda) - lgamma(1 / nu) - log(sigma)
  logc - (abs(z) / lambda)^nu
}

#' Moments for GED (excess kurtosis).
ged_moments <- function(mu, sigma, nu) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(nu) || nu <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  kurt <- gamma(5 / nu) * gamma(1 / nu) / (gamma(3 / nu)^2) - 3
  list(mean = mu, sd = sigma, skew = 0, kurtosis = kurt)
}

#' Log-density for skewed GED (Fernandez-Steel skewing).
sged_logpdf <- function(x, mu, sigma, nu, kappa) {
  if (!is.finite(kappa) || kappa <= 0) return(rep(-Inf, length(x)))
  z <- (x - mu) / sigma
  adj <- ifelse(z >= 0, z / kappa, z * kappa)
  lambda <- sqrt(gamma(1 / nu) / gamma(3 / nu))
  logc0 <- log(nu) - log(2) - log(lambda) - lgamma(1 / nu)
  logc <- log(2) - log(kappa + 1 / kappa) - log(sigma)
  logc + logc0 - (abs(adj) / lambda)^nu
}

#' Numeric moments for skewed GED (excess kurtosis).
sged_moments <- function(mu, sigma, nu, kappa) {
  dens <- function(x) exp(sged_logpdf(x, mu, sigma, nu, kappa))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for hyperbolic secant distribution.
sech_logpdf <- function(x, mu, sigma) {
  if (!is.finite(sigma) || sigma <= 0) return(rep(-Inf, length(x)))
  z <- (x - mu) / sigma
  -log(2 * sigma) + log(1 / cosh(pi * z / 2))
}

#' Normalizing integral for Champernowne distribution.
champernowne_I <- function(lambda) {
  if (!is.finite(lambda) || lambda < 0) return(NA_real_)
  if (abs(lambda - 2) < 1e-8) return(1)
  if (lambda < 2) {
    d <- sqrt(4 - lambda^2)
    return((2 / d) * (pi / 2 - atan(lambda / d)))
  }
  d <- sqrt(lambda^2 - 4)
  (1 / d) * log((lambda + d) / (lambda - d))
}

#' Log-density for Champernowne distribution.
#' Parameterization matches hyperbolic secant when lambda = 0.
champernowne_logpdf <- function(x, mu, sigma, lambda) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(lambda) || lambda < 0) {
    return(rep(-Inf, length(x)))
  }
  a <- pi / (2 * sigma)
  I <- champernowne_I(lambda)
  if (!is.finite(I) || I <= 0) return(rep(-Inf, length(x)))
  logc <- log(a) - log(I)
  ax <- a * (x - mu)
  if (lambda > 0) {
    m <- pmax(ax, -ax, log(lambda))
    log_denom <- m + log(exp(ax - m) + exp(-ax - m) + exp(log(lambda) - m))
  } else {
    m <- pmax(ax, -ax)
    log_denom <- m + log(exp(ax - m) + exp(-ax - m))
  }
  logc - log_denom
}

#' Numeric moments for Champernowne distribution (excess kurtosis).
champernowne_moments <- function(mu, sigma, lambda) {
  dens <- function(x) exp(champernowne_logpdf(x, mu, sigma, lambda))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Normal-Laplace distribution.
normal_laplace_logpdf <- function(x, mu, sigma, b) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(b) || b <= 0) {
    return(rep(-Inf, length(x)))
  }
  z <- x - mu
  logc <- -log(2 * b) + (sigma^2) / (2 * b^2)
  t1 <- z / b + pnorm(-z / sigma - sigma / b, log.p = TRUE)
  t2 <- -z / b + pnorm(z / sigma - sigma / b, log.p = TRUE)
  m <- pmax(t1, t2)
  logc + m + log(exp(t1 - m) + exp(t2 - m))
}

#' Moments for Normal-Laplace distribution (excess kurtosis).
normal_laplace_moments <- function(mu, sigma, b) {
  if (!is.finite(sigma) || sigma <= 0 || !is.finite(b) || b <= 0) {
    return(list(mean = NA_real_, sd = NA_real_, skew = NA_real_, kurtosis = NA_real_))
  }
  var_val <- sigma^2 + 2 * b^2
  mean_val <- mu
  sd_val <- sqrt(var_val)
  skew_val <- 0
  kurt_val <- if (var_val > 0) 12 * b^4 / var_val^2 else NA_real_
  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Moments for hyperbolic secant distribution (excess kurtosis).
sech_moments <- function(mu, sigma) {
  list(
    mean = mu,
    sd = sigma,
    skew = 0,
    kurtosis = 2
  )
}

#' Compute parameters for the GSH distribution from t.
gsh_params <- function(t) {
  if (!is.finite(t)) return(NULL)
  if (abs(t) < 1e-8) {
    c2 <- pi / sqrt(3)
    return(list(a = 1, c1 = c2, c2 = c2))
  }
  if (t <= -pi) return(NULL)
  if (t < 0) {
    a <- cos(t)
    c2 <- sqrt((pi^2 - t^2) / 3)
    c1 <- (sin(t) / t) * c2
    return(list(a = a, c1 = c1, c2 = c2))
  }
  a <- cosh(t)
  c2 <- sqrt((pi^2 + t^2) / 3)
  c1 <- (sinh(t) / t) * c2
  list(a = a, c1 = c1, c2 = c2)
}

#' Log-density for generalized secant hyperbolic (GSH) distribution.
gsh_logpdf <- function(x, mu, sigma, t) {
  if (!is.finite(sigma) || sigma <= 0) return(rep(-Inf, length(x)))
  pars <- gsh_params(t)
  if (is.null(pars)) return(rep(-Inf, length(x)))
  z <- (x - mu) / sigma
  u <- pars$c2 * z
  m <- pmax(2 * u, 0)
  logden <- m + log(exp(2 * u - m) + 2 * pars$a * exp(u - m) + exp(-m))
  log(pars$c1) + u - logden - log(sigma)
}

#' Log-density for skew generalized secant hyperbolic (SGSH) distribution.
sgsh_logpdf <- function(x, mu, sigma, t, kappa) {
  if (!is.finite(kappa) || kappa <= 0) return(rep(-Inf, length(x)))
  z <- (x - mu) / sigma
  w <- 2 / (kappa + 1 / kappa)
  logw <- log(w) - log(sigma)
  out <- numeric(length(z))
  neg <- z < 0
  out[neg] <- logw + gsh_logpdf(z[neg] / kappa, 0, 1, t)
  out[!neg] <- logw + gsh_logpdf(z[!neg] * kappa, 0, 1, t)
  out
}

#' Numeric moments for GSH distribution (excess kurtosis).
gsh_moments <- function(mu, sigma, t) {
  dens <- function(x) exp(gsh_logpdf(x, mu, sigma, t))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Numeric moments for SGSH distribution (excess kurtosis).
sgsh_moments <- function(mu, sigma, t, kappa) {
  dens <- function(x) exp(sgsh_logpdf(x, mu, sigma, t, kappa))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Complex log-gamma via Lanczos approximation.
loggamma_complex <- function(z) {
  g <- 7
  p <- c(
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  )
  z <- z - 1
  x <- p[1]
  for (i in 2:length(p)) {
    x <- x + p[i] / (z + i - 1)
  }
  t <- z + g + 0.5
  0.5 * log(2 * pi) + (z + 0.5) * log(t) - t + log(x)
}

#' Log-density for generalized hyperbolic secant (GHS) distribution.
ghs_logpdf <- function(x, lambda) {
  if (!is.finite(lambda) || lambda <= 0) return(rep(-Inf, length(x)))
  z <- lambda + 1i * x / 2
  logc <- (2 * lambda - 1) * log(2) - log(2 * pi) - lgamma(2 * lambda)
  logc + 2 * Re(loggamma_complex(z))
}

#' Log-density for NEF-GHS distribution.
nef_ghs_logpdf <- function(x, mu, sigma, lambda, beta) {
  if (!is.finite(sigma) || sigma <= 0) return(rep(-Inf, length(x)))
  if (!is.finite(lambda) || lambda <= 0) return(rep(-Inf, length(x)))
  if (!is.finite(beta)) return(rep(-Inf, length(x)))
  z <- (x - mu) / sigma
  theta <- atan(beta)
  log_base <- ghs_logpdf(z, lambda) - log(sigma)
  log_skew <- theta * z - lambda * log1p(beta^2)
  log_base + log_skew
}

#' Numeric moments for NEF-GHS distribution (excess kurtosis).
nef_ghs_moments <- function(mu, sigma, lambda, beta) {
  dens <- function(x) exp(nef_ghs_logpdf(x, mu, sigma, lambda, beta))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for non-central t with location and scale.
nct_logpdf <- function(x, mu, sigma, nu, ncp) {
  z <- (x - mu) / sigma
  dt(z, df = nu, ncp = ncp, log = TRUE) - log(sigma)
}

#' Numeric moments for non-central t with location and scale.
nct_moments <- function(mu, sigma, nu, ncp) {
  dens <- function(x) exp(nct_logpdf(x, mu, sigma, nu, ncp))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  sd_val <- NA_real_
  if (is.finite(mean_val)) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  skew_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3)) {
      skew_val <- m3 / sd_val^3
    }
  }

  kurt_val <- NA_real_
  if (is.finite(sd_val) && sd_val > 0) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4)) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}

#' Log-density for Jones-Faddy skew-t with location and scale.
jf_skew_t_logpdf <- function(x, mu, sigma, a, b) {
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0 || !is.finite(sigma) || sigma <= 0) {
    return(rep(-Inf, length(x)))
  }
  z <- (x - mu) / sigma
  s <- sqrt(a + b + z^2)
  logc <- -(a + b - 1) * log(2) - lbeta(a, b) - 0.5 * log(a + b)
  logc + (a + 0.5) * log1p(z / s) + (b + 0.5) * log1p(-z / s) - log(sigma)
}

#' Numeric moments for Jones-Faddy skew-t with location and scale.
jf_skew_t_moments <- function(mu, sigma, a, b) {
  df <- a + b
  mean_val <- NA_real_
  sd_val <- NA_real_
  skew_val <- NA_real_
  kurt_val <- NA_real_

  if (!is.finite(df) || df <= 1) {
    return(list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val))
  }

  dens <- function(x) exp(jf_skew_t_logpdf(x, mu, sigma, a, b))
  mean_val <- try(integrate(function(x) x * dens(x), -Inf, Inf)$value, silent = TRUE)
  if (inherits(mean_val, "try-error")) mean_val <- NA_real_

  if (is.finite(mean_val) && df > 2) {
    var_val <- try(integrate(function(x) (x - mean_val)^2 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(var_val, "try-error") && is.finite(var_val)) {
      sd_val <- sqrt(var_val)
    }
  }

  if (is.finite(sd_val) && df > 3) {
    m3 <- try(integrate(function(x) (x - mean_val)^3 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m3, "try-error") && is.finite(m3) && sd_val > 0) {
      skew_val <- m3 / sd_val^3
    }
  }

  if (is.finite(sd_val) && df > 4) {
    m4 <- try(integrate(function(x) (x - mean_val)^4 * dens(x), -Inf, Inf)$value, silent = TRUE)
    if (!inherits(m4, "try-error") && is.finite(m4) && sd_val > 0) {
      kurt_val <- m4 / sd_val^4 - 3
    }
  }

  list(mean = mean_val, sd = sd_val, skew = skew_val, kurtosis = kurt_val)
}
