skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  mu <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mu)^3) / (s^3)
}

kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) return(NA_real_)
  mu <- mean(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  mean((x - mu)^4) / (s^4)  # non-excess kurtosis
}

summ_stats <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  jb_p <- NA_real_
  dago_p <- NA_real_
  ad_p <- NA_real_

  if (n >= 3) {
    if (requireNamespace("tseries", quietly = TRUE)) {
      jb_p <- tryCatch(tseries::jarque.bera.test(x)$p.value, error = function(e) NA_real_)
    }
    if (requireNamespace("fBasics", quietly = TRUE)) {
      dago_p <- tryCatch({
        pv <- fBasics::dagoTest(x)@test$p.value
        if (length(pv) > 1) pv[[1]] else pv
      }, error = function(e) NA_real_)
    }
    if (requireNamespace("nortest", quietly = TRUE)) {
      ad_p <- tryCatch(nortest::ad.test(x)$p.value, error = function(e) NA_real_)
    }
  }

  data.frame(
    n = n,
    median = if (n > 0) median(x) else NA_real_,
    mean = if (n > 0) mean(x) else NA_real_,
    sd = if (n > 1) sd(x) else NA_real_,
    skew = skewness(x),
    kurt = kurtosis(x),
    min = if (n > 0) min(x) else NA_real_,
    max = if (n > 0) max(x) else NA_real_,
    jb_p = jb_p,
    dago_p = dago_p,
    ad_p = ad_p,
    row.names = NULL
  )
}

compute_returns <- function(p, type = c("log", "simple")) {
  type <- match.arg(type)
  p <- as.numeric(p)
  if (type == "simple") {
    r <- p[-1] / p[-length(p)] - 1
  } else {
    r <- log(p[-1] / p[-length(p)])
  }
  r
}
