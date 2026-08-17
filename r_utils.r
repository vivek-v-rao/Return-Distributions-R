# r_utils.r
#
# Shared utility functions sourced by all analysis scripts.

print_script_header <- function() {
  # Print the running script's name and current timestamp to stdout.
  f <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) cat("program:", sub("--file=", "", f[[1L]], fixed = TRUE),
                     "  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
}

read_price_file <- function(infile, max.assets = Inf) {
  # Read a CSV file containing dates in the first column and asset prices
  # in the remaining columns.
  #
  # Args:
  #   infile: Name of the CSV file to read.
  #   max.assets: Maximum number of asset price columns to read.
  #     Use Inf to read all asset price columns.
  #
  # Returns:
  #   A list with components:
  #     dates: The first column of the input file.
  #     prices: A data frame of numeric asset prices.

  dat <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(dat) < 2) {
    stop("Need a date column and at least one price column")
  }

  dates  <- dat[[1]]
  prices <- dat[-1]

  if (!is.infinite(max.assets)) {
    if (!is.numeric(max.assets) || length(max.assets) != 1 ||
        !is.finite(max.assets) || max.assets < 1) {
      stop("max.assets must be a positive number or Inf")
    }

    max.assets <- as.integer(max.assets)
    max.assets <- min(max.assets, ncol(prices))
    prices <- prices[seq_len(max.assets)]
  }

  for (j in seq_along(prices)) {
    prices[[j]] <- as.numeric(prices[[j]])
  }

  list(
    dates  = dates,
    prices = as.data.frame(prices, check.names = FALSE)
  )
}


keep_selected_symbols <- function(prices, symbols = character(0)) {
  # Keep only selected asset columns if symbols is non-empty.
  #
  # Args:
  #   prices: Data frame of asset prices.
  #   symbols: Character vector of symbols to keep. If empty, keep all columns.
  #
  # Returns:
  #   A data frame containing only the requested symbols.

  if (length(symbols) == 0) {
    return(prices)
  }

  symbols <- as.character(symbols)
  missing <- setdiff(symbols, names(prices))

  if (length(missing) > 0) {
    stop(
      "symbols not found in price file: ",
      paste(missing, collapse = ", ")
    )
  }

  prices[, symbols, drop = FALSE]
}


compute_log_returns <- function(prices) {
  # Compute log returns from asset prices.
  #
  # Args:
  #   prices: Data frame whose columns are price series.
  #
  # Returns:
  #   A data frame of log returns, with one fewer row than prices.
  #   Non-positive prices give NA returns.

  nr <- nrow(prices)
  nc <- ncol(prices)

  ret <- matrix(NA_real_, nrow = nr - 1, ncol = nc)
  colnames(ret) <- colnames(prices)

  for (j in seq_len(nc)) {
    x <- prices[[j]]
    x[!is.finite(x) | x <= 0] <- NA_real_
    ret[, j] <- log(x[-1]) - log(x[-nr])
  }

  as.data.frame(ret, check.names = FALSE)
}


ewma_vol <- function(x, lambda = 0.94, warmup = 100L, center = FALSE) {
  # Trailing RiskMetrics EWMA volatility forecast.
  #
  #   sigma^2_t = lambda * sigma^2_{t-1} + (1 - lambda) * (x_{t-1} - mu)^2
  #
  # sigma_t uses observations through t-1 only, so x_t / sigma_t is an
  # out-of-sample standardization with no look-ahead. The recursion is
  # seeded with the sample variance of the first warmup observations,
  # which themselves are returned as NA.
  #
  # Args:
  #   x: Numeric return vector, possibly containing NA.
  #   lambda: Decay factor in (0, 1). RiskMetrics uses 0.94 for daily data.
  #   warmup: Number of leading observations used to seed the recursion.
  #   center: If TRUE, subtract the sample mean before squaring.
  #     RiskMetrics assumes a zero mean, so the default is FALSE.
  #
  # Returns:
  #   Numeric vector of the same length as x, with NA in the first
  #   warmup positions.

  if (!is.finite(lambda) || lambda <= 0 || lambda >= 1) {
    stop("lambda must lie strictly between 0 and 1")
  }

  n      <- length(x)
  warmup <- as.integer(warmup)

  if (warmup < 2L) stop("warmup must be at least 2")

  sig <- rep(NA_real_, n)

  if (n <= warmup) return(sig)

  mu <- if (center) mean(x[is.finite(x)]) else 0.0
  d  <- x - mu

  seed <- d[seq_len(warmup)]
  seed <- seed[is.finite(seed)]

  if (length(seed) < 2) return(sig)

  v <- mean(seed^2)

  for (t in seq.int(warmup + 1L, n)) {
    prev <- d[t - 1L]
    # A missing return leaves the variance forecast unchanged.
    if (is.finite(prev)) v <- lambda * v + (1 - lambda) * prev^2
    sig[t] <- sqrt(v)
  }

  sig
}


standardize_by_ewma_vol <- function(ret, lambda = 0.94, warmup = 100L,
                                    center = FALSE) {
  # Divide each return series by its own trailing EWMA volatility.
  #
  # Args:
  #   ret: Data frame whose columns are return series.
  #   lambda, warmup, center: Passed to ewma_vol().
  #
  # Returns:
  #   A data frame of standardized returns with the same shape as ret.
  #   The first warmup rows of each column are NA.

  for (j in seq_along(ret)) {
    x <- as.numeric(ret[[j]])
    s <- ewma_vol(x, lambda = lambda, warmup = warmup, center = center)
    s[is.finite(s) & s <= 0] <- NA_real_
    ret[[j]] <- x / s
  }

  ret
}


safe_moments <- function(x) {
  # Compute skewness and excess kurtosis.
  #
  # Args:
  #   x: Numeric return vector.
  #
  # Returns:
  #   A named numeric vector with skew and excess_kurtosis.

  x <- x[is.finite(x)]

  if (length(x) < 2 || sd(x) <= 0) {
    return(c(skew = NA_real_, excess_kurtosis = NA_real_))
  }

  z <- (x - mean(x)) / sd(x)

  c(
    skew            = mean(z^3),
    excess_kurtosis = mean(z^4) - 3
  )
}
