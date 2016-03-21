OPT_BOUND <- 1e10
THETA_LIMIT <- 1e4

log_convolve_power <- function(log_pmf, L, desired_alpha, desired_delta) {
  if (L == 0) {
    return(c(0))
  } else if (L == 1) {
    return(log_pmf)
  }

  list[alpha, delta] <- accurate_error_bounds(L, 1 / desired_alpha, desired_delta)

  answer <- NULL
  pmf_power <- log_pmf

  while (T) {
    need_to_add <- (L %% 2) == 1
    L <- L %/% 2
    need_to_square <- L > 0
    squared <- F
    check_na(pmf_power)
    if (need_to_add) {
      if (is.null(answer)) {
        answer <- pmf_power
      } else {
        if (need_to_square) {
          list[answer, pmf_power] <- log_convolve_and_square(pmf_power,
                                                             answer,
                                                             alpha,
                                                             delta)
          squared <- T
        } else {
          answer <- log_convolve(pmf_power, answer, alpha, delta)
        }
      }
    }
    if (!need_to_square) {
      break
    }
    if (!squared) {
      pmf_power <- log_convolve_square(pmf_power, alpha, delta)
    }
  }

  answer
}

log_pvalue <- function(log_pmf, s0, L, desired_beta) {
  # convert to 0-based s0 like the python
  s0 <- s0 - 1
  total_len <- iterated_convolution_lengths(length(log_pmf), L)[[1]]

  if (s0 >= total_len) {
    return(-Inf)
  }

  theta <- compute_sisfft_theta(log_pmf, s0, L)
  if (theta < -1) {
    p <- log_pvalue(rev(log_pmf), total_len - s0 + 1, L, desired_beta)
    return(log1subexp(p))
  }

  theta <- clamp(theta, -THETA_LIMIT, THETA_LIMIT)
  list[shifted_pmf, log_mgf] <- shift(log_pmf, theta)

  alpha <- 2 / desired_beta
  log_delta <- lower_bound(log_pmf, shifted_pmf, theta, log_mgf,
                           s0, L,
                           desired_beta)
  delta <- exp(log_delta)
  conv <- log_convolve_power(shifted_pmf, L, alpha, delta)
  pvalue <- log_sum(unshift(conv, theta, log_mgf, L)[(s0 + 1):total_len])
  return(pvalue)
}

lower_bound <- function(log_pmf, shifted_pmf, theta, log_mgf, s0, L, desired_beta) {
  stopifnot(theta >= -1)

  list[log_f0, fft_len] <- power_fft(shifted_pmf, L - 1)
  f0 <- exp(log_f0)
  error_estimate <- error_threshold_factor(fft_len) * (L - 1)
  f_theta = log(pmax(f0 - error_estimate,
                     0.0))

  n <- length(log_pmf)
  n1 <- n - 1
  tail_sums <- rep.int(0, n)
  tail_sums[-1] = log_pmf[-1]
  for (i in n1:1) {
    tail_sums[i] <- logaddexp(tail_sums[i + 1], log_pmf[i])
  }

  limit <- length(f_theta)
  low0 <- max(s0 - n1, 0)
  k1_0 <- low0:(min(s0, limit) - 1)
  k1 <- k1_0 + 1
  q1 <- log_sum(f_theta[k1] + (-k1_0 * theta + (L - 1)* log_mgf)
                + tail_sums[(s0 - k1_0) + 1])

  k2_0 <- forward_seq(min(s0, limit), limit - 1)
  k2 <- k2_0 + 1
  # need to handle empty k2_0 with - operator
  q2 <- log_sum(f_theta[k2] + ((0-k2_0) * theta + (L - 1) * log_mgf))
  q <- logaddexp(q1, q2)
  factor <- L * n1 + 1 - s0

  frac <- if (abs(theta) < 1e-16) {
    -log(factor)
  } else if (theta > 0) {
    log1subexp(-theta) - log1subexp(-factor * theta)
  } else {
    logsubexp(-theta, 0) - logsubexp(-factor * theta, 0)
  }

  gamma <- q + (theta * s0 - L * log_mgf) + frac + log(desired_beta / 2)
  gamma
}

compute_sisfft_theta <- function(log_pmf, s0, L) {
  s0L <- s0 / L
  f <- function(theta) {
    log_mgf(log_pmf, theta) - theta * s0L
  }
  optimize(f, lower = -OPT_BOUND, upper = OPT_BOUND)$minimum
}

accurate_error_bounds <- function(L, beta, gamma) {
  beta2 <- (1 + beta)^(1 / (L - 1)) - 1

  count <- floor(log2(L))
  accum <- 0
  ijs <- rep.int(0, count + 1)
  Ljs <- rep.int(0, count + 1)
  popcnt <- 0
  for (i in 0:count) {
    this <- 2**i
    if (bitwAnd(L, this)) {
      accum <- accum + this
      popcnt <- popcnt + 1
      ijs[popcnt] <- i
      Ljs[popcnt] <- accum
    }
  }
  ijs <- ijs[1:popcnt]
  Ljs <- Ljs[1:popcnt]

  r <- function(i) { (1 + beta2)**(2**i - 1) - 1 }
  rbar <- function(j) { (1 + beta2)**(Ljs[j] - 1) - 1 }

  d <- function(is) {
    sapply(is, function(i) {
      k <- arange(i)
      sum(2**(i - k) * (1 + r(k)))
    })
  }
  dbar <- function(j) {
    k <- forward_seq(2, j)
    sum(d(ijs[1:j])) + sum(2 + rbar(k) + r(k))
  }
  delta <- gamma / dbar(popcnt)
  list(1/beta2, delta)
}
