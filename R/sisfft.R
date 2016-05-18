OPT_BOUND <- 1e10
THETA_LIMIT <- 1e4


log_convolve_power <- function(log_pmf, L, alpha, delta) {
  if (L == 0) {
    return(c(0))
  } else if (L == 1) {
    return(log_pmf)
  }

  list[alpha_, delta_] <- accurate_error_bounds(L, 1 / alpha, delta)

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
                                                             alpha_,
                                                             delta_)
          squared <- T
        } else {
          answer <- log_convolve(pmf_power, answer, alpha_, delta_)
        }
      }
    }
    if (!need_to_square) {
      break
    }
    if (!squared) {
      pmf_power <- log_convolve_square(pmf_power, alpha_, delta_)
    }
  }

  answer
}

log_pvalue <- function(log_pmf, s0, L, beta) {
  # convert to 0-based s0 like the python
  s0 <- s0 - 1
  total_len <- iterated_convolution_lengths(length(log_pmf), L)[[1]]

  if (s0 >= total_len) {
    return(-Inf)
  }

  list[ignored, p_lower_preshift, p_upper_preshift] <- bounds(log_pmf, log_pmf, 0, 0.0,
                                                              s0, L, beta)
  list[sfft_good_preshift, sfft_pval_preshift] <- check_sfft_pvalue(p_lower_preshift,
                                                                    p_upper_preshift,
                                                                    beta)

  if (sfft_good_preshift) {
    return(sfft_pval_preshift)
  }

  theta <- compute_sisfft_theta(log_pmf, s0, L)

  theta <- clamp(theta, 0, THETA_LIMIT)
  list[shifted_pmf, log_mgf] <- shift(log_pmf, theta)

  list[log_delta, p_lower, p_upper] <- bounds(log_pmf, shifted_pmf, theta, log_mgf,
                                              s0, L, beta)
  list[sfft_good, sfft_pval] <- check_sfft_pvalue(p_lower,
                                                  p_upper,
                                                  beta)
  if (sfft_good) {
    return(sfft_pval)
  }
  delta <- exp(log_delta)
  conv <- log_convolve_power(shifted_pmf, L, alpha, delta)
  pvalue <- log_sum(unshift(conv, theta, log_mgf, L)[(s0 + 1):total_len])
  return(pvalue)
}

check_sfft_pvalue <- function(p_lower, p_upper, desired_beta) {
  added <- logaddexp(p_upper, p_lower)
  sfft_pval <- log(2) + p_lower + p_upper - added
  sfft_accuracy <- logsubexp(p_upper, p_lower) - added
  list(sfft_accuracy < log(desired_beta), sfft_pval)
}

bounds <- function(log_pmf, shifted_pmf, theta, log_mgf, s0, L, desired_beta) {
  stopifnot(theta >= -1)

  list[log_f0, fft_len] <- power_fft(shifted_pmf, L - 1)
  f0 <- exp(log_f0)
  error_estimate <- error_threshold_factor(fft_len) * (L - 1)
  v_lower <- log(pmax(f0 - error_estimate,
                      0.0))
  v_upper <- log(f0 + error_estimate)

  Q <- length(log_pmf)
  Q1 <- Q - 1
  tail_sums <- rep.int(0, Q)
  tail_sums[-1] = log_pmf[-1]
  for (i in Q1:1) {
    tail_sums[i] <- logaddexp(tail_sums[i + 1], log_pmf[i])
  }

  pval_estimate <- function(v) {
    limit <- length(v)
    low0 <- max(s0 - Q1, 0)
    mid0 <- min(s0, limit)

    if (theta != 0) { v <- v + (0:(limit - 1)) * -theta }
    if (log_mgf != 0) { v <- v + (L - 1) * log_mgf }

    k1_0 <- forward_seq(low0, mid0 - 1)
    k1 <- k1_0 + 1
    q1 <- log_sum(v[k1] + tail_sums[(s0 - k1_0) + 1])

    k2_0 <- forward_seq(mid0, limit - 1)
    k2 <- k2_0 + 1
    q2 <- log_sum(v[k2])
    q <- logaddexp(q1, q2)
    q
  }

  p_lower <- pval_estimate(v_lower)
  p_upper <- pval_estimate(v_upper)

  factor <- L * Q1 + 1 - s0

  frac <- if (abs(theta) < 1e-16) {
    -log(factor)
  } else if (theta > 0) {
    log1subexp(-theta) - log1subexp(-factor * theta)
  } else {
    logsubexp(-theta, 0) - logsubexp(-factor * theta, 0)
  }

  gamma <- p_lower + (theta * s0 - L * log_mgf) + frac + log(desired_beta / 2)
  list(gamma, p_lower, p_upper)
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
