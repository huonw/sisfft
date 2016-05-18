COST_RATIO <- 0.5
COST_RATIO_SQUARE <- COST_RATIO * 0.5

OPT_BOUND <- 1e4

# hack to distinguish these from vectors
ESTIMATE_ONE_SPLIT <- list(1)
ESTIMATE_TWO_SPLITS <- list(2)

log_convolve <- function(log_pmf1, log_pmf2, alpha, delta = 0) {
  if (delta != 0) {
    no_shift(log_pmf1, log_pmf2, alpha, delta,
             pairwise = T,
             square_1 = F)[[1]]
  } else {
    no_lower_bound(log_pmf1, log_pmf2, alpha)
  }
}

log_convolve_square <- function(log_pmf, alpha, delta = 0) {
  if (delta == 0) {
    log_convolve(log_pmf, log_pmf, alpha, delta)
  } else {
    no_shift(log_pmf, c(), alpha, delta,
             pairwise = F,
             square_1 = T)[[2]]
  }
}

log_convolve_and_square <- function(log_pmf1, log_pmf2, alpha, delta = 0) {
  if (delta == 0) {
    # this is a suboptimal approach
    co <- convolve(log_pmf1, log_pmf2, alpha, 0)
    sq <- convolve_square(log_pmf1, alpha, 0)
    list(co, sq)
  } else {
    no_shift(log_pmf1, log_pmf2, alpha, delta,
             pairwise = T,
             square_1 = T)
  }
}

no_lower_bound <- function(log_pmf1, log_pmf2, alpha) {
  list[true_conv_len, fft_conv_len] <- pairwise_convolution_lengths(length(log_pmf1),
                                                                    length(log_pmf2))
  pmf1 <- pad_to_length(exp(log_pmf1), fft_conv_len)

  fft1 <- fft(pmf1)
  list[direct, bad_places] <- direct_fft_conv(log_pmf1, pmf1, fft1,
                                              log_pmf2,
                                              true_conv_len, fft_conv_len,
                                              alpha, 0)
  used_nc <- use_nc_if_better(log_pmf1, ESTIMATE_ONE_SPLIT,
                              log_pmf2, ESTIMATE_TWO_SPLITS,
                              direct, bad_places,
                              COST_RATIO)
  if (!is.null(used_nc)) {
    return(used_nc);
  }

  theta <- compute_theta(log_pmf1, log_pmf2)
  list[s1, log_mgf1] <- shift(log_pmf1, theta)
  list[s2, log_mgf2] <- shift(log_pmf2, theta)
  convolved <- no_shift(s1, s2, alpha, 0,
                        pairwise = T,
                        square_1 = F)[[1]]
  unshift(convolved, theta, c(log_mgf1, log_mgf2), c(1, 1))
}

no_shift <- function(log_pmf1, log_pmf2, alpha, delta,
                     pairwise, square_1) {
  stopifnot(pairwise || square_1)

  len1 <- length(log_pmf1)
  len2 <- length(log_pmf2)
  list[true_conv_len, fft_conv_len] <- pairwise_convolution_lengths(len1,
                                                                    len2)
  list[true_conv_len_sq, fft_conv_len_sq] <- pairwise_convolution_lengths(len1,
                                                                          len1)
  can_reuse_pairwise <- pairwise && fft_conv_len == fft_conv_len_sq

  answer <- NULL
  answer_sq <- NULL

  raw_pmf1 <- exp(log_pmf1)
  if (pairwise) {
    pmf1 <- pad_to_length(raw_pmf1, fft_conv_len)
    fft1 <- fft(pmf1)

    list[direct, bad_places] <- direct_fft_conv(log_pmf1, pmf1, fft1,
                                                log_pmf2,
                                                true_conv_len, fft_conv_len,
                                                alpha, delta)
  }

  if (square_1) {
    list[direct_sq, bad_places_sq] <- if (can_reuse_pairwise) {
      direct_fft_conv(log_pmf1, raw_pmf1, fft1,
                      NULL,
                      true_conv_len_sq, fft_conv_len_sq,
                      alpha, delta)
    } else {
      fft1_sq <- fft(pad_to_length(raw_pmf1, fft_conv_len_sq))
      direct_fft_conv(log_pmf1, raw_pmf1, fft1_sq,
                      NULL,
                      true_conv_len_sq, fft_conv_len_sq,
                      alpha, delta)
    }
  }

  if (pairwise) {
    used_nc <- use_nc_if_better(log_pmf1, ESTIMATE_ONE_SPLIT,
                                log_pmf2, ESTIMATE_TWO_SPLITS,
                                direct, bad_places,
                                COST_RATIO)
    if (!is.null(used_nc)) {
      answer <- used_nc;
    }
  }
  if (square_1) {
    used_nc <- use_nc_if_better(log_pmf1, ESTIMATE_TWO_SPLITS,
                                log_pmf1, ESTIMATE_TWO_SPLITS,
                                direct_sq, bad_places_sq,
                                COST_RATIO_SQUARE)
    if (!is.null(used_nc)) {
      answer_sq <- used_nc;
    }
  }

  need_to_pairwise <- is.null(answer) && pairwise
  need_to_square <- is.null(answer_sq) && square_1
  can_reuse_pairwise <- can_reuse_pairwise && need_to_pairwise

  if (need_to_pairwise) {
    limit <- split_limit(len1, len2, NULL,
                         length(bad_places),
                         COST_RATIO)
    maxima1 <- split_maxima(log_pmf1, fft_conv_len,
                            alpha, delta,
                            limit)
    if (!is.null(maxima1)) {
      limit <- split_limit(len1, len2,
                           length(maxima1),
                           length(bad_places),
                           COST_RATIO)
      maxima2 <- split_maxima(log_pmf2, fft_conv_len,
                              alpha, delta,
                              limit)
    } else {
      maxima2 <- list(NULL, NULL)
    }
  }

  if (need_to_square) {
    if (can_reuse_pairwise) {
      maxima1_sq <- maxima1
    } else {
      limit <- split_limit(len1, NULL, NULL,
                           length(bad_places_sq),
                           COST_RATIO_SQUARE)
      maxima1_sq <- split_maxima(log_pmf1, fft_conv_len_sq,
                                 alpha, delta,
                                 limit)
    }
  }

  if (need_to_pairwise) {
    used_nc <- use_nc_if_better(log_pmf1, maxima1,
                                log_pmf2, maxima2,
                                direct, bad_places,
                                COST_RATIO)
    if (!is.null(used_nc)) {
      answer <- used_nc
    }
  }
  if (need_to_square) {
    used_nc <- use_nc_if_better(log_pmf1, maxima1_sq,
                                log_pmf1, maxima1_sq,
                                direct_sq, bad_places_sq,
                                COST_RATIO_SQUARE)
    if (!is.null(used_nc)) {
      answer_sq <- direct_sq
    }
  }

  need_to_pairwise <- need_to_pairwise && is.null(answer)
  need_to_square <- need_to_square && is.null(answer_sq)
  can_reuse_pairwise <- can_reuse_pairwise && need_to_pairwise

  if (need_to_pairwise) {
    splits1 <- splits_from_maxima(log_pmf1, maxima1, delta,
                                  fft_conv_len)
    splits2 <- splits_from_maxima(log_pmf2, maxima2, delta,
                                  fft_conv_len)
  }
  if (need_to_square) {
    if (can_reuse_pairwise) {
      splits1_sq <- splits1
    } else {
      splits1_sq <- splits_from_maxima(log_pmf1, maxima1_sq, delta,
                                       fft_conv_len_sq)
    }
  }

  if (need_to_pairwise) {
    ffts1 <- mvfft(exp(splits1))
    ffts2 <- mvfft(exp(splits2))
  }
  if (need_to_square) {
    if (can_reuse_pairwise) {
      ffts1_sq <- ffts1
    } else {
      ffts1_sq <- mvfft(exp(splits1_sq))
    }
  }
  if (need_to_pairwise) {
    accum <- rep.int(-Inf, true_conv_len)
    for (i in 1:length(maxima1)) {
      normaliser1 <- maxima1[i]
      fft1 <- ffts1[,i]
      for (j in 1:length(maxima2)) {
        normaliser2 <- maxima2[j]
        fft2 <- ffts2[,j]

        conv <- filtered_mult_ifft(fft1, normaliser1,
                                   fft2, normaliser2,
                                   true_conv_len,
                                   fft_conv_len)
        accum <- logaddexp(accum, conv)
      }
    }
    answer <- accum
  }
  if (need_to_square) {
    accum_sq <- rep.int(-Inf, true_conv_len_sq)

    for (i in 1:length(maxima1_sq)) {
      normaliser1 <- maxima1_sq[i]
      fft1 <- ffts1_sq[,i]
      conv_self <- filtered_mult_ifft(fft1, normaliser1,
                                      fft1, normaliser1,
                                      true_conv_len_sq,
                                      fft_conv_len_sq)
      accum_sq <- logaddexp(accum_sq, conv_self)
      for (j in forward_seq(i + 1, length(maxima1_sq))) {
        normaliser2 <- maxima1_sq[j]
        fft2 <- ffts1_sq[, j]
        conv <- filtered_mult_ifft(fft1, normaliser1,
                                   fft2, normaliser2,
                                   true_conv_len_sq,
                                   fft_conv_len_sq)
        accum_sq <- logaddexp(accum_sq, conv + log(2))
      }
    }
    answer_sq <- accum_sq
  }

  if (pairwise) stopifnot(length(answer) == true_conv_len)
  if (square_1) stopifnot(length(answer_sq) == true_conv_len_sq)
  list(answer, answer_sq)
}

direct_fft_conv <- function(log_pmf1, pmf1, fft1, log_pmf2, true_conv_len, fft_conv_len,
                            alpha, delta) {
  if (is.null(log_pmf2)) {
    norms <- norm(pmf1, type = "2")^2
    fft_conv <- fft1^2
  } else {
    raw_pmf2 <- exp(log_pmf2)
    pmf2 <- pad_to_length(raw_pmf2, fft_conv_len)
    fft2 <- fft(pmf2)
    norms <- norm(pmf1, type = "2") * norm(raw_pmf2, type = "2")
    fft_conv <- fft1 * fft2
  }

  raw_conv <- abs(fft(fft_conv, inverse = T)[1:true_conv_len]) / fft_conv_len
  error_level <- error_threshold_factor(fft_conv_len) * norms
  threshold <- error_level * (alpha + 1)
  places_of_interest <- raw_conv <= threshold

  if (delta > error_level) {
    ignore_level <- delta - error_level

    conv_to_ignore <- raw_conv < ignore_level
    raw_conv[conv_to_ignore] <- 0
    places_of_interest <- places_of_interest & !conv_to_ignore
  }

  log_conv <- log(raw_conv)
  if (!any(places_of_interest)) {
    return(list(log_conv, c()))
  }

  Q <- true_conv_len
  support1 <- log_pmf1 > -Inf
  if (is.null(log_pmf2)) {
    support2 <- c()
  } else {
    support2 <- log_pmf2 > -Inf
  }

  if (Q <= 2**42 && (!all(support1) || !all(support2))) {
    padded_support1 <- pad_to_length(support1)
    fft_support1 <- fft(padded_support1)
    if (is.null(log_pmf2)) {
      fft_support <- fft_support1^2
    } else {
      support2 <- log_pmf > -Inf
      padded_support2 <- pad_to_length(support2)
      fft_support2 <- fft(padded_support2)
      fft_support <- fft_support1 * fft_support2
    }

    raw_support <- abs(fft(fft_support, inverse = T)[1:true_conv_len]) / fft_conv_len
    threshold <- error_threshold_factor(fft_conv_len) * 2 * Q
    zeros <- raw_support <= threshold
    log_conv[zeros] <- -Inf
    bad_places <- which(!zeros & places_of_interest)
  } else {
    bad_places <- which(places_of_interest)
  }

  return(list(log_conv, bad_places))
}

split_maxima <- function(log_pmf, fft_conv_len, alpha, delta,
                         split_limit) {
  sort_idx <- order(log_pmf, decreasing = T)
  raw_threshold <- error_threshold_factor(fft_conv_len) * alpha
  log_threshold <- -log(raw_threshold) / 2
  log_delta <- log(delta)
  log_current_max <- log_pmf[sort_idx[1]]
  log_current_norm_sq <- -Inf
  maxes <- expandingVector()
  maxes$add(log_current_max)

  for (i in 1:length(sort_idx)) {
    idx <- sort_idx[i]
    log_val <- log_pmf[idx]
    if (log_val < log_delta) {
      break
    }
    log_next_norm_sq <- logaddexp(log_current_norm_sq, log_val * 2)
    r <- log_next_norm_sq / 2 - log_val
    if (r < log_threshold) {
      log_current_norm_sq <- log_next_norm_sq
    } else {
      maxes$add(log_val)
      log_current_norm_sq <- log_val * 2
      log_current_max <- log_val
      if (maxes$length() > split_limit) {
        return(NULL)
      }
    }
  }

  maxes$as.vector()
}

splits_from_maxima <- function(log_pmf, split_maxima, delta,
                               fft_conv_len) {
  log_delta <- log(delta)
  pmf_len <- length(log_pmf)
  num_splits <- length(split_maxima)
  splits <- matrix(-Inf, fft_conv_len, num_splits)

  for (i in 1:num_splits) {
    lo <- if (i == num_splits) {
      log_delta * (1 - .Machine$double.eps / 2)
    } else {
      split_maxima[i + 1]
    }
    hi <- split_maxima[i]

    elems <- which((lo < log_pmf) & (log_pmf <= hi))
    splits[elems,i] <- log_pmf[elems] - hi
  }
  splits
}

split_limit <- function(len1, len2, maxima1, len_bad_places, cost_ratio) {
  Q <- if (!is.null(len2)) {
    max(len1, len2)
  } else {
    len1
  }
  nc_cost <- Q * len_bad_places
  split_ratio <- nc_cost / (Q * log2(Q) * cost_ratio)
  lim <- if (is.null(len2)) {
    sqrt(split_ratio)
  } else {
    other_split <- if (is.null(maxima1)) { 1 } else { maxima1 }
    split_ratio / other_split
  }
  lim
}

maxima_to_length <- function(maxima) {
  if (is.list(maxima)) {
    maxima[[1]]
  } else {
    length(maxima)
  }
}

is_nc_faster <- function(len1, maxima1,
                         len2, maxima2,
                         len_bad_places,
                         cost_ratio) {
  if (len_bad_places == 0) {
    T
  } else if (is.null(maxima1) || is.null(maxima2)) {
    T
  } else {
    Q <- max(len1, len2)
    nc_cost <- min(Q * len_bad_places, len1 * len2)
    len1 <- maxima_to_length(maxima1)
    len2 <- maxima_to_length(maxima2)
    psfft_cost <- len1 * len2 * Q * log2(Q)
    scaled_cost <- psfft_cost * cost_ratio

    scaled_cost > nc_cost
  }
}

use_nc_if_better <- function(log_pmf1, maxima1,
                             log_pmf2, maxima2,
                             direct, bad_places,
                             cost_ratio) {
  nc_is_better <- is_nc_faster(length(log_pmf1), maxima1,
                               length(log_pmf2), maxima2,
                               length(bad_places),
                               cost_ratio)

  if (nc_is_better) {
    convolve_naive_into(direct, bad_places, log_pmf1, log_pmf2)
  } else {
    NULL
  }
}

filtered_mult_ifft <- function(fft1, normaliser1, fft2, normaliser2,
                               true_conv_len, fft_conv_len) {
  norm1 <- norm(fft1, type="2")
  norm2 <- norm(fft2, type="2")

  threshold <- error_threshold_factor(fft_conv_len) * norm1 * norm2 / fft_conv_len
  entire_conv <- fft(fft1 * fft2, inverse = T)[1:true_conv_len] / fft_conv_len
  filtered <- ifelse(abs(entire_conv) > threshold,
                     Re(entire_conv),
                     0)

  log(filtered) + (normaliser1 + normaliser2)
}

compute_theta <- function(log_pmf1, log_pmf2) {
  f <- function(theta) {
    r1 <- log_dynamic_range_shifted(log_pmf1, theta)
    r2 <- log_dynamic_range_shifted(log_pmf2, theta)
    return(r1 * r2)
  }
  optimise(f, lower = -OPT_BOUND, upper = OPT_BOUND)$minimum
}
