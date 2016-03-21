log_convolve_naive <- function(log_u, log_v) {
  nu <- length(log_u)
  nv <- length(log_v)
  nc <- nu + nv - 1
  log_c <- rep.int(-Inf, nc)
  convolve_naive_into(log_c, 1:nc, log_u, log_v)
}

convolve_naive_into <- function(log_c, locations, log_u, log_v) {
  nu <- length(log_u)
  nv <- length(log_v)
  nc <- nu + nv - 1
  stopifnot(length(log_c) == nc)

  for (k in locations) {
    k0 <- k - 1
    low_j0 <- max(0, k0 - nv + 1)
    hi_j0 <- min(k0 + 1, nu)
    low_j <- low_j0 + 1
    hi_j <- hi_j0 + 1

    slice_u <- log_u[low_j:(hi_j - 1)]
    slice_v <- log_v[(k0 - low_j0 + 1):(k0 - hi_j0 + 2)]
    log_c[k] <- log_sum(slice_u + slice_v)
  }
  return(log_c)
}

power_fft <- function(log_u, L) {
  list[true_len, fft_len] <- iterated_convolution_lengths(length(log_u), L)
  fft_ <- fft(pad_to_length(exp(log_u), fft_len))
  conv <- fft(fft_^L, inverse = T)[1:true_len] / fft_len
  list(log(abs(conv)), fft_len)
}
