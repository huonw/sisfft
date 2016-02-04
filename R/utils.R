# allow destructuring a list: list[var1, var2] <- func()
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

forward_seq <- function(lo, hi) {
  if (hi < lo) { c() } else { lo:hi }
}

arange <- function(n) {
  forward_seq(0, n - 1)
}

EPS <- .Machine$double.eps / 2

clamp <- function(x, lo, hi) {
  pmax.int(pmin.int(x, hi), lo)
}

log_min_pos <- function(log_pmf) {
  min(log_pmf[log_pmf > -Inf])
}

log_dynamic_range <- function(log_pmf) {
  hi <- log_sum(2 * log_pmf) / 2
  lo <- log_min_pos(log_pmf)
  hi - lo
}

shift <- function(log_pmf, theta) {
  shifted <- log_pmf + theta * arange(length(log_pmf))
  log_mgf <- log_sum(shifted)
  list(shifted - log_mgf, log_mgf)
}

log_dynamc_range_shifted <- function(log_pmf, theta) {
  shifted <- log_pmf + arange(length(log_pmf)) * theta
  lo <- log_min_pos(shifted)
  hi <- log_sum(shifted * 2) / 2
  hi - lo
}

unshift <- function(convolved, theta, mgfs, multiplicities) {
  c <- convolved - theta * arange(length(convolved))
  for (i in 1:length(mgfs)) {
    c <- c + multiplicities[i] * mgfs[i]
  }
  c
}

logaddexp <- function(log_x, log_y) {
  ifelse(log_x == log_y, log_x + log(2),
         {
          lo <- pmin(log_x, log_y)
          hi <- pmax(log_x, log_y)
          hi + log1p(exp(lo - hi))
         })
}

logsubexp <- function(log_x, log_y) {
  stopifnot(log_x >= log_y)
  log_x + log1p(-exp(log_y - log_x))
}

log1subexp <- function(log_y) {
  stopifnot(log_y <= 0)
  log1p(-exp(log_y))
}

log_sum <- function(log_u) {
  if (length(log_u) == 0) {
    return -Inf
  }

  maxi <- which.max(log_u)
  max <- log_u[maxi]
  if (max == -Inf) {
    -Inf
  } else {
    e <- exp(log_u - max)
    max + log1p(sum(e[forward_seq(1, (maxi - 1))]) +
              sum(e[forward_seq(maxi + 1, length(e))]))
  }
}

log_mgf <- function(log_pmf, theta) {
  log_sum(log_pmf + theta * arange(length(log_pmf)))
}

error_threshold_factor <- function(conv_len) {
  c = if (conv_len > 2**5) {
    13.5
  } else {
    16
  }
  EPS * c * log(conv_len, 2)
}

pairwise_convolution_lengths <- function(a, b) {
  true <- a + b - 1
  list(true, nextn(true, factors = 2))
}
iterated_convolution_lengths <- function(a, L) {
  true <- (a - 1) * L + 1
  list(true, nextn(true, factors = 2))
}

# need this to get O(1) appends
# http://stackoverflow.com/a/32870310/1256624
expandingVector <- function(capacity = 10) {
  buffer <- vector('logical', capacity)
  length <- 0

  methods <- list()

  methods$double.size <- function() {
    buffer <<- c(buffer, vector('logical', capacity))
    capacity <<- capacity * 2
  }

  methods$add <- function(val) {
    if(length == capacity) {
      methods$double.size()
    }

    length <<- length + 1
    buffer[[length]] <<- val
  }

  methods$length <- function() {
    length
  }

  methods$as.vector <- function() {
    b <- buffer[0:length]
    return(b)
  }

  methods
}
