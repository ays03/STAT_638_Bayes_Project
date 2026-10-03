# Bayesian harmonic regression with AR(p), heavy-tailed errors, and missing days.
#
# MODEL
#   Response on a chosen scale: z_t = g(y_t), with g = log or identity.
#   Defined on a gap-free daily grid t = 1, ..., T, so the error process exists
#   on every calendar day whether or not y_t was observed.
#
#       z_t = x_t' beta + eps_t
#       eps_t = sum_{k=1}^{p} phi_k eps_{t-k} + eta_t
#       eta_t ~ N(0, sigma^2 / lambda_t),   lambda_t ~ Gamma(nu/2, nu/2)
#
#   Marginalising lambda_t gives eta_t ~ t_nu(0, sigma^2). Setting heavy=FALSE
#   fixes lambda_t = 1 (Gaussian errors); setting p = 0 gives independent errors.
#
#   x_t contains: intercept, centred linear trend in years, J harmonic pairs at
#   period 365 days, and six day-of-week dummies (Monday is the baseline).
#
#   (eps_1, ..., eps_p) is drawn from the exact stationary distribution of the
#   AR(p), so the chain starts stationary and no burn-in of the error process
#   is needed.
#
# PRIORS
#   beta        ~ N(0, diag(tau^2)), tau set by response scale (weakly informative)
#   sigma       ~ half-t(3, 0, A), via the standard inverse-gamma auxiliary variable
#   phi         ~ Uniform over the stationary region
#   nu - 2      ~ Exponential(mean 10)
#   z_t missing ~ implied by the AR process (imputed, not dropped)
#
# SAMPLING
#   Gibbs for beta (conjugate normal), sigma^2 and its auxiliary xi
#   (inverse-gamma), lambda (gamma), and the missing z_t (normal).
#   Random-walk Metropolis for phi (rejecting outside the stationary region)
#   and for log(nu - 2).
#
# R translation of src/model.py. Day indices in this file are 1-based.

PRIOR_TAU <- list(log = c(5.0, 2.0), identity = c(100.0, 40.0))
PRIOR_SIGMA_A <- list(log = 1.0, identity = 20.0)

# ----------------------------------------------------------------------------
# design matrix
# ----------------------------------------------------------------------------

build_design <- function(t, dow, J, t_center = NULL) {
  # t: 1-based integer day index. dow: 0 = Monday. Returns list(X, names).
  t <- as.numeric(t)
  if (is.null(t_center)) t_center <- mean(t)
  cols <- list(rep(1, length(t)), (t - t_center) / 365.0)
  xnames <- c("intercept", "trend_per_year")
  if (J >= 1) {
    for (j in seq_len(J)) {
      cols[[length(cols) + 1]] <- sin(2 * pi * j * t / 365.0)
      cols[[length(cols) + 1]] <- cos(2 * pi * j * t / 365.0)
      xnames <- c(xnames, paste0("sin_h", j), paste0("cos_h", j))
    }
  }
  daynames <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  for (k in 1:6) {
    cols[[length(cols) + 1]] <- as.numeric(as.integer(dow) == k)
    xnames <- c(xnames, paste0("dow_", daynames[k + 1]))
  }
  list(X = do.call(cbind, cols), names = xnames)
}

# ----------------------------------------------------------------------------
# AR(p) helpers
# ----------------------------------------------------------------------------

is_stationary <- function(phi) {
  # Roots of 1 - phi_1 B - ... - phi_p B^p outside the unit circle.
  p <- length(phi)
  if (p == 0) return(TRUE)
  companion <- matrix(0, p, p)
  companion[1, ] <- phi
  if (p > 1) companion[2:p, 1:(p - 1)] <- diag(p - 1)
  max(Mod(eigen(companion, symmetric = FALSE, only.values = TRUE)$values)) < 1 - 1e-8
}

stationary_cov <- function(phi, sigma2) {
  # Covariance of (eps_1, ..., eps_p) under the stationary AR(p).
  p <- length(phi)
  if (p == 0) return(matrix(numeric(0), 0, 0))
  if (p == 1) return(matrix(sigma2 / (1 - phi[1]^2), 1, 1))
  A <- matrix(0, p, p)
  A[1, ] <- phi
  if (p > 1) A[2:p, 1:(p - 1)] <- diag(p - 1)
  Q <- matrix(0, p, p)
  Q[1, 1] <- sigma2
  # vec(V) = (I - A kron A)^{-1} vec(Q), column-major vec.
  M <- diag(p * p) - kronecker(A, A)
  V <- matrix(solve(M, as.vector(Q)), p, p)
  0.5 * (V + t(V))
}

ar_resid <- function(eps, phi) {
  # eta_t for t = p+1, ..., T: eps_t - sum phi_k eps_{t-k}.
  p <- length(phi)
  if (p == 0) return(eps)
  n <- length(eps)
  out <- eps[(p + 1):n]
  for (k in seq_len(p)) {
    out <- out - phi[k] * eps[(p - k + 1):(n - k)]
  }
  out
}

quasi_difference <- function(M, phi) {
  # Apply (1 - phi_1 B - ... - phi_p B^p) down the rows of M.
  if (is.null(dim(M))) M <- cbind(as.numeric(M))
  p <- length(phi)
  if (p == 0) return(M)
  n <- nrow(M)
  out <- M[(p + 1):n, , drop = FALSE]
  for (k in seq_len(p)) {
    out <- out - phi[k] * M[(p - k + 1):(n - k), , drop = FALSE]
  }
  out
}

marginal_sd_ratio <- function(phi) {
  # sd of the stationary AR(p) divided by the one-step innovation sd.
  if (length(phi) == 0) return(1)
  sqrt(stationary_cov(as.numeric(phi), 1)[1, 1])
}

log_t <- function(e, sigma, nu, heavy) {
  # Log density of a scaled t (or normal if heavy is FALSE). sigma and nu may
  # be scalars or vectors recycling against e.
  if (!isTRUE(heavy) || !all(is.finite(nu))) {
    return(-0.5 * log(2 * pi) - log(sigma) - 0.5 * (e / sigma)^2)
  }
  lgamma(0.5 * (nu + 1)) - lgamma(0.5 * nu) - 0.5 * log(pi * nu) - log(sigma) -
    0.5 * (nu + 1) * log1p((e / sigma)^2 / nu)
}

# ----------------------------------------------------------------------------
# sampler
# ----------------------------------------------------------------------------

fit <- function(
  y,
  observed,
  t,
  dow,
  J = 1,
  p = 1,
  heavy = TRUE,
  scale = "log",
  n_iter = 20000,
  burn = 5000,
  thin = 5,
  seed = 1,
  ll_index = NULL,
  n_paths = 400,
  verbose = FALSE
) {
  # Run one chain. y and observed are on the full daily grid (length T);
  # y may contain NA wherever observed is FALSE.
  #
  # ll_index fixes the days (1-based) at which the pointwise log-likelihood is
  # recorded. Pass the same index set to every model so WAIC compares like with
  # like; if NULL, the model's own maximal usable set is used.
  #
  # Returns a list of posterior draws (rows = saved iterations).
  set.seed(seed)
  y <- as.numeric(y)
  observed <- as.logical(observed)
  nT <- length(y)
  des <- build_design(t, dow, J)
  X <- des$X
  xnames <- des$names
  k <- ncol(X)

  z <- rep(NA_real_, nT)
  if (scale == "log") z[observed] <- log(y[observed]) else z[observed] <- y[observed]
  z[!observed] <- mean(z[observed])
  miss_idx <- which(!observed)

  tau_int <- PRIOR_TAU[[scale]][1]
  tau_oth <- PRIOR_TAU[[scale]][2]
  prior_prec <- diag(c(1 / tau_int^2, rep(1 / tau_oth^2, k - 1)))
  sigma_A <- PRIOR_SIGMA_A[[scale]]
  nu_sigma <- 3.0
  nu_exp_mean <- 10.0

  beta <- qr.solve(X[observed, , drop = FALSE], z[observed])
  eps <- as.numeric(z - X %*% beta)
  ev <- eps[observed]
  sigma2 <- sum((ev - mean(ev))^2) / (length(ev) - k)
  xi <- 1.0
  phi <- if (p > 0) rep(0.3 / max(p, 1), p) else numeric(0)
  lam <- rep(1, nT)
  nu <- 8.0

  step_phi <- if (p > 0) rep(0.05, p) else numeric(0)
  step_nu <- 0.3
  acc_phi <- 0
  acc_nu <- 0
  n_prop <- 0L

  n_keep <- length(seq(from = burn, to = n_iter - 1L, by = thin))
  out <- list(
    beta = matrix(0, n_keep, k),
    sigma = numeric(n_keep),
    phi = matrix(0, n_keep, p),
    nu = numeric(n_keep),
    eps_last = matrix(0, n_keep, max(p, 1L)),
    z_imputed = matrix(0, n_keep, length(miss_idx)),
    loglik = NULL,
    ll_idx = integer(0),
    names = xnames
  )
  # Pointwise log-likelihood for WAIC: days that are observed and whose p
  # preceding calendar days are also observed.
  if (!is.null(ll_index)) {
    ll_idx <- as.integer(ll_index)
  } else if (p > 0) {
    cand <- seq.int(p + 1L, nT)
    keep <- vapply(cand, function(i) observed[i] && all(observed[(i - p):(i - 1L)]), logical(1))
    ll_idx <- cand[keep]
  } else {
    ll_idx <- which(observed)
  }
  out$loglik <- matrix(0, n_keep, length(ll_idx))
  out$ll_idx <- ll_idx

  path_every <- max(1L, n_keep %/% max(as.integer(n_paths), 1L))
  path_slots <- length(seq(from = 0L, to = n_keep - 1L, by = path_every))
  out$eps_path <- matrix(0, path_slots, nT)
  out$lam_path <- matrix(0, path_slots, nT)
  out$path_pars <- matrix(0, path_slots, 3L + p)
  sv <- 0L
  psv <- 0L

  log_prior_nu <- function(v) if (v > 2) -(v - 2) / nu_exp_mean else -Inf

  sum_log_lo <- 0
  sum_lo <- 0
  m_lo <- 0
  ll_nu <- function(v) {
    a <- 0.5 * v
    m_lo * (a * log(a) - lgamma(a)) + (a - 1) * sum_log_lo - a * sum_lo
  }

  log_post_phi <- function(ph) {
    if (!is_stationary(ph)) return(-Inf)
    e <- ar_resid(eps, ph)
    Vunit <- stationary_cov(ph, 1)
    logdet <- as.numeric(determinant(Vunit, logarithm = TRUE)$modulus)
    e0 <- eps[seq_len(length(ph))]
    q0 <- sum(e0 * solve(Vunit, e0))
    -0.5 * sum(lam[(length(ph) + 1L):length(lam)] * e^2) / sigma2 -
      0.5 * logdet - 0.5 * q0 / sigma2
  }

  for (it in 0:(n_iter - 1L)) {
    # ---------------- beta | phi, sigma2, lambda, z ----------------
    Xt <- quasi_difference(X, phi)
    zt <- as.numeric(quasi_difference(z, phi))
    w <- lam[(p + 1L):nT] / sigma2
    Xt_w <- Xt * rep(w, times = ncol(Xt))
    prec <- prior_prec + crossprod(Xt_w, Xt)
    rhs <- as.numeric(crossprod(Xt_w, zt))
    if (p > 0) {
      V0 <- stationary_cov(phi, sigma2)
      X0 <- X[seq_len(p), , drop = FALSE]
      z0 <- z[seq_len(p)]
      P0X <- solve(V0, X0)
      prec <- prec + crossprod(X0, P0X)
      rhs <- rhs + as.numeric(crossprod(X0, solve(V0, z0)))
    }
    prec <- (prec + t(prec)) / 2
    U <- chol(prec)
    beta <- backsolve(U, forwardsolve(t(U), rhs)) + backsolve(U, rnorm(k))
    eps <- as.numeric(z - X %*% beta)

    # ---------------- sigma2 | rest ----------------
    eta <- ar_resid(eps, phi)
    ss <- sum(lam[(p + 1L):nT] * eta^2)
    n_eff <- nT - p
    if (p > 0) {
      # The initial block contributes eps_{1:p}' (sigma2 V)^{-1} eps_{1:p},
      # and V is proportional to sigma2, so it scales out cleanly.
      Vunit <- stationary_cov(phi, 1)
      e0 <- eps[seq_len(p)]
      ss <- ss + sum(e0 * solve(Vunit, e0))
      n_eff <- n_eff + p
    }
    shape <- 0.5 * (nu_sigma + n_eff)
    rate <- nu_sigma / (2 * xi) + 0.5 * ss
    sigma2 <- rate / rgamma(1, shape = shape, scale = 1)
    xi <- (1 / sigma_A^2 + nu_sigma / sigma2) / rgamma(1, shape = 0.5 * (nu_sigma + 1), scale = 1)

    # ---------------- lambda | rest ----------------
    if (heavy) {
      eta <- ar_resid(eps, phi)
      shape_l <- 0.5 * (nu + 1)
      rate_l <- 0.5 * (nu + eta^2 / sigma2)
      lam[(p + 1L):nT] <- rgamma(nT - p, shape = shape_l, rate = rate_l)
      if (p > 0) lam[seq_len(p)] <- 1
    } else {
      lam[] <- 1
    }

    # ---------------- nu | lambda  (Metropolis on log(nu - 2)) ----------------
    if (heavy) {
      n_prop <- n_prop + 1L
      lo <- lam[(p + 1L):nT]
      prop <- 2 + exp(log(nu - 2) + step_nu * rnorm(1))
      sum_log_lo <- sum(log(lo))
      sum_lo <- sum(lo)
      m_lo <- length(lo)
      log_r <- ll_nu(prop) + log_prior_nu(prop) + log(prop - 2) -
        ll_nu(nu) - log_prior_nu(nu) - log(nu - 2)
      if (is.finite(log_r) && log(runif(1)) < log_r) {
        nu <- prop
        acc_nu <- acc_nu + 1
      }
    }

    # ---------------- phi | rest  (Metropolis) ----------------
    if (p > 0) {
      cur <- log_post_phi(phi)
      prop <- phi + step_phi * rnorm(p)
      prop_lp <- log_post_phi(prop)
      log_r <- prop_lp - cur
      if (is.finite(log_r) && log(runif(1)) < log_r) {
        phi <- prop
        acc_phi <- acc_phi + 1
      } else if (is.infinite(log_r) && log_r > 0) {
        phi <- prop
        acc_phi <- acc_phi + 1
      }
    }

    # ---------------- missing z_t | rest ----------------
    if (length(miss_idx) > 0) {
      for (i in miss_idx) {
        d_terms <- numeric(0)
        s_terms <- integer(0)
        if (i > p) {
          d_terms <- c(d_terms, 1)
          s_terms <- c(s_terms, i)
        }
        if (p > 0) {
          for (kk in seq_len(p)) {
            s <- i + kk
            if (s <= nT && s > p) {
              d_terms <- c(d_terms, -phi[kk])
              s_terms <- c(s_terms, s)
            }
          }
        }
        if (length(d_terms) == 0) next
        e_all <- ar_resid(eps, phi)
        grad <- 0
        curv <- 0
        for (j in seq_along(d_terms)) {
          s <- s_terms[j]
          dcoef <- d_terms[j]
          eta_s <- e_all[s - p]
          grad <- grad + lam[s] * eta_s * dcoef
          curv <- curv + lam[s] * dcoef * dcoef
        }
        if (i <= p && p > 0) {
          Vunit <- stationary_cov(phi, 1)
          Pi <- solve(Vunit)
          grad <- grad + sum(Pi[i, ] * eps[seq_len(p)])
          curv <- curv + Pi[i, i]
        }
        if (curv <= 0) next
        cond_mean_eps <- eps[i] - grad / curv
        eps[i] <- cond_mean_eps + sqrt(sigma2 / curv) * rnorm(1)
        z[i] <- sum(X[i, ] * beta) + eps[i]
      }
    }

    # ---------------- adapt during burn-in ----------------
    # n_prop increments only for heavy-tailed models, so the phi step is left
    # at its initial value when heavy = FALSE. That matches the Python sampler.
    if (it < burn && (it + 1L) %% 200L == 0L && n_prop > 0L) {
      if (p > 0) {
        r <- acc_phi / 200
        step_phi <- step_phi * exp((r - 0.3) * 1)
        step_phi <- pmin(pmax(step_phi, 1e-4), 1)
        acc_phi <- 0
      }
      if (heavy) {
        r <- acc_nu / 200
        step_nu <- step_nu * exp((r - 0.3) * 1)
        step_nu <- min(max(step_nu, 1e-3), 3)
        acc_nu <- 0
      }
      n_prop <- 0L
    }

    # ---------------- save ----------------
    if (it >= burn && (it - burn) %% thin == 0L && sv < n_keep) {
      sv <- sv + 1L
      out$beta[sv, ] <- beta
      out$sigma[sv] <- sqrt(sigma2)
      if (p > 0) {
        out$phi[sv, ] <- phi
        out$eps_last[sv, seq_len(p)] <- rev(eps[(nT - p + 1L):nT])
      } else {
        out$eps_last[sv, 1] <- 0
      }
      out$nu[sv] <- if (heavy) nu else Inf
      if (length(miss_idx) > 0) out$z_imputed[sv, ] <- z[miss_idx]
      e <- if (p > 0) ar_resid(eps, phi) else eps
      ee <- if (p > 0) e[ll_idx - p] else e[ll_idx]
      ll <- log_t(ee, sqrt(sigma2), nu, heavy)
      if (scale == "log") ll <- ll - z[ll_idx]
      out$loglik[sv, ] <- ll
      if ((sv - 1L) %% path_every == 0L && psv < path_slots) {
        psv <- psv + 1L
        out$eps_path[psv, ] <- eps
        out$lam_path[psv, ] <- lam
        out$path_pars[psv, 1] <- sqrt(sigma2)
        out$path_pars[psv, 2] <- if (heavy) nu else Inf
        if (p > 0) out$path_pars[psv, 4:(3L + p)] <- phi
      }
    }

    if (verbose && (it + 1L) %% 5000L == 0L) {
      message(sprintf(
        "    iter %d/%d  sigma=%.3f phi=%s nu=%.1f",
        it + 1L, n_iter, sqrt(sigma2), paste(round(phi, 3), collapse = " "), nu
      ))
    }
  }

  out$miss_idx <- miss_idx
  out$X <- X
  out$scale <- scale
  out$J <- J
  out$p <- p
  out$heavy <- heavy
  out
}

# ----------------------------------------------------------------------------
# multi-chain driver
# ----------------------------------------------------------------------------

stack_along_chains <- function(xs) {
  x1 <- xs[[1]]
  d <- dim(x1)
  nC <- length(xs)
  if (is.null(d)) {
    out <- matrix(0, nC, length(x1))
    for (i in seq_len(nC)) out[i, ] <- xs[[i]]
    return(out)
  }
  out <- array(0, c(nC, d[1], d[2]))
  if (d[2] == 0) return(out)
  for (i in seq_len(nC)) out[i, , ] <- xs[[i]]
  out
}

fit_chains <- function(y, observed, t, dow, n_chains = 4, seed = 1, ...) {
  chains <- lapply(seq_len(n_chains) - 1L, function(ch) {
    fit(y, observed, t, dow, seed = seed + 100 * ch, ...)
  })
  merged <- list()
  for (key in c("beta", "sigma", "phi", "nu", "eps_last", "z_imputed", "loglik")) {
    merged[[key]] <- stack_along_chains(lapply(chains, `[[`, key))
  }
  for (key in c("eps_path", "lam_path", "path_pars")) {
    merged[[key]] <- do.call(rbind, lapply(chains, `[[`, key))
  }
  for (key in c("names", "ll_idx", "miss_idx", "X", "scale", "J", "p", "heavy")) {
    merged[[key]] <- chains[[1]][[key]]
  }
  merged
}

common_ll_index <- function(observed, max_lag = 2) {
  # Days that are observed and whose max_lag preceding days are also observed.
  # Use one shared set across models so WAIC is computed on identical terms.
  # Indices are 1-based.
  observed <- as.logical(observed)
  nT <- length(observed)
  if (max_lag <= 0) return(which(observed))
  idx <- seq.int(max_lag + 1L, nT)
  keep <- vapply(idx, function(i) {
    observed[i] && all(observed[(i - max_lag):(i - 1L)])
  }, logical(1))
  idx[keep]
}

one_step_diagnostics <- function(merged, z_obs = NULL, observed) {
  # From the saved error paths, compute per-day one-step-ahead standardised
  # residuals and posterior predictive p-values under the fitted model.
  # Returns a list of length-T vectors (NA on days that are not evaluable).
  eps <- merged$eps_path
  pars <- merged$path_pars
  sigma <- pars[, 1]
  nu <- pars[, 2]
  p <- merged$p
  nT <- ncol(eps)
  S <- nrow(eps)
  phi <- if (p > 0) pars[, 4:(3L + p), drop = FALSE] else matrix(0, S, 0)

  pred <- matrix(0, S, nT)
  if (p > 0) {
    for (k in seq_len(p)) {
      pred[, (k + 1L):nT] <- pred[, (k + 1L):nT] + phi[, k] * eps[, seq_len(nT - k)]
    }
  }
  e <- eps - pred
  std <- e / sigma

  if (isTRUE(merged$heavy)) {
    cdf <- std
    cdf[] <- pt(as.numeric(std), df = rep(nu, times = nT))
  } else {
    cdf <- std
    cdf[] <- pnorm(as.numeric(std))
  }

  ok <- as.logical(observed)
  if (p > 0) {
    obs <- as.logical(observed)
    for (k in seq_len(p)) {
      ok[(k + 1L):nT] <- ok[(k + 1L):nT] & obs[seq_len(nT - k)]
    }
    ok[seq_len(p)] <- FALSE
  }

  list(
    std_resid_mean = ifelse(ok, colMeans(std), NA_real_),
    std_resid_sd = ifelse(ok, apply(std, 2, sd), NA_real_),
    ppp = ifelse(ok, colMeans(cdf), NA_real_),
    lam_mean = colMeans(merged$lam_path),
    evaluable = ok
  )
}

flatten <- function(merged, key) {
  # Collapse the (chain, draw, ...) axes into one. An empty trailing dimension
  # (p = 0 leaves phi empty) becomes a matrix with zero columns.
  a <- merged[[key]]
  d <- dim(a)
  if (length(d) == 2) {
    as.vector(t(a))
  } else if (length(d) == 3) {
    nC <- d[1]
    nD <- d[2]
    k <- d[3]
    if (k == 0) return(matrix(0, nC * nD, 0))
    out <- matrix(0, nC * nD, k)
    r <- 1L
    for (ch in seq_len(nC)) {
      sl <- array(a[ch, , , drop = FALSE], c(nD, k))
      out[r:(r + nD - 1L), ] <- sl
      r <- r + nD
    }
    out
  } else {
    stop("flatten() expected a matrix or a 3-way array")
  }
}

# ----------------------------------------------------------------------------
# convergence diagnostics (rank-normalised split-Rhat and bulk ESS)
# ----------------------------------------------------------------------------

rank_normalise <- function(x) {
  out <- x
  flat <- as.numeric(x)
  out[] <- qnorm((rank(flat, ties.method = "first") - 0.375) / (length(flat) + 0.25))
  out
}

split_rhat <- function(draws) {
  # draws: matrix with chains in rows and iterations in columns.
  n <- ncol(draws)
  half <- n %/% 2
  s <- rbind(draws[, seq_len(half), drop = FALSE], draws[, (half + 1L):(2L * half), drop = FALSE])
  n2 <- ncol(s)
  means <- rowMeans(s)
  W <- mean(apply(s, 1, var))
  B <- n2 * var(means)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  var_hat <- ((n2 - 1) * W + B) / n2
  sqrt(var_hat / W)
}

ess <- function(draws) {
  # Bulk effective sample size via Geyer's initial positive sequence.
  m <- nrow(draws)
  n <- ncol(draws)
  if (isTRUE(all.equal(as.numeric(draws), rep(draws[1], length(draws)), tolerance = 1e-5))) {
    return(as.numeric(m * n))
  }
  acov <- matrix(0, m, n)
  for (ch in seq_len(m)) {
    x <- draws[ch, ] - mean(draws[ch, ])
    f <- fft(c(x, rep(0, n)))
    ac <- Re(fft(f * Conj(f), inverse = TRUE)) / (2 * n)
    acov[ch, ] <- ac[seq_len(n)] / n
  }
  chain_mean <- rowMeans(draws)
  mean_var <- mean(acov[, 1]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n + (if (m > 1) var(chain_mean) else 0)
  rho_hat_even <- 1
  rho_hat_odd <- 1 - (mean_var - mean(acov[, 2])) / var_plus
  tau <- rho_hat_even + 2 * rho_hat_odd
  kk <- 1
  while (kk + 2 < n - 2 && (rho_hat_even + rho_hat_odd) > 0) {
    rho_hat_even <- 1 - (mean_var - mean(acov[, kk + 2])) / var_plus
    rho_hat_odd <- 1 - (mean_var - mean(acov[, kk + 3])) / var_plus
    if (rho_hat_even + rho_hat_odd > 0) tau <- tau + 2 * (rho_hat_even + rho_hat_odd)
    kk <- kk + 2
  }
  tau <- max(tau, 1 / log10(max(m * n, 11)))
  (m * n) / tau
}

diagnose <- function(draws) {
  # draws: matrix chain x iteration. Rank-normalised split-Rhat and bulk ESS.
  rn <- rank_normalise(draws)
  list(rhat = split_rhat(rn), ess = ess(rn))
}

# ----------------------------------------------------------------------------
# WAIC
# ----------------------------------------------------------------------------

waic <- function(loglik) {
  # loglik: array chain x draw x n_obs of pointwise log densities.
  d <- dim(loglik)
  nC <- d[1]
  nD <- d[2]
  n <- d[3]
  ll <- matrix(0, nC * nD, n)
  r <- 1L
  for (ch in seq_len(nC)) {
    ll[r:(r + nD - 1L), ] <- array(loglik[ch, , , drop = FALSE], c(nD, n))
    r <- r + nD
  }
  m <- apply(ll, 2, max)
  lppd_i <- log(colMeans(exp(sweep(ll, 2, m, "-")))) + m
  pwaic_i <- apply(ll, 2, var)
  elpd_i <- lppd_i - pwaic_i
  list(
    elpd_waic = sum(elpd_i),
    p_waic = sum(pwaic_i),
    waic = -2 * sum(elpd_i),
    se = sqrt(n * var(elpd_i)),
    n = n,
    elpd_i = elpd_i
  )
}

# ----------------------------------------------------------------------------
# forecasting
# ----------------------------------------------------------------------------

forecast <- function(
  merged, t_future, dow_future, t_center, eps_init = NULL,
  n_draws = NULL, draw_idx = NULL, seed = 7
) {
  # Simulate the posterior predictive distribution of y on future days.
  #
  # eps_init, if given, is (N, p) with columns (eps_origin, eps_origin-1, ...),
  # overriding the error state saved at the end of the fitted series. This is
  # what makes rolling-origin forecasts possible without refitting: given beta,
  # the error at any day with an observed response is just z_t - x_t'beta.
  #
  # Returns a list with
  #   y      (N, H) predictive draws on the kWh scale
  #   mu     (N, H) systematic component only, kWh scale
  #   loc    (N, H) conditional mean of the response on the MODELLING scale,
  #                 before the final innovation is added. Together with
  #                 (sigma, nu) this gives an exact mixture representation of
  #                 the predictive density, so log scores need no smoothing.
  #   sigma, nu, scale
  set.seed(seed)
  beta <- flatten(merged, "beta")
  sigma <- flatten(merged, "sigma")
  phi <- flatten(merged, "phi")
  nu <- flatten(merged, "nu")
  eps_last <- flatten(merged, "eps_last")
  if (is.null(dim(beta))) beta <- matrix(beta, nrow = 1)
  if (is.null(dim(phi))) phi <- matrix(phi, nrow = length(sigma))
  if (is.null(dim(eps_last))) eps_last <- matrix(eps_last, nrow = length(sigma))
  N <- nrow(beta)
  if (!is.null(draw_idx)) {
    sel <- as.integer(draw_idx)
  } else if (!is.null(n_draws) && n_draws < N) {
    sel <- sample.int(N, n_draws)
  } else {
    sel <- seq_len(N)
  }
  if (length(sel) != N) {
    beta <- beta[sel, , drop = FALSE]
    sigma <- sigma[sel]
    phi <- if (ncol(phi) == 0) phi[sel, , drop = FALSE] else phi[sel, , drop = FALSE]
    nu <- nu[sel]
    eps_last <- eps_last[sel, , drop = FALSE]
    N <- length(sel)
  }
  if (!is.null(eps_init) && nrow(as.matrix(eps_init)) != N) {
    stop(sprintf(
      "eps_init has %d rows but %d draws are in use; build it from the same draw_idx passed here.",
      nrow(as.matrix(eps_init)), N
    ))
  }

  p <- merged$p
  Xf <- build_design(t_future, dow_future, merged$J, t_center = t_center)$X
  H <- length(t_future)
  mu <- beta %*% t(Xf)

  if (p > 0) {
    hist <- if (is.null(eps_init)) {
      eps_last[, seq_len(p), drop = FALSE]
    } else {
      as.matrix(eps_init)
    }
  } else {
    hist <- matrix(0, N, 1)
  }

  z <- matrix(0, N, H)
  loc <- matrix(0, N, H)
  heavy <- merged$heavy
  for (h in seq_len(H)) {
    pred_eps <- if (p > 0) rowSums(phi * hist) else rep(0, N)
    loc[, h] <- mu[, h] + pred_eps
    lam <- if (heavy) rgamma(N, shape = nu / 2, scale = 2 / nu) else rep(1, N)
    eta <- rnorm(N) * sigma / sqrt(lam)
    e_new <- pred_eps + eta
    if (p > 0) {
      hist <- if (p == 1) cbind(e_new) else cbind(e_new, hist[, seq_len(p - 1L), drop = FALSE])
    }
    z[, h] <- mu[, h] + e_new
  }

  out <- list(sigma = sigma, nu = nu, scale = merged$scale, loc = loc, heavy = heavy, sel = sel)
  if (identical(merged$scale, "log")) {
    out$y <- exp(z)
    out$mu <- exp(mu)
  } else {
    out$y <- z
    out$mu <- mu
  }
  out
}

log_pred_density <- function(fc, y_actual) {
  # Exact-mixture log predictive density of the observed kWh values.
  # The final innovation is integrated analytically (t_nu or normal given the
  # draw), so this is a Rao-Blackwellised estimate rather than a kernel smooth
  # of the predictive sample. Returns a vector of length H.
  y_actual <- as.numeric(y_actual)
  H <- ncol(fc$loc)
  out <- rep(NA_real_, H)
  z_act <- if (identical(fc$scale, "log")) log(y_actual) else y_actual
  for (h in seq_len(H)) {
    if (!is.finite(z_act[h])) next
    ll <- log_t(z_act[h] - fc$loc[, h], fc$sigma, fc$nu, fc$heavy)
    if (identical(fc$scale, "log")) ll <- ll - z_act[h]
    m <- max(ll)
    out[h] <- m + log(mean(exp(ll - m)))
  }
  out
}

pit <- function(fc, y_actual) {
  # Probability integral transform: the predictive CDF at the actual value,
  # computed from the mixture rather than from sample ranks. Well-calibrated
  # forecasts give a uniform histogram.
  y_actual <- as.numeric(y_actual)
  z_act <- if (identical(fc$scale, "log")) log(y_actual) else y_actual
  H <- ncol(fc$loc)
  out <- rep(NA_real_, H)
  for (h in seq_len(H)) {
    if (!is.finite(z_act[h])) next
    s <- (z_act[h] - fc$loc[, h]) / fc$sigma
    out[h] <- if (isTRUE(fc$heavy)) mean(pt(s, df = fc$nu)) else mean(pnorm(s))
  }
  out
}

crps <- function(samples, y_actual) {
  # Exact CRPS of an ensemble forecast, one value per column of samples (N, H).
  y_actual <- as.numeric(y_actual)
  N <- nrow(samples)
  H <- ncol(samples)
  out <- rep(NA_real_, H)
  i <- 0:(N - 1L)
  for (h in seq_len(H)) {
    if (!is.finite(y_actual[h])) next
    x <- sort(samples[, h])
    term1 <- mean(abs(x - y_actual[h]))
    term2 <- 2 * sum((2 * i - N + 1) * x) / (N * N)
    out[h] <- term1 - 0.5 * term2
  }
  out
}
