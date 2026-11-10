#!/usr/bin/env Rscript
# Step 8: the frequentist alternative for Report section 7. Fits the final
# model of src/06_fit_sqrt.R (sqrt kWh, trend, two harmonics, day of week,
# AR(1) errors with Student-t innovations) by maximum likelihood, then repeats
# the section 5 summaries and the holdout forecasts with the MLE. Takes about
# a minute.
#
#   Rscript src/08_frequentist.R
#
# Writes results/08_freq/summary.rds and table_freq_vs_bayes.csv, and
# R_figures/10_freq_vs_bayes.png. The Bayesian side of the comparison is read
# from results/06_sqrt/summary.rds.
#
# The likelihood conditions on the first observed day of each run, so the 23
# missing days and the first day after each of the 8 gaps contribute no term
# (1,408 terms on the full data). Wald intervals use the observed information
# on the scale (beta, log sigma, atanh phi, log(nu - 2)), which keeps them
# inside the parameter space. kWh quantities are nonlinear in beta, so their
# intervals come from draws of the asymptotic normal distribution of the MLE
# (Krinsky and Robb, 1986), passed through the same functions that
# src/07_sqrt_results.R applies to posterior draws.

library(here)
source(here("src", "07_sqrt_results.R"))

OUT <- here("results", "08_freq")
N_SIM <- 20000L
TRAIN_END <- as.Date("2009-11-30")
ROLL_H <- 60L
ROLL_STEP <- 3L
GREEN <- "#148f77"

# ---------------------------------------------------------------------------
# maximum likelihood
# ---------------------------------------------------------------------------

fit_mle <- function(kwh, observed, t, dow, J, t_center = mean(t)) {
  des <- build_design(t, dow, J, t_center = t_center)
  X <- des$X
  k <- ncol(X)
  z <- ifelse(observed, sqrt(kwh), NA_real_)
  idx <- common_ll_index(observed, max_lag = 1)
  Xd <- X[idx, , drop = FALSE]
  Xl <- X[idx - 1L, , drop = FALSE]
  zd <- z[idx]
  zl <- z[idx - 1L]

  unpack <- function(th) {
    list(beta = th[1:k], sigma = exp(th[k + 1]), phi = tanh(th[k + 2]), nu = 2 + exp(th[k + 3]))
  }
  parts <- function(th) {
    p <- unpack(th)
    el <- zl - as.numeric(Xl %*% p$beta)
    e <- (zd - as.numeric(Xd %*% p$beta)) - p$phi * el
    c(p, list(e = e, el = el))
  }
  nll <- function(th) {
    q <- parts(th)
    -sum(log_t(q$e, q$sigma, q$nu, TRUE))
  }
  grad <- function(th) {
    q <- parts(th)
    s2 <- q$sigma^2
    u <- q$e^2 / (q$nu * s2)
    w <- (q$nu + 1) / (q$nu + q$e^2 / s2)
    Xt <- Xd - q$phi * Xl
    d_beta <- colSums(Xt * (w * q$e / s2))
    d_phi <- sum(w * q$e * q$el / s2) * (1 - q$phi^2)
    d_lsig <- sum(-1 + w * q$e^2 / s2)
    d_nu <- sum(0.5 * digamma(0.5 * (q$nu + 1)) - 0.5 * digamma(0.5 * q$nu) - 0.5 / q$nu -
                  0.5 * log1p(u) + 0.5 * (q$nu + 1) * u / (q$nu * (1 + u)))
    -c(d_beta, d_lsig, d_phi, d_nu * (q$nu - 2))
  }

  ok <- observed
  b0 <- qr.solve(X[ok, , drop = FALSE], z[ok])
  th <- c(b0, log(sd(z[ok] - X[ok, ] %*% b0)), atanh(0.4), log(6))
  for (i in 1:3) {
    o <- optim(th, nll, grad, method = "BFGS", control = list(maxit = 5000, reltol = 1e-15))
    th <- o$par
  }
  H <- optimHess(th, nll, grad)
  V <- solve(H)
  list(
    theta = th, V = V, se = sqrt(diag(V)), loglik = -o$value, n = length(idx),
    n_par = length(th), convergence = o$convergence, max_grad = max(abs(grad(th))),
    names = c(des$names, "log_sigma", "atanh_phi", "log_nu_minus_2"),
    beta_names = des$names, k = k, X = X, z = z, J = J, t_center = t_center, unpack = unpack,
    nll = nll, grad = grad
  )
}

# rows are draws of theta from N(theta_hat, V)
sim_theta <- function(m, n, seed) {
  set.seed(seed)
  Z <- matrix(rnorm(n * m$n_par), n, m$n_par)
  sweep(Z %*% chol(m$V), 2, m$theta, `+`)
}

wald <- function(m, name, f = identity) {
  i <- match(name, m$names)
  f(m$theta[i] + c(0, -1, 1) * qnorm(0.975) * m$se[i])
}

# ---------------------------------------------------------------------------
# kWh quantities, per row of b (same definitions as src/07_sqrt_results.R)
# ---------------------------------------------------------------------------

derived <- function(b, nm) {
  yr <- seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "day")
  win <- seq(as.Date("2009-07-01"), as.Date("2010-06-30"), by = "day")
  M2 <- kwh_curve(b, nm, yr, 2)
  annual <- rowMeans(M2)
  peak <- apply(kwh_curve(b, nm, win, 2), 1, max)
  trough <- apply(M2, 1, min)
  add2 <- apply(abs(M2 - kwh_curve(b, nm, yr, 1)), 1, max)
  dl <- dow_full(b, nm)
  A <- b[, match("intercept", nm)] + seas_matrix(b, nm, yr, 2)
  slope <- 2 * b[, match("trend_per_year", nm)] * (rowMeans(A) + rowMeans(dl))
  a1 <- rowMeans(A)
  a2 <- rowMeans(A^2)
  E <- a2 + 2 * dl * a1 + dl^2
  P <- E - rowMeans(E)
  wkend <- rowMeans(P[, 6:7, drop = FALSE]) - rowMeans(P[, 1:5, drop = FALSE])
  cbind(annual = annual, slope_kwh = slope, slope_pct = 100 * slope / annual, peak_kwh = peak,
        trough_kwh = trough, swing = peak - trough, add2 = add2, wkend = wkend, thu = P[, 4])
}

# ---------------------------------------------------------------------------
# forecasts: wrap parameter rows so that forecast() in model.R can use them
# ---------------------------------------------------------------------------

as_draws <- function(th, m, eps_last) {
  N <- nrow(th)
  k <- m$k
  list(
    beta = array(th[, 1:k], c(1, N, k)), sigma = matrix(exp(th[, k + 1]), 1),
    phi = array(tanh(th[, k + 2]), c(1, N, 1)), nu = matrix(2 + exp(th[, k + 3]), 1),
    eps_last = array(eps_last, c(1, N, 1)), p = 1L, J = m$J, heavy = TRUE, scale = "sqrt"
  )
}

holdout <- function(d, n_single = 4000L, n_roll = 1500L) {
  train <- d[d$date <= TRAIN_END, ]
  test <- d[d$date > TRAIN_END, ]
  rownames(test) <- NULL
  t_center <- mean(train$t)
  m <- fit_mle(train$kwh, train$observed, train$t, train$dow, J = 2, t_center = t_center)
  last <- nrow(train)
  stopifnot(train$observed[last])
  Xtest <- build_design(test$t, test$dow, 2, t_center = t_center)$X
  z_test <- sqrt(test$kwh)

  # plug-in: every path uses the MLE; "with estimation error": each path
  # first draws theta from the asymptotic normal distribution of the MLE
  sets <- list(
    plugin = function(n) matrix(m$theta, n, m$n_par, byrow = TRUE),
    est_err = function(n) sim_theta(m, n, seed = 81)
  )
  out <- list(train_phi = m$unpack(m$theta)$phi, n_train = sum(train$observed))
  for (nm in names(sets)) {
    th <- sets[[nm]](n_single)
    eps_last <- m$z[last] - as.numeric(th[, 1:m$k] %*% m$X[last, ])
    fc <- forecast(as_draws(th, m, eps_last), test$t, test$dow, t_center, seed = 99)
    q <- apply(fc$y, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    sg <- data.frame(kwh = test$kwh, lo95 = q[1, ], lo50 = q[2, ], med = q[3, ], hi50 = q[4, ],
                     hi95 = q[5, ], crps = crps(fc$y, test$kwh))[test$observed, ]
    single <- c(c95 = mean(sg$kwh >= sg$lo95 & sg$kwh <= sg$hi95),
                c50 = mean(sg$kwh >= sg$lo50 & sg$kwh <= sg$hi50),
                crps = mean(sg$crps), mae = mean(abs(sg$kwh - sg$med)), w95 = mean(sg$hi95 - sg$lo95))

    # rolling origins: the same rule as src/06_fit_sqrt.R, so the days match
    th <- sets[[nm]](n_roll)
    eps_test <- matrix(z_test, n_roll, length(z_test), byrow = TRUE) - th[, 1:m$k] %*% t(Xtest)
    origins <- seq(2L, nrow(test) - 1L, by = ROLL_STEP)
    origins <- origins[test$observed[origins] & test$observed[origins - 1L]]
    rows <- vector("list", length(origins))
    for (oi in seq_along(origins)) {
      o <- origins[oi]
      idx <- (o + 1L):min(o + ROLL_H, nrow(test))
      fcr <- forecast(as_draws(th, m, eps_test[, o]), test$t[idx], test$dow[idx], t_center,
                      eps_init = eps_test[, o, drop = FALSE], seed = 7L * oi)
      qq <- apply(fcr$y, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975))
      rows[[oi]] <- data.frame(
        origin = o, h = seq_along(idx), observed = test$observed[idx], kwh = test$kwh[idx],
        lo95 = qq[1, ], lo50 = qq[2, ], hi50 = qq[3, ], hi95 = qq[4, ], crps = crps(fcr$y, test$kwh[idx])
      )
    }
    r <- do.call(rbind, rows)
    r <- r[r$observed, ]
    r$w95 <- r$hi95 - r$lo95
    r$c95 <- r$kwh >= r$lo95 & r$kwh <= r$hi95
    r$c50 <- r$kwh >= r$lo50 & r$kwh <= r$hi50
    out[[nm]] <- list(single = single, byh = aggregate(cbind(w95, c95, c50, crps) ~ h, data = r, FUN = mean),
                      n_origins = length(origins))
  }
  out
}

# ---------------------------------------------------------------------------

main <- function() {
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  d <- read_daily(here("data", "daily.csv"))
  B <- readRDS(file.path(IN, "summary.rds"))
  S <- list()

  # ---------------- full-data fits: J = 2 (final) and J = 1 ----------------
  m2 <- fit_mle(d$kwh, d$observed, d$t, d$dow, J = 2)
  m1 <- fit_mle(d$kwh, d$observed, d$t, d$dow, J = 1)
  stopifnot(m2$convergence == 0, m1$convergence == 0, m2$n == m1$n)
  lr <- 2 * (m2$loglik - m1$loglik)
  S$fit <- list(
    n_terms = m2$n, loglik = c(J2 = m2$loglik, J1 = m1$loglik), max_grad = m2$max_grad,
    theta = setNames(m2$theta, m2$names), se = setNames(m2$se, m2$names),
    lrt_J2_vs_J1 = c(stat = lr, df = 2, p = pchisq(lr, 2, lower.tail = FALSE)),
    aic = c(J2 = 2 * m2$n_par - 2 * m2$loglik, J1 = 2 * m1$n_par - 2 * m1$loglik)
  )

  # ---------------- parameters and kWh quantities ----------------
  b_hat <- matrix(m2$theta[1:m2$k], 1, dimnames = list(NULL, m2$beta_names))
  point <- derived(b_hat, m2$beta_names)[1, ]
  th <- sim_theta(m2, N_SIM, seed = 8)
  sims <- derived(th[, 1:m2$k, drop = FALSE], m2$beta_names)
  ci <- apply(sims, 2, q95)
  est <- lapply(colnames(sims), function(q) c(est = point[[q]], lo = ci[1, q], hi = ci[2, q]))
  names(est) <- colnames(sims)
  est$phi <- wald(m2, "atanh_phi", tanh)
  est$sigma <- wald(m2, "log_sigma", exp)
  est$nu <- wald(m2, "log_nu_minus_2", function(x) 2 + exp(x))
  zt <- m2$theta[match("trend_per_year", m2$names)] / m2$se[match("trend_per_year", m2$names)]
  S$est <- est
  # one-sided p is for H0: beta1 >= 0 against a decline
  S$trend_test <- c(z = zt, p_two_sided = 2 * pnorm(-abs(zt)), p_one_sided = pnorm(zt))
  S$sim_vs_wald_beta_sd <- max(abs(apply(th[, 1:m2$k], 2, sd) / m2$se[1:m2$k] - 1))

  # ---------------- holdout forecasts ----------------
  S$holdout <- holdout(d)
  bfc <- B$fc

  # ---------------- comparison table ----------------
  bp <- B$post
  f1 <- function(s, dp = 1) sprintf(paste0("%.", dp, "f [%.", dp, "f, %.", dp, "f]"), s[1], s[2], s[3])
  pct <- function(x) sprintf("%.1f%%", 100 * x)
  hp <- S$holdout$plugin
  roll_cov <- function(b) paste(pct(mean(b$c95)), "/", pct(mean(b$c50)))
  widths <- function(b) sprintf("%.1f / %.1f", b$w95[1], mean(b$w95[b$h >= 3]))
  tab <- data.frame(
    Quantity = c(
      "Typical-day consumption, annual mean (kWh/day)",
      "Change in annual mean (kWh/day per year)",
      "Evidence about the trend",
      "Winter peak (kWh/day)",
      "Summer trough (kWh/day)",
      "Largest change from adding the second harmonic (kWh/day)",
      "Second harmonic: model comparison",
      "Weekend minus weekday (kWh/day)",
      "AR(1) coefficient $\\phi$",
      "Degrees of freedom $\\nu$",
      "Innovation scale $\\sigma$ ($\\sqrt{\\text{kWh}}$)",
      "Rolling forecasts, 1-60 days ahead: 95% / 50% interval coverage",
      "12-month forecast: mean absolute error (kWh/day)",
      "Mean 95% interval width, 1 day / 3-60 days ahead (kWh)"
    ),
    Bayesian = c(
      f1(bp$annual), f1(bp$slope_kwh, 2),
      sprintf("$P(\\text{decline} \\mid \\text{data}) = %.2f$", bp$p_decline),
      f1(bp$peak_kwh), f1(bp$trough_kwh), f1(bp$add2),
      sprintf("WAIC gain %.1f (SE %.1f)", B$waic$J2_vs_J1[["diff"]], B$waic$J2_vs_J1[["se"]]),
      f1(bp$wkend), f1(bp$phi, 2), f1(bp$nu), f1(bp$sigma, 2),
      roll_cov(bfc$byh), sprintf("%.1f", bfc$single[["mae"]]), widths(bfc$byh)
    ),
    Frequentist = c(
      f1(est$annual), f1(est$slope_kwh, 2),
      sprintf("$p = %.2f$ two-sided, $%.2f$ one-sided (decline)",
              S$trend_test[["p_two_sided"]], S$trend_test[["p_one_sided"]]),
      f1(est$peak_kwh), f1(est$trough_kwh), f1(est$add2),
      sprintf("LRT $\\chi^2_2 = %.1f$, $p %s$; AIC lower by %.1f", lr,
              if (S$fit$lrt_J2_vs_J1[["p"]] < 0.001) "< 0.001" else sprintf("= %.3f", S$fit$lrt_J2_vs_J1[["p"]]),
              S$fit$aic[["J1"]] - S$fit$aic[["J2"]]),
      f1(est$wkend), f1(est$phi, 2), f1(est$nu), f1(est$sigma, 2),
      roll_cov(hp$byh), sprintf("%.1f", hp$single[["mae"]]), widths(hp$byh)
    ),
    check.names = FALSE
  )
  names(tab)[2:3] <- c("Bayesian: posterior mean [95% CrI]", "Frequentist: MLE [95% CI]")
  write.csv(tab, file.path(OUT, "table_freq_vs_bayes.csv"), row.names = FALSE)

  # ---------------- figure ----------------
  keys <- c("annual", "peak_kwh", "trough_kwh", "swing", "wkend", "add2")
  labs <- c("Annual mean", "Winter peak", "Summer trough", "Peak-to-trough swing",
            "Weekend premium", "Second-harmonic effect")
  png(here("R_figures", "10_freq_vs_bayes.png"), width = 12, height = 4.2, units = "in", res = 140)
  layout(matrix(1:3, nrow = 1), widths = c(1.35, 1, 1))
  par(mar = c(4.5, 10.5, 2.5, 1))
  yy <- rev(seq_along(keys))
  bb <- t(sapply(keys, function(q) bp[[q]]))
  ff <- t(sapply(keys, function(q) est[[q]]))
  plot(NA, xlim = range(c(bb, ff)), ylim = c(0.5, length(keys) + 0.5), yaxt = "n", ylab = "",
       xlab = "kWh/day", main = "(a) Estimates with 95% intervals", las = 1)
  axis(2, at = yy, labels = labs, las = 1, cex.axis = 0.9)
  segments(bb[, 2], yy + 0.12, bb[, 3], yy + 0.12, col = NAV, lwd = 2)
  points(bb[, 1], yy + 0.12, pch = 16, col = NAV)
  segments(ff[, 2], yy - 0.12, ff[, 3], yy - 0.12, col = RED, lwd = 2)
  points(ff[, 1], yy - 0.12, pch = 17, col = RED)
  legend("bottomright", legend = c("Bayesian (posterior mean, CrI)", "Frequentist (MLE, CI)"),
         col = c(NAV, RED), pch = c(16, 17), lwd = 2, bty = "n", cex = 0.85)
  par(mar = c(4.5, 4.5, 2.5, 1))
  he <- S$holdout$est_err
  plot(bfc$byh$h, bfc$byh$w95, type = "l", lwd = 2, col = NAV, ylim = c(0, 1.1 * max(bfc$byh$w95, hp$byh$w95)),
       xlab = "forecast horizon (days ahead)", ylab = "mean width of 95% interval (kWh)",
       main = "(b) Predictive uncertainty by horizon", las = 1)
  lines(hp$byh$h, hp$byh$w95, lwd = 2, col = RED, lty = 2)
  lines(he$byh$h, he$byh$w95, lwd = 1.5, col = GREEN, lty = 3)
  legend("bottomright", legend = c("Bayesian predictive", "frequentist plug-in", "plug-in + estimation error"),
         col = c(NAV, RED, GREEN), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.85)
  bins <- cut(bfc$byh$h, c(0, 1, 7, 14, 30, 60), labels = c("1", "2-7", "8-14", "15-30", "31-60"))
  cov_bin <- function(b, key) 100 * tapply(b[[key]], bins, mean)
  xb <- seq_along(levels(bins))
  plot(NA, xlim = c(0.5, length(xb) + 0.5), ylim = c(35, 100), xaxt = "n", las = 1,
       xlab = "forecast horizon (days ahead)", ylab = "empirical coverage (%)",
       main = "(c) Calibration by horizon")
  axis(1, at = xb, labels = levels(bins))
  abline(h = c(95, 50), lty = 3)
  points(xb - 0.1, cov_bin(bfc$byh, "c95"), pch = 16, col = NAV)
  points(xb + 0.1, cov_bin(hp$byh, "c95"), pch = 17, col = RED)
  points(xb - 0.1, cov_bin(bfc$byh, "c50"), pch = 1, col = NAV)
  points(xb + 0.1, cov_bin(hp$byh, "c50"), pch = 2, col = RED)
  text(length(xb) + 0.45, c(95, 50), c("95%", "50%"), pos = 3, cex = 0.8)
  legend("right", legend = c("Bayesian predictive", "frequentist plug-in", "filled: 95% interval",
                             "open: 50% interval"),
         col = c(NAV, RED, "grey30", "grey30"), pch = c(16, 17, 16, 1), bty = "n", cex = 0.85)
  dev.off()

  S$table <- tab
  saveRDS(S, file.path(OUT, "summary.rds"))
  invisible(S)
}

if (sys.nframe() == 0L) {
  t0 <- proc.time()
  S <- main()
  print(S$table, right = FALSE)
  str(S[c("fit", "trend_test", "sim_vs_wald_beta_sd")], give.attr = FALSE)
  cat(sprintf("done in %.0f s\n", (proc.time() - t0)[["elapsed"]]))
}
