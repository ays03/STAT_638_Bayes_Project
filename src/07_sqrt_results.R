#!/usr/bin/env Rscript
# Step 7: tables, figures and numbers for Report sections 4 and 5, computed
# from the fits saved by src/06_fit_sqrt.R. Takes under a minute.
#
#   Rscript src/07_sqrt_results.R
#
# Writes results/06_sqrt/summary.rds (every number quoted in the text),
# results/06_sqrt/table_convergence.csv, table_posterior.csv, episodes.csv,
# and R_figures/06_trace.png, 07_season_week.png, 08_forecast.png, 09_anomalies.png.
#
# All kWh quantities are posterior draws transformed draw by draw. On the
# square-root scale the AR(1) + t error is symmetric about zero, so
# (x_t' beta)^2 is the median daily consumption given the parameters; we call
# it the typical-day consumption.

library(here)
source(here("src", "model.R"))

IN <- here("results", "06_sqrt")
FIG <- here("R_figures")
ORIGIN <- as.Date("2006-12-17")   # t = 1
DAYNAMES <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
CHAIN_COLS <- c("#1f4e79", "#b03a2e", "#148f77", "#b9770e")
NAV <- "#1f4e79"
RED <- "#b03a2e"
WIN <- 14L   # window, in days, of the block-level predictive check
read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}
load_fit <- function(name) readRDS(file.path(IN, paste0(name, ".rds")))
q95 <- function(v) unname(quantile(v, c(0.025, 0.975)))
summ <- function(v) c(mean = mean(v), lo = q95(v)[1], hi = q95(v)[2])

# ---------------------------------------------------------------------------
# convergence
# ---------------------------------------------------------------------------

param_mats <- function(f) {
  out <- list()
  for (i in seq_along(f$names)) out[[f$names[i]]] <- f$beta[, , i]
  out$sigma <- f$sigma
  for (i in seq_len(f$p)) out[[paste0("phi_", i)]] <- f$phi[, , i]
  out$nu <- f$nu
  out
}

diag_group <- function(mats) {
  g <- lapply(mats, diagnose)
  c(rhat = max(vapply(g, function(x) max(x$rhat, x$rhat_folded), numeric(1))),
    ess = min(vapply(g, `[[`, numeric(1), "ess")),
    ess_tail = min(vapply(g, `[[`, numeric(1), "ess_tail")))
}

worst_of <- function(f) {
  m <- param_mats(f)
  m <- c(m, lapply(seq_len(dim(f$z_imputed)[3]), function(j) f$z_imputed[, , j]))
  diag_group(m)
}

# ---------------------------------------------------------------------------
# kWh-scale summaries of the regression part
# ---------------------------------------------------------------------------

seas_matrix <- function(b, nm, dates, J) {
  tt <- as.numeric(dates - ORIGIN) + 1
  S <- matrix(0, nrow(b), length(dates))
  for (j in seq_len(J)) {
    S <- S + outer(b[, match(paste0("sin_h", j), nm)], sin(2 * pi * j * tt / 365)) +
      outer(b[, match(paste0("cos_h", j), nm)], cos(2 * pi * j * tt / 365))
  }
  S
}

dow_full <- function(b, nm) cbind(0, b[, match(paste0("dow_", DAYNAMES[-1]), nm), drop = FALSE])

# typical-day kWh on each date at the record midpoint, averaged over weekdays
kwh_curve <- function(b, nm, dates, J) {
  dl <- dow_full(b, nm)
  A <- b[, match("intercept", nm)] + seas_matrix(b, nm, dates, J)
  A^2 + 2 * A * rowMeans(dl) + rowMeans(dl^2)
}

main <- function() {
  d <- read_daily(here("data", "daily.csv"))
  f <- load_fit("main")
  f1 <- load_fit("J1")
  f2 <- load_fit("ar2")
  fk <- load_fit("kwh")
  ho <- load_fit("holdout")
  ho2 <- load_fit("holdout_ar2")
  S <- list()

  # ---------------- Section 4: convergence table ----------------
  n_chain <- dim(f$beta)[1]
  n_draw <- dim(f$beta)[2]
  pm <- param_mats(f)
  acc <- setNames(colMeans(f$acc_rate), c("phi", "nu"))
  imp <- lapply(seq_len(dim(f$z_imputed)[3]), function(j) f$z_imputed[, , j])
  groups <- list(
    list("Regression coefficients (12)", "Gibbs: multivariate normal", pm[f$names], NA),
    list("Innovation scale $\\sigma$", "Gibbs: inverse gamma", pm["sigma"], NA),
    list("Autoregressive coefficient $\\phi$", "Random-walk Metropolis", pm["phi_1"], acc[["phi"]]),
    list("Degrees of freedom $\\nu$", "Random-walk Metropolis on $\\log(\\nu - 2)$", pm["nu"], acc[["nu"]]),
    list(sprintf("Missing days (%d)", length(imp)), "Gibbs: normal", imp, NA)
  )
  conv <- do.call(rbind, lapply(groups, function(g) {
    dg <- diag_group(g[[3]])
    data.frame(
      Parameters = g[[1]], Update = g[[2]],
      Acceptance = if (is.na(g[[4]])) "exact draw" else sprintf("%.0f%%", 100 * g[[4]]),
      `Worst R-hat` = sprintf("%.3f", dg[["rhat"]]),
      `Min bulk ESS` = sprintf("%.0f", dg[["ess"]]),
      `Min tail ESS` = sprintf("%.0f", dg[["ess_tail"]]),
      check.names = FALSE
    )
  }))
  write.csv(conv, file.path(IN, "table_convergence.csv"), row.names = FALSE)
  all_main <- worst_of(f)
  S$conv <- list(
    n_chain = n_chain, n_draw_per_chain = n_draw, n_total = n_chain * n_draw,
    mcmc = f$mcmc, acc = acc, worst = all_main,
    mcse_max_frac = 1 / sqrt(all_main[["ess"]]),
    others = rbind(J1 = worst_of(f1), ar2 = worst_of(f2), kwh = worst_of(fk)),
    holdout = rbind(ar1 = ho$conv, ar2 = ho2$conv)
  )

  # ---------------- trace figure (warm-up shaded) ----------------
  show <- list(
    list(key = "trend_per_year", lab = expression(beta[1] ~ "(trend)")),
    list(key = "sigma", lab = expression(sigma)),
    list(key = "phi_1", lab = expression(phi)),
    list(key = "nu", lab = expression(nu))
  )
  warm_names <- c("sigma", "nu", paste0("phi_", seq_len(f$p)), f$names)
  png(file.path(FIG, "06_trace.png"), width = 12, height = 8, units = "in", res = 140)
  layout(matrix(1:8, ncol = 2, byrow = TRUE), widths = c(3, 1))
  par(mar = c(3, 4.5, 1.5, 1), oma = c(1.5, 0, 0, 0))
  for (s in show) {
    warm <- f$warmup[, , match(s$key, warm_names)]
    post <- pm[[s$key]]
    both <- cbind(warm, post)
    xs <- (seq_len(ncol(both)) - 1) * f$mcmc$thin
    plot(NA, xlim = range(xs), ylim = range(both), xlab = "", ylab = s$lab, las = 1)
    rect(0, par("usr")[3], f$mcmc$burn, par("usr")[4], col = "grey92", border = NA)
    for (ch in seq_len(n_chain)) lines(xs, both[ch, ], col = adjustcolor(CHAIN_COLS[ch], 0.6), lwd = 0.4)
    if (s$key == show[[1]]$key) {
      text(f$mcmc$burn / 2, par("usr")[4], "burn-in (discarded)", pos = 1, cex = 0.8)
      legend("topright", legend = paste("chain", 1:n_chain), col = CHAIN_COLS, lwd = 1.5,
             horiz = TRUE, bty = "n", cex = 0.75)
    }
    dens <- lapply(seq_len(n_chain), function(ch) density(post[ch, ]))
    plot(NA, xlim = range(post), ylim = c(0, max(vapply(dens, function(z) max(z$y), 0))),
         xlab = "", ylab = "", yaxt = "n")
    for (ch in seq_len(n_chain)) lines(dens[[ch]], col = CHAIN_COLS[ch], lwd = 1.2)
  }
  mtext("iteration", side = 1, outer = TRUE, adj = 0.37, line = 0.3, cex = 0.85)
  mtext("kept draws, by chain", side = 1, outer = TRUE, adj = 0.93, line = 0.3, cex = 0.85)
  dev.off()

  # ---------------- Section 5: posterior quantities ----------------
  b <- flatten(f, "beta")
  nm <- f$names
  yr <- seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "day")
  win <- seq(as.Date("2009-07-01"), as.Date("2010-06-30"), by = "day")
  M2 <- kwh_curve(b, nm, yr, 2)
  annual <- rowMeans(M2)
  Mw <- kwh_curve(b, nm, win, 2)
  peak_kwh <- apply(Mw, 1, max)
  peak_day <- win[apply(Mw, 1, which.max)]
  trough_kwh <- apply(M2, 1, min)
  trough_day <- yr[apply(M2, 1, which.min)]
  pk_q <- sort(peak_day)[round(c(0.025, 0.5, 0.975) * length(peak_day))]
  tr_q <- sort(trough_day)[round(c(0.025, 0.5, 0.975) * length(trough_day))]

  # what the second harmonic adds, inside the J = 2 posterior
  M1_in <- kwh_curve(b, nm, yr, 1)
  add2 <- apply(abs(M2 - M1_in), 1, max)
  amp <- function(j) sqrt(b[, match(paste0("sin_h", j), nm)]^2 + b[, match(paste0("cos_h", j), nm)]^2)
  ratio21 <- amp(2) / amp(1)
  b1 <- flatten(f1, "beta")
  M1_model <- kwh_curve(b1, f1$names, yr, 1)

  paired <- function(a, z) {
    di <- a$waic$elpd_i - z$waic$elpd_i
    c(diff = sum(di), se = sqrt(length(di) * var(di)))
  }
  S$waic <- list(
    J2_vs_J1 = paired(f, f1), ar2_vs_ar1 = paired(f2, f), sqrt_vs_kwh = paired(f, fk),
    p_waic = c(main = f$waic$p_waic, J1 = f1$waic$p_waic, ar2 = f2$waic$p_waic, kwh = fk$waic$p_waic),
    n = f$waic$n
  )

  # trend: the change per year in the annual-mean typical day is 2 * beta1 * mean(x'beta)
  dl <- dow_full(b, nm)
  A <- b[, match("intercept", nm)] + seas_matrix(b, nm, yr, 2)
  b_tr <- b[, match("trend_per_year", nm)]
  slope_kwh <- 2 * b_tr * (rowMeans(A) + rowMeans(dl))
  slope_pct <- 100 * slope_kwh / annual
  years <- (nrow(d) - 1) / 365

  # weekly profile in kWh, averaged over the year
  a1 <- rowMeans(A)
  a2 <- rowMeans(A^2)
  E <- a2 + 2 * dl * a1 + dl^2
  P <- E - rowMeans(E)
  wkend <- rowMeans(P[, 6:7]) - rowMeans(P[, 1:5])

  phi <- as.numeric(f$phi)
  zi <- flatten(f, "z_imputed")
  ki <- pmax(zi, 0)^2

  S$post <- list(
    annual = summ(annual),
    slope_kwh = summ(slope_kwh), slope_pct = summ(slope_pct),
    total_change = summ(years * slope_kwh), years = years,
    p_decline = mean(b_tr < 0), p_rope1 = mean(abs(slope_pct) < 1),
    peak_kwh = summ(peak_kwh), peak_day = pk_q, trough_kwh = summ(trough_kwh), trough_day = tr_q,
    swing = summ(peak_kwh - trough_kwh), peak_ratio = summ(peak_kwh / trough_kwh),
    amp_ratio = summ(ratio21), add2 = summ(add2), p_add2_gt1 = mean(add2 > 1),
    J1_vs_J2_curve_maxdiff = max(abs(colMeans(M2) - colMeans(M1_model))),
    J1_vs_J2_where = yr[which.max(abs(colMeans(M2) - colMeans(M1_model)))],
    dow = lapply(setNames(seq_along(DAYNAMES), DAYNAMES), function(k) summ(P[, k])),
    wkend = summ(wkend), p_wkend = mean(wkend > 0),
    p_thu_lowest = mean(apply(P, 1, which.min) == 4),
    sigma = summ(as.numeric(f$sigma)), phi = summ(phi), nu = summ(as.numeric(f$nu)),
    sd_ratio = summ(1 / sqrt(1 - phi^2)),
    imputed_sd = range(apply(ki, 2, sd)), imputed_width = range(apply(ki, 2, function(v) diff(q95(v))))
  )

  fmt3 <- function(s, dp = 1) sprintf(paste0("%.", dp, "f [%.", dp, "f, %.", dp, "f]"), s[1], s[2], s[3])
  pd <- S$post
  dfmt <- function(x) sub("^0", "", format(x, "%d %b"))
  post_tab <- data.frame(
    Question = c("Level", "Long-term trend", "", "Annual cycle", "", "", "",
                 "Weekly pattern", "", "Day-to-day memory", "Tail weight"),
    Quantity = c(
      "Typical-day consumption, annual mean (kWh/day)",
      "Change in annual mean (kWh/day per year)",
      "Same, as % per year",
      sprintf("Winter peak (kWh/day), median date %s", dfmt(pd$peak_day[2])),
      sprintf("Summer trough (kWh/day), median date %s", dfmt(pd$trough_day[2])),
      "Peak-to-trough swing (kWh/day)",
      "Largest change from adding the second harmonic (kWh/day)",
      "Weekend minus weekday (kWh/day)",
      "Thursday versus the weekly mean (kWh/day)",
      "AR(1) coefficient $\\phi$",
      "Degrees of freedom $\\nu$"
    ),
    `Posterior mean [95% CrI]` = c(
      fmt3(pd$annual), fmt3(pd$slope_kwh, 2), fmt3(pd$slope_pct), fmt3(pd$peak_kwh),
      fmt3(pd$trough_kwh), fmt3(pd$swing), fmt3(pd$add2), fmt3(pd$wkend), fmt3(pd$dow$Thu),
      fmt3(pd$phi, 2), fmt3(pd$nu)
    ),
    check.names = FALSE
  )
  write.csv(post_tab, file.path(IN, "table_posterior.csv"), row.names = FALSE)

  # ---------------- Figure 7: annual cycle and weekly pattern ----------------
  mq <- apply(M2, 2, quantile, c(0.025, 0.5, 0.975))
  m1q <- apply(M1_model, 2, quantile, c(0.025, 0.5, 0.975))
  obs <- d[d$observed, ]
  mon <- as.integer(format(obs$date, "%m"))
  mon_med <- tapply(obs$kwh, mon, median)
  mids <- as.Date(sprintf("2009-%02d-15", 1:12))
  pq <- apply(P, 2, quantile, c(0.025, 0.5, 0.975))
  png(file.path(FIG, "07_season_week.png"), width = 12, height = 4.4, units = "in", res = 140)
  layout(matrix(1:2, nrow = 1), widths = c(1.7, 1))
  par(mar = c(3, 4.5, 2.5, 1))
  plot(yr, mq[2, ], type = "n", ylim = range(c(mq, m1q, mon_med)), xlab = "", ylab = "typical-day kWh",
       main = "(a) Annual cycle", las = 1, xaxt = "n")
  axis(1, at = mids, labels = substr(month.abb, 1, 1))
  polygon(c(yr, rev(yr)), c(mq[1, ], rev(mq[3, ])), col = adjustcolor(NAV, 0.22), border = NA)
  polygon(c(yr, rev(yr)), c(m1q[1, ], rev(m1q[3, ])), col = adjustcolor(RED, 0.15), border = NA)
  lines(yr, mq[2, ], col = NAV, lwd = 2)
  lines(yr, m1q[2, ], col = RED, lwd = 1.6, lty = 2)
  points(mids, mon_med, pch = 16, cex = 0.8)
  legend("bottomleft", legend = c("two harmonics (final model)", "one harmonic", "observed monthly median"),
         col = c(NAV, RED, "black"), lwd = c(2, 1.6, NA), lty = c(1, 2, NA), pch = c(NA, NA, 16),
         bty = "n", cex = 0.8)
  plot(1:7, pq[2, ], ylim = range(pq), pch = 16, col = NAV, xaxt = "n", xlab = "",
       ylab = "kWh/day vs weekly mean", main = "(b) Day of week", las = 1)
  arrows(1:7, pq[1, ], 1:7, pq[3, ], angle = 90, code = 3, length = 0.04, col = NAV)
  axis(1, at = 1:7, labels = DAYNAMES)
  abline(h = 0, lty = 3)
  dev.off()

  # ---------------- forecasting (holdout) ----------------
  r <- ho$roll[ho$roll$observed, ]
  r$w95 <- r$hi95 - r$lo95
  r$c95 <- r$kwh >= r$lo95 & r$kwh <= r$hi95
  r$c50 <- r$kwh >= r$lo50 & r$kwh <= r$hi50
  byh <- aggregate(cbind(w95, c95, c50, crps) ~ h, data = r, FUN = mean)
  sg <- ho$single[ho$single$observed, ]
  S$fc <- list(
    n_train = ho$n_train, n_test = ho$n_test, n_origins = ho$n_origins, train_end = ho$train_end,
    byh = byh,
    single = c(c95 = mean(sg$kwh >= sg$lo95 & sg$kwh <= sg$hi95),
               c50 = mean(sg$kwh >= sg$lo50 & sg$kwh <= sg$hi50),
               crps = mean(sg$crps), mae = mean(abs(sg$kwh - sg$med)), w95 = mean(sg$hi95 - sg$lo95)),
    phi_holdout = mean(ho$phi), theory_ratio = 1 / sqrt(1 - mean(ho$phi)^2)
  )

  # AR(2) versus AR(1) on identical held-out forecasts (positive favours AR(2))
  key <- function(z) paste(z$origin, z$h)
  r2 <- ho2$roll[ho2$roll$observed, ]
  stopifnot(identical(key(r), key(r2)))
  # forecasts from one origin overlap, so average within origin before the SE
  per_origin <- function(v) {
    m <- tapply(v, r$origin, mean)
    c(diff = mean(m), se = sd(m) / sqrt(length(m)))
  }
  s2 <- ho2$single[ho2$single$observed, ]
  S$ar_order <- list(
    roll_lpd = per_origin(r2$lpd - r$lpd), roll_crps = per_origin(r$crps - r2$crps),
    single_lpd = mean(s2$lpd - sg$lpd), single_crps = mean(sg$crps - s2$crps),
    phi_ar2 = colMeans(ho2$phi)
  )
  write.csv(byh, file.path(IN, "horizon.csv"), row.names = FALSE)

  png(file.path(FIG, "08_forecast.png"), width = 12, height = 7.5, units = "in", res = 140)
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))
  par(mar = c(4.5, 4.5, 2.5, 1))
  s1 <- ho$single
  tr <- ho$train_tail
  yl <- range(c(tr$kwh, s1$kwh, s1$hi95, s1$lo95), na.rm = TRUE)
  plot(tr$date, tr$kwh, type = "l", col = "grey35", lwd = 0.8, xlim = c(min(tr$date), max(s1$date)),
       ylim = c(yl[1], yl[2] * 1.15), xlab = "", ylab = "kWh/day",
       main = sprintf("(a) Forecast of the final year from %s (nothing after the dashed line was used in fitting)",
                      format(ho$train_end, "%d %b %Y")), las = 1)
  polygon(c(s1$date, rev(s1$date)), c(s1$lo95, rev(s1$hi95)), col = adjustcolor(NAV, 0.2), border = NA)
  polygon(c(s1$date, rev(s1$date)), c(s1$lo50, rev(s1$hi50)), col = adjustcolor(NAV, 0.35), border = NA)
  lines(s1$date, s1$med, col = NAV, lwd = 1.4)
  points(s1$date, s1$kwh, pch = 16, cex = 0.35, col = RED)
  abline(v = ho$train_end, lty = 2)
  legend("topright", legend = c("training data", "50% / 95% predictive interval", "predictive median", "held-out actual"),
         col = c("grey35", adjustcolor(NAV, 0.5), NAV, RED), lty = c(1, NA, 1, NA), pch = c(NA, 15, NA, 16),
         pt.cex = c(NA, 1.6, NA, 0.8), bty = "n", cex = 0.8, ncol = 4)
  plot(byh$h, byh$w95, type = "b", pch = 16, cex = 0.5, col = NAV, ylim = c(0, max(byh$w95) * 1.08),
       xlab = "forecast horizon (days ahead)", ylab = "mean width of 95% interval (kWh)",
       main = "(b) Predictive uncertainty by horizon", las = 1)
  plot(byh$h, 100 * byh$c95, type = "b", pch = 16, cex = 0.5, col = NAV, ylim = c(30, 100),
       xlab = "forecast horizon (days ahead)", ylab = "empirical coverage (%)",
       main = "(c) Calibration by horizon", las = 1)
  lines(byh$h, 100 * byh$c50, type = "b", pch = 1, cex = 0.5, col = RED)
  abline(h = c(95, 50), lty = 2, col = c(NAV, RED))
  legend("right", legend = c("95% interval", "50% interval"), col = c(NAV, RED), pch = c(16, 1), bty = "n", cex = 0.8)
  dev.off()

  # ---------------- Q6: unusual days and periods ----------------
  od <- one_step_diagnostics(f, NULL, d$observed)
  lam <- od$lam_mean
  low_w <- d$observed & lam < 0.5
  sgn <- sign(od$std_resid_mean)
  absence <- rep(FALSE, nrow(d))
  rl <- rle(d$observed & d$kwh < 12)
  ends <- cumsum(rl$lengths)
  for (i in which(rl$values & rl$lengths >= 5)) absence[(ends[i] - rl$lengths[i] + 1):ends[i]] <- TRUE
  ok <- is.finite(od$ppp)

  # block-level check: each 14-day mean against replicated series from the
  # fitted model (systematic part plus a stationary AR(1) + t error), one
  # series per posterior draw. p-values are compared with 0.01 and 0.99, so
  # the check needs thousands of series to be stable.
  set.seed(3)
  rows_b <- seq_len(nrow(b))
  n_paths <- length(rows_b)
  mu_s <- b[rows_b, , drop = FALSE] %*% t(f$X)
  nT <- ncol(mu_s)
  sig_s <- flatten(f, "sigma")[rows_b]
  nu_s <- flatten(f, "nu")[rows_b]
  phi_s <- as.matrix(flatten(f, "phi"))[rows_b, 1]
  lam_s <- matrix(rgamma(n_paths * nT, shape = rep(nu_s / 2, nT), rate = rep(nu_s / 2, nT)), n_paths, nT)
  eta_s <- matrix(rnorm(n_paths * nT), n_paths, nT) * sig_s / sqrt(lam_s)
  eps_s <- matrix(0, n_paths, nT)
  eps_s[, 1] <- eta_s[, 1] / sqrt(1 - phi_s^2)
  for (i in 2:nT) eps_s[, i] <- phi_s * eps_s[, i - 1] + eta_s[, i]
  y_rep <- pmax(mu_s + eps_s, 0)^2
  y_rep[, !d$observed] <- NA
  roll <- function(v) {
    k <- stats::filter(as.numeric(!is.na(v)), rep(1, WIN), sides = 1)
    s <- stats::filter(ifelse(is.na(v), 0, v), rep(1, WIN), sides = 1)
    ifelse(k >= 0.7 * WIN, s / pmax(k, 1), NA_real_)
  }
  obs_roll <- roll(d$kwh)
  rep_roll <- t(apply(y_rep, 1, roll))
  p_blk <- colMeans(rep_roll <= rep(obs_roll, each = n_paths), na.rm = TRUE)
  p_blk[!is.finite(obs_roll)] <- NA
  flag <- is.finite(p_blk) & (p_blk < 0.01 | p_blk > 0.99)
  dir_hi <- p_blk > 0.5
  fitted_typ <- colMeans(pmax(b %*% t(f$X), 0)^2)
  eps_rows <- list()
  i <- 1L
  while (i <= nT) {
    if (isTRUE(flag[i])) {
      j <- i
      while (j + 1L <= nT && isTRUE(flag[j + 1L]) && dir_hi[j + 1L] == dir_hi[i]) j <- j + 1L
      lo <- max(1L, i - WIN + 1L)
      eps_rows[[length(eps_rows) + 1]] <- c(lo, j, as.integer(dir_hi[i]))
      i <- j + 1L
    } else i <- i + 1L
  }
  merged <- list()
  for (e in eps_rows) {
    k <- length(merged)
    if (k && merged[[k]][3] == e[3] && e[1] <= merged[[k]][2] + 1) merged[[k]][2] <- max(merged[[k]][2], e[2])
    else merged[[k + 1]] <- e
  }
  episodes <- do.call(rbind, lapply(merged, function(e) {
    ix <- e[1]:e[2]
    data.frame(
      from = d$date[e[1]], to = d$date[e[2]], days = length(ix),
      direction = if (e[3] == 1) "high" else "low",
      observed_kwh = mean(d$kwh[ix], na.rm = TRUE), typical_kwh = mean(fitted_typ[ix]),
      in_absence_block = mean(absence[ix])
    )
  }))
  episodes$ratio <- episodes$observed_kwh / episodes$typical_kwh
  write.csv(episodes, file.path(IN, "episodes.csv"), row.names = FALSE)

  S$q6 <- list(
    n_eval = sum(ok), ppp_lo = sum(od$ppp[ok] < 0.01), ppp_hi = sum(od$ppp[ok] > 0.99),
    n_low_weight = sum(low_w), n_low_weight_high = sum(low_w & sgn > 0, na.rm = TRUE),
    n_low_weight_low = sum(low_w & sgn < 0, na.rm = TRUE),
    absence_days = sum(absence), absence_low_weight = sum(low_w & absence),
    lam_absence_mean = mean(lam[absence]), lam_other_mean = mean(lam[d$observed & !absence]),
    n_windows = sum(is.finite(p_blk)), n_flag_lo = sum(p_blk < 0.01, na.rm = TRUE),
    n_flag_hi = sum(p_blk > 0.99, na.rm = TRUE), episodes = episodes
  )

  png(file.path(FIG, "09_anomalies.png"), width = 12, height = 4.4, units = "in", res = 140)
  par(mar = c(3, 4.5, 2.5, 1))
  plot(d$date, d$kwh, type = "n", xlab = "", ylab = "kWh/day", las = 1,
       ylim = c(0, 1.12 * max(d$kwh, na.rm = TRUE)),
       main = "Periods the fitted model does not explain (14-day predictive check)")
  for (k in seq_len(nrow(episodes))) {
    rect(episodes$from[k], par("usr")[3], episodes$to[k], par("usr")[4], border = NA,
         col = adjustcolor(if (episodes$direction[k] == "high") RED else NAV, 0.18))
  }
  lines(d$date, d$kwh, col = "grey45", lwd = 0.5)
  lines(d$date, stats::filter(fitted_typ, rep(1 / 7, 7), sides = 2), col = "black", lwd = 1.6)
  points(d$date[low_w], d$kwh[low_w], pch = 1, cex = 0.55, col = RED)
  legend("topright", legend = c("observed", "typical day from the model (7-day mean)", "weight below 1/2",
                                "unusually high period", "unusually low period"),
         col = c("grey45", "black", RED, adjustcolor(RED, 0.4), adjustcolor(NAV, 0.4)),
         lty = c(1, 1, NA, NA, NA), pch = c(NA, NA, 1, 15, 15), pt.cex = c(NA, NA, 0.8, 1.8, 1.8),
         bty = "n", cex = 0.75, ncol = 5)
  dev.off()

  saveRDS(S, file.path(IN, "summary.rds"))
  invisible(S)
}

if (sys.nframe() == 0L) {
  S <- main()
  str(S, max.level = 2, give.attr = FALSE)
}
