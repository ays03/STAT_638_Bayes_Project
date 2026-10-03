#!/usr/bin/env Rscript
# Step 5: turn posterior draws into the quantities the scientific questions ask
# about, and run the posterior predictive check for unusual days.
#
#   Rscript src/05_derived.R
# Expects results/03_fits/M3.rds and M7.rds from src/03_fit_models.R.

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

ROOT <- script_root()
source(file.path(ROOT, "src", "model.R"))

DAYNAMES <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

load_fit <- function(mid) readRDS(file.path(ROOT, "results", "03_fits", paste0(mid, ".rds")))

ci <- function(v, lo = 0.025, hi = 0.975) quantile(v, c(lo, hi))

fmt <- function(v, unit = "", dp = 3) {
  qs <- ci(v)
  sprintf(paste0("%.", dp, "f [%.", dp, "f, %.", dp, "f]%s"), mean(v), qs[1], qs[2], unit)
}

dmon <- function(x) sub("^0", "", format(as.Date(x), "%d %b"))
dmony <- function(x) sub("^0", "", format(as.Date(x), "%d %b %Y"))

roll_mean_vec <- function(a, w) {
  finite <- is.finite(a)
  csum <- c(0, cumsum(ifelse(finite, a, 0)))
  nsum <- c(0, cumsum(as.numeric(finite)))
  n <- length(a)
  s <- csum[(w + 1):(n + 1)] - csum[seq_len(n + 1 - w)]
  k <- nsum[(w + 1):(n + 1)] - nsum[seq_len(n + 1 - w)]
  ifelse(k >= 0.7 * w, s / pmax(k, 1), NA_real_)
}

roll_mean_mat <- function(a, w) {
  finite <- is.finite(a)
  a0 <- ifelse(finite, a, 0)
  csum <- t(apply(a0, 1, cumsum))
  nsum <- t(apply(finite, 1, function(z) cumsum(as.numeric(z))))
  csum <- cbind(0, csum)
  nsum <- cbind(0, nsum)
  nT <- ncol(a)
  s <- csum[, (w + 1):(nT + 1), drop = FALSE] - csum[, seq_len(nT + 1 - w), drop = FALSE]
  k <- nsum[, (w + 1):(nT + 1), drop = FALSE] - nsum[, seq_len(nT + 1 - w), drop = FALSE]
  ifelse(k >= 0.7 * w, s / pmax(k, 1), NA_real_)
}

main <- function() {
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  dir.create(file.path(ROOT, "results", "tables"), showWarnings = FALSE, recursive = TRUE)
  d <- read_daily(file.path(ROOT, "data", "daily.csv"))
  t <- as.numeric(d$t)
  date0 <- d$date[1]
  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("DERIVED POSTERIOR QUANTITIES")
  add("M7 = J=2, AR(1) + t errors, kWh scale  (best out-of-sample forecaster)")
  add("M3 = J=2, AR(1) + t errors, log scale  (best by WAIC among J=2 models)")
  add("Intervals are 95% central posterior credible intervals.")
  add("")

  for (mid in c("M7", "M3")) {
    f <- load_fit(mid)
    nm <- f$names
    b <- flatten(f, "beta")
    sig <- flatten(f, "sigma")
    phi <- flatten(f, "phi")
    nu <- flatten(f, "nu")
    scale <- f$scale
    unit <- if (scale == "identity") "kWh" else "log kWh"
    add(sprintf("================ %s  (%s scale) ================", mid, unit))
    i0 <- match("intercept", nm)
    dow_cols <- match(paste0("dow_", DAYNAMES[-1]), nm)
    dow_full <- cbind(0, b[, dow_cols, drop = FALSE])
    level <- b[, i0] + rowMeans(dow_full)
    add(sprintf("  average-day level at record midpoint : %s", fmt(level, paste0(" ", unit), 3)))
    i1 <- match("trend_per_year", nm)
    tr <- b[, i1]
    add("")
    add("  Q2 LONG-TERM TREND")
    add(sprintf("    slope                              : %s", fmt(tr, paste0(" ", unit, "/year"), 4)))
    if (scale == "identity") {
      pct <- 100 * tr / level
      add(sprintf("    as a percentage of the mean level   : %s", fmt(pct, " %/year", 2)))
      add(sprintf("    total change over the 3.95-y record: %s", fmt(3.95 * tr, " kWh/day", 3)))
    } else {
      add(sprintf("    as a percentage                    : %s", fmt(100 * (exp(tr) - 1), " %/year", 2)))
    }
    add(sprintf("    P(slope < 0 | data)                : %.3f", mean(tr < 0)))
    add(sprintf("    P(|slope| < 1%% of level per year)  : %.3f", mean(abs(tr) < 0.01 * abs(level))))

    add("")
    add("  Q3 ANNUAL SEASONALITY")
    amps <- list()
    for (j in 1:2) {
      if (!paste0("sin_h", j) %in% nm) next
      a_ <- b[, match(paste0("sin_h", j), nm)]
      c_ <- b[, match(paste0("cos_h", j), nm)]
      R <- sqrt(a_^2 + c_^2)
      amps[[length(amps) + 1]] <- R
      psi <- atan2(a_, c_)
      period <- 365 / j
      peak_t <- (psi %% (2 * pi)) * period / (2 * pi)
      add(sprintf("    harmonic %d: amplitude            : %s", j, fmt(R, paste0(" ", unit), 3)))
      add(sprintf("                peak-to-trough swing : %s", fmt(2 * R, paste0(" ", unit), 3)))
      pk <- date0 + as.numeric(quantile(peak_t, c(0.025, 0.5, 0.975)))
      add(sprintf("                peaks at day-of-cycle: %.1f of %.1f  (calendar %s, 95%% CI %s to %s)",
                  median(peak_t), period, format(pk[2], "%d %b"), format(pk[1], "%d %b"), format(pk[3], "%d %b")))
    }
    if (length(amps) == 2) {
      ratio <- amps[[2]] / amps[[1]]
      add(sprintf("    amplitude ratio  A2 / A1           : %s", fmt(ratio, "", 3)))
      add(sprintf("    P(A2 > 0.15 * A1 | data)           : %.3f", mean(ratio > 0.15)))
      add(sprintf("    P(A2 > 0.25 * A1 | data)           : %.3f", mean(ratio > 0.25)))
      add("    (a second harmonic this large makes the annual curve visibly")
      add("     asymmetric: a sharper winter peak and a flatter summer floor)")
    }
    tt <- 1:365
    seas <- matrix(0, nrow(b), length(tt))
    for (j in 1:2) {
      sn <- paste0("sin_h", j)
      cn <- paste0("cos_h", j)
      if (sn %in% nm) {
        wave_s <- sin(2 * pi * j * tt / 365)
        wave_c <- cos(2 * pi * j * tt / 365)
        seas <- seas + outer(b[, match(sn, nm)], wave_s) + outer(b[, match(cn, nm)], wave_c)
      }
    }
    rng_season <- apply(seas, 1, max) - apply(seas, 1, min)
    add(sprintf("    full seasonal range (max - min)    : %s", fmt(rng_season, paste0(" ", unit), 3)))
    add("")
    add("  Q4 WEEKLY PATTERN (deviation from the weekly mean)")
    prof <- dow_full - rowMeans(dow_full)
    for (k in seq_along(DAYNAMES)) add(sprintf("    %-34s : %s", DAYNAMES[k], fmt(prof[, k], paste0(" ", unit), 3)))
    wkend <- rowMeans(prof[, 6:7, drop = FALSE]) - rowMeans(prof[, 1:5, drop = FALSE])
    add(sprintf("    weekend minus weekday              : %s", fmt(wkend, paste0(" ", unit), 3)))
    add(sprintf("    P(weekend > weekday | data)        : %.4f", mean(wkend > 0)))
    add("")
    add("  ERROR PROCESS")
    add(sprintf("    sigma (one-step innovation sd)     : %s", fmt(sig, paste0(" ", unit), 3)))
    add(sprintf("    phi_1                              : %s", fmt(phi[, 1], "", 3)))
    add(sprintf("    nu (t degrees of freedom)          : %s", fmt(nu, "", 2)))
    ratio_sd <- 1 / sqrt(1 - phi[, 1]^2)
    add(sprintf("    marginal sd / one-step sd          : %s", fmt(ratio_sd, "", 4)))
    add("    -> this ratio is the ceiling on how much wider a long-horizon")
    add("       predictive interval is than a one-day-ahead one")
    add("")
  }

  f <- load_fit("M7")
  diag <- one_step_diagnostics(f, NULL, d$observed)
  d$std_resid <- diag$std_resid_mean
  d$ppp <- diag$ppp
  d$lam <- diag$lam_mean
  add("================ Q6 UNUSUAL PERIODS (model M7) ================")
  add("One-step-ahead posterior predictive p-values p_t = P(y_t^rep <= y_t | y).")
  add("Under a well-specified model these are uniform, so about 1% of days fall")
  add("below 0.01 and 1% above 0.99 by chance alone.")
  ok <- is.finite(d$ppp)
  pv <- d$ppp[ok]
  n <- sum(ok)
  add("")
  add(sprintf("  evaluable days                    : %d", n))
  add(sprintf("  p_t < 0.01 (unusually low)        : %d  (expected %.1f if calibrated)", sum(pv < 0.01), 0.01 * n))
  add(sprintf("  p_t > 0.99 (unusually high)       : %d  (expected %.1f if calibrated)", sum(pv > 0.99), 0.01 * n))
  add(sprintf("  p_t < 0.05                        : %d  (expected %.1f)", sum(pv < 0.05), 0.05 * n))
  add(sprintf("  p_t > 0.95                        : %d  (expected %.1f)", sum(pv > 0.95), 0.05 * n))
  add("")
  add("  Note these are ONE-STEP-AHEAD checks, so a sustained absence is flagged")
  add("  only on its first day or two: once the AR(1) term has absorbed the drop,")
  add("  the model expects the low level to continue. The blocks are therefore")
  add("  better seen in the lambda_t weights below.")
  add("")

  Xf <- f$X
  mu_mean <- colMeans(flatten(f, "beta") %*% t(Xf))
  d$fitted <- mu_mean
  flag <- d[ok & (d$ppp < 0.01 | d$ppp > 0.99), ]
  flag$extremity <- pmin(flag$ppp, 1 - flag$ppp)
  flag <- flag[order(flag$extremity), ]
  add(sprintf("  FLAGGED DAYS (%d total), most extreme first", nrow(flag)))
  add(sprintf("    %12s %4s %7s %7s %7s %7s", "date", "dow", "kWh", "fitted", "z", "p_t"))
  if (nrow(flag) > 0) {
    for (i in seq_len(nrow(flag))) {
      r <- flag[i, ]
      add(sprintf("    %12s %4s %7.2f %7.2f %+7.2f %7.4f",
                  r$date, DAYNAMES[as.integer(r$dow) + 1L], r$kwh, r$fitted, r$std_resid, r$ppp))
    }
  }
  add("")

  WIN <- 14L
  set.seed(3)
  pars <- f$path_pars
  sig_s <- pars[, 1]
  nu_s <- pars[, 2]
  phi_s <- pars[, 4]
  S <- nrow(pars)
  beta_s <- flatten(f, "beta")
  sel <- sample.int(nrow(beta_s), S)
  mu_s <- beta_s[sel, , drop = FALSE] %*% t(Xf)
  nT <- nrow(Xf)
  lam_s <- matrix(rgamma(S * nT, shape = rep(nu_s / 2, each = nT), scale = rep(2 / nu_s, each = nT)),
                  nrow = S, byrow = TRUE)
  eta_s <- matrix(rnorm(S * nT), nrow = S, byrow = TRUE) * sig_s / sqrt(lam_s)
  eps_s <- matrix(0, S, nT)
  eps_s[, 1] <- eta_s[, 1] / sqrt(1 - phi_s^2)
  for (i in 2:nT) eps_s[, i] <- phi_s * eps_s[, i - 1] + eta_s[, i]
  y_rep <- mu_s + eps_s

  obs_roll <- roll_mean_vec(d$kwh, WIN)
  rep_roll <- roll_mean_mat(y_rep, WIN)
  cmp <- rep_roll <= rep(obs_roll, each = S)
  cmp[is.na(rep_roll)] <- FALSE
  p_block <- colMeans(cmp)
  p_block[!is.finite(obs_roll)] <- NA_real_
  centre <- seq_along(p_block) + (WIN %/% 2)

  add(sprintf("  BLOCK-LEVEL CHECK: %d-DAY MEAN CONSUMPTION", WIN))
  add(sprintf("  Test statistic is the mean of each %d-day window, referred to its", WIN))
  add("  posterior predictive distribution under the stationary AR(1)+t process.")
  add("  This is the diagnostic that can see sustained anomalies.")
  okb <- is.finite(p_block)
  add("")
  add(sprintf("    windows evaluated                : %d", sum(okb)))
  add(sprintf("    p < 0.01                         : %d", sum(p_block[okb] < 0.01)))
  add(sprintf("    p > 0.99                         : %d", sum(p_block[okb] > 0.99)))
  add("")
  add("  SUSTAINED ANOMALOUS PERIODS")
  add("  Runs of consecutive flagged windows, each expanded to the days it covers.")
  add("  Because the windows overlap, two runs of the same direction whose covered")
  add("  days overlap describe one episode and are merged.")
  flagb <- okb & (p_block < 0.01 | p_block > 0.99)
  high <- p_block > 0.99
  runs <- list()
  i <- 1L
  while (i <= length(flagb)) {
    if (isTRUE(flagb[i])) {
      j <- i
      while (j + 1L <= length(flagb) && isTRUE(flagb[j + 1L]) && isTRUE(high[j + 1L] == high[i])) j <- j + 1L
      lo <- max(1L, centre[i] - WIN %/% 2)
      hi <- min(nrow(d), centre[j] + WIN %/% 2)
      pv <- if (isTRUE(high[i])) max(p_block[i:j]) else min(p_block[i:j])
      runs[[length(runs) + 1]] <- c(lo, hi, as.integer(isTRUE(high[i])), pv)
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  merged <- list()
  for (r in runs) {
    if (length(merged) && merged[[length(merged)]][3] == r[3] && r[1] <= merged[[length(merged)]][2] + 1) {
      merged[[length(merged)]][2] <- max(merged[[length(merged)]][2], r[2])
      if (r[3] == 1) merged[[length(merged)]][4] <- max(merged[[length(merged)]][4], r[4])
      else merged[[length(merged)]][4] <- min(merged[[length(merged)]][4], r[4])
    } else {
      merged[[length(merged) + 1]] <- r
    }
  }
  rows <- list()
  add("")
  add(sprintf("    %12s %12s %5s %8s %8s %6s %8s", "from", "to", "days", "obs kWh", "exp kWh", "ratio", "p"))
  for (r in merged) {
    lo <- r[1]; hi <- r[2]; pv <- r[4]
    sub <- d[lo:hi, ]
    o <- mean(sub$kwh, na.rm = TRUE)
    e <- mean(sub$fitted, na.rm = TRUE)
    add(sprintf("    %12s %12s %5d %8.2f %8.2f %6.2f %8.4f", sub$date[1], sub$date[nrow(sub)], hi - lo + 1, o, e, o / e, pv))
    rows[[length(rows) + 1]] <- list(a0 = sub$date[1], a1 = sub$date[nrow(sub)], nd = hi - lo + 1, o = o, e = e, rt = o / e, pv = pv)
  }
  con <- file(file.path(ROOT, "results", "tables", "anomaly_periods.tex"), "w")
  writeLines("\\begin{tabular}{@{}llrrrrr@{}}", con)
  writeLines("\\toprule", con)
  writeLines("From & To & Days & Observed & Expected & Ratio & $p$ \\\\", con)
  writeLines("\\midrule", con)
  for (r in rows) {
    bold <- r$rt < 0.6 || r$rt > 1.3
    w <- if (bold) function(s) paste0("\\textbf{", s, "}") else function(s) s
    pv_s <- if (r$pv < 0.001) "$<$0.001" else sprintf("%.3f", min(r$pv, 1 - r$pv))
    writeLines(sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      w(dmony(r$a0)), w(dmony(r$a1)), w(as.character(r$nd)),
      w(sprintf("%.2f", r$o)), w(sprintf("%.2f", r$e)), w(sprintf("%.2f", r$rt)), w(pv_s)
    ), con)
  }
  writeLines("\\bottomrule", con)
  writeLines("\\end{tabular}", con)
  close(con)
  add("")
  add(sprintf("  (%d episodes; LaTeX version written to results/tables/anomaly_periods.tex)", length(rows)))
  add("")
  add("    These are the load-bearing anomalies: multi-week absences (the French")
  add("    August holiday is the clearest) and unusually cold or mild stretches.")
  add("    They are behavioural and meteorological, and no calendar-based model")
  add("    can predict them, which is exactly why the t likelihood is needed.")
  d$p_block <- NA_real_
  d$p_block[centre] <- p_block
  write.csv(d, file.path(ROOT, "results", "05_daily_with_diagnostics.csv"), row.names = FALSE)
  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "05_derived.txt"))

  f7 <- load_fit("M7")
  b <- flatten(f7, "beta")
  nm <- f7$names
  tt <- 1:365
  comp <- list()
  for (j in 1:2) {
    comp[[j]] <- outer(b[, match(paste0("sin_h", j), nm)], sin(2 * pi * j * tt / 365)) +
      outer(b[, match(paste0("cos_h", j), nm)], cos(2 * pi * j * tt / 365))
  }
  tot <- comp[[1]] + comp[[2]]
  cal <- date0 + (tt - 1)
  xn <- as.numeric(cal)
  png(file.path(ROOT, "figures", "05_seasonal.png"), width = 13, height = 4.2, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  q1 <- apply(comp[[1]], 2, quantile, probs = c(0.025, 0.5, 0.975))
  qt <- apply(tot, 2, quantile, probs = c(0.025, 0.5, 0.975))
  ylim <- range(c(q1, qt))
  plot(cal, qt[2, ], type = "n", ylim = ylim, xlab = "", ylab = "deviation from annual mean (kWh/day)",
       main = "Annual seasonal component, model M7\n(bands are 95% credible)")
  polygon(c(xn, rev(xn)), c(q1[1, ], rev(q1[3, ])), col = adjustcolor("#b03a2e", 0.18), border = NA)
  polygon(c(xn, rev(xn)), c(qt[1, ], rev(qt[3, ])), col = adjustcolor("#1f4e79", 0.18), border = NA)
  lines(cal, q1[2, ], col = "#b03a2e", lwd = 1.8)
  lines(cal, qt[2, ], col = "#1f4e79", lwd = 1.8)
  abline(h = 0)
  legend("topright", legend = c("first harmonic only", "first + second harmonic"),
         col = c("#b03a2e", "#1f4e79"), lwd = 1.8, bty = "n", cex = 0.8)
  q2 <- apply(comp[[2]], 2, quantile, probs = c(0.025, 0.5, 0.975))
  plot(cal, q2[2, ], type = "n", ylim = range(q2), xlab = "", ylab = "kWh/day",
       main = "Second harmonic alone: what it adds\n(sharpens the winter peak, flattens the summer floor)")
  polygon(c(xn, rev(xn)), c(q2[1, ], rev(q2[3, ])), col = adjustcolor("#2e7d32", 0.20), border = NA)
  lines(cal, q2[2, ], col = "#2e7d32", lwd = 1.8)
  abline(h = 0)
  dev.off()

  dow_cols <- match(paste0("dow_", DAYNAMES[-1]), nm)
  dow_full <- cbind(0, b[, dow_cols, drop = FALSE])
  prof <- dow_full - rowMeans(dow_full)
  png(file.path(ROOT, "figures", "05_weekly.png"), width = 12, height = 3.8, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  q <- apply(prof, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  plot(0:6, q[3, ], ylim = range(q), pch = 16, col = "#1f4e79", xaxt = "n",
       xlab = "", ylab = "kWh/day vs weekly mean", main = "Day-of-week effect, model M7")
  axis(1, at = 0:6, labels = DAYNAMES)
  arrows(0:6, q[1, ], 0:6, q[5, ], angle = 90, code = 3, length = 0.05, col = "#1f4e79")
  arrows(0:6, q[2, ], 0:6, q[4, ], angle = 90, code = 3, length = 0, lwd = 3, col = "#1f4e79")
  abline(h = 0)
  legend("topleft", legend = c("95% CI", "50% CI"), lwd = c(1, 3), col = "#1f4e79", bty = "n", cex = 0.8)
  wkend <- rowMeans(prof[, 6:7, drop = FALSE]) - rowMeans(prof[, 1:5, drop = FALSE])
  hist(wkend, breaks = 60, col = "#1f4e79", xlab = "weekend minus weekday (kWh/day)",
       main = sprintf("Posterior of the weekend effect\nP(> 0) = %.4f", mean(wkend > 0)))
  abline(v = 0)
  dev.off()

  png(file.path(ROOT, "figures", "05_anomalies.png"), width = 12, height = 6, units = "in", res = 140)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  plot(d$date, d$kwh, type = "l", lwd = 0.5, col = "#888888", xlab = "", ylab = "kWh/day",
       main = "Posterior predictive check: which days the model cannot explain")
  lines(d$date, d$fitted, col = "#1f4e79", lwd = 1)
  fl <- d[is.finite(d$ppp) & (d$ppp < 0.01 | d$ppp > 0.99), ]
  points(fl$date, fl$kwh, pch = 1, col = "#b03a2e")
  legend("topright", legend = c("observed", "fitted systematic part", "flagged, p < 0.01 or > 0.99"),
         col = c("#888888", "#1f4e79", "#b03a2e"), lty = c(1, 1, NA), pch = c(NA, NA, 1), cex = 0.75, ncol = 3, bty = "n")
  plot(d$date, d$lam, type = "l", lwd = 0.7, col = "#2e7d32", xlab = "date", ylab = "posterior mean lambda")
  abline(h = 0.5, lty = 2, col = "#b03a2e")
  legend("topright", legend = "lambda = 0.5 downweighting threshold", col = "#b03a2e", lty = 2, bty = "n", cex = 0.8)
  title("Heavy-tail weights: low lambda marks days treated as outliers")
  dev.off()
  cat("\nwrote results/05_derived.txt, 05_daily_with_diagnostics.csv, 3 figures\n")
}

if (sys.nframe() == 0L) main()
