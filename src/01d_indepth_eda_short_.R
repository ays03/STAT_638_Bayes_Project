#!/usr/bin/env Rscript
# Step 1d: exploratory analysis of daily household demand (kWh scale).
#
# Reads only data/daily.csv (no raw minute file). All paths are built with here(),
# which finds the project root from the empty marker file project/.here, so the
# script runs from any working directory inside the project (e.g. project/ or project/src/).
#
# Outputs
#   R_figures/01d_overview.png         EDA Fig 1: daily series, daily coverage, annual shape (J=1 vs J=2),
#                                      weekday vs weekend, weekly pattern by day
#   R_figures/01d_residual_checks.png  EDA Fig 2: residuals over time, residual vs fitted, normal Q-Q
#   R_figures/01d_residual.png         EDA Fig 3: residual ACF, PACF, residual SD by month
#   results/tables/table1_scale_check.csv   Table 1: residual summary on kWh and log scales
#   results/01d_eda.txt                     text summary of the numbers behind the figures
#
# WHERE THINGS ARE:
#   [HELPERS]   small functions used below (stats, harmonics, ACF/PACF, absence blocks)
#   [DATA]      read daily.csv, set up y / t / dow / month / year, find absence blocks
#   [SUMMARY]   numbers written to results/01d_eda.txt
#   [RESID]     OLS fit of trend + 2 harmonics + day of week -> residuals used everywhere below
#   [TABLE 1]   residual skew / kurtosis / lag-1 ACF on kWh vs log scale
#   [FIGURE 1]  data overview
#   [FIGURE 2]  residual checks
#   [FIGURE 3]  residual ACF, PACF, SD by month

ORIGIN <- as.Date("2006-12-17")   # t = 1
MONTHS <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
DOW_NAMES <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
MIN_COVERAGE <- 0.95              # 95% retention rule used in 01_aggregate.R (for Figure 1b)

library(here)   # install.packages("here") once; project root = folder containing .here

# ============================================================================
# [HELPERS]  functions used in main()
# ============================================================================
read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

skew_adj <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  xc <- x - mean(x)
  m2 <- sum(xc^2) / n
  m3 <- sum(xc^3) / n
  (m3 / m2^1.5) * sqrt(n * (n - 1)) / (n - 2)
}

# plain (not excess) kurtosis: a normal distribution gives 3; > 3 means heavy tails (used in Table 1)
kurt_raw <- function(x) { xc <- x - mean(x); mean(xc^4) / mean(xc^2)^2 }

ols <- function(X, y) {
  beta <- qr.solve(X, y)
  list(beta = beta, resid = as.numeric(y - X %*% beta))
}

clim_day <- function(dates) {
  dates <- as.Date(dates)
  doy <- as.integer(format(dates, "%j"))
  y <- as.integer(format(dates, "%Y"))
  m <- as.integer(format(dates, "%m"))
  day <- as.integer(format(dates, "%d"))
  leap <- (y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0))
  d <- as.numeric(doy)
  d[leap & m > 2] <- d[leap & m > 2] - 1
  d[m == 2 & day == 29] <- NA_real_
  d
}

smooth_on_clim <- function(clim, values, window = 15, circular = TRUE, min_frac = 0.6) {
  half <- window %/% 2
  s <- numeric(365)
  cnt <- numeric(365)
  for (i in seq_along(clim)) {
    if (is.na(clim[i]) || is.na(values[i])) next
    s[as.integer(clim[i])] <- s[as.integer(clim[i])] + values[i]
    cnt[as.integer(clim[i])] <- cnt[as.integer(clim[i])] + 1
  }
  out <- rep(NA_real_, 365)
  for (i in 0:364) {
    if (circular) idx <- ((i + (-half:half)) %% 365) + 1L
    else {
      ks <- (-half:half)
      ks <- ks[i + ks >= 0 & i + ks < 365]
      idx <- i + ks + 1L
    }
    cc <- sum(cnt[idx])
    if (cc >= window * min_frac) out[i + 1L] <- sum(s[idx]) / cc
  }
  out
}

# runs of >= `run` consecutive retained days below `threshold` kWh (absence blocks)
low_use_mask <- function(kwh, run = 5, threshold = 12) {
  low <- kwh < threshold
  block <- rep(FALSE, length(kwh))
  i <- 1L
  n <- length(kwh)
  while (i <= n) {
    if (isTRUE(low[i])) {
      j <- i
      while (j + 1L <= n && isTRUE(low[j + 1L])) j <- j + 1L
      if (j - i + 1L >= run) block[i:j] <- TRUE
      i <- j + 1L
    } else i <- i + 1L
  }
  block
}

dow_dummies <- function(dow, drop = 0) {
  levels <- setdiff(0:6, drop)
  out <- vapply(levels, function(k) as.numeric(dow == k), numeric(length(dow)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  list(X = out, levels = levels)
}

harmonic_cols <- function(t, J) {
  if (J <= 0) return(matrix(0, length(t), 0))
  cols <- vector("list", 2 * J)
  k <- 1L
  for (j in seq_len(J)) {
    w <- 2 * pi * j * t / 365
    cols[[k]] <- sin(w); k <- k + 1L
    cols[[k]] <- cos(w); k <- k + 1L
  }
  do.call(cbind, cols)
}

harmonic_on_dates <- function(beta_harm, dates, J) {
  tt <- as.numeric(as.Date(dates) - ORIGIN) + 1
  h <- numeric(length(dates))
  for (j in seq_len(J)) {
    w <- 2 * pi * j * tt / 365
    h <- h + beta_harm[2 * j - 1] * sin(w) + beta_harm[2 * j] * cos(w)
  }
  h
}

# autocorrelation at calendar lags, so gaps from missing days are respected
acf_calendar <- function(resid, t, nlags) {
  tt <- as.integer(round(t))
  pos <- integer(max(tt))
  pos[tt] <- seq_along(tt)
  out <- numeric(nlags)
  for (k in seq_len(nlags)) {
    prev <- tt - k
    ok <- prev >= 1L & prev <= length(pos)
    j <- integer(length(tt))
    j[ok] <- pos[prev[ok]]
    ok <- ok & j > 0L
    a <- resid[ok]
    b <- resid[j[ok]]
    a <- a - mean(a)
    b <- b - mean(b)
    out[k] <- sum(a * b) / sqrt(sum(a * a) * sum(b * b))
  }
  out
}

# partial autocorrelation from the ACF via Yule-Walker
pacf_yw <- function(rho) {
  pac <- numeric(length(rho))
  for (k in seq_along(rho)) {
    r <- rho[seq_len(k)]
    if (k == 1) {
      pac[k] <- r[1]
      next
    }
    R <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(k)) {
      lag <- abs(i - j)
      R[i, j] <- if (lag == 0) 1 else rho[lag]
    }
    phi <- tryCatch(solve(R, r), error = function(e) qr.solve(R, r))
    pac[k] <- phi[k]
  }
  pac
}

dmon <- function(x) sub("^0", "", format(as.Date(x), "%d %b"))

# ============================================================================
# main(): runs everything; called on the last line of the file
# ============================================================================
main <- function() {

  # ==========================================================================
  # [DATA]  read the cleaned daily series and set up the variables
  # ==========================================================================
  dir.create(here("R_figures"), showWarnings = FALSE)
  dir.create(here("results", "tables"), showWarnings = FALSE, recursive = TRUE)

  d <- read_daily(here("data", "daily.csv"))
  obs <- d[d$observed, ]
  rownames(obs) <- NULL
  t_check <- as.numeric(obs$date - ORIGIN) + 1
  if (!isTRUE(all.equal(t_check, as.numeric(obs$t)))) stop("t index does not match days since 2006-12-17")

  y <- as.numeric(obs$kwh)
  t <- as.numeric(obs$t)
  dow <- as.integer(obs$dow)
  month <- as.integer(format(obs$date, "%m"))
  year <- as.integer(obs$year)
  n <- length(y)
  block <- low_use_mask(y)
  Dd <- dow_dummies(dow, 0)$X

  # ==========================================================================
  # [SUMMARY]  everything passed to add() is written to results/01d_eda.txt
  # ==========================================================================
  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("EDA SUMMARY  (descriptive OLS; not posterior estimates)")
  add(sprintf("Retained days: %d.  Missing days: %d.  Low-use block days: %d.", n, sum(!d$observed), sum(block)))
  add("")

  # ---------------- low-use (absence) blocks ----------------
  add("LOW-USE BLOCKS (>= 5 consecutive retained days below 12 kWh)")
  i <- 1L
  while (i <= n) {
    if (block[i]) {
      j <- i
      while (j < n && block[j + 1L]) j <- j + 1L
      add(sprintf("  %s to %s  (%2d days, mean %5.2f kWh)", obs$date[i], obs$date[j], j - i + 1L, mean(y[i:j])))
      i <- j + 1L
    } else i <- i + 1L
  }
  add("")

  # ---------------- harmonic shape (for Fig 1b) ----------------
  add("HARMONIC SHAPE ON THE kWh SCALE  (OLS with trend and day of week)")
  ref_dates <- seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "day")
  curves <- list()
  for (J in 1:2) {
    fitJ <- ols(cbind(1, t / 365, harmonic_cols(t, J), Dd), y)
    harm <- fitJ$beta[3:(2 + 2 * J)]
    raw_curve <- harmonic_on_dates(harm, ref_dates, J)
    curve_J <- raw_curve - mean(raw_curve) + mean(y[year %in% 2007:2009])
    curves[[J]] <- curve_J
    peak <- which.max(curve_J)
    trough <- which.min(curve_J)
    add(sprintf("  J=%d  peak %s (%.2f kWh)  trough %s (%.2f kWh)  range %.2f  resid sd %.2f",
                J, dmon(ref_dates[peak]), curve_J[peak], dmon(ref_dates[trough]), curve_J[trough],
                curve_J[peak] - curve_J[trough], sd(fitJ$resid)))
  }
  month_means <- vapply(1:12, function(m) mean(y[month == m]), numeric(1))
  add("")

  # ==========================================================================
  # [RESID]  OLS pre-fit of the planned mean structure (Ch 9)
  #   X_j2d    design matrix: intercept, trend (per year), 2 harmonic pairs, 6 day-of-week dummies
  #   resid_j2 kWh residuals  -> used in Table 1, Figure 2 and Figure 3
  #   fit_j2   fitted values  -> used in Figure 2(b)
  # ==========================================================================
  X_j2d <- cbind(1, t / 365, harmonic_cols(t, 2), Dd)
  resid_j2 <- ols(X_j2d, y)$resid
  fit_j2 <- y - resid_j2

  # ==========================================================================
  # [TABLE 1]  residual summary on the kWh and log scales (same mean structure X_j2d)
  #   Purpose: decides kWh vs log, and supports t errors.
  #   Columns: Residual skew (normal = 0), Kurtosis (normal = 3), Lag-1 ACF (calendar-consecutive days)
  #   Output:  results/tables/table1_scale_check.csv (read by the Rmd) + a copy in 01d_eda.txt
  #   Uses:    X_j2d from [RESID]; kurt_raw(), skew_adj(), acf_calendar(), ols() from [HELPERS]
  # ==========================================================================
  tab1 <- do.call(rbind, lapply(c("kWh", "log"), function(sc) {
    yy <- if (sc == "kWh") y else log(y)          # response on this scale
    rr <- ols(X_j2d, yy)$resid                    # residuals on this scale
    data.frame(Scale = sc,
               `Residual skew` = round(skew_adj(rr), 3),
               Kurtosis        = round(kurt_raw(rr), 2),
               `Lag-1 ACF`     = round(acf_calendar(rr, t, 1), 3),
               check.names = FALSE)
  }))
  add("TABLE 1: RESIDUAL SUMMARY ON BOTH SCALES  (normal errors: skew 0, kurtosis 3)")
  add(paste(capture.output(print(tab1, row.names = FALSE)), collapse = "\n"))
  write.csv(tab1, here("results", "tables", "table1_scale_check.csv"), row.names = FALSE)   # for the Rmd
  add("")

  # ---------------- residual dependence ----------------
  rho <- acf_calendar(resid_j2, t, 21)
  rho_ex <- acf_calendar(resid_j2[!block], t[!block], 21)
  pac <- pacf_yw(rho)
  pac_ex <- pacf_yw(rho_ex)
  add("RESIDUAL ACF (calendar lags)")
  add(sprintf("%4s %10s %10s %10s", "lag", "all days", "AR1 pred", "no blocks"))
  for (k in c(1, 2, 3, 4, 5, 6, 7, 14, 21)) add(sprintf("%4d %10.3f %10.3f %10.3f", k, rho[k], rho[1]^k, rho_ex[k]))
  add("PACF")
  for (k in 1:7) add(sprintf("  lag %d: all days %+.3f    no blocks %+.3f", k, pac[k], pac_ex[k]))
  ss <- resid_j2^2
  add(sprintf("Low-use blocks: %.1f%% of days, %.1f%% of the residual sum of squares",
              100 * mean(block), 100 * sum(ss[block]) / sum(ss)))
  add(sprintf("Residual sd all days %.2f, excluding blocks %.2f", sd(resid_j2), sd(resid_j2[!block])))
  add("")

  # ---------------- residual spread by month ----------------
  add("RESIDUAL SD BY MONTH")
  add(sprintf("%4s %8s %12s", "mo", "sd", "sd ex-block"))
  sd_m <- numeric(12)
  sd_m_ex <- numeric(12)
  for (m in 1:12) {
    sel <- month == m
    sd_m[m] <- sd(resid_j2[sel])
    sd_m_ex[m] <- if (any(sel & !block)) sd(resid_j2[sel & !block]) else NA_real_
    add(sprintf("%4d %8.2f %12.2f", m, sd_m[m], sd_m_ex[m]))
  }

  # ==========================================================================
  # [FIGURE 1]  data overview  ->  R_figures/01d_overview.png
  #   row 1: (a) daily series            (from 01_series_and_coverage, top panel)
  #   row 2: (b) daily data coverage     (from 01_series_and_coverage, bottom panel)
  #   row 3: (c) annual shape by year    (from 01d_seasonal)
  #          (d) weekday vs weekend      (Mon-Fri vs Sat-Sun): the headline
  #          (e) weekly pattern by day   (from 01_exploratory_season_week): the detail
  # ==========================================================================
  year_colors <- c("2007" = "#6c3483", "2008" = "#1f4e79", "2009" = "#148f77", "2010" = "#b9770e")
  clim <- clim_day(obs$date)
  smooths <- lapply(names(year_colors), function(yr) {
    sel <- year == as.integer(yr)
    smooth_on_clim(clim[sel], y[sel], window = 15, circular = yr != "2010")
  })
  mids <- integer(12)
  for (m in 1:12) {
    days <- ref_dates[as.integer(format(ref_dates, "%m")) == m]
    mids[m] <- as.integer(format(days[length(days) %/% 2 + 1L], "%j"))
  }

  png(here("R_figures", "01d_overview.png"), width = 14, height = 10.5, units = "in", res = 140)
  layout(matrix(c(1, 1, 1,
                  2, 2, 2,
                  3, 4, 5), nrow = 3, byrow = TRUE),
         heights = c(1, 0.6, 1.2), widths = c(1.5, 0.6, 1.1))
  par(mar = c(4, 4, 3, 1))

  # (a) daily series on the full calendar (missing days appear as gaps)
  plot(d$date, d$kwh, type = "l", lwd = 0.5, col = "#1f4e79", xlab = "", ylab = "kWh / day",
       main = "(a) Daily household electricity consumption, Dec 2006 - Nov 2010")
  rug(d$date[!d$observed], col = "#b03a2e", lwd = 1.5)
  legend("topright", legend = c("daily kWh", "missing day"), col = c("#1f4e79", "#b03a2e"),
         lwd = c(1, 1.5), cex = 0.75, bty = "n")

  # (b) daily data coverage: fraction of the 1,440 minutes with a reading; days below the
  #     dashed 95% line are treated as missing (d = full calendar from daily.csv)
  plot(d$date, d$coverage, type = "h", lwd = 0.6, col = "#b03a2e",
       ylab = "fraction of 1440 min observed", xlab = "", main = "(b) Daily data coverage", ylim = c(0, 1))
  abline(h = MIN_COVERAGE, lty = 2)
  legend("bottomright", legend = sprintf("retention threshold %.2f", MIN_COVERAGE), lty = 2, bty = "n", cex = 0.8)

  # (c) annual shape by year, with one- and two-harmonic fits and monthly means
  plot(NA, xlim = c(1, 365), ylim = range(c(unlist(smooths), unlist(curves), month_means), na.rm = TRUE),
       xaxt = "n", xlab = "", ylab = "kWh / day", main = "(c) Annual shape by year")
  for (i in seq_along(smooths)) lines(1:365, smooths[[i]], col = year_colors[[i]], lwd = 1.15)
  lines(1:365, curves[[1]], col = "#2c3e50", lwd = 1.15, lty = 2)
  lines(1:365, curves[[2]], col = "#b03a2e", lwd = 2)
  points(mids, month_means, pch = 16, cex = 0.6)
  axis(1, at = mids, labels = substr(MONTHS, 1, 1))
  legend("bottomleft", legend = c(names(year_colors), "one harmonic", "two harmonics", "month mean"),
         col = c(year_colors, "#2c3e50", "#b03a2e", "black"), lty = c(1, 1, 1, 1, 2, 1, NA),
         pch = c(NA, NA, NA, NA, NA, NA, 16), lwd = c(1, 1, 1, 1, 1, 2, NA), cex = 0.65, ncol = 2, bty = "n")

  # (d) weekday vs weekend (dow 0-4 = Mon-Fri, 5-6 = Sat-Sun); same y-axis as (e) for comparison
  boxplot(split(y, factor(ifelse(dow >= 5, "Weekend", "Weekday"), levels = c("Weekday", "Weekend"))),
          col = "#d6eaf8", ylim = range(y), ylab = "kWh / day", main = "(d) Weekday vs weekend")

  # (e) weekly pattern by day: do individual days differ beyond the weekday/weekend split?
  boxplot(split(y, factor(dow, levels = 0:6)), names = DOW_NAMES, col = "#d6eaf8",
          ylim = range(y), ylab = "kWh / day", main = "(e) Weekly pattern by day")
  dev.off()

  # ==========================================================================
  # [FIGURE 2]  residual checks on the kWh scale  ->  R_figures/01d_residual_checks.png
  # resid_j2 = kWh residuals after trend + 2 harmonics + day of week (defined above);
  # fit_j2 = y - resid_j2 (also defined above)
  # ==========================================================================
  png(here("R_figures", "01d_residual_checks.png"), width = 12, height = 7.5, units = "in", res = 140)
  layout(matrix(c(1, 1,
                  2, 3), nrow = 2, byrow = TRUE))
  par(mar = c(4, 4, 3, 1))

  # (a) residuals over time
  plot(obs$date, resid_j2, type = "h", col = ifelse(block, "#b03a2e", "grey55"),
       xlab = "", ylab = "residual (kWh)",
       main = "(a) Residuals over time after trend, two harmonics, day of week")
  abline(h = 0)
  sm <- lowess(as.numeric(obs$date), resid_j2, f = 0.1)     # should stay near 0
  lines(as.Date(sm$x, origin = "1970-01-01"), sm$y, col = "#1f4e79", lwd = 2)
  legend("topright", legend = c("residual", "low-use block day", "smoothed residual"),
         col = c("grey55", "#b03a2e", "#1f4e79"), lwd = c(1, 1, 2), cex = 0.75, bty = "n")

  # (b) residual vs fitted
  plot(fit_j2, resid_j2, pch = 16, cex = 0.4, col = adjustcolor("#1f4e79", 0.4),
       xlab = "fitted (kWh)", ylab = "residual (kWh)", main = "(b) Residual vs fitted")
  abline(h = 0)
  lines(lowess(fit_j2, resid_j2), col = "#b03a2e", lwd = 1.5)

  # (c) normal Q-Q
  qqnorm(resid_j2 / sd(resid_j2), pch = 16, cex = 0.4, col = adjustcolor("#1f4e79", 0.5),
         main = sprintf("(c) Normal Q-Q (skew %.2f)", skew_adj(resid_j2)))
  qqline(resid_j2 / sd(resid_j2))
  dev.off()

  # ==========================================================================
  # [FIGURE 3]  residual dependence and spread  ->  R_figures/01d_residual.png
  #   (a) ACF, (b) PACF: all days vs absence blocks removed; (c) residual SD by month
  #   Uses rho, rho_ex, pac, pac_ex, sd_m, sd_m_ex computed in [SUMMARY]
  # ==========================================================================
  phi <- rho[1]
  lags <- 1:21
  ci <- 1.96 / sqrt(n)
  png(here("R_figures", "01d_residual.png"), width = 15, height = 4, units = "in", res = 140)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  bp <- barplot(rho, names.arg = lags, col = "#1f4e79", xlab = "lag (days)", ylab = "residual autocorrelation",
                main = "(a) Residual ACF", ylim = range(c(rho, rho_ex, phi^lags, ci, -ci)))
  lines(bp, rho_ex, type = "b", pch = 16, cex = 0.5, col = "#b03a2e")
  lines(bp, phi^lags, lty = 2, col = "#7f8c8d")
  abline(h = c(ci, -ci), col = "gray", lty = 3)
  legend("topright", legend = c("all retained days", "low-use blocks removed", expression(paste("AR(1) from ", hat(rho)[1]))),
         col = c("#1f4e79", "#b03a2e", "#7f8c8d"), lty = c(NA, 1, 2), pch = c(15, 16, NA), cex = 0.7, bty = "n")

  bp <- barplot(pac, names.arg = lags, col = "#1f4e79", xlab = "lag (days)", ylab = "partial autocorrelation",
                main = "(b) Residual PACF", ylim = range(c(pac, pac_ex, ci, -ci)))
  lines(bp, pac_ex, type = "b", pch = 16, cex = 0.5, col = "#b03a2e")
  abline(h = c(ci, -ci), col = "gray", lty = 3)
  legend("topright", legend = c("all retained days", "low-use blocks removed"),
         col = c("#1f4e79", "#b03a2e"), lty = c(NA, 1), pch = c(15, 16), cex = 0.7, bty = "n")

  bp <- barplot(sd_m, names.arg = substr(MONTHS, 1, 1), col = "#1f4e79", ylab = "residual SD (kWh)",
                main = "(c) Residual SD by month", ylim = c(0, 1.05 * max(c(sd_m, sd_m_ex), na.rm = TRUE)))
  points(bp, sd_m_ex, pch = 16, col = "#b03a2e")
  legend("topright", legend = c("all retained days", "low-use blocks removed"),
         col = c("#1f4e79", "#b03a2e"), pch = c(15, 16), cex = 0.75, bty = "n")
  dev.off()

  txt <- paste(lines, collapse = "\n")
  out <- here("results", "01d_eda.txt")
  writeLines(txt, out)
  cat(txt, "\n", sep = "")
  cat(sprintf("\nwrote %s, three figures in R_figures/, Table 1 in results/tables/\n", out))
}

main()
