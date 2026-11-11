#!/usr/bin/env Rscript
# Step 9: sensitivity analysis and model checking (Report section 6).
#
# A. Prior sensitivity. The final model (J = 2, AR(1) + t errors, sqrt kWh) is
#    refitted on all data, changing one prior at a time:
#      coef_half    coefficient prior SDs halved:  intercept N(0, 10^2), others N(0, 2.5^2)
#      coef_double  coefficient prior SDs doubled: intercept N(0, 40^2), others N(0, 10^2)
#      nu_heavy     nu - 2 ~ Exponential(mean 3), favouring heavy tails
#      nu_light     nu - 2 ~ Exponential(mean 30), favouring nearly Gaussian tails
#    Seed and likelihood terms match `main` in src/06_fit_sqrt.R, so every
#    difference from the reported posterior comes from the prior. MCMC settings
#    match too, except that nu_light runs twice as long to reach Rhat < 1.01.
# B. Posterior predictive checks of the final model: whole series replicated
#    from the posterior and compared with the data on six statistics.
# C. Coverage of the held-out predictive intervals by season.
#
#   Rscript src/09_sensitivity.R    # missing fits (about 2 minutes each), then summaries
#
# Needs results/06_sqrt/main.rds and holdout.rds (src/06_fit_sqrt.R).
# Writes results/09_sensitivity/<fit>.rds (draws, not tracked), summary.rds,
# table_prior.csv, table_ppc.csv, table_season.csv and R_figures/11_ppc.png.

library(here)
source(here("src", "model.R"))
h <- new.env()
sys.source(here("src", "07_sqrt_results.R"), envir = h)   # kwh_curve(), diag_group(), ...

OUT <- here("results", "09_sensitivity")
IN6 <- here("results", "06_sqrt")
FIG <- here("R_figures")
MCMC <- list(n_iter = 40000, burn = 10000, thin = 5, n_chains = 4)
SEED_MAIN <- 6100   # seed of `main` in src/06_fit_sqrt.R
TRAIN_END <- as.Date("2009-11-30")
LOW_KWH <- 12       # absence-block threshold used in the EDA
N_REP <- 2000
NAV <- "#1f4e79"
RED <- "#b03a2e"

PRIORS <- list(
  coef_half = list(prior_tau = c(10, 2.5), label = "Coefficient SDs halved"),
  coef_double = list(prior_tau = c(40, 10), label = "Coefficient SDs doubled"),
  nu_heavy = list(nu_exp_mean = 3, label = "nu - 2 ~ Exp(mean 3)"),
  # nu mixes more slowly when its prior allows large values, so this chain is twice as long
  nu_light = list(nu_exp_mean = 30, label = "nu - 2 ~ Exp(mean 30)",
                  mcmc = list(n_iter = 80000, burn = 20000, thin = 10, n_chains = 4))
)

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}
summ <- function(v) c(mean = mean(v), lo = unname(quantile(v, 0.025)), hi = unname(quantile(v, 0.975)))

run_prior <- function(name, d) {
  pr <- PRIORS[[name]]
  mc <- if (is.null(pr$mcmc)) MCMC else pr$mcmc
  f <- fit_chains(
    d$kwh, d$observed, d$t, d$dow,
    n_chains = mc$n_chains, seed = SEED_MAIN,
    J = 2, p = 1, heavy = TRUE, scale = "sqrt",
    n_iter = mc$n_iter, burn = mc$burn, thin = mc$thin,
    ll_index = common_ll_index(d$observed, max_lag = 2),
    n_paths = 1, disperse_init = TRUE,
    prior_tau = if (is.null(pr$prior_tau)) PRIOR_TAU$sqrt else pr$prior_tau,
    nu_exp_mean = if (is.null(pr$nu_exp_mean)) 10 else pr$nu_exp_mean
  )
  f$waic <- waic(f$loglik)
  f$conv <- h$diag_group(h$param_mats(f))
  f[c("loglik", "eps_path", "lam_path", "warmup", "z_imputed", "eps_last")] <- NULL
  f
}

# The Section 5 conclusions, computed from one fit exactly as src/07_sqrt_results.R does.
conclusions <- function(f) {
  b <- flatten(f, "beta")
  nm <- f$names
  yr <- seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "day")
  win <- seq(as.Date("2009-07-01"), as.Date("2010-06-30"), by = "day")
  M2 <- h$kwh_curve(b, nm, yr, 2)
  add2 <- apply(abs(M2 - h$kwh_curve(b, nm, yr, 1)), 1, max)
  dl <- h$dow_full(b, nm)
  A <- b[, match("intercept", nm)] + h$seas_matrix(b, nm, yr, 2)
  b_tr <- b[, match("trend_per_year", nm)]
  E <- rowMeans(A^2) + 2 * dl * rowMeans(A) + dl^2
  P <- E - rowMeans(E)
  wkend <- rowMeans(P[, 6:7]) - rowMeans(P[, 1:5])
  list(
    annual = summ(rowMeans(M2)),
    slope = summ(2 * b_tr * (rowMeans(A) + rowMeans(dl))),
    p_decline = mean(b_tr < 0),
    peak = summ(apply(h$kwh_curve(b, nm, win, 2), 1, max)),
    trough = summ(apply(M2, 1, min)),
    add2 = summ(add2), p_add2_gt1 = mean(add2 > 1),
    wkend = summ(wkend), p_wkend = mean(wkend > 0),
    phi = summ(as.numeric(f$phi)), sigma = summ(as.numeric(f$sigma)), nu = summ(as.numeric(f$nu))
  )
}

# ---------------------------------------------------------------------------
# B. posterior predictive statistics (y in kWh on the full calendar, NA = missing)
# ---------------------------------------------------------------------------

longest_run <- function(x) {
  r <- rle(!is.na(x) & x < LOW_KWH)
  if (any(r$values)) max(r$lengths[r$values]) else 0L
}
ppc_stats <- function(y, mon) {
  ok <- !is.na(y)
  v <- y[ok]
  pair <- ok[-1] & ok[-length(y)]
  winter <- ok & mon %in% c(12, 1, 2)
  summer <- ok & mon %in% c(6, 7, 8)
  c(
    lag1 = cor(y[-1][pair], y[-length(y)][pair]),
    sd_ratio = sd(y[winter]) / sd(y[summer]),
    skew = mean((v - mean(v))^3) / sd(v)^3,
    max = max(v),
    n_low = sum(v < LOW_KWH),
    run_low = longest_run(y)
  )
}
PPC_LABELS <- c(
  lag1 = "Lag-1 autocorrelation of daily kWh",
  sd_ratio = "SD in Dec-Feb / SD in Jun-Aug",
  skew = "Skewness of daily kWh",
  max = "Highest day (kWh)",
  n_low = sprintf("Days below %d kWh", LOW_KWH),
  run_low = sprintf("Longest run of days below %d kWh", LOW_KWH)
)

replicate_series <- function(f, observed, n_rep, seed = 9) {
  # Whole series from the fitted model: systematic part plus a stationary
  # AR(1) error with t innovations, one posterior draw per replicate.
  set.seed(seed)
  b <- flatten(f, "beta")
  sel <- sample.int(nrow(b), n_rep)
  sig <- as.numeric(f$sigma)[sel]
  nu <- as.numeric(f$nu)[sel]
  phi <- as.numeric(f$phi)[sel]
  mu <- b[sel, , drop = FALSE] %*% t(f$X)
  nT <- ncol(mu)
  lam <- matrix(rgamma(n_rep * nT, shape = nu / 2, rate = nu / 2), n_rep, nT)
  eta <- matrix(rnorm(n_rep * nT), n_rep, nT) * sig / sqrt(lam)
  eps <- matrix(0, n_rep, nT)
  eps[, 1] <- eta[, 1] / sqrt(1 - phi^2)
  for (i in 2:nT) eps[, i] <- phi * eps[, i - 1] + eta[, i]
  y_rep <- pmax(mu + eps, 0)^2
  y_rep[, !observed] <- NA
  y_rep
}

main <- function() {
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  d <- read_daily(here("data", "daily.csv"))
  f0 <- readRDS(file.path(IN6, "main.rds"))

  # ---------------- A. prior sensitivity ----------------
  for (name in names(PRIORS)) {
    path <- file.path(OUT, paste0(name, ".rds"))
    if (file.exists(path)) next
    t0 <- proc.time()
    saveRDS(run_prior(name, d), path)
    cat(sprintf("%s done in %.0f s\n", name, (proc.time() - t0)[["elapsed"]]))
  }
  fits <- c(list(main = f0), lapply(setNames(names(PRIORS), names(PRIORS)),
                                    function(n) readRDS(file.path(OUT, paste0(n, ".rds")))))
  res <- lapply(fits, conclusions)
  conv <- t(vapply(names(PRIORS), function(n) fits[[n]]$conv, numeric(3)))
  d_elpd <- vapply(fits, function(f) {
    di <- f$waic$elpd_i - f0$waic$elpd_i
    c(diff = sum(di), se = sqrt(length(di) * var(di)))
  }, numeric(2))

  fmt <- function(s, dp = 1) sprintf(paste0("%.", dp, "f [%.", dp, "f, %.", dp, "f]"), s[1], s[2], s[3])
  col <- function(r, n) c(
    fmt(r$annual), fmt(r$slope, 2), sprintf("%.2f", r$p_decline), fmt(r$peak), fmt(r$trough),
    fmt(r$add2), fmt(r$wkend), fmt(r$phi, 2), fmt(r$nu),
    if (n == "main") "0 (reference)" else sprintf("%+.1f", d_elpd["diff", n])
  )
  tab <- data.frame(
    Quantity = c("Typical day, annual mean (kWh/day)", "Change in annual mean (kWh/day per year)",
                 "P(decline)", "Winter peak (kWh/day)", "Summer trough (kWh/day)",
                 "Largest change from second harmonic (kWh/day)", "Weekend minus weekday (kWh/day)",
                 "AR(1) coefficient $\\phi$", "Degrees of freedom $\\nu$", "Change in WAIC elpd"),
    Reported = col(res$main, "main"),
    check.names = FALSE
  )
  heads <- c(coef_half = "SDs halved", coef_double = "SDs doubled",
             nu_heavy = "$\\nu - 2$, mean 3", nu_light = "$\\nu - 2$, mean 30")
  for (n in names(PRIORS)) tab[[heads[[n]]]] <- col(res[[n]], n)
  write.csv(tab, file.path(OUT, "table_prior.csv"), row.names = FALSE)

  # largest shift of a kWh-scale posterior mean, relative to its posterior SD-ish (half the CrI width)
  kwh_keys <- c("annual", "slope", "peak", "trough", "add2", "wkend")
  shift <- sapply(names(PRIORS), function(n) sapply(kwh_keys, function(k) {
    abs(res[[n]][[k]][["mean"]] - res$main[[k]][["mean"]]) /
      ((res$main[[k]][["hi"]] - res$main[[k]][["lo"]]) / 3.92)
  }))

  # ---------------- B. posterior predictive checks ----------------
  mon <- as.integer(format(d$date, "%m"))
  y_rep <- replicate_series(f0, d$observed, N_REP)
  obs_stat <- ppc_stats(d$kwh, mon)
  rep_stat <- t(apply(y_rep, 1, ppc_stats, mon = mon))
  ppp <- colMeans(rep_stat >= rep(obs_stat, each = N_REP))
  rq <- apply(rep_stat, 2, quantile, c(0.025, 0.5, 0.975))
  dp <- c(lag1 = 2, sd_ratio = 2, skew = 2, max = 1, n_low = 0, run_low = 0)
  ppc_tab <- data.frame(
    Statistic = unname(PPC_LABELS[names(obs_stat)]),
    Observed = vapply(names(obs_stat), function(k) formatC(obs_stat[[k]], format = "f", digits = dp[[k]]), ""),
    `Replicated median [95% interval]` = vapply(names(obs_stat), function(k) {
      sprintf(paste0("%.", dp[[k]], "f [%.", dp[[k]], "f, %.", dp[[k]], "f]"), rq[2, k], rq[1, k], rq[3, k])
    }, ""),
    `P(replicated $\\geq$ observed)` = sprintf("%.3f", ppp),
    check.names = FALSE
  )
  write.csv(ppc_tab, file.path(OUT, "table_ppc.csv"), row.names = FALSE)

  png(file.path(FIG, "11_ppc.png"), width = 12, height = 6.5, units = "in", res = 140)
  par(mfrow = c(2, 3), mar = c(4.2, 4, 3, 1))
  for (k in names(obs_stat)) {
    v <- rep_stat[, k]
    xl <- range(c(quantile(v, c(0.005, 0.995)), obs_stat[[k]]))
    br <- if (k == "run_low") seq(min(v, obs_stat[[k]]) - 0.5, max(v, obs_stat[[k]]) + 0.5, by = 1) else 60
    hist(v, breaks = br, col = adjustcolor(NAV, 0.35), border = "white", xlim = xl,
         main = PPC_LABELS[[k]], xlab = sprintf("p = %.3f", ppp[[k]]), ylab = "replicated series", cex.main = 0.95)
    abline(v = obs_stat[[k]], col = RED, lwd = 2.5)
  }
  dev.off()

  # ---------------- C. held-out interval coverage by season ----------------
  ho <- readRDS(file.path(IN6, "holdout.rds"))
  test <- d[d$date > TRAIN_END, ]
  r <- ho$roll[ho$roll$observed, ]
  r$date <- test$date[r$origin + r$h]
  m <- as.integer(format(r$date, "%m"))
  r$season <- factor(ifelse(m %in% c(12, 1, 2), "Winter (Dec-Feb)",
                     ifelse(m %in% 3:5, "Spring (Mar-May)",
                     ifelse(m %in% 6:8, "Summer (Jun-Aug)", "Autumn (Sep-Nov)"))),
                     levels = c("Winter (Dec-Feb)", "Spring (Mar-May)", "Summer (Jun-Aug)", "Autumn (Sep-Nov)"))
  low <- d$observed & d$kwh < LOW_KWH
  rl <- rle(low)
  ends <- cumsum(rl$lengths)
  absence <- rep(FALSE, nrow(d))
  for (i in which(rl$values & rl$lengths >= 5)) absence[(ends[i] - rl$lengths[i] + 1):ends[i]] <- TRUE
  r$absence <- absence[match(r$date, d$date)]
  cover <- function(x) c(
    n = nrow(x), c95 = mean(x$kwh >= x$lo95 & x$kwh <= x$hi95),
    c50 = mean(x$kwh >= x$lo50 & x$kwh <= x$hi50), w95 = mean(x$hi95 - x$lo95)
  )
  by_season <- t(sapply(split(r, r$season), cover))
  by_season_ex <- t(sapply(split(r[!r$absence, ], r$season[!r$absence]), cover))
  season_tab <- data.frame(
    Season = rownames(by_season),
    `Forecasts scored` = sprintf("%.0f", by_season[, "n"]),
    `95% coverage` = sprintf("%.1f%%", 100 * by_season[, "c95"]),
    `95% coverage, absence days removed` = sprintf("%.1f%%", 100 * by_season_ex[, "c95"]),
    `50% coverage` = sprintf("%.1f%%", 100 * by_season[, "c50"]),
    `Mean 95% width (kWh)` = sprintf("%.1f", by_season[, "w95"]),
    check.names = FALSE
  )
  write.csv(season_tab, file.path(OUT, "table_season.csv"), row.names = FALSE)

  S <- list(
    priors = PRIORS, res = res, conv = conv, d_elpd = d_elpd, max_shift = max(shift), shift = shift,
    ppc = list(obs = obs_stat, rep_q = rq, ppp = ppp, n_rep = N_REP),
    season = by_season, season_ex = by_season_ex,
    absence_holdout_days = sum(absence[d$date > TRAIN_END])
  )
  saveRDS(S, file.path(OUT, "summary.rds"))
  invisible(S)
}

if (sys.nframe() == 0L) {
  S <- main()
  print(S$res$main$nu); print(S$conv); print(S$d_elpd); print(round(S$shift, 2))
  print(S$ppc$obs); print(S$ppc$rep_q); print(S$ppc$ppp); print(S$season); print(S$season_ex)
}
