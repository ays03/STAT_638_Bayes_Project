#!/usr/bin/env Rscript
# Step 4: how accurately can future demand be predicted, and how does predictive
# uncertainty grow with the forecast horizon?
#
# Models are refitted on data through 30 Nov 2009 and never see the final year.
#
#   A. SINGLE ORIGIN, 12-month forecast from 30 Nov 2009.
#   B. ROLLING ORIGIN, horizons 1-60 from every 3rd day of the holdout year.
#
#   Rscript src/04_holdout.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

ROOT <- script_root()
source(file.path(ROOT, "src", "model.R"))

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

SPECS <- list(
  list(id = "M1", J = 1, p = 0, heavy = FALSE, scale = "log", label = "J=1, iid Gaussian, log"),
  list(id = "M2", J = 1, p = 1, heavy = TRUE, scale = "log", label = "J=1, AR(1)+t, log"),
  list(id = "M3", J = 2, p = 1, heavy = TRUE, scale = "log", label = "J=2, AR(1)+t, log"),
  list(id = "M4", J = 3, p = 1, heavy = TRUE, scale = "log", label = "J=3, AR(1)+t, log"),
  list(id = "M5", J = 2, p = 2, heavy = TRUE, scale = "log", label = "J=2, AR(2)+t, log"),
  list(id = "M6", J = 2, p = 1, heavy = FALSE, scale = "log", label = "J=2, AR(1)+Gaussian, log"),
  list(id = "M7", J = 2, p = 1, heavy = TRUE, scale = "identity", label = "J=2, AR(1)+t, kWh")
)
MCMC <- list(n_iter = 30000, burn = 8000, thin = 5, n_chains = 4)
BUCKETS <- list(c(1, 7), c(8, 30), c(31, 90), c(91, 365))
ROLL_H <- 60L
ROLL_STEP <- 3L
TRAIN_END <- as.Date("2009-11-30")

ar_psi <- function(phi, h) {
  p <- length(phi)
  psi <- numeric(max(h, 1L))
  psi[1] <- 1
  if (length(psi) >= 2) {
    for (j in 2:length(psi)) {
      kk <- seq_len(min(p, j - 1L))
      psi[j] <- sum(phi[kk] * psi[j - kk])
    }
  }
  psi
}

climatology_baseline <- function(train, test) {
  tr <- train[train$observed, ]
  preds <- numeric(nrow(test))
  sds <- numeric(nrow(test))
  for (i in seq_len(nrow(test))) {
    doy <- test$doy[i]
    wk <- test$dow[i] >= 5
    dd <- pmin(abs(tr$doy - doy), 365 - abs(tr$doy - doy))
    sel <- tr[dd <= 7 & ((tr$dow >= 5) == wk), ]
    if (nrow(sel) < 5) sel <- tr[dd <= 14, ]
    lz <- log(sel$kwh)
    preds[i] <- mean(lz)
    sds[i] <- sd(lz)
  }
  list(mu = preds, sd = sds)
}

main <- function() {
  OUT <- file.path(ROOT, "results", "04_holdout")
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  d <- read_daily(file.path(ROOT, "data", "daily.csv"))
  train <- d[d$date <= TRAIN_END, ]
  test <- d[d$date > TRAIN_END, ]
  rownames(train) <- NULL
  rownames(test) <- NULL
  t_center <- mean(train$t)
  cat(sprintf("train: %s to %s (%d days, %d observed)\n",
              min(train$date), max(train$date), nrow(train), sum(train$observed)))
  cat(sprintf("test : %s to %s (%d days, %d observed)\n",
              min(test$date), max(test$date), nrow(test), sum(test$observed)))

  y_test <- test$kwh
  obs_test <- test$observed
  h_arr <- seq_len(nrow(test))
  MAX_P <- max(vapply(SPECS, function(s) s$p, numeric(1)))
  origins <- integer(0)
  for (o in seq(1L, nrow(test) - 1L, by = ROLL_STEP)) {
    if (o >= MAX_P && all(obs_test[o - (0:(MAX_P - 1L))])) origins <- c(origins, o)
  }
  cat(sprintf("%d rolling origins shared by all models (every %drd day, horizons 1-%d)\n",
              length(origins), ROLL_STEP, ROLL_H))

  rows_single <- list()
  roll_store <- list()
  fits <- list()

  for (spec in SPECS) {
    t0 <- proc.time()
    mid <- spec$id
    cat(sprintf("\nfitting %s on training data ...\n", mid))
    fitres <- fit_chains(
      train$kwh, train$observed, train$t, train$dow,
      seed = 5000 + 7 * as.integer(substring(mid, 2)),
      n_paths = 1,
      J = spec$J, p = spec$p, heavy = spec$heavy, scale = spec$scale,
      n_iter = MCMC$n_iter, burn = MCMC$burn, thin = MCMC$thin, n_chains = MCMC$n_chains
    )
    fits[[mid]] <- fitres

    fc <- forecast(fitres, test$t, test$dow, t_center, n_draws = 4000, seed = 99 + as.integer(substring(mid, 2)))
    lpd <- log_pred_density(fc, y_test)
    cr <- crps(fc$y, y_test)
    pit_vals <- pit(fc, y_test)
    q <- apply(fc$y, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975))
    cov95 <- (y_test >= q[1, ]) & (y_test <= q[4, ])
    cov50 <- (y_test >= q[2, ]) & (y_test <= q[3, ])
    width95 <- q[4, ] - q[1, ]
    single <- list(
      lpd = lpd, crps = cr, q = q, cov95 = cov95, cov50 = cov50, width95 = width95, pit = pit_vals,
      ypred_mean = colMeans(fc$y), ypred_med = apply(fc$y, 2, median)
    )
    saveRDS(single, file.path(OUT, paste0("single_", mid, ".rds")))
    rows_single[[mid]] <- list(
      id = mid, label = spec$label, lpd = lpd, crps = cr,
      cov95 = cov95, cov50 = cov50, width95 = width95, obs = obs_test, q = q,
      ypred_med = single$ypred_med, pit = pit_vals
    )

    beta <- flatten(fitres, "beta")
    p <- spec$p
    Xtest <- build_design(test$t, test$dow, spec$J, t_center = t_center)$X
    z_test <- if (spec$scale == "log") log(y_test) else y_test
    mu_te <- beta %*% t(Xtest)
    eps_test <- mu_te
    eps_test[] <- rep(z_test, each = nrow(mu_te)) - as.numeric(mu_te)

    acc <- lapply(seq_len(ROLL_H), function(h) {
      list(lpd = numeric(0), crps = numeric(0), c95 = logical(0), c50 = logical(0),
           w95 = numeric(0), w50 = numeric(0), day = integer(0))
    })
    set.seed(as.integer(substring(mid, 2)))
    n_sub <- 1500L
    sub <- sample.int(nrow(beta), min(n_sub, nrow(beta)))
    for (oi in seq_along(origins)) {
      o <- origins[oi]
      H <- min(ROLL_H, nrow(test) - o)
      if (H <= 0) next
      idx <- (o + 1L):(o + H)
      e_init <- if (p > 0) eps_test[sub, o - (0:(p - 1L)), drop = FALSE] else NULL
      fcr <- forecast(
        fitres, test$t[idx], test$dow[idx], t_center,
        eps_init = e_init, draw_idx = sub, seed = 7 * (oi - 1L) + 3L
      )
      lp <- log_pred_density(fcr, y_test[idx])
      cc <- crps(fcr$y, y_test[idx])
      qq <- apply(fcr$y, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975))
      for (j in seq_len(H)) {
        if (!obs_test[idx[j]]) next
        h <- j
        acc[[h]]$lpd <- c(acc[[h]]$lpd, lp[j])
        acc[[h]]$crps <- c(acc[[h]]$crps, cc[j])
        acc[[h]]$c95 <- c(acc[[h]]$c95, qq[1, j] <= y_test[idx[j]] && y_test[idx[j]] <= qq[4, j])
        acc[[h]]$c50 <- c(acc[[h]]$c50, qq[2, j] <= y_test[idx[j]] && y_test[idx[j]] <= qq[3, j])
        acc[[h]]$w95 <- c(acc[[h]]$w95, qq[4, j] - qq[1, j])
        acc[[h]]$w50 <- c(acc[[h]]$w50, qq[3, j] - qq[2, j])
        acc[[h]]$day <- c(acc[[h]]$day, idx[j])
      }
    }
    roll_store[[mid]] <- acc
    if (mid != SPECS[[1]]$id) {
      first <- roll_store[[SPECS[[1]]$id]]
      for (h in seq_len(ROLL_H)) {
        if (!identical(as.integer(acc[[h]]$day), as.integer(first[[h]]$day))) {
          stop(sprintf("%s scored different days than %s at h=%d", mid, SPECS[[1]]$id, h))
        }
      }
    }
    saveRDS(acc, file.path(OUT, paste0("roll_", mid, ".rds")))
    cat(sprintf("   %s done in %.1fs  (%d rolling origins, n per horizon ~%d)\n",
                mid, (proc.time() - t0)[["elapsed"]], length(origins), length(acc[[1]]$lpd)))
  }

  cl <- climatology_baseline(train, test)
  z_act <- log(ifelse(obs_test, y_test, NA_real_))
  cl_lpd <- dnorm(z_act, cl$mu, cl$sd, log = TRUE) - z_act
  qs <- qnorm(c(0.025, 0.25, 0.75, 0.975))
  cl_q <- outer(qs, cl$sd, `*`)
  cl_q[] <- cl_q + rep(cl$mu, each = length(qs))
  cl_q <- exp(cl_q)
  cl_cov95 <- (y_test >= cl_q[1, ]) & (y_test <= cl_q[4, ])
  set.seed(0)
  cl_samp <- matrix(rnorm(4000 * nrow(test)), nrow = 4000)
  cl_samp[] <- exp(as.numeric(cl_samp) * rep(cl$sd, each = 4000) + rep(cl$mu, each = 4000))
  cl_crps <- crps(cl_samp, ifelse(obs_test, y_test, NA_real_))

  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("OUT-OF-SAMPLE FORECAST EVALUATION")
  add(sprintf("Trained on %d observed days through %s;", sum(train$observed), TRAIN_END))
  add(sprintf("scored on %d observed days from %s to %s.", sum(obs_test), min(test$date), max(test$date)))
  add("Scores are per day on the kWh scale. lpd: higher is better (log predictive")
  add("density). CRPS: lower is better, in kWh. Coverage should match the nominal level.")
  add("")
  add("=== DESIGN A: single origin, 12-month-ahead forecast ===")
  add("")
  add(sprintf("%4s %28s %9s %7s %7s %7s %13s", "id", "specification", "mean lpd", "CRPS", "cov95", "cov50", "mean width95"))
  for (r in rows_single) {
    m <- r$obs
    add(sprintf("%4s %28s %9.3f %7.3f %6.1f%% %6.1f%% %12.2f",
                r$id, r$label, mean(r$lpd, na.rm = TRUE), mean(r$crps, na.rm = TRUE),
                100 * mean(r$cov95[m]), 100 * mean(r$cov50[m]), mean(r$width95[m])))
  }
  add(sprintf("%4s %28s %9.3f %7.3f %6.1f%% %7s %12.2f",
              "--", "day-of-year climatology", mean(cl_lpd, na.rm = TRUE), mean(cl_crps, na.rm = TRUE),
              100 * mean(cl_cov95[obs_test]), "", mean((cl_q[4, ] - cl_q[1, ])[obs_test])))
  add("")
  add("Scores by horizon bucket (CRPS in kWh, coverage of the 95% interval):")
  hdr <- sprintf("%4s", "id")
  for (b in BUCKETS) hdr <- paste0(hdr, sprintf(" | %18s", sprintf("%d-%dd", b[1], b[2])))
  add(hdr)
  add(paste0(sprintf("%4s", ""), paste(rep(sprintf(" | %18s", "CRPS  lpd  cov"), length(BUCKETS)), collapse = "")))
  for (r in rows_single) {
    line <- sprintf("%4s", r$id)
    for (b in BUCKETS) {
      m <- r$obs & h_arr >= b[1] & h_arr <= b[2]
      if (!any(m)) line <- paste0(line, sprintf(" | %18s", "--"))
      else line <- paste0(line, sprintf(" | %5.2f %5.2f %5.1f%%",
                                        mean(r$crps[m], na.rm = TRUE), mean(r$lpd[m], na.rm = TRUE),
                                        100 * mean(r$cov95[m])))
    }
    add(line)
  }
  add("")
  add("  Caveat: the 361 horizons come from a single realisation of the holdout")
  add("  year, so these bucket scores are correlated. Design B is the replicated")
  add("  version and is the one to quote for the horizon question.")
  add("")
  add("=== DESIGN B: rolling origins, horizons 1-60 ===")
  add("")
  hs_show <- c(1, 2, 3, 5, 7, 14, 21, 30, 45, 60)
  add(paste0(sprintf("%4s ", "id"), paste(sprintf("%7s", paste0("h=", hs_show)), collapse = " ")))
  blocks <- list(
    list(key = "crps", title = "mean CRPS (kWh, lower better)", pct = FALSE, fmt = "%7.3f"),
    list(key = "lpd", title = "mean log predictive density (higher better)", pct = FALSE, fmt = "%7.3f"),
    list(key = "w95", title = "mean 95% interval width (kWh)", pct = FALSE, fmt = "%7.2f"),
    list(key = "c95", title = "coverage of 95% interval", pct = TRUE, fmt = "%6.1f%%"),
    list(key = "c50", title = "coverage of 50% interval", pct = TRUE, fmt = "%6.1f%%")
  )
  for (blk in blocks) {
    add("")
    add(paste0("  ", blk$title))
    for (spec in SPECS) {
      vals <- vapply(hs_show, function(h) {
        v <- roll_store[[spec$id]][[h]][[blk$key]]
        if (blk$pct) 100 * mean(v) else mean(v, na.rm = TRUE)
      }, numeric(1))
      add(paste0(sprintf("%4s ", spec$id), paste(sprintf(blk$fmt, vals), collapse = " ")))
    }
  }
  add("")
  add("HOW UNCERTAINTY GROWS WITH HORIZON (best model by design-A lpd)")
  best_id <- names(rows_single)[which.max(vapply(rows_single, function(r) mean(r$lpd, na.rm = TRUE), numeric(1)))]
  rs <- roll_store[[best_id]]
  w1 <- mean(rs[[1]]$w95)
  add(sprintf("  model %s", best_id))
  add(sprintf("%4s %9s %13s %7s %7s", "h", "width95", "ratio to h=1", "CRPS", "cov95"))
  for (h in c(1, 2, 3, 4, 5, 7, 10, 14, 21, 30, 45, 60)) {
    w <- mean(rs[[h]]$w95)
    add(sprintf("%4d %9.2f %13.3f %7.3f %6.1f%%", h, w, w / w1, mean(rs[[h]]$crps, na.rm = TRUE), 100 * mean(rs[[h]]$c95)))
  }
  add("")
  add("  The saturation point is the answer to the horizon question: once the")
  add("  ratio stops rising, recent consumption has stopped being informative and")
  add("  the forecast is seasonal-plus-weekly climatology.")
  add("")
  add("MODELLING DECISIONS RE-TESTED OUT OF SAMPLE (paired differences in lpd)")
  add("  positive favours the first model; se is the paired standard error")
  paired <- function(a, b, what, design = "A") {
    if (design == "A") {
      m <- rows_single[[a]]$obs
      dv <- rows_single[[a]]$lpd[m] - rows_single[[b]]$lpd[m]
    } else {
      dv <- unlist(lapply(seq_len(ROLL_H), function(h) roll_store[[a]][[h]]$lpd - roll_store[[b]][[h]]$lpd))
    }
    dlt <- mean(dv, na.rm = TRUE)
    se <- sd(dv, na.rm = TRUE) / sqrt(sum(is.finite(dv)))
    verdict <- if (dlt > 2 * se) "supported" else if (dlt < -2 * se) "no support" else "inconclusive"
    add(sprintf("  [%s] %-42s %s-%s = %+.4f +/- %.4f  %s", design, what, a, b, dlt, se, verdict))
  }
  for (design in c("A", "B")) {
    paired("M2", "M1", "AR(1)+t vs iid Gaussian", design)
    paired("M3", "M2", "second harmonic (J=2 vs J=1)", design)
    paired("M4", "M3", "third harmonic (J=3 vs J=2)", design)
    paired("M5", "M3", "AR(2) vs AR(1)", design)
    paired("M3", "M6", "t errors vs Gaussian", design)
    paired("M3", "M7", "log scale vs kWh scale", design)
    add("")
  }
  txt <- paste(lines, collapse = "\n")
  cat("\n", txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "04_holdout.txt"))

  hs <- seq_len(ROLL_H)
  png(file.path(ROOT, "figures", "04_horizon.png"), width = 13, height = 7.5, units = "in", res = 140)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  cols <- c("#1f4e79", "#b03a2e", "#1e8449", "#b9770e", "#6c3483", "#1a5276", "#7f8c8d")
  plot(NA, xlim = c(1, ROLL_H), ylim = range(unlist(lapply(SPECS, function(s) vapply(hs, function(h) mean(roll_store[[s$id]][[h]]$w95), numeric(1))))),
       xlab = "forecast horizon h (days)", ylab = "mean width of 95% interval (kWh)",
       main = "(a) Predictive uncertainty vs horizon (rolling origins)")
  for (i in seq_along(SPECS)) lines(hs, vapply(hs, function(h) mean(roll_store[[SPECS[[i]]$id]][[h]]$w95), numeric(1)), col = cols[i], lwd = 1.3)
  legend("topleft", legend = vapply(SPECS, function(s) s$id, character(1)), col = cols, lwd = 1.3, cex = 0.7, ncol = 2, bty = "n")

  w <- vapply(hs, function(h) mean(rs[[h]]$w95), numeric(1))
  plot(hs, w / w[1], type = "b", pch = 16, cex = 0.5, col = "#1f4e79",
       xlab = "forecast horizon h (days)", ylab = "width relative to h = 1",
       main = "(b) Uncertainty saturates within about a week")
  phi_bar <- colMeans(flatten(fits[[best_id]], "phi"))
  if (length(phi_bar) > 0 && any(is.finite(phi_bar))) {
    theo <- vapply(hs, function(h) sqrt(sum(ar_psi(phi_bar, h)^2)), numeric(1))
    lines(hs, theo / theo[1], lty = 2, col = "#b03a2e")
    abline(h = marginal_sd_ratio(phi_bar), lty = 3, col = "gray")
    legend("bottomright", legend = c(paste(best_id, "empirical"), "AR(p) theory", "stationary limit"),
           col = c("#1f4e79", "#b03a2e", "gray"), lty = c(1, 2, 3), pch = c(16, NA, NA), cex = 0.75, bty = "n")
  }
  plot(NA, xlim = c(1, ROLL_H),
       ylim = range(c(mean(cl_crps, na.rm = TRUE), unlist(lapply(SPECS, function(s) vapply(hs, function(h) mean(roll_store[[s$id]][[h]]$crps, na.rm = TRUE), numeric(1)))))),
       xlab = "forecast horizon h (days)", ylab = "mean CRPS (kWh)", main = "(c) Forecast accuracy vs horizon")
  for (i in seq_along(SPECS)) lines(hs, vapply(hs, function(h) mean(roll_store[[SPECS[[i]]$id]][[h]]$crps, na.rm = TRUE), numeric(1)), col = cols[i])
  abline(h = mean(cl_crps, na.rm = TRUE), lty = 2)
  legend("topleft", legend = c(vapply(SPECS, function(s) s$id, character(1)), "climatology"),
         col = c(cols, "black"), lty = c(rep(1, length(SPECS)), 2), cex = 0.65, ncol = 2, bty = "n")
  plot(NA, xlim = c(1, ROLL_H), ylim = c(60, 102), xlab = "forecast horizon h (days)",
       ylab = "empirical coverage of 95% interval (%)", main = "(d) Calibration vs horizon")
  for (i in seq_along(SPECS)) lines(hs, vapply(hs, function(h) 100 * mean(roll_store[[SPECS[[i]]$id]][[h]]$c95), numeric(1)), col = cols[i])
  abline(h = 95, lty = 2)
  legend("bottomleft", legend = vapply(SPECS, function(s) s$id, character(1)), col = cols, lwd = 1, cex = 0.65, ncol = 2, bty = "n")
  dev.off()

  z <- rows_single[[best_id]]
  png(file.path(ROOT, "figures", "04_holdout_fan.png"), width = 12, height = 4.2, units = "in", res = 140)
  par(mar = c(4, 4, 3, 1))
  tr_i <- max(1, nrow(train) - 199):nrow(train)
  ylim <- range(c(train$kwh[tr_i], y_test, z$q), na.rm = TRUE)
  plot(train$date[tr_i], train$kwh[tr_i], type = "l", lwd = 0.8, col = "#444444",
       xlim = c(train$date[tr_i[1]], max(test$date)), ylim = ylim,
       xlab = "", ylab = "kWh / day",
       main = sprintf("12-month-ahead forecast from %s, model %s\n(nothing after the dashed line was used in fitting)", TRAIN_END, best_id))
  polygon(c(test$date, rev(test$date)), c(z$q[1, ], rev(z$q[4, ])), col = adjustcolor("#1f4e79", 0.20), border = NA)
  polygon(c(test$date, rev(test$date)), c(z$q[2, ], rev(z$q[3, ])), col = adjustcolor("#1f4e79", 0.40), border = NA)
  lines(test$date, z$ypred_med, col = "#1f4e79", lwd = 1.3)
  points(test$date, y_test, pch = 16, cex = 0.35, col = "#b03a2e")
  abline(v = TRAIN_END, lty = 2)
  legend("topright", legend = c("training data", "95% predictive", "50% predictive", "predictive median", "actual (held out)"),
         col = c("#444444", adjustcolor("#1f4e79", 0.35), adjustcolor("#1f4e79", 0.7), "#1f4e79", "#b03a2e"),
         lty = c(1, NA, NA, 1, NA), pch = c(NA, 15, 15, NA, 16), pt.cex = 1.2, cex = 0.7, ncol = 3, bty = "n")
  dev.off()

  show <- unique(c(best_id, c("M1", "M6", "M7")))
  show <- show[seq_len(min(4, length(show)))]
  png(file.path(ROOT, "figures", "04_pit.png"), width = 3.4 * length(show), height = 3.4, units = "in", res = 140)
  par(mfrow = c(1, length(show)), mar = c(4, 4, 3, 1))
  for (mid in show) {
    v <- rows_single[[mid]]$pit
    v <- v[is.finite(v)]
    hist(v, breaks = seq(0, 1, length.out = 21), col = "#1f4e79", xlim = c(0, 1),
         xlab = "predictive CDF at the actual", main = sprintf("%s", mid))
    abline(h = length(v) / 20, lty = 2)
  }
  dev.off()
  cat("\nwrote results/04_holdout.txt and 3 figures\n")
}

if (sys.nframe() == 0L) main()
