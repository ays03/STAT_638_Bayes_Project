#!/usr/bin/env Rscript
# Step 1b: justify the modelling scale and the error structure using OLS pre-fits.
#
# Compares raw-kWh, sqrt and log responses on the same mean structure
# (intercept + linear trend + J harmonics + day-of-week), and reports residual
# skewness, kurtosis, heteroscedasticity and the residual autocorrelation
# function. This is the evidence used to justify the Bayesian model in step 2.
#
#   Rscript src/01b_transform_check.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

design <- function(t, dow, J) {
  cols <- list(rep(1, length(t)), t / 365)
  xnames <- c("intercept", "trend_per_year")
  if (J >= 1) {
    for (j in seq_len(J)) {
      cols[[length(cols) + 1]] <- sin(2 * pi * j * t / 365)
      cols[[length(cols) + 1]] <- cos(2 * pi * j * t / 365)
      xnames <- c(xnames, paste0("sin", j), paste0("cos", j))
    }
  }
  for (k in 1:6) {
    cols[[length(cols) + 1]] <- as.numeric(dow == k)
    xnames <- c(xnames, paste0("dow_", k))
  }
  list(X = do.call(cbind, cols), names = xnames)
}

acf_np <- function(x, nlags) {
  x <- x - mean(x)
  denom <- sum(x * x)
  vapply(seq_len(nlags), function(k) sum(x[seq_len(length(x) - k)] * x[(k + 1):length(x)]) / denom, numeric(1))
}

skew_raw <- function(x) {
  xc <- x - mean(x)
  m2 <- mean(xc^2)
  mean(xc^3) / m2^1.5
}

kurt_raw <- function(x) {
  xc <- x - mean(x)
  m2 <- mean(xc^2)
  mean(xc^4) / m2^2
}

main <- function() {
  ROOT <- script_root()
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  dir.create(file.path(ROOT, "results"), showWarnings = FALSE)
  d <- read_daily(file.path(ROOT, "data", "daily.csv"))
  obs <- d[d$observed, ]
  obs <- obs[order(obs$t), ]
  t <- as.numeric(obs$t)
  dow <- obs$dow
  des <- design(t, dow, J = 2)
  X <- des$X

  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("RESPONSE-SCALE COMPARISON (OLS, intercept + trend + 2 harmonics + day-of-week)")
  add(sprintf("%8s %7s %11s %11s %10s %7s", "scale", "R2", "resid skew", "resid kurt", "BP het. p", "rho1"))
  resids <- list()
  ys <- list(kwh = obs$kwh, sqrt = sqrt(obs$kwh), log = log(obs$kwh))
  for (label in names(ys)) {
    y <- ys[[label]]
    b <- qr.solve(X, y)
    fitv <- as.numeric(X %*% b)
    r <- y - fitv
    # np.var is the population variance; the ratio matches the sample-variance ratio.
    r2 <- 1 - var(r) / var(y)
    Z <- cbind(1, fitv)
    g <- qr.solve(Z, r^2)
    ss_tot <- sum((r^2 - mean(r^2))^2)
    ss_res <- sum((r^2 - as.numeric(Z %*% g))^2)
    lm_stat <- length(r) * (1 - ss_res / ss_tot)
    p_het <- pchisq(lm_stat, 1, lower.tail = FALSE)
    adj <- diff(t) == 1
    rho1 <- cor(r[-length(r)][adj], r[-1][adj])
    add(sprintf("%8s %7.3f %11.3f %11.2f %10.4f %7.3f", label, r2, skew_raw(r), kurt_raw(r), p_het, rho1))
    resids[[label]] <- list(r = r, fit = fitv)
  }
  add("")
  add("  Interpretation: the log scale should show the flattest variance-vs-mean")
  add("  relationship (largest BP p-value). Gaussian kurtosis is 3; values well")
  add("  above 3 motivate a t likelihood.")
  add("")

  add("HARMONIC ORDER, OLS on the log scale")
  add(sprintf("%3s %7s %7s %8s %9s %10s", "J", "params", "R2", "adj R2", "resid sd", "BIC"))
  ylog <- log(obs$kwh)
  for (J in 0:4) {
    XJ <- design(t, dow, J)$X
    b <- qr.solve(XJ, ylog)
    r <- ylog - as.numeric(XJ %*% b)
    n <- length(r)
    p <- ncol(XJ)
    r2 <- 1 - var(r) / var(ylog)
    adjr <- 1 - (1 - r2) * (n - 1) / (n - p)
    bic <- n * log(sum(r^2) / n) + p * log(n)
    add(sprintf("%3d %7d %7.3f %8.3f %9.3f %10.1f", J, p, r2, adjr, sd(r), bic))
  }
  add("")
  add("  Note: these OLS criteria ignore the strong residual autocorrelation and so")
  add("  overstate the evidence for extra harmonics. They are exploratory only; the")
  add("  formal comparison in step 3 uses WAIC and held-out forecast scores.")
  add("")

  r_log <- resids$log$r
  add("RESIDUAL AUTOCORRELATION, log scale, J = 2")
  a <- acf_np(r_log, 21)
  for (k in c(1, 2, 3, 4, 5, 6, 7, 14, 21)) add(sprintf("  lag %2d  %+0.3f", k, a[k]))
  add("")
  add("  Geometric decay consistent with AR(1); an AR(1) with rho1 as fitted predicts:")
  for (k in c(2, 3, 7)) add(sprintf("  lag %2d  %+0.3f (AR(1) prediction) vs %+0.3f (observed)", k, a[1]^k, a[k]))
  add("")
  add("LEFT-TAIL DAYS ON THE LOG SCALE (standardised OLS residual < -2.5)")
  z <- r_log / sd(r_log)
  for (i in which(z < -2.5)) {
    add(sprintf("  %s  kwh=%6.2f  z=%+0.2f", obs$date[i], obs$kwh[i], z[i]))
  }

  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "01b_transform_check.txt"))

  png(file.path(ROOT, "figures", "01b_scale_diagnostics.png"), width = 13.5, height = 7, units = "in", res = 140)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (label in c("kwh", "sqrt", "log")) {
    r <- resids[[label]]$r
    fitv <- resids[[label]]$fit
    plot(fitv, r, pch = 16, cex = 0.35, col = adjustcolor("#1f4e79", 0.4),
         xlab = "fitted", ylab = "residual", main = sprintf("%s: residual vs fitted", label))
    abline(h = 0)
  }
  for (label in c("kwh", "sqrt", "log")) {
    r <- resids[[label]]$r
    qqnorm(r / sd(r), pch = 16, cex = 0.35, col = adjustcolor("#1f4e79", 0.5),
           main = sprintf("%s: normal Q-Q", label))
    qqline(r / sd(r))
  }
  dev.off()

  png(file.path(ROOT, "figures", "01b_residual_acf.png"), width = 7, height = 3.4, units = "in", res = 140)
  par(mar = c(4, 4, 3, 1))
  lags <- 1:21
  bp <- barplot(a, names.arg = lags, col = "#1f4e79", xlab = "lag (days)", ylab = "residual autocorrelation",
                main = "Residual ACF after trend + 2 harmonics + day-of-week (log scale)")
  lines(bp, a[1]^lags, type = "b", pch = 16, cex = 0.6, col = "#b03a2e")
  ci <- 1.96 / sqrt(length(r_log))
  abline(h = c(ci, -ci), col = "gray", lty = 2)
  legend("topright", legend = expression(paste("AR(1) fit: ", rho[1]^k)), col = "#b03a2e", pch = 16, lty = 1, bty = "n")
  dev.off()
}

if (sys.nframe() == 0L) main()
