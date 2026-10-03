#!/usr/bin/env Rscript
# Step 1c: how much memory does the error process have?
#
# The step-1b ACF decays far more slowly than AR(1) implies, so this script fits
# AR(p) error models by conditional least squares on gap-free runs of observed
# days and compares orders. It also checks whether the slow decay is an artefact
# of the multi-week absence blocks by refitting with those days removed.
#
#   Rscript src/01c_ar_order_check.R

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
  if (J >= 1) {
    for (j in seq_len(J)) {
      cols[[length(cols) + 1]] <- sin(2 * pi * j * t / 365)
      cols[[length(cols) + 1]] <- cos(2 * pi * j * t / 365)
    }
  }
  for (k in 1:6) cols[[length(cols) + 1]] <- as.numeric(dow == k)
  do.call(cbind, cols)
}

ar_fit <- function(e, t, p) {
  # Conditional LS for an AR(p). Only rows whose p lags are the immediately
  # preceding calendar days are used. No intercept: residuals are mean zero.
  tt <- as.integer(round(t))
  if (p == 0) return(list(phi = numeric(0), sigma = sd(e), n = length(e), resid = e))
  pos <- integer(max(tt))
  pos[tt] <- seq_along(tt)
  ok <- rep(TRUE, length(tt))
  cols <- vector("list", p)
  for (k in seq_len(p)) {
    prev <- tt - k
    valid <- prev >= 1L & prev <= length(pos)
    j <- integer(length(tt))
    j[valid] <- pos[prev[valid]]
    valid <- valid & j > 0L
    ok <- ok & valid
    cols[[k]] <- j
  }
  idx <- which(ok)
  y <- e[idx]
  X <- vapply(cols, function(j) e[j[idx]], numeric(length(idx)))
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  phi <- as.numeric(qr.solve(X, y))
  r <- as.numeric(y - X %*% phi)
  list(phi = phi, sigma = sd(r), n = length(y), resid = r)
}

report <- function(e, t, label, add) {
  add(sprintf("\n%s   (n = %d days)", label, length(e)))
  add(sprintf("%3s %6s %10s %10s %10s  coefficients", "p", "n_eff", "sigma_eta", "AIC", "BIC"))
  for (p in 0:5) {
    fit <- ar_fit(e, t, p)
    aic <- fit$n * log(fit$sigma^2) + 2 * (p + 1)
    bic <- fit$n * log(fit$sigma^2) + (p + 1) * log(fit$n)
    co <- paste(sprintf("%+.3f", fit$phi), collapse = "  ")
    add(sprintf("%3d %6d %10.4f %10.1f %10.1f  %s", p, fit$n, fit$sigma, aic, bic, co))
  }
  fit2 <- ar_fit(e, t, 2)
  if (length(fit2$phi) == 2) {
    p1 <- fit2$phi[1]
    p2 <- fit2$phi[2]
    stat <- (p1 + p2 < 1) && (p2 - p1 < 1) && (abs(p2) < 1)
    rho1 <- p1 / (1 - p2)
    marg <- fit2$sigma^2 * (1 - p2) / ((1 + p2) * ((1 - p2)^2 - p1^2))
    add(sprintf("    AR(2): stationary=%s, implied rho1=%.3f, marginal sd=%.3f",
                if (stat) "True" else "False", rho1, sqrt(marg)))
    add(sprintf("    ratio marginal sd / one-step sd = %.2f  (how much wider a long-horizon interval is than a one-day-ahead one)",
                sqrt(marg) / fit2$sigma))
  }
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
  y <- log(obs$kwh)
  X <- design(t, dow, J = 2)
  e <- as.numeric(y - X %*% qr.solve(X, y))

  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("AR ORDER SELECTION FOR THE ERROR PROCESS")
  add("Residuals from OLS on log kWh with intercept + trend + 2 harmonics + day-of-week.")
  report(e, t, "ALL OBSERVED DAYS", add)

  low <- obs$kwh < 12
  block <- rep(FALSE, nrow(obs))
  i <- 1L
  while (i <= nrow(obs)) {
    if (isTRUE(low[i])) {
      j <- i
      while (j + 1L <= nrow(obs) && isTRUE(low[j + 1L])) j <- j + 1L
      if (j - i + 1L >= 5L) block[i:j] <- TRUE
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  add("")
  add(sprintf("SUSTAINED LOW-USE BLOCKS (>=5 consecutive days under 12 kWh): %d days", sum(block)))
  i <- 1L
  while (i <= nrow(obs)) {
    if (isTRUE(block[i])) {
      j <- i
      while (j + 1L <= nrow(obs) && isTRUE(block[j + 1L])) j <- j + 1L
      add(sprintf("  %s to %s  (%2d days, mean %5.2f kWh)",
                  obs$date[i], obs$date[j], j - i + 1L, mean(obs$kwh[i:j])))
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }

  keep <- !block
  Xk <- X[keep, , drop = FALSE]
  yk <- y[keep]
  ek <- as.numeric(yk - Xk %*% qr.solve(Xk, yk))
  report(ek, t[keep], "EXCLUDING SUSTAINED LOW-USE BLOCKS", add)
  add("")
  add("  If the preferred order drops sharply once the absence blocks are removed, the")
  add("  long memory is largely a regime effect and a heavy-tailed AR(2) is the right")
  add("  compromise: AR(2) for genuine weather persistence, t errors for the blocks.")

  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "01c_ar_order.txt"))

  png(file.path(ROOT, "figures", "01c_pacf.png"), width = 11, height = 3.6, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  panels <- list(
    list(e = e, t = t, lab = "all days"),
    list(e = ek, t = t[keep], lab = "absence blocks removed")
  )
  for (pan in panels) {
    pac <- vapply(1:10, function(p) tail(ar_fit(pan$e, pan$t, p)$phi, 1), numeric(1))
    ci <- 1.96 / sqrt(length(pan$e))
    barplot(pac, names.arg = 1:10, col = "#1f4e79", ylim = range(c(pac, ci, -ci)),
            main = sprintf("Partial autocorrelation: %s", pan$lab), xlab = "lag (days)")
    abline(h = c(ci, -ci), col = "gray", lty = 2)
  }
  dev.off()
}

if (sys.nframe() == 0L) main()
