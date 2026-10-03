#!/usr/bin/env Rscript
# Step 2: validate the sampler on simulated data with known truth.
#
# Simulates from the model itself (known beta, phi, sigma, nu), blanks out days in
# contiguous blocks to mimic the real missingness pattern, then checks that 95%
# credible intervals cover the true values and that the posterior means are close.
#
#   Rscript src/02_validate_sampler.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

ROOT <- script_root()
source(file.path(ROOT, "src", "model.R"))

simulate <- function(Tlen = 1440, J = 2, p = 1, seed = 0) {
  set.seed(seed)
  t <- seq_len(Tlen)
  dow <- (t + 5) %% 7
  des <- build_design(t, dow, J)
  X <- des$X
  beta_true <- c(3.20, -0.04, 0.05, 0.42, 0.02, -0.11, 0.06, 0.07, -0.05, 0.03, 0.18, 0.17)
  stopifnot(length(beta_true) == ncol(X))
  phi_true <- 0.45
  if (p != 1) phi_true <- rep(0.45, p)[seq_len(p)]
  phi_true <- as.numeric(phi_true)[seq_len(p)]
  sigma_true <- 0.24
  nu_true <- 5.0

  eps <- numeric(Tlen)
  if (p > 0) {
    V0 <- stationary_cov(phi_true, sigma_true^2)
    eps[seq_len(p)] <- as.numeric(t(chol(V0)) %*% rnorm(p))
  }
  lam <- rgamma(Tlen, shape = nu_true / 2, scale = 2 / nu_true)
  eta <- rnorm(Tlen) * sigma_true / sqrt(lam)
  if (p > 0) {
    for (i in (p + 1):Tlen) eps[i] <- sum(phi_true * eps[(i - 1):(i - p)]) + eta[i]
  } else {
    eps <- eta
  }
  z <- as.numeric(X %*% beta_true + eps)
  y <- exp(z)

  observed <- rep(TRUE, Tlen)
  for (start in sample.int(Tlen - 10, 6)) {
    observed[start:(start + sample.int(4, 1) - 1L)] <- FALSE
  }
  y <- ifelse(observed, y, NA_real_)
  list(
    y = y, observed = observed, t = t, dow = dow, names = des$names,
    beta = beta_true, phi = phi_true, sigma = sigma_true, nu = nu_true, J = J, p = p
  )
}

main <- function() {
  sim <- simulate()
  cat(sprintf("simulated T=%d, observed=%d, missing=%d\n",
              length(sim$y), sum(sim$observed), sum(!sim$observed)))
  cat("fitting 4 chains x 12000 iterations ...\n")
  fitres <- fit_chains(
    sim$y, sim$observed, sim$t, sim$dow,
    J = sim$J, p = sim$p, heavy = TRUE, scale = "log",
    n_iter = 12000, burn = 4000, thin = 4, n_chains = 4, seed = 11
  )

  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  n_cov <- 0L
  n_tot <- 0L
  add("SAMPLER VALIDATION ON SIMULATED DATA")
  add("True values are known by construction; 95% intervals should cover them.")
  add("")
  add(sprintf("%16s %8s %10s %8s %8s %6s %6s %7s",
              "parameter", "truth", "post mean", "2.5%", "97.5%", "cover", "Rhat", "ESS"))

  row <- function(label, truth, draws2d) {
    flat <- as.numeric(draws2d)
    qs <- quantile(flat, c(0.025, 0.975))
    dg <- diagnose(draws2d)
    ok <- qs[1] <= truth && truth <= qs[2]
    n_cov <<- n_cov + as.integer(ok)
    n_tot <<- n_tot + 1L
    add(sprintf("%16s %8.3f %10.3f %8.3f %8.3f %6s %6.3f %7.0f",
                label, truth, mean(flat), qs[1], qs[2], if (ok) "yes" else "NO", dg$rhat, dg$ess))
  }

  for (i in seq_along(sim$names)) row(sim$names[i], sim$beta[i], fitres$beta[, , i])
  row("sigma", sim$sigma, fitres$sigma)
  if (sim$p > 0) {
    for (i in seq_len(sim$p)) row(sprintf("phi_%d", i), sim$phi[i], fitres$phi[, , i])
  }
  row("nu", sim$nu, fitres$nu)

  add("")
  add(sprintf("coverage of 95%% intervals: %d/%d parameters", n_cov, n_tot))
  add("  (with ~15 parameters, 14 or 15 of 15 is the expected outcome)")
  rhats <- c(
    vapply(seq_len(dim(fitres$beta)[3]), function(i) diagnose(fitres$beta[, , i])$rhat, numeric(1)),
    diagnose(fitres$sigma)$rhat,
    diagnose(fitres$nu)$rhat
  )
  esss <- c(
    vapply(seq_len(dim(fitres$beta)[3]), function(i) diagnose(fitres$beta[, , i])$ess, numeric(1)),
    diagnose(fitres$sigma)$ess,
    diagnose(fitres$nu)$ess
  )
  add(sprintf("worst Rhat = %.4f (target < 1.01),  min bulk ESS = %.0f (target > 400)", max(rhats), min(esss)))
  add("")
  add("MISSING-DAY IMPUTATION")
  zi <- flatten(fitres, "z_imputed")
  sds <- apply(zi, 2, sd)
  add(sprintf("  %d imputed days; posterior sd of imputed log kWh ranges %.3f to %.3f",
              ncol(zi), min(sds), max(sds)))
  add("  (should be near the marginal error sd, larger for days inside long gaps)")

  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n", sep = "")
  dir.create(file.path(ROOT, "results"), showWarnings = FALSE)
  writeLines(txt, file.path(ROOT, "results", "02_sampler_validation.txt"))
}

if (sys.nframe() == 0L) main()
