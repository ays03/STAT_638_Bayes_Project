#!/usr/bin/env Rscript
# Step 3: fit the model family to the full daily series and compare.
#
# Seven specifications share the same priors, the same data, and the same set of
# pointwise likelihood terms (days whose two preceding days are also observed).
# Every log-likelihood is a density on the kWh scale.
#
#   M1  J=1  iid Gaussian errors, log scale
#   M2  J=1  AR(1) + t errors, log scale
#   M3  J=2  AR(1) + t errors, log scale
#   M4  J=3  AR(1) + t errors, log scale
#   M5  J=2  AR(2) + t errors, log scale
#   M6  J=2  AR(1) + Gaussian errors, log scale
#   M7  J=2  AR(1) + t errors, kWh scale
#
#   Rscript src/03_fit_models.R

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
MCMC <- list(n_iter = 40000, burn = 10000, thin = 5, n_chains = 4)

param_table <- function(fitres) {
  out <- list()
  for (i in seq_along(fitres$names)) out[[length(out) + 1]] <- list(name = fitres$names[i], draws = fitres$beta[, , i])
  out[[length(out) + 1]] <- list(name = "sigma", draws = fitres$sigma)
  if (fitres$p > 0) {
    for (i in seq_len(fitres$p)) out[[length(out) + 1]] <- list(name = sprintf("phi_%d", i), draws = fitres$phi[, , i])
  }
  if (isTRUE(fitres$heavy)) out[[length(out) + 1]] <- list(name = "nu", draws = fitres$nu)
  out
}

kde_curve <- function(v) {
  v <- as.numeric(v)
  s <- sqrt(mean((v - mean(v))^2))
  bw <- 0.9 * s * length(v)^(-0.2)
  xs <- seq(min(v), max(v), length.out = 120)
  if (!is.finite(bw) || bw < 1e-12) return(list(x = xs, y = rep(0, length(xs))))
  kde <- rowSums(exp(-0.5 * (outer(xs, v, "-") / bw)^2))
  list(x = xs, y = kde / max(kde))
}

main <- function() {
  FITDIR <- file.path(ROOT, "results", "03_fits")
  dir.create(FITDIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  d <- read_daily(file.path(ROOT, "data", "daily.csv"))
  y <- d$kwh
  observed <- d$observed
  t <- d$t
  dow <- d$dow
  ll_index <- common_ll_index(observed, max_lag = 2)
  cat(sprintf("T = %d days, observed = %d, shared likelihood terms = %d\n",
              length(y), sum(observed), length(ll_index)))

  fits <- list()
  conv_lines <- c(
    "CONVERGENCE DIAGNOSTICS",
    sprintf("%d chains x %d iterations, %d discarded, thin %d -> %d draws",
            MCMC$n_chains, MCMC$n_iter, MCMC$burn, MCMC$thin,
            MCMC$n_chains * (MCMC$n_iter - MCMC$burn) %/% MCMC$thin),
    "Rhat is rank-normalised split-Rhat; ESS is bulk effective sample size."
  )
  waic_rows <- list()

  for (spec in SPECS) {
    t0 <- proc.time()
    cat(sprintf("\nfitting %s: %s ...\n", spec$id, spec$label))
    fitres <- fit_chains(
      y, observed, t, dow, ll_index = ll_index,
      seed = 1000 + 7 * as.integer(substring(spec$id, 2)),
      J = spec$J, p = spec$p, heavy = spec$heavy, scale = spec$scale,
      n_iter = MCMC$n_iter, burn = MCMC$burn, thin = MCMC$thin, n_chains = MCMC$n_chains
    )
    fits[[spec$id]] <- fitres
    el <- (proc.time() - t0)[["elapsed"]]
    pars <- param_table(fitres)
    rhats <- numeric(length(pars))
    esss <- numeric(length(pars))
    conv_lines <- c(conv_lines, "", sprintf("--- %s: %s   (%.1f s) ---", spec$id, spec$label, el),
                    sprintf("%16s %10s %9s %10s %10s %7s %8s",
                            "parameter", "mean", "sd", "2.5%", "97.5%", "Rhat", "ESS"))
    for (i in seq_along(pars)) {
      flat <- as.numeric(pars[[i]]$draws)
      dg <- diagnose(pars[[i]]$draws)
      rhats[i] <- dg$rhat
      esss[i] <- dg$ess
      qs <- quantile(flat, c(0.025, 0.975))
      conv_lines <- c(conv_lines, sprintf("%16s %10.4f %9.4f %10.4f %10.4f %7.4f %8.0f",
                                          pars[[i]]$name, mean(flat), sd(flat), qs[1], qs[2], dg$rhat, dg$ess))
    }
    conv_lines <- c(conv_lines, sprintf("  worst Rhat %.4f (target < 1.01), min ESS %.0f (target > 400)",
                                        max(rhats), min(esss)))
    cat(sprintf("   done in %.1fs  worst Rhat %.4f  min ESS %.0f\n", el, max(rhats), min(esss)))
    w <- waic(fitres$loglik)
    waic_rows[[spec$id]] <- list(
      id = spec$id, label = spec$label, elpd = w$elpd_waic, p_waic = w$p_waic,
      waic = w$waic, se = w$se, n = w$n, elpd_i = w$elpd_i,
      worst_rhat = max(rhats), min_ess = min(esss), secs = el
    )
    saveRDS(fitres, file.path(FITDIR, paste0(spec$id, ".rds")))
  }

  writeLines(conv_lines, file.path(ROOT, "results", "03_convergence.txt"))

  best <- waic_rows[[which.max(vapply(waic_rows, function(r) r$elpd, numeric(1)))]]
  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("WAIC MODEL COMPARISON")
  add(sprintf("All log-likelihoods are densities of daily kWh on %d shared days", best$n))
  add("(log-scale models carry the log-Jacobian), so all rows are comparable.")
  add("elpd_waic: higher is better. d_elpd is relative to the best model;")
  add("se_diff is the standard error of that paired difference.")
  add("")
  add(sprintf("%4s %28s %10s %8s %9s %8s %7s %7s",
              "id", "specification", "elpd_waic", "p_waic", "d_elpd", "se_diff", "Rhat", "ESS"))
  ord <- order(vapply(waic_rows, function(r) -r$elpd, numeric(1)))
  for (r in waic_rows[ord]) {
    dif <- r$elpd - best$elpd
    se_d <- if (r$id == best$id) 0 else {
      di <- r$elpd_i - best$elpd_i
      sqrt(length(di) * var(di))
    }
    add(sprintf("%4s %28s %10.1f %8.1f %9.1f %8.1f %7.4f %7.0f",
                r$id, r$label, r$elpd, r$p_waic, dif, se_d, r$worst_rhat, r$min_ess))
  }
  add("")
  add("TARGETED COMPARISONS (each isolates one modelling decision)")
  cmp <- function(a, b, what) {
    di <- waic_rows[[a]]$elpd_i - waic_rows[[b]]$elpd_i
    dlt <- sum(di)
    se <- sqrt(length(di) * var(di))
    verdict <- if (dlt > 2 * se) "supported" else if (dlt < -2 * se) "no support" else "inconclusive"
    add(sprintf("  %-44s %s - %s = %+8.1f +/- %5.1f   %s", what, a, b, dlt, se, verdict))
  }
  cmp("M2", "M1", "AR(1) + t errors vs iid Gaussian")
  cmp("M3", "M2", "second harmonic (J=2 vs J=1)")
  cmp("M4", "M3", "third harmonic (J=3 vs J=2)")
  cmp("M5", "M3", "AR(2) vs AR(1)")
  cmp("M3", "M6", "t errors vs Gaussian, given AR(1)")
  cmp("M3", "M7", "log scale vs kWh scale")
  add("")
  add("  'supported' means the elpd difference exceeds two standard errors.")
  add("  WAIC on a conditional (one-step-ahead) factorisation understates how much")
  add("  a term matters for long-horizon forecasts, so step 4 re-tests each of")
  add("  these decisions on genuinely held-out data.")
  txt <- paste(lines, collapse = "\n")
  cat("\n", txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "03_waic.txt"))

  for (mid in c("M3", "M5")) {
    pars <- param_table(fits[[mid]])
    n <- length(pars)
    png(file.path(ROOT, "figures", sprintf("03_trace_%s.png", mid)),
        width = 11, height = 1.5 * n, units = "in", res = 120)
    layout(matrix(seq_len(2 * n), ncol = 2, byrow = TRUE), widths = c(2.2, 1))
    par(mar = c(2, 4, 2, 1))
    chain_cols <- c("#1f4e79", "#b03a2e", "#1e8449", "#b9770e")
    for (i in seq_along(pars)) {
      dr <- pars[[i]]$draws
      matplot(t(dr), type = "l", lty = 1, lwd = 0.4, col = chain_cols, ylab = pars[[i]]$name, xlab = "")
      if (i == 1) title(sprintf("%s: trace (4 chains)", mid), cex.main = 0.9)
      for (ch in seq_len(nrow(dr))) {
        kd <- kde_curve(dr[ch, ])
        if (ch == 1) {
          plot(kd$x, kd$y, type = "l", lwd = 0.8, col = chain_cols[ch], yaxt = "n", xlab = "", ylab = "")
          if (i == 1) title("marginal posterior by chain", cex.main = 0.9)
        } else {
          lines(kd$x, kd$y, lwd = 0.8, col = chain_cols[ch])
        }
      }
    }
    dev.off()
  }
  cat(sprintf("\nwrote %d fits to results/03_fits/, convergence + WAIC tables, trace plots\n", length(SPECS)))
}

if (sys.nframe() == 0L) main()
