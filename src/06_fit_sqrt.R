#!/usr/bin/env Rscript
# Step 6: fit the final model on the square-root scale (Report sections 4 and 5).
#
#   main     J=2, AR(1) + t errors, sqrt(kWh), all data: the reported posterior
#   J1       J=1, otherwise identical: does the second harmonic matter?
#   ar2      AR(2), otherwise identical: is one day of memory enough?
#   kwh      J=2, AR(1) + t errors on the untransformed kWh scale
#   holdout  main model refitted on data through 2009-11-30, then forecast the
#            final year (single 12-month origin and rolling origins, h = 1-60)
#   holdout_ar2  the same with AR(2) errors, scored on the same days
#
#   Rscript src/06_fit_sqrt.R               # everything, about 12 minutes
#   Rscript src/06_fit_sqrt.R main holdout  # a subset
#
# Writes results/06_sqrt/<name>.rds. The pointwise log-likelihood is reduced to
# its WAIC summary before saving, which keeps each file small.

library(here)
source(here("src", "model.R"))

OUT <- here("results", "06_sqrt")
MCMC <- list(n_iter = 40000, burn = 10000, thin = 5, n_chains = 4)
TRAIN_END <- as.Date("2009-11-30")
ROLL_H <- 60L
ROLL_STEP <- 3L

FULL_FITS <- list(
  main = list(J = 2, p = 1, scale = "sqrt", seed = 6100, n_paths = 100),
  J1 = list(J = 1, p = 1, scale = "sqrt", seed = 6200, n_paths = 1),
  ar2 = list(J = 2, p = 2, scale = "sqrt", seed = 6300, n_paths = 1),
  kwh = list(J = 2, p = 1, scale = "identity", seed = 6400, n_paths = 1)
)

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

run_full <- function(name, d) {
  s <- FULL_FITS[[name]]
  f <- fit_chains(
    d$kwh, d$observed, d$t, d$dow,
    n_chains = MCMC$n_chains, seed = s$seed,
    J = s$J, p = s$p, heavy = TRUE, scale = s$scale,
    n_iter = MCMC$n_iter, burn = MCMC$burn, thin = MCMC$thin,
    # identical likelihood terms for every model, so WAIC differences are paired
    ll_index = common_ll_index(d$observed, max_lag = 2),
    n_paths = s$n_paths, disperse_init = TRUE
  )
  f$waic <- waic(f$loglik)
  f$loglik <- NULL
  f$t_center <- mean(d$t)
  f$mcmc <- MCMC
  f
}

run_holdout <- function(d, p = 1, seed = 6500) {
  train <- d[d$date <= TRAIN_END, ]
  test <- d[d$date > TRAIN_END, ]
  rownames(test) <- NULL
  t_center <- mean(train$t)
  f <- fit_chains(
    train$kwh, train$observed, train$t, train$dow,
    n_chains = MCMC$n_chains, seed = seed,
    J = 2, p = p, heavy = TRUE, scale = "sqrt",
    n_iter = MCMC$n_iter, burn = MCMC$burn, thin = MCMC$thin,
    n_paths = 1, disperse_init = TRUE
  )
  f$loglik <- NULL
  scalars <- c(list(f$sigma, f$nu), lapply(seq_len(p), function(i) f$phi[, , i]),
               lapply(seq_along(f$names), function(i) f$beta[, , i]))
  dg <- lapply(scalars, diagnose)
  conv <- c(
    rhat = max(vapply(dg, function(g) max(g$rhat, g$rhat_folded), numeric(1))),
    ess = min(vapply(dg, `[[`, numeric(1), "ess")),
    ess_tail = min(vapply(dg, `[[`, numeric(1), "ess_tail"))
  )

  # A. one 12-month forecast from the last training day
  fc <- forecast(f, test$t, test$dow, t_center, n_draws = 4000, seed = 99)
  q <- apply(fc$y, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  single <- data.frame(
    date = test$date, kwh = test$kwh, observed = test$observed,
    lo95 = q[1, ], lo50 = q[2, ], med = q[3, ], hi50 = q[4, ], hi95 = q[5, ],
    crps = crps(fc$y, test$kwh), lpd = log_pred_density(fc, test$kwh), pit = pit(fc, test$kwh)
  )

  # B. rolling origins: from every 3rd holdout day whose day and previous day
  #    are observed, forecast 1-60 days ahead. Given beta, the AR state at an
  #    observed origin is z_o - x_o'beta, so no refitting is needed. The origin
  #    rule does not depend on p, so AR(1) and AR(2) are scored on the same days.
  beta <- flatten(f, "beta")
  set.seed(64)
  sub <- sample.int(nrow(beta), 1500)
  Xtest <- build_design(test$t, test$dow, 2, t_center = t_center)$X
  z_test <- to_model_scale(test$kwh, "sqrt")
  eps_test <- matrix(z_test, nrow(beta), length(z_test), byrow = TRUE) - beta %*% t(Xtest)
  origins <- seq(2L, nrow(test) - 1L, by = ROLL_STEP)
  origins <- origins[test$observed[origins] & test$observed[origins - 1L]]
  rows <- vector("list", length(origins))
  for (oi in seq_along(origins)) {
    o <- origins[oi]
    idx <- (o + 1L):min(o + ROLL_H, nrow(test))
    fcr <- forecast(f, test$t[idx], test$dow[idx], t_center,
                    eps_init = eps_test[sub, o - (seq_len(p) - 1L), drop = FALSE],
                    draw_idx = sub, seed = 7L * oi)
    qq <- apply(fcr$y, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975))
    rows[[oi]] <- data.frame(
      origin = o, h = seq_along(idx), observed = test$observed[idx], kwh = test$kwh[idx],
      lo95 = qq[1, ], lo50 = qq[2, ], hi50 = qq[3, ], hi95 = qq[4, ],
      crps = crps(fcr$y, test$kwh[idx]), lpd = log_pred_density(fcr, test$kwh[idx])
    )
  }
  list(
    single = single, roll = do.call(rbind, rows), n_origins = length(origins),
    n_train = sum(train$observed), n_test = sum(test$observed),
    train_end = TRAIN_END, phi = matrix(f$phi, ncol = p), conv = conv,
    train_tail = train[train$date > TRAIN_END - 200, c("date", "kwh")]
  )
}

main <- function(which) {
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  d <- read_daily(here("data", "daily.csv"))
  for (name in which) {
    t0 <- proc.time()
    res <- switch(name,
      holdout = run_holdout(d, p = 1, seed = 6500),
      holdout_ar2 = run_holdout(d, p = 2, seed = 6600),
      run_full(name, d)
    )
    saveRDS(res, file.path(OUT, paste0(name, ".rds")))
    cat(sprintf("%s done in %.0f s\n", name, (proc.time() - t0)[["elapsed"]]))
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  main(if (length(args)) args else c(names(FULL_FITS), "holdout", "holdout_ar2"))
}
