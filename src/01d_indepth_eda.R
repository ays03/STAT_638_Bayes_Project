#!/usr/bin/env Rscript
# Step 1d: in-depth exploratory analysis of daily household demand.
#
# Descriptive, pre-model work. Dependence diagnostics are on the kWh scale.
#
#   Rscript src/01d_indepth_eda.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

ORIGIN <- as.Date("2006-12-17")
MONTHS <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
DOW_NAMES <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
ENDUSE <- list(
  list("remainder", "Unmetered", "#1f4e79"),
  list("sub3", "Water heater and air conditioner", "#e67e22"),
  list("sub2", "Laundry", "#1e8449"),
  list("sub1", "Kitchen", "#7d3c98")
)
SEASON_MONTHS <- list(winter = c(12, 1, 2), spring = c(3, 4, 5), summer = c(6, 7, 8), autumn = c(9, 10, 11))

read_daily <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d$observed <- tolower(as.character(d$observed)) %in% c("true", "t", "1")
  d
}

r2_of <- function(y, resid) 1 - sum(resid^2) / sum((y - mean(y))^2)

skew_adj <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  xc <- x - mean(x)
  m2 <- sum(xc^2) / n
  m3 <- sum(xc^3) / n
  (m3 / m2^1.5) * sqrt(n * (n - 1)) / (n - 2)
}

kurt_adj <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  xc <- x - mean(x)
  g2 <- mean(xc^4) / mean(xc^2)^2 - 3
  ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * g2 + 6)
}

nw_cov <- function(X, resid, L = 14) {
  xe <- X * rep(resid, times = ncol(X))
  meat <- crossprod(xe)
  n <- nrow(xe)
  for (lag in seq_len(L)) {
    w <- 1 - lag / (L + 1)
    g <- crossprod(xe[(lag + 1):n, , drop = FALSE], xe[1:(n - lag), , drop = FALSE])
    meat <- meat + w * (g + t(g))
  }
  inv <- solve(crossprod(X))
  inv %*% meat %*% inv
}

ols <- function(X, y) {
  beta <- qr.solve(X, y)
  list(beta = beta, resid = as.numeric(y - X %*% beta))
}

easter <- function(year) {
  a <- year %% 19
  b <- year %/% 100
  c <- year %% 100
  d <- b %/% 4
  e <- b %% 4
  f <- (b + 8) %/% 25
  g <- (b - f + 1) %/% 3
  h <- (19 * a + b - d - g + 15) %% 30
  i <- c %/% 4
  k <- c %% 4
  ell <- (32 + 2 * e + 2 * i - h - k) %% 7
  m <- (a + 11 * h + 22 * ell) %/% 451
  month <- (h + ell - 7 * m + 114) %/% 31
  day <- ((h + ell - 7 * m + 114) %% 31) + 1
  as.Date(sprintf("%04d-%02d-%02d", year, month, day))
}

french_holidays <- function(years) {
  out <- list()
  for (year in years) {
    eas <- easter(year)
    fixed <- list(
      list(as.Date(sprintf("%d-01-01", year)), "New Year"),
      list(eas + 1, "Easter Monday"),
      list(as.Date(sprintf("%d-05-01", year)), "Labour Day"),
      list(as.Date(sprintf("%d-05-08", year)), "Victory Day"),
      list(eas + 39, "Ascension"),
      list(eas + 50, "Whit Monday"),
      list(as.Date(sprintf("%d-07-14", year)), "Bastille Day"),
      list(as.Date(sprintf("%d-08-15", year)), "Assumption"),
      list(as.Date(sprintf("%d-11-01", year)), "All Saints"),
      list(as.Date(sprintf("%d-11-11", year)), "Armistice"),
      list(as.Date(sprintf("%d-12-25", year)), "Christmas")
    )
    seen <- as.Date(character())
    for (item in fixed) {
      if (item[[1]] %in% seen) next
      seen <- c(seen, item[[1]])
      out[[length(out) + 1]] <- item
    }
  }
  data.frame(date = as.Date(vapply(out, function(z) as.character(z[[1]]), character(1))),
             name = vapply(out, function(z) z[[2]], character(1)), stringsAsFactors = FALSE)
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

month_dummies <- function(month, drop = 1) {
  levels <- which(vapply(1:12, function(m) m != drop && any(month == m), logical(1)))
  if (length(levels) == 0) return(matrix(0, length(month), 0))
  out <- vapply(levels, function(m) as.numeric(month == m), numeric(length(month)))
  if (is.null(dim(out))) matrix(out, ncol = 1) else out
}

dow_dummies <- function(dow, drop = 0) {
  levels <- setdiff(0:6, drop)
  out <- vapply(levels, function(k) as.numeric(dow == k), numeric(length(dow)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  list(X = out, levels = levels)
}

year_dummies <- function(year) {
  lv <- sort(unique(year))
  if (length(lv) <= 1) return(matrix(0, length(year), 0))
  out <- vapply(lv[-1], function(v) as.numeric(year == v), numeric(length(year)))
  if (is.null(dim(out))) matrix(out, ncol = 1) else out
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

season_of <- function(month) {
  out <- rep(NA_character_, length(month))
  for (nm in names(SEASON_MONTHS)) out[month %in% SEASON_MONTHS[[nm]]] <- nm
  out
}

write_tabular <- function(path, colspec, header, rows) {
  lines <- c(sprintf("\\begin{tabular}{%s}", colspec), "\\toprule", header, "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, path)
}

load_minutes <- function(path) {
  df <- read.csv(path, sep = ";", na.strings = "?", stringsAsFactors = FALSE,
                 colClasses = c("character", "character", "numeric", "NULL", "NULL", "NULL",
                                "numeric", "numeric", "numeric"))
  names(df) <- c("Date", "Time", "Global_active_power", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
  df$date <- as.Date(df$Date, format = "%d/%m/%Y")
  df$hour <- as.integer(substr(df$Time, 1, 2))
  df
}

dmon <- function(x) sub("^0", "", format(as.Date(x), "%d %b"))

main <- function() {
  ROOT <- script_root()
  RAW <- file.path(ROOT, "data", "household_power_consumption.txt")
  dir.create(file.path(ROOT, "figures"), showWarnings = FALSE)
  dir.create(file.path(ROOT, "results", "tables"), showWarnings = FALSE, recursive = TRUE)

  d <- read_daily(file.path(ROOT, "data", "daily.csv"))
  obs <- d[d$observed, ]
  rownames(obs) <- NULL
  t_check <- as.numeric(obs$date - ORIGIN) + 1
  if (!isTRUE(all.equal(t_check, as.numeric(obs$t)))) stop("t index does not match days since 2006-12-17")
  stopifnot(easter(2007) == as.Date("2007-04-08"), easter(2008) == as.Date("2008-03-23"),
            easter(2009) == as.Date("2009-04-12"), easter(2010) == as.Date("2010-04-04"),
            clim_day(as.Date("2008-03-01")) == 60, is.na(clim_day(as.Date("2008-02-29"))))

  y <- as.numeric(obs$kwh)
  t <- as.numeric(obs$t)
  dow <- as.integer(obs$dow)
  month <- as.integer(format(obs$date, "%m"))
  year <- as.integer(obs$year)
  n <- length(y)
  block <- low_use_mask(y)
  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("IN-DEPTH EDA  (descriptive OLS and direct summaries; not posterior estimates)")
  add(sprintf("Retained days: %d.  Low-use block days: %d.", n, sum(block)))
  add("")
  miss <- d[!d$observed, ]
  add("DAYS NOT RETAINED")
  add(sprintf("  n = %d", nrow(miss)))
  for (i in seq_len(nrow(miss))) {
    add(sprintf("  %s  coverage=%.3f  n_obs=%d", miss$date[i], miss$coverage[i], as.integer(miss$n_obs[i])))
  }
  add("")

  add("SEQUENTIAL R2 ON kWh  (each row adds the named columns to everything above)")
  add(sprintf("%-42s %7s %7s %9s", "specification", "R2", "dR2", "resid sd"))
  cols <- list(rep(1, n), t / 365)
  labels_seq <- c("intercept only", "+ linear trend")
  mats <- list(matrix(1, n, 1), do.call(cbind, cols))
  for (J in 1:2) {
    H <- harmonic_cols(t, J)
    cols[[length(cols) + 1]] <- H[, ncol(H) - 1]
    cols[[length(cols) + 1]] <- H[, ncol(H)]
    labels_seq <- c(labels_seq, sprintf("+ harmonic %d", J))
    mats[[length(mats) + 1]] <- do.call(cbind, cols)
  }
  Dd <- dow_dummies(dow, 0)$X
  for (j in seq_len(ncol(Dd))) cols[[length(cols) + 1]] <- Dd[, j]
  labels_seq <- c(labels_seq, "+ day of week (Mon baseline)")
  mats[[length(mats) + 1]] <- do.call(cbind, cols)
  for (J in 3:4) {
    H <- harmonic_cols(t, J)
    cols[[length(cols) + 1]] <- H[, ncol(H) - 1]
    cols[[length(cols) + 1]] <- H[, ncol(H)]
    labels_seq <- c(labels_seq, sprintf("+ harmonic %d", J))
    mats[[length(mats) + 1]] <- do.call(cbind, cols)
  }
  seq_rows <- list()
  prev <- 0
  for (i in seq_along(labels_seq)) {
    fit <- ols(mats[[i]], y)
    r2 <- r2_of(y, fit$resid)
    add(sprintf("%-42s %7.3f %7.3f %9.2f", labels_seq[i], r2, r2 - prev, sd(fit$resid)))
    seq_rows[[length(seq_rows) + 1]] <- list(label = labels_seq[i], r2 = r2, dr2 = r2 - prev, sd = sd(fit$resid))
    prev <- r2
  }
  X_month <- cbind(1, t / 365, month_dummies(month))
  fit_m <- ols(X_month, y)
  r2_month <- r2_of(y, fit_m$resid)
  X_month_dow <- cbind(X_month, Dd)
  fit_md <- ols(X_month_dow, y)
  r2_month_dow <- r2_of(y, fit_md$resid)
  weekend <- as.numeric(dow >= 5)
  r2_month_we <- r2_of(y, ols(cbind(X_month, weekend), y)$resid)
  X_j2 <- mats[[match("+ harmonic 2", labels_seq)]]
  r2_j2_we <- r2_of(y, ols(cbind(X_j2, weekend), y)$resid)
  r2_j2_dow <- seq_rows[[match("+ day of week (Mon baseline)", labels_seq)]]$r2
  add("")
  add("SEASONAL CEILING AND WEEKEND RESTRICTION")
  add(sprintf("  trend + month dummies                         R2 = %.3f", r2_month))
  add(sprintf("  trend + month + weekend indicator             R2 = %.3f  (dR2 vs month = %.3f)", r2_month_we, r2_month_we - r2_month))
  add(sprintf("  trend + month + full day-of-week              R2 = %.3f  (dR2 vs weekend = %.3f)", r2_month_dow, r2_month_dow - r2_month_we))
  add(sprintf("  trend + 2 harmonics + weekend                 R2 = %.3f", r2_j2_we))
  add(sprintf("  trend + 2 harmonics + full day-of-week        R2 = %.3f  (dR2 vs weekend = %.3f)", r2_j2_dow, r2_j2_dow - r2_j2_we))
  add(sprintf("  month indicators minus 2 harmonics, both + trend + DOW: R2 gap = %.3f", r2_month_dow - r2_j2_dow))
  add("")
  add("HARMONIC SHAPE ON THE kWh SCALE  (OLS with trend and day of week)")
  ref_dates <- seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "day")
  curves <- list()
  for (J in 1:2) {
    fitJ <- ols(cbind(1, t / 365, harmonic_cols(t, J), Dd), y)
    harm <- fitJ$beta[3:(2 + 2 * J)]
    raw_curve <- harmonic_on_dates(harm, ref_dates, J)
    curve <- raw_curve - mean(raw_curve) + mean(y[year %in% 2007:2009])
    curves[[J]] <- curve
    peak <- which.max(curve)
    trough <- which.min(curve)
    add(sprintf("  J=%d  peak %s (%.2f kWh)  trough %s (%.2f kWh)  range %.2f  resid sd %.2f",
                J, dmon(ref_dates[peak]), curve[peak], dmon(ref_dates[trough]), curve[trough],
                curve[peak] - curve[trough], sd(fitJ$resid)))
    for (j in seq_len(J)) {
      a <- harm[2 * j - 1]; c_ <- harm[2 * j]
      add(sprintf("       harmonic %d amplitude %.2f kWh", j, sqrt(a^2 + c_^2)))
    }
  }
  add("")
  add("MONTHLY DISTRIBUTION OF DAILY kWh  (all retained days)")
  add(sprintf("%4s %5s %8s %7s %8s %8s %8s", "mo", "n", "mean", "sd", "p10", "p50", "p90"))
  month_stats <- vector("list", 12)
  for (m in 1:12) {
    v <- y[month == m]
    month_stats[[m]] <- c(m, length(v), mean(v), sd(v), quantile(v, 0.10), median(v), quantile(v, 0.90))
    row <- month_stats[[m]]
    add(sprintf("%4d %5d %8.2f %7.2f %8.2f %8.2f %8.2f", row[1], row[2], row[3], row[4], row[5], row[6], row[7]))
  }
  add("")
  add("MONTHLY MEAN BY YEAR")
  add(paste("  year", paste(sprintf("%6d", 1:12), collapse = " ")))
  for (yr in sort(unique(year))) {
    cells <- vapply(1:12, function(m) {
      v <- y[year == yr & month == m]
      if (length(v)) sprintf("%6.1f", mean(v)) else sprintf("%6s", "\u00b7")
    }, character(1))
    add(sprintf("  %d %s", yr, paste(cells, collapse = " ")))
  }
  dec_ex <- y[month == 12 & year != 2006]
  add("")
  add(sprintf("DECEMBER MEAN  all retained %.2f   excluding 2006 %.2f (n=%d)", mean(y[month == 12]), mean(dec_ex), length(dec_ex)))
  aug <- y[month == 8]
  aug_ex <- y[month == 8 & !block]
  add(sprintf("AUGUST MEAN    all retained %.2f (median %.2f)   excluding low-use blocks %.2f (median %.2f, n=%d)",
              mean(aug), median(aug), mean(aug_ex), median(aug_ex), length(aug_ex)))
  jul <- y[month == 7]
  add(sprintf("JULY MEAN      all retained %.2f (median %.2f)", mean(jul), median(jul)))
  add("")
  add("CORRELATION OF MONTHLY MEANS BETWEEN YEARS (months observed in both)")
  years <- 2007:2010
  mm <- lapply(years, function(yr) vapply(1:12, function(m) {
    v <- y[year == yr & month == m]
    if (length(v)) mean(v) else NA_real_
  }, numeric(1)))
  names(mm) <- years
  for (i in seq_len(length(years) - 1)) for (b in years[(i + 1):length(years)]) {
    okm <- is.finite(mm[[as.character(years[i])]]) & is.finite(mm[[as.character(b)]])
    add(sprintf("  %d vs %d: r = %.3f  (%d months)", years[i], b, cor(mm[[as.character(years[i])]][okm], mm[[as.character(b)]][okm]), sum(okm)))
  }
  add("")
  add("SEASONAL MEANS")
  for (nm in names(SEASON_MONTHS)) {
    v <- y[month %in% SEASON_MONTHS[[nm]]]
    add(sprintf("  %-8s n=%4d  mean=%6.2f  sd=%5.2f  skew=%5.2f", nm, length(v), mean(v), sd(v), skew_adj(v)))
  }
  add("")

  add("READING MINUTE FILE FOR SUB-METERS AND THE DIURNAL PROFILE ...")
  minutes <- load_minutes(RAW)
  miss_hour <- tapply(is.na(minutes$Global_active_power), minutes$hour, mean)
  add("FRACTION OF ROWS WITH MISSING POWER, BY HOUR")
  add(paste(" ", paste(sprintf("%02d:%.3f", 0:23, miss_hour[as.character(0:23)]), collapse = "  ")))
  gap_na <- is.na(minutes$Global_active_power)
  sub_na <- is.na(minutes$Sub_metering_1) | is.na(minutes$Sub_metering_2) | is.na(minutes$Sub_metering_3)
  add(sprintf("  minutes with GAP missing: %d", sum(gap_na)))
  add(sprintf("  minutes with any sub-meter missing: %d", sum(sub_na)))
  add(sprintf("  disagreement between those two masks: %d", sum(gap_na != sub_na)))
  okm <- minutes[!gap_na, ]
  g <- as.integer(okm$date)
  zn <- function(x) ifelse(is.na(x), 0, x)
  kwh_re <- rowsum(okm$Global_active_power / 60, g, reorder = TRUE)
  sub1 <- rowsum(zn(okm$Sub_metering_1) / 1000, g, reorder = TRUE)
  sub2 <- rowsum(zn(okm$Sub_metering_2) / 1000, g, reorder = TRUE)
  sub3 <- rowsum(zn(okm$Sub_metering_3) / 1000, g, reorder = TRUE)
  daily_sub <- data.frame(
    date = as.Date(as.integer(rownames(kwh_re)), origin = "1970-01-01"),
    kwh_re = as.numeric(kwh_re), sub1 = as.numeric(sub1), sub2 = as.numeric(sub2), sub3 = as.numeric(sub3)
  )
  ix <- match(obs$date, daily_sub$date)
  obs$kwh_re <- daily_sub$kwh_re[ix]
  obs$sub1 <- daily_sub$sub1[ix]
  obs$sub2 <- daily_sub$sub2[ix]
  obs$sub3 <- daily_sub$sub3[ix]
  add(sprintf("  max |recomputed daily kWh - daily.csv|: %.6f", max(abs(obs$kwh - obs$kwh_re), na.rm = TRUE)))
  obs$remainder <- obs$kwh - obs$sub1 - obs$sub2 - obs$sub3
  add(sprintf("  days with negative remainder: %d", sum(obs$remainder < -1e-6, na.rm = TRUE)))
  add(sprintf("  minimum remainder: %.3f kWh", min(obs$remainder, na.rm = TRUE)))
  add("")
  add("END-USE MEANS, ALL RETAINED DAYS")
  total_mean <- mean(obs$kwh)
  for (eu in ENDUSE) add(sprintf("  %-40s %6.2f kWh/day   %5.1f%%", eu[[2]], mean(obs[[eu[[1]]]]), 100 * mean(obs[[eu[[1]]]]) / total_mean))
  add("")
  add("END USE BY MONTH (mean kWh/day)")
  add(sprintf("%4s %8s %8s %8s %8s %8s", "mo", "total", "remain", "sub3", "sub2", "sub1"))
  for (m in 1:12) {
    sub <- obs[month == m, ]
    add(sprintf("%4d %8.2f %8.2f %8.2f %8.2f %8.2f", m, mean(sub$kwh), mean(sub$remainder), mean(sub$sub3), mean(sub$sub2), mean(sub$sub1)))
  }
  dec <- obs[month == 12, ]
  aug_df <- obs[month == 8, ]
  add("")
  add("DECEMBER MINUS AUGUST, BY END USE  (all retained days in those months)")
  d_tot <- mean(dec$kwh) - mean(aug_df$kwh)
  for (eu in ENDUSE) {
    delta <- mean(dec[[eu[[1]]]]) - mean(aug_df[[eu[[1]]]])
    add(sprintf("  %-40s %+6.2f kWh/day   %5.1f%% of the gap", eu[[2]], delta, 100 * delta / d_tot))
  }
  aug_in <- obs[month == 8 & !block, ]
  add("")
  add("AUGUST EXCLUDING LOW-USE BLOCKS, BY END USE")
  add(sprintf("  %-40s %6.2f", "Total", mean(aug_in$kwh)))
  for (eu in ENDUSE) add(sprintf("  %-40s %6.2f", eu[[2]], mean(aug_in[[eu[[1]]]])))
  add(sprintf("  December minus this August, total: %+.2f", mean(dec$kwh) - mean(aug_in$kwh)))
  vac <- obs[obs$date >= as.Date("2008-08-06") & obs$date <= as.Date("2008-08-30"), ]
  jul08 <- obs[year == 2008 & month == 7, ]
  aug07 <- obs[year == 2007 & month == 8, ]
  add("")
  add("VACANCY BLOCK 6-30 Aug 2008 VERSUS OCCUPIED COMPARISONS")
  add(sprintf("%-28s %4s %8s %8s %8s %8s %8s", "slice", "n", "total", "remain", "sub3", "sub2", "sub1"))
  for (sl in list(list("6-30 Aug 2008", vac), list("July 2008", jul08), list("August 2007", aug07))) {
    s <- sl[[2]]
    add(sprintf("%-28s %4d %8.2f %8.2f %8.2f %8.2f %8.2f", sl[[1]], nrow(s), mean(s$kwh), mean(s$remainder), mean(s$sub3), mean(s$sub2), mean(s$sub1)))
  }
  add(sprintf("  vacancy average power: %.3f kW", mean(vac$kwh) / 24))
  add(sprintf("  July 2008 average power: %.3f kW", mean(jul08$kwh) / 24))
  add(sprintf("  December average power: %.3f kW", mean(dec$kwh) / 24))

  add("")
  add("DAY OF WEEK AFTER MONTH EFFECTS  (deviations from the 7-day mean, Newey-West L=14)")
  Xw <- cbind(1, month_dummies(month), Dd)
  fitw <- ols(Xw, y)
  rw <- fitw$resid
  Vw <- nw_cov(Xw, rw, L = 14)
  b_dow <- tail(fitw$beta, 6)
  V_dow <- Vw[(nrow(Vw) - 5):nrow(Vw), (ncol(Vw) - 5):ncol(Vw), drop = FALSE]
  AA <- matrix(-1 / 7, 7, 6)
  for (k in 1:6) AA[k + 1, k] <- AA[k + 1, k] + 1
  dev <- as.numeric(AA %*% b_dow)
  se_dev <- sqrt(diag(AA %*% V_dow %*% t(AA)))
  for (k in 1:7) {
    add(sprintf("  %-4s %+6.2f   95%% NW [%+.2f, %+.2f]", DOW_NAMES[k], dev[k], dev[k] - 1.96 * se_dev[k], dev[k] + 1.96 * se_dev[k]))
  }
  add("")
  add("WEEKEND MINUS WEEKDAY, WITHIN SEASON, CONTROLLING FOR MONTH  (Newey-West L=14)")
  season <- season_of(month)
  we_rows <- list()
  for (name in c("winter", "spring", "summer", "autumn")) {
    sel <- season == name
    present <- which(vapply(1:12, function(m) any(month[sel] == m), logical(1)))
    X <- cbind(1, month_dummies(month[sel], drop = present[1]), weekend[sel])
    fit <- ols(X, y[sel])
    V <- nw_cov(X, fit$resid, L = 14)
    est <- fit$beta[length(fit$beta)]
    se <- sqrt(V[nrow(V), ncol(V)])
    we_rows[[length(we_rows) + 1]] <- list(name = name, n = sum(sel), est = est, se = se)
    add(sprintf("  %-8s n=%4d  weekend-weekday %+5.2f  NW se %.2f  95%% [%+.2f, %+.2f]",
                name, sum(sel), est, se, est - 1.96 * se, est + 1.96 * se))
  }
  Xwe <- cbind(1, month_dummies(month), weekend)
  fitwe <- ols(Xwe, y)
  Vwe <- nw_cov(Xwe, fitwe$resid, L = 14)
  sig2 <- sum(fitwe$resid^2) / (n - ncol(Xwe))
  ols_se <- sqrt(sig2 * solve(crossprod(Xwe))[ncol(Xwe), ncol(Xwe)])
  add(sprintf("  ALL      n=%4d  weekend-weekday %+5.2f  NW se %.2f  OLS se %.2f",
              n, fitwe$beta[length(fitwe$beta)], sqrt(Vwe[nrow(Vwe), ncol(Vwe)]), ols_se))
  for (L in c(7, 21)) {
    VL <- nw_cov(Xwe, fitwe$resid, L = L)
    add(sprintf("    sensitivity L=%d: NW se %.2f", L, sqrt(VL[nrow(VL), ncol(VL)])))
  }
  add("")
  add("END USE BY DAY OF WEEK, WINTER (DJF) AND SUMMER (JJA), RAW MEANS")
  for (item in list(list("winter", c(12, 1, 2)), list("summer", c(6, 7, 8)))) {
    add(paste(" ", item[[1]]))
    sl <- obs[month %in% item[[2]], ]
    for (k in 0:6) {
      s <- sl[sl$dow == k, ]
      add(sprintf("    %s n=%3d  total %6.2f  rem %6.2f  sub3 %6.2f  sub2 %5.2f  sub1 %5.2f",
                  DOW_NAMES[k + 1], nrow(s), mean(s$kwh), mean(s$remainder), mean(s$sub3), mean(s$sub2), mean(s$sub1)))
    }
  }

  add("")
  add("LINEAR TREND, kWh PER YEAR, WITH MONTH AND DAY-OF-WEEK CONTROLS")
  add("Newey-West L=14. The slope is the coefficient on t/365, matching the later model.")
  trend_row <- function(label, sel) {
    present <- which(vapply(1:12, function(m) any(month[sel] == m), logical(1)))
    Dsel <- dow_dummies(dow[sel], 0)$X
    X <- cbind(1, t[sel] / 365, month_dummies(month[sel], drop = present[1]), Dsel)
    fit <- ols(X, y[sel])
    V <- nw_cov(X, fit$resid, L = 14)
    est <- fit$beta[2]
    se <- sqrt(V[2, 2])
    add(sprintf("  %-40s n=%4d  %+6.3f  se %.3f  95%% [%+.3f, %+.3f]  %+.2f%%/year",
                label, sum(sel), est, se, est - 1.96 * se, est + 1.96 * se, 100 * est / mean(y[sel])))
    list(label = label, est = est, se = se, n = sum(sel), mu = mean(y[sel]))
  }
  trend_rows <- list(
    trend_row("All retained days", rep(TRUE, n)),
    trend_row("Drop 2006", year != 2006),
    trend_row("Drop low-use blocks", !block),
    trend_row("Complete years 2007-2009", year %in% 2007:2009)
  )
  stopifnot(all(vapply(trend_rows, length, integer(1)) == 5L))

  add("")
  add("PUBLIC HOLIDAYS IN THE RETAINED SAMPLE")
  hol <- french_holidays(2006:2010)
  obs$name <- hol$name[match(obs$date, hol$date)]
  obs$r_month_dow <- rw
  hol_days <- obs[!is.na(obs$name), ]
  add(sprintf("%-12s %-16s %-4s %8s %8s", "date", "name", "dow", "kwh", "resid"))
  for (i in seq_len(nrow(hol_days))) {
    r <- hol_days[i, ]
    add(sprintf("%-12s %-16s %-4s %8.2f %+8.2f", r$date, r$name, DOW_NAMES[as.integer(r$dow) + 1L], r$kwh, r$r_month_dow))
  }
  group_residual <- function(label, mask) {
    if (!any(mask)) {
      add(sprintf("  %s: no days", label))
      return()
    }
    add(sprintf("  %-28s n=%3d  mean kWh %6.2f  mean month+DOW residual %+6.2f",
                label, sum(mask), mean(y[mask]), mean(rw[mask])))
  }
  add("")
  add("HOLIDAY GROUPS  (residual is after month and day-of-week means)")
  names_h <- ifelse(is.na(obs$name), "", obs$name)
  group_residual("All public holidays", names_h != "")
  group_residual("Christmas Day", names_h == "Christmas")
  group_residual("New Year's Day", names_h == "New Year")
  group_residual("Bastille + Assumption", names_h %in% c("Bastille Day", "Assumption"))
  group_residual("May cluster + Easter/Whitsun", names_h %in% c("Easter Monday", "Labour Day", "Victory Day", "Ascension", "Whit Monday"))
  group_residual("November holidays", names_h %in% c("All Saints", "Armistice"))
  day <- as.integer(format(obs$date, "%d"))
  in_dec <- month == 12
  late_dec <- in_dec & day >= 24
  early_dec <- in_dec & !late_dec
  jan1 <- month == 1 & day == 1
  add("")
  add("YEAR-END STRETCH")
  add(sprintf("  1-23 Dec   n=%3d  mean %6.2f", sum(early_dec), mean(y[early_dec])))
  add(sprintf("  24-31 Dec  n=%3d  mean %6.2f", sum(late_dec), mean(y[late_dec])))
  add(sprintf("  1 Jan      n=%3d  mean %6.2f", sum(jan1), mean(y[jan1])))
  Xd <- cbind(1, year_dummies(year[in_dec]), as.numeric(day[in_dec] >= 24))
  fitd <- ols(Xd, y[in_dec])
  Vd <- nw_cov(Xd, fitd$resid, L = 7)
  add(sprintf("  late-minus-early December, year controls: %+.2f  NW se %.2f",
              fitd$beta[length(fitd$beta)], sqrt(Vd[nrow(Vd), ncol(Vd)])))
  add("")
  add("HIGHEST DAYS")
  for (i in order(obs$kwh, decreasing = TRUE)[1:8]) {
    r <- obs[i, ]
    add(sprintf("  %s %-3s %6.2f  remain %5.1f  sub3 %5.1f  sub2 %5.1f  sub1 %5.1f",
                r$date, DOW_NAMES[as.integer(r$dow) + 1L], r$kwh, r$remainder, r$sub3, r$sub2, r$sub1))
  }
  add("LOWEST DAYS")
  for (i in order(obs$kwh)[1:8]) {
    r <- obs[i, ]
    add(sprintf("  %s %-3s %6.2f  remain %5.1f  sub3 %5.1f  sub2 %5.1f  sub1 %5.1f",
                r$date, DOW_NAMES[as.integer(r$dow) + 1L], r$kwh, r$remainder, r$sub3, r$sub2, r$sub1))
  }

  add("")
  add("RESIDUAL DEPENDENCE ON THE kWh SCALE")
  add("Residuals from OLS: intercept + trend + 2 harmonics + day of week.")
  fitj <- ols(cbind(1, t / 365, harmonic_cols(t, 2), Dd), y)
  resid_j2 <- fitj$resid
  rho <- acf_calendar(resid_j2, t, 21)
  rho_ex <- acf_calendar(resid_j2[!block], t[!block], 21)
  add(sprintf("%4s %10s %10s %10s", "lag", "all days", "AR1 pred", "no blocks"))
  for (k in c(1, 2, 3, 4, 5, 6, 7, 14, 21)) add(sprintf("%4d %10.3f %10.3f %10.3f", k, rho[k], rho[1]^k, rho_ex[k]))
  pac <- pacf_yw(rho)
  pac_ex <- pacf_yw(rho_ex)
  add("PACF")
  for (k in 1:7) add(sprintf("  lag %d: all days %+.3f    no blocks %+.3f", k, pac[k], pac_ex[k]))
  ss <- resid_j2^2
  add("")
  add("SHARE OF J=2 RESIDUAL SUM OF SQUARES IN THE LOW-USE BLOCKS")
  add(sprintf("  days %.1f%%    residual SS %.1f%%", 100 * mean(block), 100 * sum(ss[block]) / sum(ss)))
  worst <- ss >= quantile(ss, 0.95)
  add(sprintf("  worst 5%% of days: %.1f%% of days, %.1f%% of residual SS", 100 * mean(worst), 100 * sum(ss[worst]) / sum(ss)))
  high <- worst & resid_j2 > 0
  low <- worst & resid_j2 < 0
  add(sprintf("    positive residuals among them: %d days, %.1f%% of that SS", sum(high), 100 * sum(ss[high]) / sum(ss[worst])))
  add(sprintf("    negative residuals among them: %d days, %.1f%% of that SS", sum(low), 100 * sum(ss[low]) / sum(ss[worst])))
  add(sprintf("    of the worst days, also inside a low-use block: %d", sum(worst & block)))
  top <- which.max(resid_j2)
  add(sprintf("    largest positive residual: %s  %+.1f kWh  (consumption %.1f)", obs$date[top], resid_j2[top], y[top]))
  add(sprintf("  residual skew %.2f   excess kurtosis %.2f", skew_adj(resid_j2), kurt_adj(resid_j2)))
  add(sprintf("  residual sd all days %.2f   excluding blocks %.2f", sd(resid_j2), sd(resid_j2[!block])))
  add("")
  add("RESIDUAL SD BY MONTH, AFTER TREND + 2 HARMONICS + DAY OF WEEK")
  add(sprintf("%4s %8s %12s %6s", "mo", "sd", "sd ex-block", "n ex"))
  sd_m <- numeric(12)
  sd_m_ex <- numeric(12)
  for (m in 1:12) {
    sel <- month == m
    sd_m[m] <- sd(resid_j2[sel])
    sd_m_ex[m] <- if (any(sel & !block)) sd(resid_j2[sel & !block]) else NA_real_
    add(sprintf("%4d %8.2f %12.2f %6d", m, sd_m[m], sd_m_ex[m], sum(sel & !block)))
  }
  win <- month %in% c(12, 1, 2)
  smr <- month %in% c(6, 7, 8)
  add(sprintf("  winter residual sd %.2f   summer residual sd %.2f", sd(resid_j2[win]), sd(resid_j2[smr])))
  add(sprintf("  winter ex-block %.2f   summer ex-block %.2f", sd(resid_j2[win & !block]), sd(resid_j2[smr & !block])))
  phi <- rho[1]
  add("")
  add(sprintf("DESCRIPTIVE AR(1) ON kWh RESIDUALS: rho1=%.3f", phi))
  add(sprintf("  marginal sd / one-step sd = %.3f", 1 / sqrt(1 - phi^2)))
  add(sprintf("  excluding blocks: rho1=%.3f   ratio %.3f", rho_ex[1], 1 / sqrt(1 - rho_ex[1]^2)))

  add("")
  add("DIURNAL PROFILE")
  retained <- obs$date
  ok_ret <- okm[okm$date %in% retained, ]
  id <- as.integer(ok_ret$date) * 24L + ok_ret$hour
  cnt <- rowsum(rep(1, nrow(ok_ret)), id, reorder = TRUE)
  n_before <- nrow(cnt)
  keep_id <- as.integer(rownames(cnt)[as.numeric(cnt) >= 30])
  sum_kw <- rowsum(ok_ret$Global_active_power, id, reorder = TRUE)
  sum_s1 <- rowsum(ok_ret$Sub_metering_1 * 0.06, id, reorder = TRUE)
  sum_s2 <- rowsum(ok_ret$Sub_metering_2 * 0.06, id, reorder = TRUE)
  sum_s3 <- rowsum(ok_ret$Sub_metering_3 * 0.06, id, reorder = TRUE)
  ids <- as.integer(rownames(cnt))
  dh <- data.frame(
    id = ids,
    date = as.Date(ids %/% 24L, origin = "1970-01-01"),
    hour = ids %% 24L,
    n = as.numeric(cnt),
    kw = as.numeric(sum_kw) / as.numeric(cnt),
    s1 = as.numeric(sum_s1) / as.numeric(cnt),
    s2 = as.numeric(sum_s2) / as.numeric(cnt),
    s3 = as.numeric(sum_s3) / as.numeric(cnt)
  )
  dh <- dh[dh$n >= 30, ]
  add(sprintf("  day-hours before / after requiring >= 30 observed minutes: %d / %d", n_before, nrow(dh)))
  dh$remainder_kw <- dh$kw - dh$s1 - dh$s2 - dh$s3
  info <- data.frame(date = obs$date, dow = obs$dow, month = month, season = season_of(month), weekend = obs$dow >= 5)
  dh <- merge(dh, info, by = "date", all.x = TRUE, sort = FALSE)
  profile <- function(mask) {
    sub <- dh[mask, , drop = FALSE]
    out <- matrix(NA_real_, 24, 5)
    colnames(out) <- c("kw", "remainder_kw", "s3", "s2", "s1")
    if (!nrow(sub)) return(out)
    for (h in 0:23) {
      sl <- sub[sub$hour == h, , drop = FALSE]
      if (!nrow(sl)) next
      out[h + 1, ] <- c(mean(sl$kw), mean(sl$remainder_kw), mean(sl$s3), mean(sl$s2), mean(sl$s1))
    }
    out
  }
  profiles <- list()
  for (season_name in c("winter", "summer")) {
    for (we in c(FALSE, TRUE)) {
      key <- sprintf("%s %s", season_name, if (we) "weekend" else "weekday")
      profiles[[key]] <- profile(dh$season == season_name & dh$weekend == we)
      pwr <- profiles[[key]][, "kw"]
      add(sprintf("  %-18s night 01-05 %.3f kW   evening 18-21 %.3f kW   peak hour %02d:00 (%.3f kW)   integral %.2f kWh",
                  key, mean(pwr[2:6]), mean(pwr[19:22]), which.max(pwr) - 1L, max(pwr), sum(pwr)))
    }
  }
  add("COMPOSITION AT THE MORNING AND EVENING PEAKS (kW)")
  for (key in c("winter weekday", "winter weekend", "summer weekday", "summer weekend")) {
    for (hour in c(7, 20)) {
      row <- profiles[[key]][hour + 1, ]
      add(sprintf("  %-18s %02d:00  total %.3f  rem %.3f  s3 %.3f  s2 %.3f  s1 %.3f",
                  key, hour, row["kw"], row["remainder_kw"], row["s3"], row["s2"], row["s1"]))
    }
  }
  win_wd <- obs[month %in% c(12, 1, 2) & obs$dow < 5, ]
  add(sprintf("  winter weekday mean of daily total: %.2f kWh (compare with integral above)", mean(win_wd$kwh)))

  tab <- file.path(ROOT, "results", "tables")
  pretty <- c(
    "intercept only" = "Intercept only",
    "+ linear trend" = "Linear trend",
    "+ harmonic 1" = "Add first harmonic",
    "+ harmonic 2" = "Add second harmonic",
    "+ harmonic 3" = "Add third harmonic",
    "+ harmonic 4" = "Add fourth harmonic",
    "+ day of week (Mon baseline)" = "Add day of week"
  )
  var_rows <- vapply(seq_rows, function(r) sprintf("%s & %.3f & %.3f & %.2f \\\\", pretty[[r$label]], r$r2, r$dr2, r$sd), character(1))
  var_rows <- c(var_rows, "\\addlinespace",
                sprintf("Trend + month + day of week & %.3f & --- & %.2f \\\\", r2_month_dow, sd(fit_md$resid)))
  write_tabular(file.path(tab, "eda_variance.tex"), "@{}lrrr@{}", "Specification & $R^2$ & Increment & Residual SD \\\\", var_rows)
  mon_rows <- vapply(month_stats, function(row) sprintf("%s & %d & %.2f & %.2f & %.2f & %.2f & %.2f \\\\",
                                                       MONTHS[row[1]], row[2], row[3], row[4], row[5], row[6], row[7]), character(1))
  write_tabular(file.path(tab, "eda_monthly.tex"), "@{}lrrrrrr@{}", "Month & $n$ & Mean & SD & 10\\% & Median & 90\\% \\\\", mon_rows)
  tr_rows <- vapply(trend_rows, function(r) sprintf("%s & %d & %+.3f & [%+.3f,\\ %+.3f] \\\\",
                                                    r$label, r$n, r$est, r$est - 1.96 * r$se, r$est + 1.96 * r$se), character(1))
  write_tabular(file.path(tab, "eda_trend.tex"), "@{}lrrr@{}", "Sample & $n$ & Slope (kWh/year) & 95\\% interval \\\\", tr_rows)
  dow_rows <- vapply(1:7, function(k) sprintf("%s & %+.2f & [%+.2f,\\ %+.2f] \\\\",
                                              DOW_NAMES[k], dev[k], dev[k] - 1.96 * se_dev[k], dev[k] + 1.96 * se_dev[k]), character(1))
  dow_rows <- c(dow_rows, "\\addlinespace", vapply(we_rows, function(r) {
    sprintf("%s weekend $-$ weekday & %+.2f & [%+.2f,\\ %+.2f] \\\\",
            paste0(toupper(substring(r$name, 1, 1)), substring(r$name, 2)), r$est, r$est - 1.96 * r$se, r$est + 1.96 * r$se)
  }, character(1)))
  write_tabular(file.path(tab, "eda_week.tex"), "@{}lrr@{}", "Contrast & kWh/day & 95\\% interval \\\\", dow_rows)
  end_rows <- vapply(ENDUSE, function(eu) {
    m <- mean(obs[[eu[[1]]]])
    delta <- mean(dec[[eu[[1]]]]) - mean(aug_df[[eu[[1]]]])
    sprintf("%s & %.2f & %.1f\\%% & %+.2f & %.1f\\%% \\\\", eu[[2]], m, 100 * m / total_mean, delta, 100 * delta / d_tot)
  }, character(1))
  end_rows <- c(end_rows, sprintf("Total & %.2f & 100\\%% & %+.2f & 100\\%% \\\\", total_mean, d_tot))
  write_tabular(file.path(tab, "eda_enduse.tex"), "@{}lrrrr@{}", "End use & Mean & Share & Dec$-$Aug & Share of gap \\\\", end_rows)

  year_colors <- c("2007" = "#6c3483", "2008" = "#1f4e79", "2009" = "#148f77", "2010" = "#b9770e")
  clim <- clim_day(obs$date)
  png(file.path(ROOT, "figures", "01d_seasonal.png"), width = 12.2, height = 4.4, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(NA, xlim = c(1, 365), ylim = range(c(y, unlist(curves)), na.rm = TRUE), xlab = "", ylab = "kWh / day", main = "Annual shape by year")
  for (yr in names(year_colors)) {
    sel <- year == as.integer(yr)
    sm <- smooth_on_clim(clim[sel], y[sel], window = 15, circular = yr != "2010")
    lines(1:365, sm, col = year_colors[[yr]], lwd = 1.15)
  }
  lines(1:365, curves[[1]], col = "#2c3e50", lwd = 1.15, lty = 2)
  lines(1:365, curves[[2]], col = "#b03a2e", lwd = 2)
  mids <- integer(12)
  for (m in 1:12) {
    days <- ref_dates[as.integer(format(ref_dates, "%m")) == m]
    mids[m] <- as.integer(format(days[length(days) %/% 2 + 1L], "%j"))
  }
  points(mids, vapply(month_stats, function(r) r[3], numeric(1)), pch = 16, cex = 0.6)
  axis(1, at = mids, labels = substr(MONTHS, 1, 1))
  legend("topright", legend = c(names(year_colors), "one harmonic", "two harmonics", "month mean"),
         col = c(year_colors, "#2c3e50", "#b03a2e", "black"), lty = c(1, 1, 1, 1, 2, 1, NA),
         pch = c(NA, NA, NA, NA, NA, NA, 16), lwd = c(1, 1, 1, 1, 1, 2, NA), cex = 0.65, ncol = 2, bty = "n")
  plot(NA, xlim = c(1, 12), ylim = range(vapply(month_stats, function(r) c(r[5], r[7]), numeric(2))),
       xlab = "", ylab = "kWh / day", main = "Level and spread by month", xaxt = "n")
  axis(1, at = 1:12, labels = substr(MONTHS, 1, 1))
  for (m in 1:12) {
    row <- month_stats[[m]]
    segments(m, row[5], m, row[7], col = adjustcolor("#1f4e79", 0.35), lwd = 4, lend = 1)
    points(m, row[6], pch = 16, col = "#1f4e79")
    points(m, row[3], pch = 4, col = "#b03a2e")
  }
  legend("topright", legend = c("10th to 90th percentile", "median", "mean"),
         col = c("#1f4e79", "#1f4e79", "#b03a2e"), pch = c(NA, 16, 4), lty = c(1, NA, NA), lwd = c(4, NA, NA), bty = "n", cex = 0.75)
  dev.off()

  png(file.path(ROOT, "figures", "01d_enduse.png"), width = 8.4, height = 4.3, units = "in", res = 140)
  par(mar = c(4, 4, 3, 1))
  M <- sapply(1:12, function(m) vapply(ENDUSE, function(eu) mean(obs[month == m, eu[[1]]]), numeric(1)))
  bp <- barplot(M, col = vapply(ENDUSE, function(eu) eu[[3]], character(1)), ylim = c(0, max(colSums(M)) * 1.22),
                names.arg = substr(MONTHS, 1, 1), ylab = "kWh / day", main = "Monthly mean consumption by end use")
  legend("top", legend = vapply(ENDUSE, function(eu) eu[[2]], character(1)),
         fill = vapply(ENDUSE, function(eu) eu[[3]], character(1)), ncol = 2, bty = "n", cex = 0.75)
  dev.off()

  png(file.path(ROOT, "figures", "01d_diurnal.png"), width = 12.4, height = 4.15, units = "in", res = 140)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  styles <- list("winter weekday" = "#1f4e79", "winter weekend" = "#5dade2", "summer weekday" = "#b03a2e", "summer weekend" = "#e59866")
  ylim <- range(vapply(profiles, function(p) p[, "kw"], numeric(24)), na.rm = TRUE)
  plot(0:23, profiles[["winter weekday"]][, "kw"], type = "n", ylim = ylim, xlab = "hour", ylab = "kW", main = "Mean power")
  for (key in names(styles)) lines(0:23, profiles[[key]][, "kw"], col = styles[[key]], lwd = 1.6)
  legend("topleft", legend = names(styles), col = unlist(styles), lwd = 1.6, cex = 0.65, bty = "n")
  stack_keys <- list(c("remainder_kw", "Unmetered", "#1f4e79"), c("s3", "Water heater / AC", "#e67e22"),
                     c("s2", "Laundry", "#1e8449"), c("s1", "Kitchen", "#7d3c98"))
  for (item in list(list("winter weekday", "Winter weekday"), list("summer weekday", "Summer weekday"))) {
    mat <- profiles[[item[[1]]]]
    bottom <- rep(0, 24)
    plot(0:23, mat[, "kw"], type = "n", ylim = ylim, xlab = "hour", ylab = "", main = item[[2]])
    for (sk in stack_keys) {
      vals <- mat[, sk[1]]
      polygon(c(0:23, 23:0), c(bottom + vals, rev(bottom)), col = sk[3], border = NA)
      bottom <- bottom + vals
    }
  }
  legend("topleft", legend = vapply(stack_keys, function(sk) sk[2], character(1)),
         fill = vapply(stack_keys, function(sk) sk[3], character(1)), cex = 0.65, bty = "n")
  dev.off()

  png(file.path(ROOT, "figures", "01d_residual.png"), width = 12, height = 4, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  lags <- 1:21
  bp <- barplot(rho, names.arg = lags, col = "#1f4e79", xlab = "lag (days)", ylab = "residual autocorrelation",
                main = "kWh residuals after trend, two harmonics, day of week",
                ylim = range(c(rho, rho_ex, phi^lags, 1.96 / sqrt(n), -1.96 / sqrt(n))))
  lines(bp, rho_ex, type = "b", pch = 16, cex = 0.5, col = "#b03a2e")
  lines(bp, phi^lags, lty = 2, col = "#7f8c8d")
  abline(h = c(1.96, -1.96) / sqrt(n), col = "gray", lty = 3)
  legend("topright", legend = c("all retained days", "low-use blocks removed", expression(paste("AR(1) from ", hat(rho)[1]))),
         col = c("#1f4e79", "#b03a2e", "#7f8c8d"), lty = c(NA, 1, 2), pch = c(15, 16, NA), cex = 0.7, bty = "n")
  bp <- barplot(sd_m, names.arg = substr(MONTHS, 1, 1), col = "#1f4e79", ylab = "residual SD (kWh)",
                main = "Scale of the residual by month", ylim = range(c(sd_m, sd_m_ex), na.rm = TRUE))
  points(bp, sd_m_ex, pch = 16, col = "#b03a2e")
  legend("topright", legend = c("all retained days", "low-use blocks removed"),
         col = c("#1f4e79", "#b03a2e"), pch = c(15, 16), cex = 0.75, bty = "n")
  dev.off()

  txt <- paste(lines, collapse = "\n")
  out <- file.path(ROOT, "results", "01d_eda.txt")
  writeLines(txt, out)
  cat(txt, "\n", sep = "")
  cat(sprintf("\nwrote %s, four figures, five tables\n", out))
}

if (sys.nframe() == 0L) main()
