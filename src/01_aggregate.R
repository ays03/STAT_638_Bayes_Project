#!/usr/bin/env Rscript
# Step 1: aggregate the UCI one-minute household power record to a daily kWh series.
#
# Global_active_power (GAP) is average real power (kW) over each one-minute
# interval, so the energy used in that minute is GAP/60 kWh.
#
# Cleaning: the partial first and last days are dropped, and every remaining day
# is checked to contain exactly 1,440 unique minutes. A day is retained only if at
# least MIN_COVERAGE (95%) of its minutes have a reading; other days stay on the
# calendar as missing (kwh = NA) rather than being summed over their gaps.
#
# Transformation: on retained days, missing minutes are filled with that day's
# average observed power (implemented by scaling the observed total up to 1,440
# minutes), so the daily total is
#
#   y_d = (1/60) * sum_{m = 1}^{1440} GAP_m       [kWh]
#
# Outputs
#   data/daily.csv          date, t, kwh, dow, n_obs, coverage, observed, sub_kwh, year, doy
#   figures/01_*.png        daily series + coverage, distributions, season and week
#   results/01_data_summary.txt
#
#   Rscript src/01_aggregate.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))   # no ".." anymore
}
# script_root <- function() {
#   args <- commandArgs(trailingOnly = FALSE)
#   file_arg <- grep("^--file=", args, value = TRUE)
#   if (length(file_arg) == 0) return(normalizePath("."))
#   script <- sub("^--file=", "", file_arg[[1]])
#   normalizePath(file.path(dirname(normalizePath(script)), ".."))
# }

skew_adj <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  xc <- x - mean(x)
  m2 <- sum(xc^2) / n
  m3 <- sum(xc^3) / n
  (m3 / m2^1.5) * sqrt(n * (n - 1)) / (n - 2)
}

comma <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)

load_minutes <- function(path) {
  df <- read.csv(
    path, sep = ";", na.strings = "?", stringsAsFactors = FALSE,
    colClasses = c("character", "character", "numeric", "NULL", "NULL", "NULL",
                   "numeric", "numeric", "numeric")
  )
  names(df) <- c("Date", "Time", "Global_active_power", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
  df$date <- as.Date(df$Date, format = "%d/%m/%Y")
  df
}

#group minutes by date and aggregates to days
aggregate_daily <- function(df) {
  g <- as.integer(df$date)
  zero_na <- function(x) ifelse(is.na(x), 0, x)
  n_obs <- rowsum(as.numeric(!is.na(df$Global_active_power)), g, reorder = TRUE)
  kwh_raw <- rowsum(zero_na(df$Global_active_power), g, reorder = TRUE) / 60
  sub_kwh <- rowsum(
    zero_na(df$Sub_metering_1) + zero_na(df$Sub_metering_2) + zero_na(df$Sub_metering_3),
    g, reorder = TRUE
  ) / 1000
  n_rows <- rowsum(rep(1, nrow(df)), g, reorder = TRUE)
  data.frame(
    date = as.Date(as.integer(rownames(n_obs)), origin = "1970-01-01"),
    n_obs = as.numeric(n_obs), # number of minutes that have a power reading
    kwh_raw = as.numeric(kwh_raw), # sum of power/60 (will be used later)
    sub_kwh = as.numeric(sub_kwh),
    n_rows = as.numeric(n_rows), # how many minute rows the day has (to check that every day has 1440 rows later)
    stringsAsFactors = FALSE
  )
}

main <- function() {
  ROOT <- script_root()
  RAW <- file.path(ROOT, "data", "household_power_consumption.txt")
  MIN_COVERAGE <- 0.95 # 95% rule
  MINUTES_PER_DAY <- 1440
  for (sub in c("figures", "results")) dir.create(file.path(ROOT, sub), showWarnings = FALSE)

  cat("reading raw file ...\n")
  minutes <- load_minutes(RAW)
  n_rows <- nrow(minutes)
  n_missing <- sum(is.na(minutes$Global_active_power))
  daily <- aggregate_daily(minutes)
  # Time used to verify one row per minute (no repeated time stamps)
  stopifnot(!any(duplicated(minutes[, c("Date", "Time")])))
  rm(minutes) # everything is on a daily scale now

  first <- daily$date[1]
  last <- daily$date[nrow(daily)]
  daily <- daily[-c(1, nrow(daily)), ] # removed partial first and last days
  # every full day should have exactly 1,440 minutes (rows)
  stopifnot(all(daily$n_rows == 1440))
  #build the full calendar and apply the 95% rule
  full_dates <- seq(min(daily$date), max(daily$date), by = "day")
  full <- merge(data.frame(date = full_dates), daily, by = "date", all.x = TRUE, sort = TRUE)
  full$n_obs[is.na(full$n_obs)] <- 0 # ensure missing days don't have NA, but have 0 instead to indicate no recordings
  full$coverage <- full$n_obs / MINUTES_PER_DAY
  full$coverage[is.na(full$coverage)] <- 0
  full$observed <- full$coverage >= MIN_COVERAGE # mark the days with 95% coverage
  #for retained days with few missing values, re-scale the daily total 
  # (equivalent to filling each missing minute on a retained day with that day's 
  #  average observed power then summing all 1440 minutes to get daily energy)
  full$kwh <- ifelse(full$observed, full$kwh_raw * MINUTES_PER_DAY / full$n_obs, NA_real_) # our final response var (Y_d)
  full$t <- seq_len(nrow(full)) # day counter (1:1440) for use in harmonics
  full$dow <- as.integer(format(full$date, "%u")) - 1L # day of week (0=Monday,..., 6=Sunday)
  full$year <- as.integer(format(full$date, "%Y")) 
  full$doy <- as.integer(format(full$date, "%j"))

  obs <- full[full$observed, ]
  bad <- sum(obs$sub_kwh > obs$kwh_raw, na.rm = TRUE) # verify that sub-meters < total

  out <- full[, c("date", "t", "kwh", "dow", "n_obs", "coverage", "observed", "sub_kwh", "year", "doy")]
  op <- options(digits = 15, scipen = 999)
  on.exit(options(op), add = TRUE)
  write.csv(out, file.path(ROOT, "data", "daily.csv"), row.names = FALSE, na = "")

################################################################################
# Summary and Figures
################################################################################
  lines <- character()
  add <- function(s) lines <<- c(lines, s)
  add("RAW FILE")
  add(sprintf("  rows                         %s", comma(n_rows)))
  add(sprintf("  missing Global_active_power  %s (%.3f%%)", comma(n_missing), 100 * n_missing / n_rows))
  add(sprintf("  first / last timestamped day %s / %s  (both dropped as partial)", first, last))
  add("")
  add("DAILY SERIES")
  add(sprintf("  calendar days on grid        %d", nrow(full)))
  add(sprintf("  days with coverage >= %.2f   %d", MIN_COVERAGE, sum(full$observed)))
  add(sprintf("  days treated as missing      %d", sum(!full$observed)))
  add(sprintf("  of which zero observed mins  %d", sum(full$n_obs == 0)))
  add(sprintf("  sub-meter > total violations %d  (expect 0)", bad))
  add("")
  add("DAILY kWh (observed days only)")
  d <- obs$kwh
  add(sprintf("  n        %d", length(d)))
  add(sprintf("  mean     %.2f", mean(d)))
  add(sprintf("  sd       %.2f", sd(d)))
  add(sprintf("  min      %.2f  on %s", min(d), obs$date[which.min(d)]))
  add(sprintf("  q25      %.2f", as.numeric(quantile(d, 0.25))))
  add(sprintf("  median   %.2f", median(d)))
  add(sprintf("  q75      %.2f", as.numeric(quantile(d, 0.75))))
  add(sprintf("  max      %.2f  on %s", max(d), obs$date[which.max(d)]))
  add(sprintf("  skewness %.2f   (log scale: %.2f)", skew_adj(d), skew_adj(log(d))))
  add("")
  # 2010 is partial bc excludes December (one of the highest-consumption months)
  # 2006 omitted bc only had December
  add("ANNUALISED CONSUMPTION  (mean daily kWh x 365,2006 omitted; only 15 December days)")
  for (yr in setdiff(sort(unique(obs$year)),2006)) {
    sub <- obs[obs$year == yr, ]
    note <- if (yr == 2010) "  (partial year: Jan-Nov)" else ""
    add(sprintf("  %d  n=%4d  mean=%5.2f kWh/day  -> %s kWh/yr%s",
                yr, nrow(sub), mean(sub$kwh), comma(round(365 * mean(sub$kwh))), note))
  }
  add("")
  add("MEAN DAILY kWh BY DAY OF WEEK  (0=Mon)")
  for (k in 0:6) {
    sub <- obs[obs$dow == k, ]
    add(sprintf("  %d  n=%4d  mean=%5.2f", k, nrow(sub), mean(sub$kwh)))
  }
  add("")
  add("MEAN DAILY kWh BY MONTH")
  mon <- as.integer(format(obs$date, "%m"))
  for (m in 1:12) {
    sub <- obs[mon == m, ]
    add(sprintf("  %2d  n=%4d  mean=%5.2f", m, nrow(sub), mean(sub$kwh)))
  }
   add("")
  # add("LAG-1 AUTOCORRELATION OF kWh AFTER REMOVING MONTH + DOW MEANS")
  # z <- obs$kwh #observed
  # X <- cbind( #coeff estimates for expected (build design matrix of 0/1 indicator columns to indicate day and month)
  #   model.matrix(~ 0 + factor(mon)),
  #   model.matrix(~ factor(obs$dow))[, -1, drop = FALSE]
  # )
  # resid <- as.numeric(z - X %*% qr.solve(X, z)) #calculate residuals
  # adj <- diff(obs$t) == 1 # ensures lag-1 (bc some days are missing, not all retained days are consecutive)
  # rho1 <- cor(resid[-length(resid)][adj], resid[-1][adj]) # lag-1 autocorrelation
  # add(sprintf("  rho1 = %.3f  (consecutive observed days)", rho1))

  txt <- paste(lines, collapse = "\n")
  cat(txt, "\n", sep = "")
  writeLines(txt, file.path(ROOT, "results", "01_data_summary.txt"))

  png(file.path(ROOT, "figures", "01_series_and_coverage.png"), width = 11, height = 9, units = "in", res = 140)
  par(mfrow = c(3, 1), mar = c(3, 4, 3, 1))
  plot(full$date, full$kwh, type = "l", lwd = 0.5, col = "#1f4e79",
       main = "Daily household electricity consumption, Dec 2006 - Nov 2010", ylab = "kWh / day", xlab = "")
  plot(full$date, log(full$kwh), type = "l", lwd = 0.5, col = "#1f4e79",
       main = "log daily consumption (modelling scale)", ylab = "log kWh", xlab = "")
  plot(full$date, full$coverage, type = "h", lwd = 0.6, col = "#b03a2e",
       ylab = "fraction of 1440 min observed", xlab = "", main = "Daily data coverage", ylim = c(0, 1))
  abline(h = MIN_COVERAGE, lty = 2)
  legend("bottomright", legend = sprintf("retention threshold %.2f", MIN_COVERAGE), lty = 2, bty = "n", cex = 0.8)
  dev.off()

  png(file.path(ROOT, "figures", "01_distributions.png"), width = 13, height = 3.6, units = "in", res = 140)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  #plot days with some missing minutes and mark the 95% line
  low <- full$coverage[full$coverage < 1]
  hist(low, breaks = 30, col = "#1f4e79",
       main = sprintf("Coverage of the %d days with any missing minutes", length(low)),
       xlab = "daily coverage", ylab = "number of days", xlim = c(0, 1))
  abline(v = MIN_COVERAGE, lty = 2, col = "#b03a2e")
  # h <- hist(full$coverage, breaks = 60, plot = FALSE)
  # h$counts[h$counts == 0] <- NA
  # plot(h, col = "#1f4e79", main = "Coverage histogram (log count)",
  #      xlab = "daily coverage", ylab = "count", log = "y")
  hist(obs$kwh, breaks = 45, col = "#1f4e79",
       main = sprintf("Daily kWh (skewness %.2f)", skew_adj(obs$kwh)), xlab = "kWh / day")
  hist(log(obs$kwh), breaks = 45, col = "#2e7d32",
       main = sprintf("log daily kWh (skewness %.2f)", skew_adj(log(obs$kwh))), xlab = "log kWh / day")
  dev.off()

  png(file.path(ROOT, "figures", "01_exploratory_season_week.png"), width = 12, height = 4, units = "in", res = 140)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(NA, xlim = c(1, 366), ylim = range(obs$kwh), xlab = "day of year", ylab = "kWh / day",
       main = "Annual cycle: consumption vs day of year")
  cols <- c("#1f4e79", "#b03a2e", "#2e7d32", "#b9770e", "#6c3483")
  i <- 1
  for (yr in sort(unique(obs$year))) {
    sub <- obs[obs$year == yr, ]
    points(sub$doy, sub$kwh, pch = 16, cex = 0.35, col = adjustcolor(cols[(i - 1) %% length(cols) + 1], 0.55))
    i <- i + 1
  }
  legend("topright", legend = sort(unique(obs$year)), col = cols, pch = 16, bty = "n", cex = 0.75)
  boxplot(split(obs$kwh, factor(obs$dow, levels = 0:6)),
          names = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"),
          ylab = "kWh / day", main = "Weekly pattern (raw)", col = "#d6eaf8")
  dev.off()

  cat("\nwrote data/daily.csv, results/01_data_summary.txt, 3 figures\n")
}

# if (sys.nframe() == 0L) main()
main()
