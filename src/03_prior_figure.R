#!/usr/bin/env Rscript
# Prior densities for Report section 3.
# Writes R_figures/03_priors.png. Numbers match src/model.R for the sqrt(kWh) model:
#   intercept N(0, 20^2), other coefficients N(0, 5^2),
#   sigma ~ half-t(3, scale 2), phi uniform on the stationary region,
#   nu - 2 ~ Exponential(mean 10).

library(here)

nav <- "#1f4e79"; teal <- "#148f77"; red <- "#b03a2e"; gold <- "#b9770e"

TAU_INT <- 20   # intercept prior SD          (PRIOR_TAU$sqrt[1] in model.R)
TAU_OTH <- 5    # other coefficients prior SD (PRIOR_TAU$sqrt[2])
A_SIGMA <- 2    # half-t scale for sigma      (PRIOR_SIGMA_A$sqrt)

# reference values from the data, used only to show the priors are wide
d   <- read.csv(here("data", "daily.csv"))
d$date <- as.Date(d$date)
sq  <- sqrt(d$kwh)                                     # NA on missing days
mon <- format(d$date, "%m")
gap   <- mean(sq[mon == "12"], na.rm = TRUE) - mean(sq[mon == "08"], na.rm = TRUE)
sd_sq <- sd(sq, na.rm = TRUE)

dir.create(here("R_figures"), showWarnings = FALSE)
png(here("R_figures", "03_priors.png"), width = 11, height = 8.2, units = "in", res = 140)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.2))

# (a) regression coefficients -------------------------------------------------
x <- seq(-50, 50, length.out = 800)
d_int <- dnorm(x, 0, TAU_INT)
d_oth <- dnorm(x, 0, TAU_OTH)
plot(x, d_oth, type = "l", lwd = 2, col = teal, ylim = c(0, 1.1 * max(d_oth)),
     xlab = "coefficient (sqrt kWh)", ylab = "prior density",
     main = "(a) Regression coefficients")
lines(x, d_int, lwd = 2, col = nav)
arrows(-gap / 2, 0.006, gap / 2, 0.006, code = 3, length = 0.05, lwd = 2, col = red)
text(gap / 2 + 1, 0.012, sprintf("Aug-Dec gap, %.1f sqrt kWh", gap), col = red, cex = 0.75, adj = 0)
legend("topright", bty = "n", lwd = 2, col = c(nav, teal), cex = 0.75,
       legend = c(sprintf("intercept, SD %g", TAU_INT), sprintf("other coefficients, SD %g", TAU_OTH)))

# (b) half-t for sigma --------------------------------------------------------
xs <- seq(0.01, 10, length.out = 600)
# Folded t_3 with scale A: f(s) = 2 * dt(s/A, 3) / A, s > 0.
ds <- 2 * dt(xs / A_SIGMA, df = 3) / A_SIGMA
plot(xs, ds, type = "l", lwd = 2, col = nav,
     xlab = "innovation scale (sqrt kWh)", ylab = "prior density",
     main = "(b) Innovation scale")
abline(v = sd_sq, lty = 2, lwd = 1.5, col = red)
text(sd_sq + 0.2, 0.9 * max(ds), sprintf("SD of daily sqrt kWh = %.2f", sd_sq), col = red, cex = 0.75, adj = 0)
text(6.5, 0.5 * max(ds), sprintf("half-t, 3 df\nscale %g", A_SIGMA), col = nav, cex = 0.75)

# (c) uniform on the AR(1) stationary interval --------------------------------
plot(NA, xlim = c(-1.45, 1.45), ylim = c(0, 0.75), xaxt = "n",
     xlab = expression("autoregressive coefficient" ~ phi),
     ylab = "prior density", main = "(c) AR(1) coefficient")
rect(-1.45, 0, -1, 0.75, col = "#f6e6e4", border = NA)
rect(1, 0, 1.45, 0.75, col = "#f6e6e4", border = NA)
rect(-1, 0, 1, 0.5, col = adjustcolor(nav, 0.18), border = nav, lwd = 2)
text(0, 0.58, "stationary: shocks fade", col = nav, cex = 0.75)
text(-1.22, 0.62, "excluded", col = red, cex = 0.65, srt = 90)
text(1.22, 0.62, "excluded", col = red, cex = 0.65, srt = 90)
axis(1, at = c(-1, 0, 1))

# (d) degrees of freedom ------------------------------------------------------
nu <- seq(2.02, 48, length.out = 500)
# nu - 2 ~ Exp(mean 10), so the density on nu is (1/10) exp(-(nu-2)/10).
dn <- dexp(nu - 2, rate = 0.1)
plot(nu, dn, type = "l", lwd = 2, col = nav,
     xlab = expression("degrees of freedom" ~ nu), ylab = "prior density",
     main = "(d) Tail weight")
abline(v = 2, lty = 2, lwd = 1.5, col = red)
abline(v = 12, lty = 3, lwd = 1.5, col = gold)
text(3.1, 0.085, expression(nu == 2), col = red, cex = 0.75, adj = 0)
text(13.2, 0.07, "prior mean of nu = 12", col = gold, cex = 0.75, adj = 0)
text(32, 0.04, "nearly\nGaussian", col = nav, cex = 0.75)

dev.off()