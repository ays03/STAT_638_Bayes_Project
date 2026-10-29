#!/usr/bin/env Rscript
# Prior densities for Report section 3.
# Writes R_figures/03_priors.png. Numbers match src/model.R:
#   kWh: intercept N(0, 100^2), other coefficients N(0, 40^2),
#        sigma ~ half-t(3, scale 20), phi uniform on (-1, 1),
#        nu - 2 ~ Exponential(mean 10).

# Project root is the directory that contains the empty .here marker.
# Works when sourced from Report.Rmd (working directory is that folder)
# and when run from the workspace parent.
root <- if (file.exists(".here")) {
  normalizePath(".")
} else if (file.exists(file.path(".github-sync", ".here"))) {
  normalizePath(".github-sync")
} else {
  stop("Cannot find .here. Run this from the project directory.")
}

nav <- "#1f4e79"
teal <- "#148f77"
red <- "#b03a2e"
gold <- "#b9770e"

png(file.path(root, "R_figures", "03_priors.png"), width = 11, height = 8.2, units = "in", res = 140)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.2))

# (a) regression coefficients -------------------------------------------------
x <- seq(-200, 200, length.out = 800)
d_int <- dnorm(x, 0, 100)
d_oth <- dnorm(x, 0, 40)
plot(x, d_oth, type = "l", lwd = 2, col = teal, ylim = c(0, 0.013),
     xlab = "coefficient (kWh)", ylab = "prior density",
     main = "(a) Regression coefficients")
lines(x, d_int, lwd = 2, col = nav)
# August mean 13.67, December mean 35.64: a 22 kWh contrast.
arrows(-11, 0.0016, 11, 0.0016, code = 3, length = 0.06, lwd = 2, col = red)
text(48, 0.0024, "Aug-Dec gap, 22 kWh", col = red, cex = 0.75, adj = 0)
legend("topright", bty = "n", lwd = 2, col = c(nav, teal), cex = 0.75,
       legend = c("intercept, SD 100", "other coefficients, SD 40"))

# (b) half-t for sigma --------------------------------------------------------
xs <- seq(0.05, 90, length.out = 600)
# Folded t_3 with scale 20: f(s) = 2 * dt(s/20, 3) / 20, s > 0.
ds <- 2 * dt(xs / 20, df = 3) / 20
plot(xs, ds, type = "l", lwd = 2, col = nav,
     xlab = "innovation standard deviation (kWh)", ylab = "prior density",
     main = "(b) Innovation scale")
abline(v = 10, lty = 2, lwd = 1.5, col = red)
text(12, 0.032, "SD of daily kWh = 10", col = red, cex = 0.75, adj = 0)
text(55, 0.020, "half-t, 3 df\nscale 20", col = nav, cex = 0.75)

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
