"""
Step 1b: justify the modelling scale and the error structure using OLS pre-fits.

Compares raw-kWh, sqrt and log responses on the same mean structure
(intercept + linear trend + J harmonics + day-of-week), and reports residual
skewness, kurtosis, heteroscedasticity and the residual autocorrelation
function. This is the evidence used to justify the Bayesian model in step 2.
"""
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def design(t, dow, J):
    cols = [np.ones_like(t, dtype=float), t / 365.0]
    names = ["intercept", "trend_per_year"]
    for j in range(1, J + 1):
        cols += [np.sin(2 * np.pi * j * t / 365.0), np.cos(2 * np.pi * j * t / 365.0)]
        names += [f"sin{j}", f"cos{j}"]
    for k in range(1, 7):  # Monday = baseline
        cols.append((dow == k).astype(float))
        names.append(f"dow_{k}")
    return np.column_stack(cols), names


def acf(x, nlags):
    x = x - x.mean()
    denom = np.dot(x, x)
    return np.array([np.dot(x[:-k], x[k:]) / denom for k in range(1, nlags + 1)])


def main():
    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    obs = d[d["observed"]].reset_index(drop=True)
    t, dow = obs["t"].values.astype(float), obs["dow"].values
    X, names = design(t, dow, J=2)

    lines = []
    A = lines.append
    A("RESPONSE-SCALE COMPARISON (OLS, intercept + trend + 2 harmonics + day-of-week)")
    A(f"{'scale':>8} {'R2':>7} {'resid skew':>11} {'resid kurt':>11} {'BP het. p':>10} {'rho1':>7}")
    resids = {}
    for label, y in [
        ("kwh", obs["kwh"].values),
        ("sqrt", np.sqrt(obs["kwh"].values)),
        ("log", np.log(obs["kwh"].values)),
    ]:
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        fit = X @ b
        r = y - fit
        r2 = 1 - r.var() / y.var()
        # Breusch-Pagan style test: regress squared resid on the fitted values
        Z = np.column_stack([np.ones_like(fit), fit])
        g = np.linalg.lstsq(Z, r**2, rcond=None)[0]
        ss_tot = ((r**2 - (r**2).mean()) ** 2).sum()
        ss_res = ((r**2 - Z @ g) ** 2).sum()
        lm = len(r) * (1 - ss_res / ss_tot)
        p_het = 1 - stats.chi2.cdf(lm, 1)
        adj = np.diff(t) == 1
        rho1 = np.corrcoef(r[:-1][adj], r[1:][adj])[0, 1]
        A(
            f"{label:>8} {r2:7.3f} {stats.skew(r):11.3f} {stats.kurtosis(r, fisher=False):11.2f} "
            f"{p_het:10.4f} {rho1:7.3f}"
        )
        resids[label] = (r, fit)
    A("")
    A("  Interpretation: the log scale should show the flattest variance-vs-mean")
    A("  relationship (largest BP p-value). Gaussian kurtosis is 3; values well")
    A("  above 3 motivate a t likelihood.")
    A("")

    A("HARMONIC ORDER, OLS on the log scale")
    A(f"{'J':>3} {'params':>7} {'R2':>7} {'adj R2':>8} {'resid sd':>9} {'BIC':>10}")
    for J in range(0, 5):
        XJ, _ = design(t, dow, J)
        b = np.linalg.lstsq(XJ, np.log(obs["kwh"].values), rcond=None)[0]
        r = np.log(obs["kwh"].values) - XJ @ b
        n, p = len(r), XJ.shape[1]
        r2 = 1 - r.var() / np.log(obs["kwh"].values).var()
        adj = 1 - (1 - r2) * (n - 1) / (n - p)
        bic = n * np.log((r**2).sum() / n) + p * np.log(n)
        A(f"{J:3d} {p:7d} {r2:7.3f} {adj:8.3f} {r.std(ddof=1):9.3f} {bic:10.1f}")
    A("")
    A("  Note: these OLS criteria ignore the strong residual autocorrelation and so")
    A("  overstate the evidence for extra harmonics. They are exploratory only; the")
    A("  formal comparison in step 3 uses WAIC and held-out forecast scores.")
    A("")

    r_log = resids["log"][0]
    A("RESIDUAL AUTOCORRELATION, log scale, J = 2")
    a = acf(r_log, 21)
    for k in (1, 2, 3, 4, 5, 6, 7, 14, 21):
        A(f"  lag {k:2d}  {a[k-1]:+.3f}")
    A("")
    A("  Geometric decay consistent with AR(1); an AR(1) with rho1 as fitted predicts:")
    for k in (2, 3, 7):
        A(f"  lag {k:2d}  {a[0]**k:+.3f} (AR(1) prediction) vs {a[k-1]:+.3f} (observed)")
    A("")
    A("LEFT-TAIL DAYS ON THE LOG SCALE (standardised OLS residual < -2.5)")
    z = r_log / r_log.std(ddof=1)
    for i in np.where(z < -2.5)[0]:
        A(f"  {obs['date'].iloc[i].date()}  kwh={obs['kwh'].iloc[i]:6.2f}  z={z[i]:+.2f}")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(ROOT, "results", "01b_transform_check.txt"), "w") as f:
        f.write(txt + "\n")

    fig, ax = plt.subplots(2, 3, figsize=(13.5, 7))
    for j, label in enumerate(["kwh", "sqrt", "log"]):
        r, fit = resids[label]
        ax[0, j].scatter(fit, r, s=3, alpha=0.4, color="#1f4e79")
        ax[0, j].axhline(0, color="k", lw=0.8)
        ax[0, j].set_title(f"{label}: residual vs fitted")
        stats.probplot(r / r.std(ddof=1), dist="norm", plot=ax[1, j])
        ax[1, j].set_title(f"{label}: normal Q-Q")
        ax[1, j].get_lines()[0].set(ms=2.5, alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01b_scale_diagnostics.png"), dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 3.4))
    lags = np.arange(1, 22)
    ax.bar(lags, a, color="#1f4e79")
    ax.plot(lags, a[0] ** lags, "o-", color="#b03a2e", ms=4, label=r"AR(1) fit: $\rho_1^k$")
    ci = 1.96 / np.sqrt(len(r_log))
    ax.axhline(ci, color="gray", ls="--", lw=0.8)
    ax.axhline(-ci, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("lag (days)")
    ax.set_ylabel("residual autocorrelation")
    ax.set_title("Residual ACF after trend + 2 harmonics + day-of-week (log scale)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01b_residual_acf.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
