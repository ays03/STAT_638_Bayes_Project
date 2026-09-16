"""
Step 1c: how much memory does the error process have?

The step-1b ACF decays far more slowly than AR(1) implies (lag 7 is +0.27 where
AR(1) predicts +0.02), so this script fits AR(p) error models by exact
conditional least squares on gap-free runs of observed days and compares orders.
It also checks whether the slow decay is an artefact of the multi-week
absence blocks by refitting with those days removed.
"""
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def design(t, dow, J):
    cols = [np.ones_like(t, dtype=float), t / 365.0]
    for j in range(1, J + 1):
        cols += [np.sin(2 * np.pi * j * t / 365.0), np.cos(2 * np.pi * j * t / 365.0)]
    for k in range(1, 7):
        cols.append((dow == k).astype(float))
    return np.column_stack(cols)


def ar_fit(e, t, p):
    """Conditional LS for an AR(p) on residuals e observed at day indices t.
    Only uses rows where all p lags are the immediately preceding calendar days."""
    idx = {int(tt): i for i, tt in enumerate(t)}
    rows_y, rows_X = [], []
    for i, tt in enumerate(t):
        lags = [idx.get(int(tt) - k) for k in range(1, p + 1)]
        if any(l is None for l in lags):
            continue
        rows_y.append(e[i])
        rows_X.append([e[l] for l in lags])
    y = np.array(rows_y)
    X = np.array(rows_X).reshape(len(y), p)
    if p == 0:
        return np.zeros(0), y.std(ddof=1), len(y), y
    phi = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ phi
    return phi, r.std(ddof=1), len(y), r


def report(e, t, label, A):
    A(f"\n{label}   (n = {len(e)} days)")
    A(f"{'p':>3} {'n_eff':>6} {'sigma_eta':>10} {'AIC':>10} {'BIC':>10}  coefficients")
    base_n = None
    for p in range(0, 6):
        phi, s, n, r = ar_fit(e, t, p)
        if base_n is None:
            _, _, base_n, _ = ar_fit(e, t, 5)
        # refit on a common sample size so AIC/BIC are comparable
        aic = n * np.log(s**2) + 2 * (p + 1)
        bic = n * np.log(s**2) + (p + 1) * np.log(n)
        co = "  ".join(f"{v:+.3f}" for v in phi)
        A(f"{p:3d} {n:6d} {s:10.4f} {aic:10.1f} {bic:10.1f}  {co}")
    phi2, s2, _, _ = ar_fit(e, t, 2)
    if len(phi2) == 2:
        # stationarity and implied long-run variance for an AR(2)
        p1, p2 = phi2
        stat = (p1 + p2 < 1) and (p2 - p1 < 1) and (abs(p2) < 1)
        rho1 = p1 / (1 - p2)
        marg = s2**2 * (1 - p2) / ((1 + p2) * ((1 - p2) ** 2 - p1**2))
        A(f"    AR(2): stationary={stat}, implied rho1={rho1:.3f}, marginal sd={np.sqrt(marg):.3f}")
        A(f"    ratio marginal sd / one-step sd = {np.sqrt(marg)/s2:.2f}  "
          f"(how much wider a long-horizon interval is than a one-day-ahead one)")


def main():
    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    obs = d[d["observed"]].reset_index(drop=True)
    t, dow = obs["t"].values.astype(float), obs["dow"].values
    y = np.log(obs["kwh"].values)
    X = design(t, dow, J=2)
    e = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]

    lines = []
    A = lines.append
    A("AR ORDER SELECTION FOR THE ERROR PROCESS")
    A("Residuals from OLS on log kWh with intercept + trend + 2 harmonics + day-of-week.")
    report(e, t, "ALL OBSERVED DAYS", A)

    # Identify sustained absence blocks: >= 5 consecutive days below 12 kWh.
    low = obs["kwh"].values < 12.0
    block = np.zeros(len(obs), dtype=bool)
    i = 0
    while i < len(obs):
        if low[i]:
            j = i
            while j + 1 < len(obs) and low[j + 1]:
                j += 1
            if j - i + 1 >= 5:
                block[i : j + 1] = True
            i = j + 1
        else:
            i += 1
    A("")
    A(f"SUSTAINED LOW-USE BLOCKS (>=5 consecutive days under 12 kWh): {int(block.sum())} days")
    runs = []
    i = 0
    while i < len(obs):
        if block[i]:
            j = i
            while j + 1 < len(obs) and block[j + 1]:
                j += 1
            runs.append((obs["date"].iloc[i].date(), obs["date"].iloc[j].date(), j - i + 1,
                         obs["kwh"].iloc[i:j+1].mean()))
            i = j + 1
        else:
            i += 1
    for s0, s1, n, m in runs:
        A(f"  {s0} to {s1}  ({n:2d} days, mean {m:5.2f} kWh)")

    keep = ~block
    Xk = X[keep]
    yk = y[keep]
    ek = yk - Xk @ np.linalg.lstsq(Xk, yk, rcond=None)[0]
    report(ek, t[keep], "EXCLUDING SUSTAINED LOW-USE BLOCKS", A)
    A("")
    A("  If the preferred order drops sharply once the absence blocks are removed, the")
    A("  long memory is largely a regime effect and a heavy-tailed AR(2) is the right")
    A("  compromise: AR(2) for genuine weather persistence, t errors for the blocks.")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(ROOT, "results", "01c_ar_order.txt"), "w") as f:
        f.write(txt + "\n")

    # PACF-style plot via successive AR fits
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for a_, ee, tt, lab in [(ax[0], e, t, "all days"), (ax[1], ek, t[keep], "absence blocks removed")]:
        pac = [ar_fit(ee, tt, p)[0][-1] for p in range(1, 11)]
        a_.bar(np.arange(1, 11), pac, color="#1f4e79")
        ci = 1.96 / np.sqrt(len(ee))
        a_.axhline(ci, color="gray", ls="--", lw=0.8)
        a_.axhline(-ci, color="gray", ls="--", lw=0.8)
        a_.set_title(f"Partial autocorrelation: {lab}")
        a_.set_xlabel("lag (days)")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01c_pacf.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
