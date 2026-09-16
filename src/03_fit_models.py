"""
Step 3: fit the model family to the full daily series and compare.

Seven specifications, all sharing the same priors, the same data, and the same
set of pointwise likelihood terms (days whose two preceding days are also
observed) so that WAIC compares like with like. Because the log-scale models
carry a Jacobian correction, every log-likelihood is a density on the original
kWh scale and all seven are directly comparable.

    M1  J=1  iid Gaussian errors, log scale      the naive baseline
    M2  J=1  AR(1) + t errors, log scale
    M3  J=2  AR(1) + t errors, log scale         adds the second harmonic
    M4  J=3  AR(1) + t errors, log scale         has the harmonic order gone too far?
    M5  J=2  AR(2) + t errors, log scale         does the error need longer memory?
    M6  J=2  AR(1) + Gaussian errors, log scale  isolates the value of heavy tails
    M7  J=2  AR(1) + t errors, kWh scale         isolates the value of the log scale

Outputs
    results/03_fits/<id>.npz     posterior draws
    results/03_convergence.txt   Rhat / ESS for every parameter of every model
    results/03_waic.txt          WAIC comparison table
    figures/03_trace_*.png       trace and density plots
"""
import os
import pickle
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import model as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FITDIR = os.path.join(ROOT, "results", "03_fits")
os.makedirs(FITDIR, exist_ok=True)

SPECS = [
    ("M1", dict(J=1, p=0, heavy=False, scale="log"), "J=1, iid Gaussian, log"),
    ("M2", dict(J=1, p=1, heavy=True, scale="log"), "J=1, AR(1)+t, log"),
    ("M3", dict(J=2, p=1, heavy=True, scale="log"), "J=2, AR(1)+t, log"),
    ("M4", dict(J=3, p=1, heavy=True, scale="log"), "J=3, AR(1)+t, log"),
    ("M5", dict(J=2, p=2, heavy=True, scale="log"), "J=2, AR(2)+t, log"),
    ("M6", dict(J=2, p=1, heavy=False, scale="log"), "J=2, AR(1)+Gaussian, log"),
    ("M7", dict(J=2, p=1, heavy=True, scale="identity"), "J=2, AR(1)+t, kWh"),
]

MCMC = dict(n_iter=40000, burn=10000, thin=5, n_chains=4)


def load():
    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    return d


def param_table(fitres):
    """Yield (label, draws(chain,iter)) for every scalar parameter."""
    for i, nm in enumerate(fitres["names"]):
        yield nm, fitres["beta"][:, :, i]
    yield "sigma", fitres["sigma"]
    for i in range(fitres["p"]):
        yield f"phi_{i+1}", fitres["phi"][:, :, i]
    if fitres["heavy"]:
        yield "nu", fitres["nu"]


def main():
    d = load()
    y = d["kwh"].values
    observed = d["observed"].values.astype(bool)
    t = d["t"].values
    dow = d["dow"].values
    t_center = float(t.mean())
    ll_index = M.common_ll_index(observed, max_lag=2)
    print(f"T = {len(y)} days, observed = {observed.sum()}, "
          f"shared likelihood terms = {len(ll_index)}")

    fits, conv_lines, waic_rows = {}, [], []
    conv_lines.append("CONVERGENCE DIAGNOSTICS")
    conv_lines.append(
        f"{MCMC['n_chains']} chains x {MCMC['n_iter']} iterations, "
        f"{MCMC['burn']} discarded, thin {MCMC['thin']} "
        f"-> {MCMC['n_chains']*(MCMC['n_iter']-MCMC['burn'])//MCMC['thin']} draws"
    )
    conv_lines.append("Rhat is rank-normalised split-Rhat; ESS is bulk effective sample size.")

    for mid, spec, label in SPECS:
        t0 = time.time()
        print(f"\nfitting {mid}: {label} ...")
        fitres = M.fit_chains(
            y, observed, t, dow, ll_index=ll_index, seed=1000 + 7 * int(mid[1:]), **spec, **MCMC
        )
        fits[mid] = fitres
        el = time.time() - t0

        rhats, esss = [], []
        conv_lines.append("")
        conv_lines.append(f"--- {mid}: {label}   ({el:.1f} s) ---")
        conv_lines.append(f"{'parameter':>16} {'mean':>10} {'sd':>9} {'2.5%':>10} {'97.5%':>10} {'Rhat':>7} {'ESS':>8}")
        for nm, dr in param_table(fitres):
            flat = dr.reshape(-1)
            rh, es = M.diagnose(dr)
            rhats.append(rh)
            esss.append(es)
            lo, hi = np.percentile(flat, [2.5, 97.5])
            conv_lines.append(
                f"{nm:>16} {flat.mean():10.4f} {flat.std(ddof=1):9.4f} {lo:10.4f} {hi:10.4f} {rh:7.4f} {es:8.0f}"
            )
        conv_lines.append(f"  worst Rhat {max(rhats):.4f} (target < 1.01), "
                          f"min ESS {min(esss):.0f} (target > 400)")
        print(f"   done in {el:.1f}s  worst Rhat {max(rhats):.4f}  min ESS {min(esss):.0f}")

        w = M.waic(fitres["loglik"])
        waic_rows.append(
            dict(id=mid, label=label, elpd=w["elpd_waic"], p_waic=w["p_waic"],
                 waic=w["waic"], se=w["se"], n=w["n"], elpd_i=w["elpd_i"],
                 worst_rhat=max(rhats), min_ess=min(esss), secs=el)
        )

        np.savez_compressed(
            os.path.join(FITDIR, f"{mid}.npz"),
            **{k: v for k, v in fitres.items() if isinstance(v, np.ndarray)},
        )
        with open(os.path.join(FITDIR, f"{mid}_meta.pkl"), "wb") as f:
            pickle.dump({k: v for k, v in fitres.items() if not isinstance(v, np.ndarray)}, f)

    with open(os.path.join(ROOT, "results", "03_convergence.txt"), "w") as f:
        f.write("\n".join(conv_lines) + "\n")

    # ---------------- WAIC table ----------------
    best = max(waic_rows, key=lambda r: r["elpd"])
    lines = []
    A = lines.append
    A("WAIC MODEL COMPARISON")
    A(f"All log-likelihoods are densities of daily kWh on {best['n']} shared days")
    A("(log-scale models carry the log-Jacobian), so all rows are comparable.")
    A("elpd_waic: higher is better. d_elpd is relative to the best model;")
    A("se_diff is the standard error of that paired difference.")
    A("")
    A(f"{'id':>4} {'specification':>28} {'elpd_waic':>10} {'p_waic':>8} {'d_elpd':>9} {'se_diff':>8} {'Rhat':>7} {'ESS':>7}")
    for r in sorted(waic_rows, key=lambda r: -r["elpd"]):
        dif = r["elpd"] - best["elpd"]
        if r["id"] == best["id"]:
            se_d = 0.0
        else:
            di = r["elpd_i"] - best["elpd_i"]
            se_d = np.sqrt(len(di) * di.var(ddof=1))
        A(f"{r['id']:>4} {r['label']:>28} {r['elpd']:10.1f} {r['p_waic']:8.1f} "
          f"{dif:9.1f} {se_d:8.1f} {r['worst_rhat']:7.4f} {r['min_ess']:7.0f}")
    A("")
    A("TARGETED COMPARISONS (each isolates one modelling decision)")
    byid = {r["id"]: r for r in waic_rows}

    def cmp(a, b, what):
        di = byid[a]["elpd_i"] - byid[b]["elpd_i"]
        d = di.sum()
        se = np.sqrt(len(di) * di.var(ddof=1))
        verdict = "supported" if d > 2 * se else ("no support" if d < -2 * se else "inconclusive")
        A(f"  {what:<44} {a} - {b} = {d:+8.1f} +/- {se:5.1f}   {verdict}")

    cmp("M2", "M1", "AR(1) + t errors vs iid Gaussian")
    cmp("M3", "M2", "second harmonic (J=2 vs J=1)")
    cmp("M4", "M3", "third harmonic (J=3 vs J=2)")
    cmp("M5", "M3", "AR(2) vs AR(1)")
    cmp("M3", "M6", "t errors vs Gaussian, given AR(1)")
    cmp("M3", "M7", "log scale vs kWh scale")
    A("")
    A("  'supported' means the elpd difference exceeds two standard errors.")
    A("  WAIC on a conditional (one-step-ahead) factorisation understates how much")
    A("  a term matters for long-horizon forecasts, so step 4 re-tests each of")
    A("  these decisions on genuinely held-out data.")

    txt = "\n".join(lines)
    print("\n" + txt)
    with open(os.path.join(ROOT, "results", "03_waic.txt"), "w") as f:
        f.write(txt + "\n")

    # ---------------- trace plots for the headline model ----------------
    for mid in ("M3", "M5"):
        fitres = fits[mid]
        pars = list(param_table(fitres))
        n = len(pars)
        fig, ax = plt.subplots(n, 2, figsize=(11, 1.5 * n), gridspec_kw={"width_ratios": [2.2, 1]})
        for i, (nm, dr) in enumerate(pars):
            for c in range(dr.shape[0]):
                ax[i, 0].plot(dr[c], lw=0.35, alpha=0.75)
            ax[i, 0].set_ylabel(nm, fontsize=8, rotation=0, ha="right", va="center")
            ax[i, 0].tick_params(labelsize=6)
            for c in range(dr.shape[0]):
                v = dr[c]
                xs = np.linspace(v.min(), v.max(), 120)
                kde = np.exp(-0.5 * ((xs[:, None] - v[None, :]) / (0.9 * v.std() * len(v) ** -0.2)) ** 2).sum(axis=1)
                ax[i, 1].plot(xs, kde / kde.max(), lw=0.8)
            ax[i, 1].tick_params(labelsize=6)
            ax[i, 1].set_yticks([])
        ax[0, 0].set_title(f"{mid}: trace (4 chains)", fontsize=9)
        ax[0, 1].set_title("marginal posterior by chain", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(ROOT, "figures", f"03_trace_{mid}.png"), dpi=120)
        plt.close(fig)

    print(f"\nwrote {len(SPECS)} fits to results/03_fits/, convergence + WAIC tables, trace plots")


if __name__ == "__main__":
    main()
