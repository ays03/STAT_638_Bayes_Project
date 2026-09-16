"""
Step 4: how accurately can future demand be predicted, and how does predictive
uncertainty grow with the forecast horizon?

WAIC scores one-step-ahead fit, which is the wrong target for a forecasting
question, so every modelling decision is re-tested out of sample. Models are
refitted on data through 30 Nov 2009 and never see the final twelve months.

Two evaluation designs, because they answer different questions:

  A. SINGLE ORIGIN, 12-month forecast from 30 Nov 2009.
     This is the honest "here is what you could have said a year in advance"
     test. It gives the fan chart and the coverage table, but the 361 horizons
     come from one realisation, so horizon-by-horizon scores are correlated
     and noisy.

  B. ROLLING ORIGIN, horizons 1-60 from every 3rd day of the holdout year.
     Parameters still come only from the training period, but each forecast
     conditions on actuals up to its own origin. Recovering the error state at
     a new origin needs no refitting: given beta, eps_t = z_t - x_t'beta is
     deterministic wherever y_t is observed. About 120 origins per horizon
     gives a clean decay curve with real replication.

Scores: log predictive density (exact mixture, integrating the final
innovation analytically), CRPS (exact ensemble formula), and empirical coverage
of central 50% and 95% intervals. All scores are on the kWh scale.
"""
import os
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
OUT = os.path.join(ROOT, "results", "04_holdout")
os.makedirs(OUT, exist_ok=True)

TRAIN_END = pd.Timestamp("2009-11-30")
SPECS = [
    ("M1", dict(J=1, p=0, heavy=False, scale="log"), "J=1, iid Gaussian, log"),
    ("M2", dict(J=1, p=1, heavy=True, scale="log"), "J=1, AR(1)+t, log"),
    ("M3", dict(J=2, p=1, heavy=True, scale="log"), "J=2, AR(1)+t, log"),
    ("M4", dict(J=3, p=1, heavy=True, scale="log"), "J=3, AR(1)+t, log"),
    ("M5", dict(J=2, p=2, heavy=True, scale="log"), "J=2, AR(2)+t, log"),
    ("M6", dict(J=2, p=1, heavy=False, scale="log"), "J=2, AR(1)+Gaussian, log"),
    ("M7", dict(J=2, p=1, heavy=True, scale="identity"), "J=2, AR(1)+t, kWh"),
]
MCMC = dict(n_iter=30000, burn=8000, thin=5, n_chains=4)
BUCKETS = [(1, 7), (8, 30), (31, 90), (91, 365)]
ROLL_H = 60
ROLL_STEP = 3


def bucket_of(h):
    for lo, hi in BUCKETS:
        if lo <= h <= hi:
            return f"{lo}-{hi}"
    return None


def climatology_baseline(train, test):
    """Reference forecast: the mean of log kWh for the same day-of-year window
    (+/- 7 days) and weekend/weekday status, from the training years only.
    Any model that cannot beat this has learned nothing useful."""
    tr = train[train["observed"]]
    preds, sds = [], []
    for _, r in test.iterrows():
        doy, wk = r["doy"], r["dow"] >= 5
        dd = np.minimum(np.abs(tr["doy"] - doy), 365 - np.abs(tr["doy"] - doy))
        sel = tr[(dd <= 7) & ((tr["dow"] >= 5) == wk)]
        if len(sel) < 5:
            sel = tr[dd <= 14]
        lz = np.log(sel["kwh"].values)
        preds.append(lz.mean())
        sds.append(lz.std(ddof=1))
    return np.array(preds), np.array(sds)


def main():
    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    train = d[d["date"] <= TRAIN_END].reset_index(drop=True)
    test = d[d["date"] > TRAIN_END].reset_index(drop=True)
    t_center = float(train["t"].values.mean())
    print(f"train: {train['date'].min().date()} to {train['date'].max().date()} "
          f"({len(train)} days, {train['observed'].sum()} observed)")
    print(f"test : {test['date'].min().date()} to {test['date'].max().date()} "
          f"({len(test)} days, {test['observed'].sum()} observed)")

    y_test = test["kwh"].values
    obs_test = test["observed"].values.astype(bool)
    h_arr = np.arange(1, len(test) + 1)

    # One common set of rolling origins for every model, using the strictest
    # requirement in the family (AR(2) needs the origin and the day before it
    # observed). Identical origins mean identical target days, so the paired
    # score differences below compare like with like.
    MAX_P = max(s["p"] for _, s, _ in SPECS)
    origins = [
        o
        for o in range(0, len(test) - 1, ROLL_STEP)
        if o - MAX_P + 1 >= 0 and all(obs_test[o - k] for k in range(MAX_P))
    ]
    print(f"{len(origins)} rolling origins shared by all models "
          f"(every {ROLL_STEP}rd day, horizons 1-{ROLL_H})")

    rows_single, roll_store, fits = [], {}, {}

    for mid, spec, label in SPECS:
        t0 = time.time()
        print(f"\nfitting {mid} on training data ...")
        fitres = M.fit_chains(
            train["kwh"].values,
            train["observed"].values.astype(bool),
            train["t"].values,
            train["dow"].values,
            seed=5000 + 7 * int(mid[1:]),
            n_paths=1,
            **spec,
            **MCMC,
        )
        fits[mid] = (fitres, spec, label)

        # ---------- design A: single origin, whole holdout year ----------
        fc = M.forecast(
            fitres, test["t"].values, test["dow"].values, t_center,
            n_draws=4000, seed=99 + int(mid[1:]),
        )
        lpd = M.log_pred_density(fc, y_test)
        cr = M.crps(fc["y"], y_test)
        pit_vals = M.pit(fc, y_test)
        q = np.percentile(fc["y"], [2.5, 25, 75, 97.5], axis=0)
        cov95 = (y_test >= q[0]) & (y_test <= q[3])
        cov50 = (y_test >= q[1]) & (y_test <= q[2])
        width95 = q[3] - q[0]
        np.savez_compressed(
            os.path.join(OUT, f"single_{mid}.npz"),
            lpd=lpd, crps=cr, q=q, cov95=cov95, cov50=cov50, width95=width95, pit=pit_vals,
            ypred_mean=fc["y"].mean(axis=0), ypred_med=np.median(fc["y"], axis=0),
        )
        rows_single.append(dict(id=mid, label=label, lpd=lpd, crps=cr,
                                cov95=cov95, cov50=cov50, width95=width95, obs=obs_test))

        # ---------- design B: rolling origins ----------
        beta = M.flatten(fitres, "beta")
        p = spec["p"]
        Xtest, _ = M.build_design(test["t"].values, test["dow"].values, spec["J"], t_center=t_center)
        z_test = np.log(y_test) if spec["scale"] == "log" else y_test
        # eps on test days, per posterior draw: (N, n_test)
        eps_test = z_test[None, :] - beta @ Xtest.T

        acc = {h: {"lpd": [], "crps": [], "c95": [], "c50": [], "w95": [], "w50": [], "day": []}
               for h in range(1, ROLL_H + 1)}
        rng = np.random.default_rng(int(mid[1:]))
        n_sub = 1500
        sub = rng.choice(beta.shape[0], min(n_sub, beta.shape[0]), replace=False)
        for oi, o in enumerate(origins):
            H = min(ROLL_H, len(test) - 1 - o)
            if H <= 0:
                continue
            idx = np.arange(o + 1, o + 1 + H)
            e_init = (
                np.column_stack([eps_test[sub][:, o - k] for k in range(p)]) if p > 0 else None
            )
            fcr = M.forecast(
                fitres, test["t"].values[idx], test["dow"].values[idx], t_center,
                eps_init=e_init, draw_idx=sub, seed=7 * oi + 3,
            )
            lp = M.log_pred_density(fcr, y_test[idx])
            cc = M.crps(fcr["y"], y_test[idx])
            qq = np.percentile(fcr["y"], [2.5, 25, 75, 97.5], axis=0)
            for j, h in enumerate(range(1, H + 1)):
                if not obs_test[idx[j]]:
                    continue
                acc[h]["lpd"].append(lp[j])
                acc[h]["crps"].append(cc[j])
                acc[h]["c95"].append(qq[0, j] <= y_test[idx[j]] <= qq[3, j])
                acc[h]["c50"].append(qq[1, j] <= y_test[idx[j]] <= qq[2, j])
                acc[h]["w95"].append(qq[3, j] - qq[0, j])
                acc[h]["w50"].append(qq[2, j] - qq[1, j])
                acc[h]["day"].append(idx[j])
        roll_store[mid] = {h: {k: np.array(v) for k, v in acc[h].items()} for h in acc}
        # the scored target days must be identical across models for pairing
        first = roll_store[SPECS[0][0]]
        for h in range(1, ROLL_H + 1):
            assert np.array_equal(roll_store[mid][h]["day"], first[h]["day"]), (
                f"{mid} scored different days than {SPECS[0][0]} at h={h}"
            )
        np.savez_compressed(
            os.path.join(OUT, f"roll_{mid}.npz"),
            **{f"h{h}_{k}": v for h in acc for k, v in roll_store[mid][h].items()},
        )
        print(f"   {mid} done in {time.time()-t0:.1f}s  "
              f"({len(origins)} rolling origins, n per horizon ~{len(acc[1]['lpd'])})")

    # ---------- climatology reference ----------
    cl_mu, cl_sd = climatology_baseline(train, test)
    from scipy.stats import norm

    z_act = np.log(np.where(obs_test, y_test, np.nan))
    cl_lpd = norm.logpdf(z_act, cl_mu, cl_sd) - z_act
    cl_q = np.exp(norm.ppf([[0.025], [0.25], [0.75], [0.975]]) * cl_sd + cl_mu)
    cl_cov95 = (y_test >= cl_q[0]) & (y_test <= cl_q[3])
    cl_samp = np.exp(np.random.default_rng(0).standard_normal((4000, len(test))) * cl_sd + cl_mu)
    cl_crps = M.crps(cl_samp, np.where(obs_test, y_test, np.nan))

    # ================= report =================
    lines = []
    A = lines.append
    A("OUT-OF-SAMPLE FORECAST EVALUATION")
    A(f"Trained on {train['observed'].sum()} observed days through {TRAIN_END.date()};")
    A(f"scored on {obs_test.sum()} observed days from {test['date'].min().date()} "
      f"to {test['date'].max().date()}.")
    A("Scores are per day on the kWh scale. lpd: higher is better (log predictive")
    A("density). CRPS: lower is better, in kWh. Coverage should match the nominal level.")
    A("")
    A("=== DESIGN A: single origin, 12-month-ahead forecast ===")
    A("")
    A(f"{'id':>4} {'specification':>28} {'mean lpd':>9} {'CRPS':>7} {'cov95':>7} {'cov50':>7} {'mean width95':>13}")
    for r in rows_single:
        m = r["obs"]
        A(f"{r['id']:>4} {r['label']:>28} {np.nanmean(r['lpd']):9.3f} {np.nanmean(r['crps']):7.3f} "
          f"{100*r['cov95'][m].mean():6.1f}% {100*r['cov50'][m].mean():6.1f}% "
          f"{r['width95'][m].mean():12.2f}")
    A(f"{'--':>4} {'day-of-year climatology':>28} {np.nanmean(cl_lpd):9.3f} "
      f"{np.nanmean(cl_crps):7.3f} {100*cl_cov95[obs_test].mean():6.1f}% "
      f"{'':>7} {(cl_q[3]-cl_q[0])[obs_test].mean():12.2f}")
    A("")
    A("Scores by horizon bucket (CRPS in kWh, coverage of the 95% interval):")
    hdr = f"{'id':>4}"
    for lo, hi in BUCKETS:
        hdr += f" | {f'{lo}-{hi}d':>18}"
    A(hdr)
    A(f"{'':>4}" + "".join(f" | {'CRPS  lpd  cov':>18}" for _ in BUCKETS))
    for r in rows_single:
        line = f"{r['id']:>4}"
        for lo, hi in BUCKETS:
            m = r["obs"] & (h_arr >= lo) & (h_arr <= hi)
            if m.sum() == 0:
                line += f" | {'--':>18}"
            else:
                line += (f" | {np.nanmean(r['crps'][m]):5.2f} {np.nanmean(r['lpd'][m]):5.2f} "
                         f"{100*r['cov95'][m].mean():5.1f}%")
        A(line)
    A("")
    A("  Caveat: the 361 horizons come from a single realisation of the holdout")
    A("  year, so these bucket scores are correlated. Design B is the replicated")
    A("  version and is the one to quote for the horizon question.")
    A("")

    A("=== DESIGN B: rolling origins, horizons 1-60 ===")
    A("")
    A(f"{'id':>4} " + " ".join(f"{'h='+str(h):>7}" for h in (1, 2, 3, 5, 7, 14, 21, 30, 45, 60)))
    for key, fmt, title in [
        ("crps", "{:7.3f}", "mean CRPS (kWh, lower better)"),
        ("lpd", "{:7.3f}", "mean log predictive density (higher better)"),
        ("w95", "{:7.2f}", "mean 95% interval width (kWh)"),
        ("c95", "{:6.1f}%", "coverage of 95% interval"),
        ("c50", "{:6.1f}%", "coverage of 50% interval"),
    ]:
        A("")
        A(f"  {title}")
        for mid, _, _ in SPECS:
            vals = []
            for h in (1, 2, 3, 5, 7, 14, 21, 30, 45, 60):
                v = roll_store[mid][h][key]
                # only the coverage rows are percentages; "crps" also starts
                # with a c, so name the coverage keys explicitly
                v = 100 * v.mean() if key in ("c95", "c50") else np.nanmean(v)
                vals.append(fmt.format(v))
            A(f"{mid:>4} " + " ".join(vals))
    A("")
    A("HOW UNCERTAINTY GROWS WITH HORIZON (best model by design-A lpd)")
    best_id = max(rows_single, key=lambda r: np.nanmean(r["lpd"]))["id"]
    rs = roll_store[best_id]
    w1 = rs[1]["w95"].mean()
    A(f"  model {best_id}")
    A(f"{'h':>4} {'width95':>9} {'ratio to h=1':>13} {'CRPS':>7} {'cov95':>7}")
    for h in (1, 2, 3, 4, 5, 7, 10, 14, 21, 30, 45, 60):
        w = rs[h]["w95"].mean()
        A(f"{h:>4} {w:9.2f} {w/w1:13.3f} {np.nanmean(rs[h]['crps']):7.3f} "
          f"{100*rs[h]['c95'].mean():6.1f}%")
    A("")
    A("  The saturation point is the answer to the horizon question: once the")
    A("  ratio stops rising, recent consumption has stopped being informative and")
    A("  the forecast is seasonal-plus-weekly climatology.")
    A("")
    A("MODELLING DECISIONS RE-TESTED OUT OF SAMPLE (paired differences in lpd)")
    A("  positive favours the first model; se is the paired standard error")

    def paired(a, b, what, design="A"):
        if design == "A":
            ra = next(r for r in rows_single if r["id"] == a)
            rb = next(r for r in rows_single if r["id"] == b)
            m = ra["obs"]
            dv = ra["lpd"][m] - rb["lpd"][m]
        else:
            dv = np.concatenate(
                [roll_store[a][h]["lpd"] - roll_store[b][h]["lpd"] for h in range(1, ROLL_H + 1)]
            )
        d, se = np.nanmean(dv), np.nanstd(dv, ddof=1) / np.sqrt(np.sum(np.isfinite(dv)))
        verdict = "supported" if d > 2 * se else ("no support" if d < -2 * se else "inconclusive")
        A(f"  [{design}] {what:<42} {a}-{b} = {d:+.4f} +/- {se:.4f}  {verdict}")

    for design in ("A", "B"):
        paired("M2", "M1", "AR(1)+t vs iid Gaussian", design)
        paired("M3", "M2", "second harmonic (J=2 vs J=1)", design)
        paired("M4", "M3", "third harmonic (J=3 vs J=2)", design)
        paired("M5", "M3", "AR(2) vs AR(1)", design)
        paired("M3", "M6", "t errors vs Gaussian", design)
        paired("M3", "M7", "log scale vs kWh scale", design)
        A("")

    txt = "\n".join(lines)
    print("\n" + txt)
    with open(os.path.join(ROOT, "results", "04_holdout.txt"), "w") as f:
        f.write(txt + "\n")

    # ================= figures =================
    fig, ax = plt.subplots(2, 2, figsize=(13, 7.5))
    # (a) interval width vs horizon, rolling design
    for mid, _, label in SPECS:
        hs = np.arange(1, ROLL_H + 1)
        w = [roll_store[mid][h]["w95"].mean() for h in hs]
        ax[0, 0].plot(hs, w, lw=1.3, label=mid)
    ax[0, 0].set_xlabel("forecast horizon h (days)")
    ax[0, 0].set_ylabel("mean width of 95% interval (kWh)")
    ax[0, 0].set_title("(a) Predictive uncertainty vs horizon (rolling origins)")
    ax[0, 0].legend(fontsize=7, ncol=2)
    # (b) normalised width for the best model, with the AR(1) theoretical curve
    rs = roll_store[best_id]
    hs = np.arange(1, ROLL_H + 1)
    w = np.array([rs[h]["w95"].mean() for h in hs])
    ax[0, 1].plot(hs, w / w[0], "o-", ms=3, color="#1f4e79", label=f"{best_id} empirical")
    phi_bar = M.flatten(fits[best_id][0], "phi").mean(axis=0)
    if len(phi_bar) > 0:
        theo = [np.sqrt(sum(_ar_psi(phi_bar, h)[: h] ** 2)) for h in hs]
        theo = np.array(theo) / theo[0]
        ax[0, 1].plot(hs, theo, "--", color="#b03a2e",
                      label=r"AR($p$) theory: $\sqrt{\sum_{j<h}\psi_j^2}$")
        ax[0, 1].axhline(M.marginal_sd_ratio(phi_bar), color="gray", ls=":",
                         label="stationary limit")
    ax[0, 1].set_xlabel("forecast horizon h (days)")
    ax[0, 1].set_ylabel("width relative to h = 1")
    ax[0, 1].set_title("(b) Uncertainty saturates within about a week")
    ax[0, 1].legend(fontsize=8)
    # (c) CRPS vs horizon
    for mid, _, _ in SPECS:
        ax[1, 0].plot(hs, [np.nanmean(roll_store[mid][h]["crps"]) for h in hs], lw=1.2, label=mid)
    ax[1, 0].axhline(np.nanmean(cl_crps), color="k", ls="--", lw=1, label="climatology")
    ax[1, 0].set_xlabel("forecast horizon h (days)")
    ax[1, 0].set_ylabel("mean CRPS (kWh)")
    ax[1, 0].set_title("(c) Forecast accuracy vs horizon")
    ax[1, 0].legend(fontsize=7, ncol=2)
    # (d) coverage vs horizon
    for mid, _, _ in SPECS:
        ax[1, 1].plot(hs, [100 * roll_store[mid][h]["c95"].mean() for h in hs], lw=1.2, label=mid)
    ax[1, 1].axhline(95, color="k", ls="--", lw=1)
    ax[1, 1].set_ylim(60, 102)
    ax[1, 1].set_xlabel("forecast horizon h (days)")
    ax[1, 1].set_ylabel("empirical coverage of 95% interval (%)")
    ax[1, 1].set_title("(d) Calibration vs horizon")
    ax[1, 1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "04_horizon.png"), dpi=140)
    plt.close(fig)

    # holdout fan chart for the best model
    z = np.load(os.path.join(OUT, f"single_{best_id}.npz"))
    q = z["q"]
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.plot(train["date"].values[-200:], train["kwh"].values[-200:], lw=0.8,
            color="#444", label="training data")
    ax.fill_between(test["date"], q[0], q[3], color="#1f4e79", alpha=0.20, label="95% predictive")
    ax.fill_between(test["date"], q[1], q[2], color="#1f4e79", alpha=0.40, label="50% predictive")
    ax.plot(test["date"], z["ypred_med"], color="#1f4e79", lw=1.3, label="predictive median")
    ax.plot(test["date"], y_test, ".", ms=3, color="#b03a2e", label="actual (held out)")
    ax.axvline(TRAIN_END, color="k", ls="--", lw=1)
    ax.set_ylabel("kWh / day")
    ax.set_title(f"12-month-ahead forecast from {TRAIN_END.date()}, model {best_id} "
                 f"(nothing after the dashed line was used in fitting)")
    ax.legend(fontsize=8, ncol=5, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "04_holdout_fan.png"), dpi=140)
    plt.close(fig)

    # PIT histograms: calibration of the holdout predictive distributions.
    # A U shape means intervals are too narrow; a hump means too wide; a tilt
    # means the level is biased.
    show = [best_id] + [m for m in ("M1", "M6", "M7") if m != best_id][:3]
    fig, ax = plt.subplots(1, len(show), figsize=(3.4 * len(show), 3.4))
    for a_, mid in zip(np.atleast_1d(ax), show):
        v = np.load(os.path.join(OUT, f"single_{mid}.npz"))["pit"]
        v = v[np.isfinite(v)]
        a_.hist(v, bins=20, range=(0, 1), color="#1f4e79")
        a_.axhline(len(v) / 20, color="k", ls="--", lw=1)
        lab = next(l for i, _, l in SPECS if i == mid)
        a_.set_title(f"{mid}: {lab}", fontsize=8)
        a_.set_xlabel("predictive CDF at the actual")
    fig.suptitle("PIT histograms, 12-month-ahead holdout forecasts "
                 "(flat = calibrated, U = overconfident)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "04_pit.png"), dpi=140)
    plt.close(fig)

    print("\nwrote results/04_holdout.txt and 3 figures")


def _ar_psi(phi, h):
    """MA(inf) coefficients psi_0..psi_{h-1} of the AR(p)."""
    p = len(phi)
    psi = np.zeros(max(h, 1))
    psi[0] = 1.0
    for j in range(1, len(psi)):
        psi[j] = sum(phi[k] * psi[j - 1 - k] for k in range(min(p, j)))
    return psi


if __name__ == "__main__":
    main()
