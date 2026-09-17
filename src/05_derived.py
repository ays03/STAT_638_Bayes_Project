"""
Step 5: turn posterior draws into the quantities the scientific questions ask
about, and run the posterior predictive check for unusual days.

Everything here is a deterministic function of the saved draws, so each answer
carries a full posterior distribution rather than a point estimate:

  trend        kWh/year and %/year, with P(declining | data)
  seasonality  amplitude and phase of each harmonic in kWh, peak/trough dates,
               and the amplitude ratio that says whether harmonic 2 matters
  weekly       the seven day-of-week deviations from the weekly mean, in kWh
  anomalies    one-step-ahead standardised residuals and posterior predictive
               p-values p_t = P(y_t^rep <= y_t | y), flagged at p<0.01 / >0.99

Outputs results/05_derived.txt and three figures.
"""
import os
import pickle
import sys

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import model as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DAYNAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def load_fit(mid):
    z = np.load(os.path.join(ROOT, "results", "03_fits", f"{mid}.npz"), allow_pickle=True)
    meta = pickle.load(open(os.path.join(ROOT, "results", "03_fits", f"{mid}_meta.pkl"), "rb"))
    d = {k: z[k] for k in z.keys()}
    d.update(meta)
    return d


def ci(v, lo=2.5, hi=97.5):
    return np.percentile(v, [lo, hi])


def fmt(v, unit="", dp=3):
    m = v.mean()
    l, h = ci(v)
    return f"{m:.{dp}f} [{l:.{dp}f}, {h:.{dp}f}]{unit}"


def main():
    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    t = d["t"].values.astype(float)
    t_center = t.mean()
    date0 = d["date"].iloc[0]

    lines = []
    A = lines.append
    A("DERIVED POSTERIOR QUANTITIES")
    A("M7 = J=2, AR(1) + t errors, kWh scale  (best out-of-sample forecaster)")
    A("M3 = J=2, AR(1) + t errors, log scale  (best by WAIC among J=2 models)")
    A("Intervals are 95% central posterior credible intervals.")
    A("")

    for mid in ("M7", "M3"):
        f = load_fit(mid)
        nm = f["names"]
        b = M.flatten(f, "beta")
        sig = M.flatten(f, "sigma")
        phi = M.flatten(f, "phi")
        nu = M.flatten(f, "nu")
        scale = f["scale"]
        unit = "kWh" if scale == "identity" else "log kWh"
        A(f"================ {mid}  ({unit} scale) ================")

        # ---------- level ----------
        i0 = nm.index("intercept")
        # the intercept is the Monday level at the record midpoint; add the mean
        # day-of-week offset to get an average-day level
        dow_cols = [nm.index(f"dow_{x}") for x in DAYNAMES[1:]]
        dow_full = np.column_stack([np.zeros(len(b)), b[:, dow_cols]])  # Mon = 0
        level = b[:, i0] + dow_full.mean(axis=1)
        A(f"  average-day level at record midpoint : {fmt(level, ' ' + unit, 3)}")

        # ---------- Q2: trend ----------
        i1 = nm.index("trend_per_year")
        tr = b[:, i1]
        A("")
        A("  Q2 LONG-TERM TREND")
        A(f"    slope                              : {fmt(tr, ' ' + unit + '/year', 4)}")
        if scale == "identity":
            pct = 100 * tr / level
            A(f"    as a percentage of the mean level   : {fmt(pct, ' %/year', 2)}")
            A(f"    total change over the 3.95-y record: {fmt(3.95 * tr, ' kWh/day', 3)}")
        else:
            A(f"    as a percentage                    : {fmt(100*(np.exp(tr)-1), ' %/year', 2)}")
        A(f"    P(slope < 0 | data)                : {np.mean(tr < 0):.3f}")
        A(f"    P(|slope| < 1% of level per year)  : "
          f"{np.mean(np.abs(tr) < 0.01*np.abs(level)):.3f}")

        # ---------- Q3: seasonality ----------
        A("")
        A("  Q3 ANNUAL SEASONALITY")
        amps = []
        for j in (1, 2):
            if f"sin_h{j}" not in nm:
                continue
            a_ = b[:, nm.index(f"sin_h{j}")]
            c_ = b[:, nm.index(f"cos_h{j}")]
            R = np.sqrt(a_**2 + c_**2)
            amps.append(R)
            # a sin(w t) + c cos(w t) = R cos(w t - psi), psi = atan2(a, c)
            psi = np.arctan2(a_, c_)
            period = 365.0 / j
            peak_t = (psi % (2 * np.pi)) * period / (2 * np.pi)
            A(f"    harmonic {j}: amplitude            : {fmt(R, ' ' + unit, 3)}")
            A(f"                peak-to-trough swing : {fmt(2*R, ' ' + unit, 3)}")
            pk = [date0 + pd.Timedelta(days=float(x)) for x in np.percentile(peak_t, [2.5, 50, 97.5])]
            A(f"                peaks at day-of-cycle: {np.median(peak_t):.1f} of {period:.1f}"
              f"  (calendar {pk[1].strftime('%d %b')}, 95% CI {pk[0].strftime('%d %b')}"
              f" to {pk[2].strftime('%d %b')})")
        if len(amps) == 2:
            ratio = amps[1] / amps[0]
            A(f"    amplitude ratio  A2 / A1           : {fmt(ratio, '', 3)}")
            A(f"    P(A2 > 0.15 * A1 | data)           : {np.mean(ratio > 0.15):.3f}")
            A(f"    P(A2 > 0.25 * A1 | data)           : {np.mean(ratio > 0.25):.3f}")
            A("    (a second harmonic this large makes the annual curve visibly")
            A("     asymmetric: a sharper winter peak and a flatter summer floor)")

        # combined seasonal shape over one year
        tt = np.arange(1, 366)
        seas = np.zeros((len(b), len(tt)))
        for j in (1, 2):
            if f"sin_h{j}" in nm:
                seas += (b[:, [nm.index(f"sin_h{j}")]] * np.sin(2 * np.pi * j * tt / 365.0)
                         + b[:, [nm.index(f"cos_h{j}")]] * np.cos(2 * np.pi * j * tt / 365.0))
        rng_season = seas.max(axis=1) - seas.min(axis=1)
        A(f"    full seasonal range (max - min)    : {fmt(rng_season, ' ' + unit, 3)}")

        # ---------- Q4: weekly ----------
        A("")
        A("  Q4 WEEKLY PATTERN (deviation from the weekly mean)")
        prof = dow_full - dow_full.mean(axis=1, keepdims=True)
        for k, name in enumerate(DAYNAMES):
            A(f"    {name}                                : {fmt(prof[:, k], ' ' + unit, 3)}")
        wkend = prof[:, 5:].mean(axis=1) - prof[:, :5].mean(axis=1)
        A(f"    weekend minus weekday              : {fmt(wkend, ' ' + unit, 3)}")
        A(f"    P(weekend > weekday | data)        : {np.mean(wkend > 0):.4f}")

        # ---------- error process ----------
        A("")
        A("  ERROR PROCESS")
        A(f"    sigma (one-step innovation sd)     : {fmt(sig, ' ' + unit, 3)}")
        A(f"    phi_1                              : {fmt(phi[:, 0], '', 3)}")
        A(f"    nu (t degrees of freedom)          : {fmt(nu, '', 2)}")
        ratio_sd = 1.0 / np.sqrt(1.0 - phi[:, 0] ** 2)
        A(f"    marginal sd / one-step sd          : {fmt(ratio_sd, '', 4)}")
        A("    -> this ratio is the ceiling on how much wider a long-horizon")
        A("       predictive interval is than a one-day-ahead one")
        A("")

    # ================= Q6: anomalies, from M7 =================
    f = load_fit("M7")
    diag = M.one_step_diagnostics(f, None, d["observed"].values.astype(bool))
    d["std_resid"] = diag["std_resid_mean"]
    d["ppp"] = diag["ppp"]
    d["lam"] = diag["lam_mean"]

    A("================ Q6 UNUSUAL PERIODS (model M7) ================")
    A("One-step-ahead posterior predictive p-values p_t = P(y_t^rep <= y_t | y).")
    A("Under a well-specified model these are uniform, so about 1% of days fall")
    A("below 0.01 and 1% above 0.99 by chance alone.")
    ok = np.isfinite(d["ppp"].values)
    pv = d["ppp"].values[ok]
    n = ok.sum()
    A("")
    A(f"  evaluable days                    : {n}")
    A(f"  p_t < 0.01 (unusually low)        : {int((pv<0.01).sum())}  "
      f"(expected {0.01*n:.1f} if calibrated)")
    A(f"  p_t > 0.99 (unusually high)       : {int((pv>0.99).sum())}  "
      f"(expected {0.01*n:.1f} if calibrated)")
    A(f"  p_t < 0.05                        : {int((pv<0.05).sum())}  "
      f"(expected {0.05*n:.1f})")
    A(f"  p_t > 0.95                        : {int((pv>0.95).sum())}  "
      f"(expected {0.05*n:.1f})")
    A("")
    A("  Note these are ONE-STEP-AHEAD checks, so a sustained absence is flagged")
    A("  only on its first day or two: once the AR(1) term has absorbed the drop,")
    A("  the model expects the low level to continue. The blocks are therefore")
    A("  better seen in the lambda_t weights below.")
    A("")

    flag = d[ok & ((d["ppp"] < 0.01) | (d["ppp"] > 0.99))].copy()
    A(f"  FLAGGED DAYS ({len(flag)} total), most extreme first")
    A(f"    {'date':>12} {'dow':>4} {'kWh':>7} {'fitted':>7} {'z':>7} {'p_t':>7}")
    Xf = f["X"]
    mu_mean = (M.flatten(f, "beta") @ Xf.T).mean(axis=0)
    d["fitted"] = mu_mean
    flag["fitted"] = d.loc[flag.index, "fitted"]
    flag["extremity"] = np.minimum(flag["ppp"], 1 - flag["ppp"])
    for _, r in flag.sort_values("extremity").iterrows():
        A(f"    {r['date'].date()!s:>12} {DAYNAMES[int(r['dow'])]:>4} {r['kwh']:7.2f} "
          f"{r['fitted']:7.2f} {r['std_resid']:+7.2f} {r['ppp']:7.4f}")
    A("")

    # ---- block-level posterior predictive check ----
    # A one-step-ahead check cannot detect a sustained anomaly: once the AR(1)
    # term has absorbed a drop, the model expects the low level to continue, so
    # innovations return to normal size and lambda_t returns to 1. The question
    # asks about PERIODS, so the test statistic has to be a block statistic.
    # Here it is the mean consumption over a rolling window, referred to its
    # posterior predictive distribution under the stationary error process.
    WIN = 14
    rng = np.random.default_rng(3)
    pars = f["path_pars"]
    sig_s, nu_s, phi_s = pars[:, 0], pars[:, 1], pars[:, 3]
    S = len(pars)
    beta_s = M.flatten(f, "beta")
    sel = rng.choice(len(beta_s), S, replace=False)
    mu_s = beta_s[sel] @ Xf.T  # (S, T)
    T = Xf.shape[0]

    # simulate replicate series from the stationary AR(1) + t error process
    lam_s = rng.gamma(nu_s[:, None] / 2, 2.0 / nu_s[:, None], size=(S, T))
    eta_s = rng.standard_normal((S, T)) * sig_s[:, None] / np.sqrt(lam_s)
    eps_s = np.empty((S, T))
    eps_s[:, 0] = eta_s[:, 0] / np.sqrt(1 - phi_s**2)
    for i in range(1, T):
        eps_s[:, i] = phi_s * eps_s[:, i - 1] + eta_s[:, i]
    y_rep = mu_s + eps_s

    def roll_mean(a, w):
        c = np.cumsum(np.where(np.isfinite(a), a, 0.0), axis=-1)
        n = np.cumsum(np.isfinite(a).astype(float), axis=-1)
        c = np.concatenate([np.zeros(a.shape[:-1] + (1,)), c], axis=-1)
        n = np.concatenate([np.zeros(a.shape[:-1] + (1,)), n], axis=-1)
        s = c[..., w:] - c[..., :-w]
        k = n[..., w:] - n[..., :-w]
        return np.where(k >= 0.7 * w, s / np.maximum(k, 1), np.nan)

    obs_roll = roll_mean(d["kwh"].values, WIN)
    rep_roll = roll_mean(y_rep, WIN)
    p_block = np.nanmean(rep_roll <= obs_roll[None, :], axis=0)
    p_block = np.where(np.isfinite(obs_roll), p_block, np.nan)
    centre = np.arange(len(p_block)) + WIN // 2

    A(f"  BLOCK-LEVEL CHECK: {WIN}-DAY MEAN CONSUMPTION")
    A(f"  Test statistic is the mean of each {WIN}-day window, referred to its")
    A("  posterior predictive distribution under the stationary AR(1)+t process.")
    A("  This is the diagnostic that can see sustained anomalies.")
    okb = np.isfinite(p_block)
    A("")
    A(f"    windows evaluated                : {int(okb.sum())}")
    A(f"    p < 0.01                         : {int((p_block[okb]<0.01).sum())}")
    A(f"    p > 0.99                         : {int((p_block[okb]>0.99).sum())}")
    A("")
    A("  SUSTAINED ANOMALOUS PERIODS")
    A("  Runs of consecutive flagged windows, each expanded to the days it covers.")
    A("  Because the windows overlap, two runs of the same direction whose covered")
    A("  days overlap describe one episode and are merged.")
    flagb = okb & ((p_block < 0.01) | (p_block > 0.99))
    high = p_block > 0.99
    runs, i = [], 0
    while i < len(flagb):
        if flagb[i]:
            j = i
            while j + 1 < len(flagb) and flagb[j + 1] and high[j + 1] == high[i]:
                j += 1
            lo = max(0, centre[i] - WIN // 2)
            hi = min(len(d) - 1, centre[j] + WIN // 2)
            runs.append([lo, hi, bool(high[i]),
                         float(p_block[i : j + 1].max() if high[i] else p_block[i : j + 1].min())])
            i = j + 1
        else:
            i += 1
    # merge same-direction episodes whose covered day ranges touch or overlap
    merged = []
    for r in runs:
        if merged and merged[-1][2] == r[2] and r[0] <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], r[1])
            merged[-1][3] = max(merged[-1][3], r[3]) if r[2] else min(merged[-1][3], r[3])
        else:
            merged.append(list(r))

    rows = []
    A("")
    A(f"    {'from':>12} {'to':>12} {'days':>5} {'obs kWh':>8} {'exp kWh':>8} "
      f"{'ratio':>6} {'p':>8}")
    for lo, hi, is_high, pv in merged:
        sub = d.iloc[lo : hi + 1]
        o, e = sub["kwh"].mean(), sub["fitted"].mean()
        A(f"    {sub['date'].iloc[0].date()!s:>12} {sub['date'].iloc[-1].date()!s:>12} "
          f"{hi-lo+1:5d} {o:8.2f} {e:8.2f} {o/e:6.2f} {pv:8.4f}")
        rows.append((sub["date"].iloc[0], sub["date"].iloc[-1], hi - lo + 1, o, e, o / e, pv))

    # emit the table as LaTeX so the report never transcribes these by hand
    tabdir = os.path.join(ROOT, "results", "tables")
    os.makedirs(tabdir, exist_ok=True)
    # emit the complete tabular: \input of bare rows breaks booktabs rules
    with open(os.path.join(tabdir, "anomaly_periods.tex"), "w") as th:
        th.write("\\begin{tabular}{@{}llrrrrr@{}}\n\\toprule\n")
        th.write("From & To & Days & Observed & Expected & Ratio & $p$ \\\\\n\\midrule\n")
        for a0, a1, nd, o, e, rt, pv in rows:
            bold = rt < 0.6 or rt > 1.3
            w = (lambda s: r"\textbf{" + s + "}") if bold else (lambda s: s)
            th.write(
                f"{w(a0.strftime('%-d %b %Y'))} & {w(a1.strftime('%-d %b %Y'))} & "
                f"{w(str(nd))} & {w(f'{o:.2f}')} & {w(f'{e:.2f}')} & {w(f'{rt:.2f}')} & "
                f"{w('$<$0.001' if pv < 0.001 else f'{min(pv,1-pv):.3f}')} \\\\\n"
            )
        th.write("\\bottomrule\n\\end{tabular}\n")
    A("")
    A(f"  ({len(rows)} episodes; LaTeX version written to results/tables/anomaly_periods.tex)")
    A("")
    A("    These are the load-bearing anomalies: multi-week absences (the French")
    A("    August holiday is the clearest) and unusually cold or mild stretches.")
    A("    They are behavioural and meteorological, and no calendar-based model")
    A("    can predict them, which is exactly why the t likelihood is needed.")
    d["p_block"] = np.nan
    d.loc[centre, "p_block"] = p_block

    d.to_csv(os.path.join(ROOT, "results", "05_daily_with_diagnostics.csv"), index=False)
    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(ROOT, "results", "05_derived.txt"), "w") as fh:
        fh.write(txt + "\n")

    # ================= figures =================
    f7 = load_fit("M7")
    b = M.flatten(f7, "beta")
    nm = f7["names"]
    tt = np.arange(1, 366)
    comp = {}
    for j in (1, 2):
        comp[j] = (b[:, [nm.index(f"sin_h{j}")]] * np.sin(2 * np.pi * j * tt / 365.0)
                   + b[:, [nm.index(f"cos_h{j}")]] * np.cos(2 * np.pi * j * tt / 365.0))
    tot = comp[1] + comp[2]
    cal = [date0 + pd.Timedelta(days=int(x)) for x in tt - 1]

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    for arr, lab, col in [(comp[1], "first harmonic only", "#b03a2e"),
                          (tot, "first + second harmonic", "#1f4e79")]:
        q = np.percentile(arr, [2.5, 50, 97.5], axis=0)
        ax[0].fill_between(cal, q[0], q[2], color=col, alpha=0.18)
        ax[0].plot(cal, q[1], color=col, lw=1.8, label=lab)
    ax[0].axhline(0, color="k", lw=0.7)
    ax[0].set_ylabel("deviation from annual mean (kWh/day)")
    ax[0].set_title("Annual seasonal component, model M7\n(bands are 95% credible)", fontsize=10)
    ax[0].legend(fontsize=8)
    import matplotlib.dates as mdates

    ax[0].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    q2 = np.percentile(comp[2], [2.5, 50, 97.5], axis=0)
    ax[1].fill_between(cal, q2[0], q2[2], color="#2e7d32", alpha=0.20)
    ax[1].plot(cal, q2[1], color="#2e7d32", lw=1.8)
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].set_title("Second harmonic alone: what it adds\n"
                    "(sharpens the winter peak, flattens the summer floor)", fontsize=10)
    ax[1].set_ylabel("kWh/day")
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "05_seasonal.png"), dpi=140)
    plt.close(fig)

    dow_cols = [nm.index(f"dow_{x}") for x in DAYNAMES[1:]]
    dow_full = np.column_stack([np.zeros(len(b)), b[:, dow_cols]])
    prof = dow_full - dow_full.mean(axis=1, keepdims=True)
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    q = np.percentile(prof, [2.5, 25, 50, 75, 97.5], axis=0)
    ax[0].errorbar(range(7), q[2], yerr=[q[2] - q[0], q[4] - q[2]], fmt="o",
                   color="#1f4e79", capsize=4, label="95% CI")
    ax[0].errorbar(range(7), q[2], yerr=[q[2] - q[1], q[3] - q[2]], fmt="o",
                   color="#1f4e79", lw=3, capsize=0, label="50% CI")
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_xticks(range(7))
    ax[0].set_xticklabels(DAYNAMES)
    ax[0].set_ylabel("kWh/day vs weekly mean")
    ax[0].set_title("Day-of-week effect, model M7", fontsize=10)
    ax[0].legend(fontsize=8)
    wkend = prof[:, 5:].mean(axis=1) - prof[:, :5].mean(axis=1)
    ax[1].hist(wkend, bins=60, color="#1f4e79")
    ax[1].axvline(0, color="k", lw=1)
    ax[1].set_xlabel("weekend minus weekday (kWh/day)")
    ax[1].set_title(f"Posterior of the weekend effect\nP(> 0) = {np.mean(wkend>0):.4f}", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "05_weekly.png"), dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    ax[0].plot(d["date"], d["kwh"], lw=0.5, color="#888", label="observed")
    ax[0].plot(d["date"], d["fitted"], lw=1.0, color="#1f4e79", label="fitted systematic part")
    fl = d[np.isfinite(d["ppp"]) & ((d["ppp"] < 0.01) | (d["ppp"] > 0.99))]
    ax[0].plot(fl["date"], fl["kwh"], "o", ms=4, mfc="none", color="#b03a2e",
               label="flagged, $p_t<0.01$ or $>0.99$")
    ax[0].set_ylabel("kWh/day")
    ax[0].legend(fontsize=8, ncol=3)
    ax[0].set_title("Posterior predictive check: which days the model cannot explain", fontsize=10)
    ax[1].plot(d["date"], d["lam"], lw=0.7, color="#2e7d32")
    ax[1].axhline(0.5, color="#b03a2e", ls="--", lw=1, label=r"$\lambda_t=0.5$ downweighting threshold")
    ax[1].set_ylabel(r"posterior mean $\lambda_t$")
    ax[1].set_xlabel("date")
    ax[1].legend(fontsize=8)
    ax[1].set_title(r"Heavy-tail weights: low $\lambda_t$ marks days treated as outliers", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "05_anomalies.png"), dpi=140)
    plt.close(fig)

    print("\nwrote results/05_derived.txt, 05_daily_with_diagnostics.csv, 3 figures")


if __name__ == "__main__":
    main()
