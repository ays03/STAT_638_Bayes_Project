"""
Step 1: aggregate the UCI one-minute household power record to a daily kWh series.

Global_active_power is average real power (kW) over each one-minute interval, so the
energy used in that minute is GAP/60 kWh and the daily total is

    y_d = (1/60) * sum_{m in day d} GAP_m       [kWh]

A day is retained only if at least MIN_COVERAGE of its 1440 minutes are observed;
partial days are recorded as missing rather than summed, because summing over gaps
manufactures artificially low days that would later be flagged as anomalies.

Outputs
    data/daily.csv          date, t, kwh, dow, n_obs, coverage, observed
    figures/01_*.png        raw series, coverage histogram, missingness map
    results/01_data_summary.txt
"""
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(ROOT, "data", "household_power_consumption.txt")
MIN_COVERAGE = 0.95
MINUTES_PER_DAY = 1440

for sub in ("figures", "results"):
    os.makedirs(os.path.join(ROOT, sub), exist_ok=True)


def load_minutes():
    """Read the raw semicolon-delimited file. '?' is the missing-value code."""
    df = pd.read_csv(
        RAW,
        sep=";",
        na_values="?",
        low_memory=False,
        usecols=[
            "Date",
            "Time",
            "Global_active_power",
            "Sub_metering_1",
            "Sub_metering_2",
            "Sub_metering_3",
        ],
        dtype={
            "Global_active_power": "float64",
            "Sub_metering_1": "float64",
            "Sub_metering_2": "float64",
            "Sub_metering_3": "float64",
        },
    )
    # dayfirst: the file uses dd/mm/yyyy
    df["date"] = pd.to_datetime(df["Date"], format="%d/%m/%Y")
    return df


def aggregate(df):
    g = df.groupby("date")
    daily = pd.DataFrame(
        {
            "n_obs": g["Global_active_power"].count(),
            "kwh_raw": g["Global_active_power"].sum() / 60.0,
            "sub_kwh": (
                g["Sub_metering_1"].sum() + g["Sub_metering_2"].sum() + g["Sub_metering_3"].sum()
            )
            / 1000.0,
            "n_rows": g["Global_active_power"].size(),
        }
    )
    daily["coverage"] = daily["n_obs"] / MINUTES_PER_DAY
    return daily.reset_index()


def main():
    print("reading raw file ...")
    minutes = load_minutes()
    n_rows = len(minutes)
    n_missing = int(minutes["Global_active_power"].isna().sum())

    daily = aggregate(minutes)

    # Drop the two partial boundary days outright: the record begins 17:24 on
    # 2006-12-16 and ends 21:02 on 2010-11-26.
    first, last = daily["date"].iloc[0], daily["date"].iloc[-1]
    daily = daily.iloc[1:-1].reset_index(drop=True)

    # Re-index onto a gap-free calendar so the AR(1) error process is defined on
    # every day, whether or not the outcome is observed.
    full = pd.DataFrame(
        {"date": pd.date_range(daily["date"].min(), daily["date"].max(), freq="D")}
    ).merge(daily, on="date", how="left")
    full["n_obs"] = full["n_obs"].fillna(0).astype(int)
    full["coverage"] = full["coverage"].fillna(0.0)

    full["observed"] = full["coverage"] >= MIN_COVERAGE
    full["kwh"] = np.where(full["observed"], full["kwh_raw"], np.nan)
    full["t"] = np.arange(1, len(full) + 1)
    full["dow"] = full["date"].dt.dayofweek  # 0 = Monday
    full["year"] = full["date"].dt.year
    full["doy"] = full["date"].dt.dayofyear

    # Sanity check: sub-metered energy must not exceed the whole-house total.
    obs = full[full["observed"]]
    bad = int((obs["sub_kwh"] > obs["kwh"]).sum())

    out = full[
        ["date", "t", "kwh", "dow", "n_obs", "coverage", "observed", "sub_kwh", "year", "doy"]
    ]
    out.to_csv(os.path.join(ROOT, "data", "daily.csv"), index=False)

    lines = []
    A = lines.append
    A("RAW FILE")
    A(f"  rows                         {n_rows:,}")
    A(f"  missing Global_active_power  {n_missing:,} ({100*n_missing/n_rows:.3f}%)")
    A(f"  first / last timestamped day {first.date()} / {last.date()}  (both dropped as partial)")
    A("")
    A("DAILY SERIES")
    A(f"  calendar days on grid        {len(full)}")
    A(f"  days with coverage >= {MIN_COVERAGE:.2f}   {int(full['observed'].sum())}")
    A(f"  days treated as missing      {int((~full['observed']).sum())}")
    A(f"  of which zero observed mins  {int((full['n_obs'] == 0).sum())}")
    A(f"  sub-meter > total violations {bad}  (expect 0)")
    A("")
    A("DAILY kWh (observed days only)")
    d = obs["kwh"]
    A(f"  n        {len(d)}")
    A(f"  mean     {d.mean():.2f}")
    A(f"  sd       {d.std(ddof=1):.2f}")
    A(f"  min      {d.min():.2f}  on {obs.loc[d.idxmin(), 'date'].date()}")
    A(f"  q25      {d.quantile(.25):.2f}")
    A(f"  median   {d.median():.2f}")
    A(f"  q75      {d.quantile(.75):.2f}")
    A(f"  max      {d.max():.2f}  on {obs.loc[d.idxmax(), 'date'].date()}")
    A(f"  skewness {d.skew():.2f}   (log scale: {np.log(d).skew():.2f})")
    A("")
    A("ANNUALISED CONSUMPTION  (mean daily kWh x 365, by calendar year)")
    for yr, sub in obs.groupby("year"):
        A(f"  {yr}  n={len(sub):4d}  mean={sub['kwh'].mean():5.2f} kWh/day  -> {365*sub['kwh'].mean():,.0f} kWh/yr")
    A("")
    A("MEAN DAILY kWh BY DAY OF WEEK  (0=Mon)")
    for k, sub in obs.groupby("dow"):
        A(f"  {k}  n={len(sub):4d}  mean={sub['kwh'].mean():5.2f}")
    A("")
    A("MEAN DAILY kWh BY MONTH")
    for m, sub in obs.groupby(obs["date"].dt.month):
        A(f"  {m:2d}  n={len(sub):4d}  mean={sub['kwh'].mean():5.2f}")
    A("")
    A("LAG-1 AUTOCORRELATION OF log(kwh) AFTER REMOVING MONTH + DOW MEANS")
    z = np.log(obs["kwh"].values)
    dmat = pd.get_dummies(obs["date"].dt.month.astype(str) + "_m").values.astype(float)
    wmat = pd.get_dummies(obs["dow"].astype(str) + "_d").values.astype(float)
    X = np.column_stack([dmat, wmat[:, 1:]])
    resid = z - X @ np.linalg.lstsq(X, z, rcond=None)[0]
    tt = obs["t"].values
    adj = np.diff(tt) == 1
    A(f"  rho1 = {np.corrcoef(resid[:-1][adj], resid[1:][adj])[0,1]:.3f}  (consecutive observed days)")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(ROOT, "results", "01_data_summary.txt"), "w") as f:
        f.write(txt + "\n")

    # ---------------- figures ----------------
    fig, ax = plt.subplots(3, 1, figsize=(11, 9))
    ax[0].plot(full["date"], full["kwh"], lw=0.5, color="#1f4e79")
    ax[0].set_title("Daily household electricity consumption, Dec 2006 - Nov 2010")
    ax[0].set_ylabel("kWh / day")
    ax[1].plot(full["date"], np.log(full["kwh"]), lw=0.5, color="#1f4e79")
    ax[1].set_title("log daily consumption (modelling scale)")
    ax[1].set_ylabel("log kWh")
    ax[2].vlines(full["date"], 0, full["coverage"], lw=0.6, color="#b03a2e")
    ax[2].axhline(MIN_COVERAGE, color="k", ls="--", lw=1, label=f"retention threshold {MIN_COVERAGE}")
    ax[2].set_ylabel("fraction of 1440 min observed")
    ax[2].set_title("Daily data coverage")
    ax[2].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01_series_and_coverage.png"), dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].hist(full["coverage"], bins=60, color="#1f4e79")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("daily coverage")
    ax[0].set_title("Coverage histogram (log count)")
    ax[1].hist(obs["kwh"], bins=45, color="#1f4e79")
    ax[1].set_xlabel("kWh / day")
    ax[1].set_title(f"Daily kWh, right-skewed (skew {obs['kwh'].skew():.2f})")
    ax[2].hist(np.log(obs["kwh"]), bins=45, color="#2e7d32")
    ax[2].set_xlabel("log kWh / day")
    ax[2].set_title(f"log kWh, near-symmetric (skew {np.log(obs['kwh']).skew():.2f})")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01_distributions.png"), dpi=140)
    plt.close(fig)

    # seasonal / weekly exploratory views
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for yr, sub in obs.groupby("year"):
        ax[0].plot(sub["doy"], sub["kwh"], ".", ms=2.5, alpha=0.55, label=str(yr))
    ax[0].set_xlabel("day of year")
    ax[0].set_ylabel("kWh / day")
    ax[0].set_title("Annual cycle: consumption vs day of year")
    ax[0].legend(fontsize=8, markerscale=3)
    names = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    ax[1].boxplot([obs.loc[obs["dow"] == k, "kwh"].values for k in range(7)], labels=names)
    ax[1].set_ylabel("kWh / day")
    ax[1].set_title("Weekly pattern (raw)")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01_exploratory_season_week.png"), dpi=140)
    plt.close(fig)

    print("\nwrote data/daily.csv, results/01_data_summary.txt, 3 figures")


if __name__ == "__main__":
    main()
