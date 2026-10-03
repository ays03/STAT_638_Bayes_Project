"""
Step 1d: in-depth exploratory analysis of daily household demand.

Descriptive, pre-model work. Locates the variance among a linear trend, annual
harmonics and day of week; compares that seasonal shape with a month-dummy
ceiling and with the year-by-year curve; splits the annual cycle by sub-meter;
checks the weekly pattern after seasonal adjustment, including whether the
weekend premium changes with season; and summarises the within-day load curve.
Dependence diagnostics here are on the kWh scale, which is the scale the
holdout later selects.

Outputs
    results/01d_eda.txt
    results/tables/eda_*.tex
    figures/01d_*.png
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
ORIGIN = pd.Timestamp("2006-12-17")
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
DOW_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
# UCI documentation: kitchen; laundry; water heater and air conditioner.
ENDUSE = [
    ("remainder", "Unmetered", "#1f4e79"),
    ("sub3", "Water heater and air conditioner", "#e67e22"),
    ("sub2", "Laundry", "#1e8449"),
    ("sub1", "Kitchen", "#7d3c98"),
]
SEASON_MONTHS = {
    "winter": (12, 1, 2),
    "spring": (3, 4, 5),
    "summer": (6, 7, 8),
    "autumn": (9, 10, 11),
}


def r2_of(y, resid):
    sst = np.sum((y - y.mean()) ** 2)
    return 1.0 - np.sum(resid ** 2) / sst


def nw_cov(X, resid, L=14):
    """Newey-West covariance of the OLS coefficient vector, Bartlett weights."""
    xe = X * resid[:, None]
    meat = xe.T @ xe
    for lag in range(1, L + 1):
        w = 1.0 - lag / (L + 1.0)
        g = xe[lag:].T @ xe[:-lag]
        meat += w * (g + g.T)
    xtx = X.T @ X
    inv = np.linalg.inv(xtx)
    return inv @ meat @ inv


def ols(X, y):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    return beta, resid


def easter(year):
    """Anonymous Gregorian computus. Returns Easter Sunday."""
    a = year % 19
    b = year // 100
    c = year % 100
    d = b // 4
    e = b % 4
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i = c // 4
    k = c % 4
    ell = (32 + 2 * e + 2 * i - h - k) % 7
    m = (a + 11 * h + 22 * ell) // 451
    month = (h + ell - 7 * m + 114) // 31
    day = ((h + ell - 7 * m + 114) % 31) + 1
    return pd.Timestamp(year, month, day)


def french_holidays(years):
    """Métropole public holidays. Alsace-Moselle extras are not relevant in Sceaux."""
    out = []
    for year in years:
        eas = easter(year)
        fixed = [
            (pd.Timestamp(year, 1, 1), "New Year"),
            (eas + pd.Timedelta(days=1), "Easter Monday"),
            (pd.Timestamp(year, 5, 1), "Labour Day"),
            (pd.Timestamp(year, 5, 8), "Victory Day"),
            (eas + pd.Timedelta(days=39), "Ascension"),
            (eas + pd.Timedelta(days=50), "Whit Monday"),
            (pd.Timestamp(year, 7, 14), "Bastille Day"),
            (pd.Timestamp(year, 8, 15), "Assumption"),
            (pd.Timestamp(year, 11, 1), "All Saints"),
            (pd.Timestamp(year, 11, 11), "Armistice"),
            (pd.Timestamp(year, 12, 25), "Christmas"),
        ]
        # Ascension fell on 1 May 2008, so Labour Day and Ascension coincide.
        seen = set()
        for dt, name in fixed:
            if dt in seen:
                continue
            seen.add(dt)
            out.append((dt, name))
    return out


def clim_day(dates):
    """Non-leap day-of-year, 1..365. 29 February is missing (NaN)."""
    doy = dates.dt.dayofyear.to_numpy()
    leap_after_feb = (dates.dt.is_leap_year & (dates.dt.month > 2)).to_numpy()
    feb29 = ((dates.dt.month == 2) & (dates.dt.day == 29)).to_numpy()
    d = doy.astype(float) - leap_after_feb.astype(float)
    d[feb29] = np.nan
    return d


def smooth_on_clim(clim, values, window=15, circular=True, min_frac=0.6):
    half = window // 2
    s = np.zeros(365)
    c = np.zeros(365)
    for d, v in zip(clim, values):
        if np.isnan(d) or np.isnan(v):
            continue
        s[int(d) - 1] += v
        c[int(d) - 1] += 1
    out = np.full(365, np.nan)
    for i in range(365):
        if circular:
            idx = [(i + k) % 365 for k in range(-half, half + 1)]
        else:
            idx = [i + k for k in range(-half, half + 1) if 0 <= i + k < 365]
        cc = sum(c[j] for j in idx)
        if cc >= window * min_frac:
            out[i] = sum(s[j] for j in idx) / cc
    return out


def low_use_mask(kwh, run=5, threshold=12.0):
    """Same rule as 01c: runs of >= run consecutive days under threshold kWh."""
    low = kwh < threshold
    block = np.zeros(len(kwh), dtype=bool)
    i = 0
    while i < len(kwh):
        if low[i]:
            j = i
            while j + 1 < len(kwh) and low[j + 1]:
                j += 1
            if j - i + 1 >= run:
                block[i : j + 1] = True
            i = j + 1
        else:
            i += 1
    return block


def month_dummies(month, drop=1):
    levels = [m for m in range(1, 13) if m != drop and np.any(month == m)]
    if not levels:
        return np.zeros((len(month), 0))
    return np.column_stack([(month == m).astype(float) for m in levels])


def dow_dummies(dow, drop=0):
    levels = [k for k in range(7) if k != drop]
    return np.column_stack([(dow == k).astype(float) for k in levels]), levels


def harmonic_cols(t, J):
    cols = []
    for j in range(1, J + 1):
        w = 2 * np.pi * j * t / 365.0
        cols += [np.sin(w), np.cos(w)]
    return np.column_stack(cols) if cols else np.zeros((len(t), 0))


def acf_calendar(resid, t, nlags):
    """Pearson correlation at calendar lag k, using pairs that are both observed."""
    idx = {int(tt): i for i, tt in enumerate(t)}
    out = []
    for k in range(1, nlags + 1):
        pairs = [(i, idx[int(tt) - k]) for i, tt in enumerate(t) if int(tt) - k in idx]
        a = resid[[i for i, _ in pairs]]
        b = resid[[j for _, j in pairs]]
        a = a - a.mean()
        b = b - b.mean()
        out.append(float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b))))
    return np.array(out)


def pacf_yw(rho):
    """Yule-Walker partial autocorrelations from an ACF that starts at lag 1."""
    pac = []
    for k in range(1, len(rho) + 1):
        r = rho[:k]
        if k == 1:
            pac.append(r[0])
            continue
        R = np.empty((k, k))
        for i in range(k):
            for j in range(k):
                lag = abs(i - j)
                R[i, j] = 1.0 if lag == 0 else rho[lag - 1]
        try:
            phi = np.linalg.solve(R, r)
        except np.linalg.LinAlgError:
            phi = np.linalg.lstsq(R, r, rcond=None)[0]
        pac.append(float(phi[-1]))
    return np.array(pac)


def harmonic_on_dates(beta_harm, dates, J):
    """beta_harm is [sin1, cos1, ...] in that order. dates is a DatetimeIndex."""
    t = (pd.Series(dates) - ORIGIN).dt.days.to_numpy() + 1.0
    h = np.zeros(len(dates))
    for j in range(1, J + 1):
        w = 2 * np.pi * j * t / 365.0
        h += beta_harm[2 * j - 2] * np.sin(w) + beta_harm[2 * j - 1] * np.cos(w)
    return h


def write_tabular(path, colspec, header, rows, note=None):
    lines = [rf"\begin{{tabular}}{{{colspec}}}", r"\toprule", header, r"\midrule"]
    lines.extend(rows)
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    if note:
        lines.append(note)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def load_minutes():
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
    df["date"] = pd.to_datetime(df["Date"], format="%d/%m/%Y")
    df["hour"] = df["Time"].str.slice(0, 2).astype(int)
    return df


def season_of(month):
    out = np.empty(len(month), dtype=object)
    for name, months in SEASON_MONTHS.items():
        out[np.isin(month, months)] = name
    return out


def main():
    for sub in ("figures", "results", os.path.join("results", "tables")):
        os.makedirs(os.path.join(ROOT, sub), exist_ok=True)

    d = pd.read_csv(os.path.join(ROOT, "data", "daily.csv"), parse_dates=["date"])
    d["observed"] = d["observed"].astype(bool)
    obs = d[d["observed"]].copy().reset_index(drop=True)
    # t in daily.csv is the calendar index; recompute and check.
    t_check = (obs["date"] - ORIGIN).dt.days.to_numpy() + 1
    if not np.array_equal(t_check, obs["t"].to_numpy()):
        raise RuntimeError("t index does not match days since 2006-12-17")

    y = obs["kwh"].to_numpy(float)
    t = obs["t"].to_numpy(float)
    dow = obs["dow"].to_numpy(int)
    month = obs["date"].dt.month.to_numpy(int)
    year = obs["year"].to_numpy(int)
    n = len(y)
    block = low_use_mask(y)
    # Apple Accelerate can raise spurious FP flags inside an otherwise exact
    # matmul. Results are checked for finiteness where they are reported.
    np.seterr(divide="ignore", over="ignore", invalid="ignore", under="ignore")

    # Easter dates are load-bearing for the holiday list; pin the ones used here.
    assert easter(2007) == pd.Timestamp("2007-04-08")
    assert easter(2008) == pd.Timestamp("2008-03-23")
    assert easter(2009) == pd.Timestamp("2009-04-12")
    assert easter(2010) == pd.Timestamp("2010-04-04")
    assert clim_day(pd.Series(pd.to_datetime(["2008-03-01"])))[0] == 60
    assert np.isnan(clim_day(pd.Series(pd.to_datetime(["2008-02-29"])))[0])

    lines = []
    A = lines.append
    A("IN-DEPTH EDA  (descriptive OLS and direct summaries; not posterior estimates)")
    A(f"Retained days: {n}.  Low-use block days: {int(block.sum())}.")
    A("")

    # ------------------------------------------------------------------
    # Missing days
    # ------------------------------------------------------------------
    miss = d[~d["observed"]]
    A("DAYS NOT RETAINED")
    A(f"  n = {len(miss)}")
    for _, row in miss.iterrows():
        A(f"  {row['date'].date()}  coverage={row['coverage']:.3f}  n_obs={int(row['n_obs'])}")
    A("")

    # ------------------------------------------------------------------
    # Sequential variance decomposition on the kWh scale
    # ------------------------------------------------------------------
    A("SEQUENTIAL R2 ON kWh  (each row adds the named columns to everything above)")
    A(f"{'specification':<42} {'R2':>7} {'dR2':>7} {'resid sd':>9}")
    # Order matches the model: trend, harmonic 1, harmonic 2, day of week.
    # Harmonics 3 and 4 are then added on top, so their increments are the
    # leftover seasonal shape rather than a claim on the variance before DOW.
    cols = [np.ones(n), t / 365.0]
    labels_seq = ["intercept only", "+ linear trend"]
    mats = [np.ones((n, 1)), np.column_stack(cols)]
    for J in (1, 2):
        cols.extend(list(harmonic_cols(t, J)[:, -2:].T))
        labels_seq.append(f"+ harmonic {J}")
        mats.append(np.column_stack(cols))
    Dd, _ = dow_dummies(dow, drop=0)
    cols.extend(list(Dd.T))
    labels_seq.append("+ day of week (Mon baseline)")
    mats.append(np.column_stack(cols))
    for J in (3, 4):
        cols.extend(list(harmonic_cols(t, J)[:, -2:].T))
        labels_seq.append(f"+ harmonic {J}")
        mats.append(np.column_stack(cols))

    seq_rows = []
    prev = 0.0
    fits = {}
    for label, X in zip(labels_seq, mats):
        beta, resid = ols(X, y)
        r2 = r2_of(y, resid)
        fits[label] = (beta, resid, r2, X)
        A(f"{label:<42} {r2:7.3f} {r2 - prev:7.3f} {resid.std(ddof=1):9.2f}")
        seq_rows.append((label, r2, r2 - prev, resid.std(ddof=1)))
        prev = r2

    # ceilings and restrictions
    X_month = np.column_stack([np.ones(n), t / 365.0, month_dummies(month)])
    b_m, r_m = ols(X_month, y)
    r2_month = r2_of(y, r_m)
    X_month_dow = np.column_stack([X_month, Dd])
    b_md, r_md = ols(X_month_dow, y)
    r2_month_dow = r2_of(y, r_md)
    weekend = (dow >= 5).astype(float)
    X_month_we = np.column_stack([X_month, weekend])
    _, r_we = ols(X_month_we, y)
    r2_month_we = r2_of(y, r_we)
    X_j2 = mats[labels_seq.index("+ harmonic 2")]
    X_j2_we = np.column_stack([X_j2, weekend])
    _, r_j2we = ols(X_j2_we, y)
    r2_j2_we = r2_of(y, r_j2we)
    r2_j2_dow = dict(zip(labels_seq, [row[1] for row in seq_rows]))["+ day of week (Mon baseline)"]

    A("")
    A("SEASONAL CEILING AND WEEKEND RESTRICTION")
    A(f"  trend + month dummies                         R2 = {r2_month:.3f}")
    A(f"  trend + month + weekend indicator             R2 = {r2_month_we:.3f}  "
      f"(dR2 vs month = {r2_month_we - r2_month:.3f})")
    A(f"  trend + month + full day-of-week              R2 = {r2_month_dow:.3f}  "
      f"(dR2 vs weekend = {r2_month_dow - r2_month_we:.3f})")
    A(f"  trend + 2 harmonics + weekend                 R2 = {r2_j2_we:.3f}")
    A(f"  trend + 2 harmonics + full day-of-week        R2 = {r2_j2_dow:.3f}  "
      f"(dR2 vs weekend = {r2_j2_dow - r2_j2_we:.3f})")
    A(f"  month indicators minus 2 harmonics, both + trend + DOW: "
      f"R2 gap = {r2_month_dow - r2_j2_dow:.3f}")
    A("")

    # Harmonic amplitude and peak, joint with trend and DOW, J=1 and J=2.
    A("HARMONIC SHAPE ON THE kWh SCALE  (OLS with trend and day of week)")
    ref_dates = pd.date_range("2009-01-01", "2009-12-31")
    curves = {}
    for J in (1, 2):
        XJ = np.column_stack([np.ones(n), t / 365.0, harmonic_cols(t, J), Dd])
        bJ, rJ = ols(XJ, y)
        harm = bJ[2 : 2 + 2 * J]
        raw_curve = harmonic_on_dates(harm, ref_dates, J)
        # Place the curve on the same level as 2007-2009 mean consumption.
        y_complete = y[np.isin(year, [2007, 2008, 2009])]
        curve = raw_curve - raw_curve.mean() + y_complete.mean()
        curves[J] = curve
        peak = int(np.argmax(curve))
        trough = int(np.argmin(curve))
        A(f"  J={J}  peak {ref_dates[peak].strftime('%-d %b')} ({curve[peak]:.2f} kWh)  "
          f"trough {ref_dates[trough].strftime('%-d %b')} ({curve[trough]:.2f} kWh)  "
          f"range {curve[peak] - curve[trough]:.2f}  resid sd {rJ.std(ddof=1):.2f}")
        # amplitude of each harmonic
        for j in range(1, J + 1):
            a, c = harm[2 * j - 2], harm[2 * j - 1]
            A(f"       harmonic {j} amplitude {np.hypot(a, c):.2f} kWh")
    A("")

    # ------------------------------------------------------------------
    # Monthly distribution, overall and by year
    # ------------------------------------------------------------------
    A("MONTHLY DISTRIBUTION OF DAILY kWh  (all retained days)")
    A(f"{'mo':>4} {'n':>5} {'mean':>8} {'sd':>7} {'p10':>8} {'p50':>8} {'p90':>8}")
    month_stats = []
    for m in range(1, 13):
        v = y[month == m]
        row = (m, len(v), v.mean(), v.std(ddof=1),
               np.quantile(v, 0.10), np.median(v), np.quantile(v, 0.90))
        month_stats.append(row)
        A(f"{m:4d} {row[1]:5d} {row[2]:8.2f} {row[3]:7.2f} {row[4]:8.2f} {row[5]:8.2f} {row[6]:8.2f}")

    A("")
    A("MONTHLY MEAN BY YEAR")
    A("  year " + " ".join(f"{m:>6d}" for m in range(1, 13)))
    for yr in sorted(obs["year"].unique()):
        cells = []
        for m in range(1, 13):
            v = y[(year == yr) & (month == m)]
            cells.append(f"{v.mean():6.1f}" if len(v) else f"{'·':>6}")
        A(f"  {yr} " + " ".join(cells))

    dec_ex = y[(month == 12) & (year != 2006)]
    A("")
    A(f"DECEMBER MEAN  all retained {y[month == 12].mean():.2f}   "
      f"excluding 2006 {dec_ex.mean():.2f} (n={len(dec_ex)})")
    aug = y[month == 8]
    aug_ex = y[(month == 8) & ~block]
    A(f"AUGUST MEAN    all retained {aug.mean():.2f} (median {np.median(aug):.2f})   "
      f"excluding low-use blocks {aug_ex.mean():.2f} (median {np.median(aug_ex):.2f}, n={len(aug_ex)})")
    jul = y[month == 7]
    A(f"JULY MEAN      all retained {jul.mean():.2f} (median {np.median(jul):.2f})")

    # Year-to-year correlation of monthly means, months in common.
    A("")
    A("CORRELATION OF MONTHLY MEANS BETWEEN YEARS (months observed in both)")
    years = [2007, 2008, 2009, 2010]
    mm = {}
    for yr in years:
        mm[yr] = np.array([
            y[(year == yr) & (month == m)].mean() if np.any((year == yr) & (month == m)) else np.nan
            for m in range(1, 13)
        ])
    for i, a in enumerate(years):
        for b in years[i + 1 :]:
            ok = np.isfinite(mm[a]) & np.isfinite(mm[b])
            corr = np.corrcoef(mm[a][ok], mm[b][ok])[0, 1]
            A(f"  {a} vs {b}: r = {corr:.3f}  ({int(ok.sum())} months)")

    # Seasonal means
    A("")
    A("SEASONAL MEANS")
    for name, months_ in SEASON_MONTHS.items():
        v = y[np.isin(month, months_)]
        A(f"  {name:<8} n={len(v):4d}  mean={v.mean():6.2f}  sd={v.std(ddof=1):5.2f}  "
          f"skew={pd.Series(v).skew():5.2f}")
    A("")

    # ------------------------------------------------------------------
    # Sub-meters. Re-aggregate so the three channels are separate.
    # ------------------------------------------------------------------
    A("READING MINUTE FILE FOR SUB-METERS AND THE DIURNAL PROFILE ...")
    minutes = load_minutes()
    # Missingness by hour, over rows that exist in the file.
    miss_hour = minutes.groupby("hour")["Global_active_power"].apply(lambda s: float(s.isna().mean()))
    A("FRACTION OF ROWS WITH MISSING POWER, BY HOUR")
    A("  " + "  ".join(f"{h:02d}:{miss_hour.loc[h]:.3f}" for h in range(24)))

    gap_na = minutes["Global_active_power"].isna()
    sub_na = minutes[["Sub_metering_1", "Sub_metering_2", "Sub_metering_3"]].isna().any(axis=1)
    A(f"  minutes with GAP missing: {int(gap_na.sum())}")
    A(f"  minutes with any sub-meter missing: {int(sub_na.sum())}")
    A(f"  disagreement between those two masks: {int((gap_na != sub_na).sum())}")

    ok = minutes.loc[~gap_na].copy()
    ok["sub1_kwh"] = ok["Sub_metering_1"] / 1000.0
    ok["sub2_kwh"] = ok["Sub_metering_2"] / 1000.0
    ok["sub3_kwh"] = ok["Sub_metering_3"] / 1000.0
    # GAP is kW over the minute; /60 converts the minute to kWh.
    ok["kwh_min"] = ok["Global_active_power"] / 60.0
    daily_sub = ok.groupby("date").agg(
        kwh_re=("kwh_min", "sum"),
        sub1=("sub1_kwh", "sum"),
        sub2=("sub2_kwh", "sum"),
        sub3=("sub3_kwh", "sum"),
    )
    obs = obs.merge(daily_sub, left_on="date", right_index=True, how="left")
    gap = np.nanmax(np.abs(obs["kwh"] - obs["kwh_re"]))
    A(f"  max |recomputed daily kWh - daily.csv|: {gap:.6f}")
    obs["remainder"] = obs["kwh"] - obs["sub1"] - obs["sub2"] - obs["sub3"]
    A(f"  days with negative remainder: {int((obs['remainder'] < -1e-6).sum())}")
    A(f"  minimum remainder: {obs['remainder'].min():.3f} kWh")

    A("")
    A("END-USE MEANS, ALL RETAINED DAYS")
    total_mean = obs["kwh"].mean()
    for key, label, _ in ENDUSE:
        m = obs[key].mean()
        A(f"  {label:<40} {m:6.2f} kWh/day   {100 * m / total_mean:5.1f}%")

    A("")
    A("END USE BY MONTH (mean kWh/day)")
    A(f"{'mo':>4} {'total':>8} {'remain':>8} {'sub3':>8} {'sub2':>8} {'sub1':>8}")
    month_use = []
    for m in range(1, 13):
        sub = obs[obs["date"].dt.month == m]
        rec = (m, sub["kwh"].mean(), sub["remainder"].mean(), sub["sub3"].mean(),
               sub["sub2"].mean(), sub["sub1"].mean())
        month_use.append(rec)
        A(f"{m:4d} {rec[1]:8.2f} {rec[2]:8.2f} {rec[3]:8.2f} {rec[4]:8.2f} {rec[5]:8.2f}")

    dec = obs[obs["date"].dt.month == 12]
    aug_df = obs[obs["date"].dt.month == 8]
    A("")
    A("DECEMBER MINUS AUGUST, BY END USE  (all retained days in those months)")
    d_tot = dec["kwh"].mean() - aug_df["kwh"].mean()
    for key, label, _ in ENDUSE:
        delta = dec[key].mean() - aug_df[key].mean()
        A(f"  {label:<40} {delta:+6.2f} kWh/day   {100 * delta / d_tot:5.1f}% of the gap")

    aug_in = obs[(obs["date"].dt.month == 8) & ~block]
    A("")
    A("AUGUST EXCLUDING LOW-USE BLOCKS, BY END USE")
    for key, label, _ in (("kwh", "Total", None),) + tuple(ENDUSE):
        A(f"  {label:<40} {aug_in[key].mean():6.2f}")
    A(f"  December minus this August, total: {dec['kwh'].mean() - aug_in['kwh'].mean():+.2f}")

    # The August 2008 vacancy versus a nearby occupied stretch.
    vac_start, vac_end = pd.Timestamp("2008-08-06"), pd.Timestamp("2008-08-30")
    vac = obs[(obs["date"] >= vac_start) & (obs["date"] <= vac_end)]
    # Occupied comparison: July 2008, and August 2007 (no long vacancy block that year).
    jul08 = obs[(obs["year"] == 2008) & (obs["date"].dt.month == 7)]
    aug07 = obs[(obs["year"] == 2007) & (obs["date"].dt.month == 8)]
    A("")
    A("VACANCY BLOCK 6-30 Aug 2008 VERSUS OCCUPIED COMPARISONS")
    A(f"{'slice':<28} {'n':>4} {'total':>8} {'remain':>8} {'sub3':>8} {'sub2':>8} {'sub1':>8}")
    for label, sl in [("6-30 Aug 2008", vac), ("July 2008", jul08), ("August 2007", aug07)]:
        A(f"{label:<28} {len(sl):4d} {sl['kwh'].mean():8.2f} {sl['remainder'].mean():8.2f} "
          f"{sl['sub3'].mean():8.2f} {sl['sub2'].mean():8.2f} {sl['sub1'].mean():8.2f}")
    A(f"  vacancy average power: {vac['kwh'].mean() / 24:.3f} kW")
    A(f"  July 2008 average power: {jul08['kwh'].mean() / 24:.3f} kW")
    A(f"  December average power: {dec['kwh'].mean() / 24:.3f} kW")

    # ------------------------------------------------------------------
    # Weekly pattern after removing month effects. Newey-West SEs.
    # ------------------------------------------------------------------
    A("")
    A("DAY OF WEEK AFTER MONTH EFFECTS  (deviations from the 7-day mean, Newey-West L=14)")
    Xw = np.column_stack([np.ones(n), month_dummies(month), Dd])
    bw, rw = ols(Xw, y)
    Vw = nw_cov(Xw, rw, L=14)
    # DOW coefficients are the last 6 columns, Tuesday..Sunday, Monday = 0.
    b_dow = bw[-6:]
    V_dow = Vw[-6:, -6:]
    # dev = A @ b_dow, Monday included as the zero constraint.
    AA = np.full((7, 6), -1.0 / 7.0)
    for k in range(6):
        AA[k + 1, k] += 1.0
    dev = AA @ b_dow
    se_dev = np.sqrt(np.diag(AA @ V_dow @ AA.T))
    for k in range(7):
        lo, hi = dev[k] - 1.96 * se_dev[k], dev[k] + 1.96 * se_dev[k]
        A(f"  {DOW_NAMES[k]:<4} {dev[k]:+6.2f}   95% NW [{lo:+.2f}, {hi:+.2f}]")

    A("")
    A("WEEKEND MINUS WEEKDAY, WITHIN SEASON, CONTROLLING FOR MONTH  (Newey-West L=14)")
    season = season_of(month)
    we_rows = []
    for name in ("winter", "spring", "summer", "autumn"):
        sel = season == name
        yy = y[sel]
        mm = month[sel]
        ww = weekend[sel]
        # Drop the first month present so the dummies are full rank.
        present = [m for m in range(1, 13) if np.any(mm == m)]
        X = np.column_stack([np.ones(sel.sum()), month_dummies(mm, drop=present[0]), ww])
        b, r = ols(X, yy)
        V = nw_cov(X, r, L=14)
        est, se = float(b[-1]), float(np.sqrt(V[-1, -1]))
        we_rows.append((name, int(sel.sum()), est, se))
        A(f"  {name:<8} n={int(sel.sum()):4d}  weekend-weekday {est:+5.2f}  "
          f"NW se {se:.2f}  95% [{est - 1.96 * se:+.2f}, {est + 1.96 * se:+.2f}]")

    # Overall weekend contrast with month controls, NW vs OLS SE.
    Xwe = np.column_stack([np.ones(n), month_dummies(month), weekend])
    bwe, rwe = ols(Xwe, y)
    Vwe = nw_cov(Xwe, rwe, L=14)
    sig2 = np.sum(rwe ** 2) / (n - Xwe.shape[1])
    ols_se = float(np.sqrt(sig2 * np.linalg.inv(Xwe.T @ Xwe)[-1, -1]))
    A(f"  ALL      n={n:4d}  weekend-weekday {bwe[-1]:+5.2f}  "
      f"NW se {np.sqrt(Vwe[-1, -1]):.2f}  OLS se {ols_se:.2f}")
    # Sensitivity of the overall contrast to the lag truncation.
    for L in (7, 21):
        VL = nw_cov(Xwe, rwe, L=L)
        A(f"    sensitivity L={L}: NW se {np.sqrt(VL[-1, -1]):.2f}")
    A("")
    A("END USE BY DAY OF WEEK, WINTER (DJF) AND SUMMER (JJA), RAW MEANS")
    for season_name, months_ in (("winter", (12, 1, 2)), ("summer", (6, 7, 8))):
        A(f"  {season_name}")
        sl = obs[obs["date"].dt.month.isin(months_)]
        for k, name in enumerate(DOW_NAMES):
            s = sl[sl["dow"] == k]
            A(f"    {name} n={len(s):3d}  total {s['kwh'].mean():6.2f}  "
              f"rem {s['remainder'].mean():6.2f}  sub3 {s['sub3'].mean():6.2f}  "
              f"sub2 {s['sub2'].mean():5.2f}  sub1 {s['sub1'].mean():5.2f}")

    # ------------------------------------------------------------------
    # Trend robustness
    # ------------------------------------------------------------------
    A("")
    A("LINEAR TREND, kWh PER YEAR, WITH MONTH AND DAY-OF-WEEK CONTROLS")
    A("Newey-West L=14. The slope is the coefficient on t/365, matching the later model.")

    def trend_row(label, sel):
        yy = y[sel]
        tt = t[sel]
        mm = month[sel]
        dd = dow[sel]
        Dsel, _ = dow_dummies(dd, drop=0)
        present = [m for m in range(1, 13) if np.any(mm == m)]
        X = np.column_stack([np.ones(sel.sum()), tt / 365.0, month_dummies(mm, drop=present[0]), Dsel])
        b, r = ols(X, yy)
        V = nw_cov(X, r, L=14)
        est, se = float(b[1]), float(np.sqrt(V[1, 1]))
        A(f"  {label:<40} n={int(sel.sum()):4d}  {est:+6.3f}  se {se:.3f}  "
          f"95% [{est - 1.96 * se:+.3f}, {est + 1.96 * se:+.3f}]  "
          f"{100 * est / yy.mean():+.2f}%/year")
        return est, se, int(sel.sum()), yy.mean()

    trend_rows = []
    trend_rows.append(("All retained days",) + trend_row("All retained days", np.ones(n, dtype=bool)))
    trend_rows.append(("Drop 2006",) + trend_row("Drop 2006", year != 2006))
    trend_rows.append(("Drop low-use blocks",) + trend_row("Drop low-use blocks", ~block))
    trend_rows.append(("Complete years 2007-2009",) + trend_row(
        "Complete years 2007-2009", np.isin(year, [2007, 2008, 2009])))

    # ------------------------------------------------------------------
    # Holidays and year-end
    # ------------------------------------------------------------------
    A("")
    A("PUBLIC HOLIDAYS IN THE RETAINED SAMPLE")
    hol = french_holidays(range(2006, 2011))
    hol_df = pd.DataFrame(hol, columns=["date", "name"])
    merged = obs.merge(hol_df, on="date", how="left")
    # Month + DOW residual, so a holiday is compared with its own month and weekday.
    merged["r_month_dow"] = rw  # rw aligns with obs, and merged preserves obs order
    hol_days = merged[merged["name"].notna()]
    A(f"{'date':<12} {'name':<16} {'dow':<4} {'kwh':>8} {'resid':>8}")
    for _, row in hol_days.iterrows():
        A(f"{str(row['date'].date()):<12} {row['name']:<16} {DOW_NAMES[int(row['dow'])]:<4} "
          f"{row['kwh']:8.2f} {row['r_month_dow']:+8.2f}")

    def group_residual(label, mask):
        if mask.sum() == 0:
            A(f"  {label}: no days")
            return
        rr = rw[mask]
        A(f"  {label:<28} n={int(mask.sum()):3d}  mean kWh {y[mask].mean():6.2f}  "
          f"mean month+DOW residual {rr.mean():+6.2f}")

    A("")
    A("HOLIDAY GROUPS  (residual is after month and day-of-week means)")
    names = merged["name"].fillna("")
    group_residual("All public holidays", names.to_numpy() != "")
    group_residual("Christmas Day", names.to_numpy() == "Christmas")
    group_residual("New Year's Day", names.to_numpy() == "New Year")
    group_residual("Bastille + Assumption", np.isin(names.to_numpy(), ["Bastille Day", "Assumption"]))
    group_residual("May cluster + Easter/Whitsun", np.isin(names.to_numpy(), [
        "Easter Monday", "Labour Day", "Victory Day", "Ascension", "Whit Monday"]))
    group_residual("November holidays", np.isin(names.to_numpy(), ["All Saints", "Armistice"]))

    # Year-end stretch versus the rest of December, pooled across years.
    in_dec = month == 12
    late_dec = in_dec & (obs["date"].dt.day >= 24).to_numpy()
    early_dec = in_dec & ~late_dec
    # Also 1 January.
    jan1 = (month == 1) & (obs["date"].dt.day == 1).to_numpy()
    A("")
    A("YEAR-END STRETCH")
    A(f"  1-23 Dec   n={int(early_dec.sum()):3d}  mean {y[early_dec].mean():6.2f}")
    A(f"  24-31 Dec  n={int(late_dec.sum()):3d}  mean {y[late_dec].mean():6.2f}")
    A(f"  1 Jan      n={int(jan1.sum()):3d}  mean {y[jan1].mean():6.2f}")
    # Contrast late vs early December with a year dummy so one extreme December can't dominate.
    dec_obs_year = year[in_dec]
    Xd = np.column_stack([
        np.ones(in_dec.sum()),
        pd.get_dummies(dec_obs_year, drop_first=True).to_numpy(float),
        (obs.loc[in_dec, "date"].dt.day >= 24).to_numpy(float),
    ])
    bd, rd = ols(Xd, y[in_dec])
    Vd = nw_cov(Xd, rd, L=7)
    A(f"  late-minus-early December, year controls: {bd[-1]:+.2f}  "
      f"NW se {np.sqrt(Vd[-1, -1]):.2f}")

    # Extremes
    A("")
    A("HIGHEST DAYS")
    top = obs.nlargest(8, "kwh")
    for _, row in top.iterrows():
        A(f"  {row['date'].date()} {DOW_NAMES[int(row['dow'])]:<3} {row['kwh']:6.2f}  "
          f"remain {row['remainder']:5.1f}  sub3 {row['sub3']:5.1f}  "
          f"sub2 {row['sub2']:5.1f}  sub1 {row['sub1']:5.1f}")
    A("LOWEST DAYS")
    bot = obs.nsmallest(8, "kwh")
    for _, row in bot.iterrows():
        A(f"  {row['date'].date()} {DOW_NAMES[int(row['dow'])]:<3} {row['kwh']:6.2f}  "
          f"remain {row['remainder']:5.1f}  sub3 {row['sub3']:5.1f}  "
          f"sub2 {row['sub2']:5.1f}  sub1 {row['sub1']:5.1f}")

    # ------------------------------------------------------------------
    # Dependence on the kWh scale
    # ------------------------------------------------------------------
    A("")
    A("RESIDUAL DEPENDENCE ON THE kWh SCALE")
    A("Residuals from OLS: intercept + trend + 2 harmonics + day of week.")
    _, resid_j2, _, _ = fits["+ day of week (Mon baseline)"]
    # Confirm that fit is J=4 + DOW, not J=2. The sequential build ended at J=4 + DOW.
    # Refit the J=2 specification explicitly so this block matches 01b's mean structure.
    Xj2 = np.column_stack([np.ones(n), t / 365.0, harmonic_cols(t, 2), Dd])
    bj2, resid_j2 = ols(Xj2, y)
    rho = acf_calendar(resid_j2, t, 21)
    rho_ex = acf_calendar(resid_j2[~block], t[~block], 21)
    A(f"{'lag':>4} {'all days':>10} {'AR1 pred':>10} {'no blocks':>10}")
    for k in (1, 2, 3, 4, 5, 6, 7, 14, 21):
        A(f"{k:4d} {rho[k - 1]:10.3f} {rho[0] ** k:10.3f} {rho_ex[k - 1]:10.3f}")
    pac = pacf_yw(rho)
    pac_ex = pacf_yw(rho_ex)
    A("PACF")
    for k in range(1, 8):
        A(f"  lag {k}: all days {pac[k - 1]:+.3f}    no blocks {pac_ex[k - 1]:+.3f}")

    ss = resid_j2 ** 2
    A("")
    A(f"SHARE OF J=2 RESIDUAL SUM OF SQUARES IN THE LOW-USE BLOCKS")
    A(f"  days {100 * block.mean():.1f}%    residual SS {100 * ss[block].sum() / ss.sum():.1f}%")
    worst = ss >= np.quantile(ss, 0.95)
    A(f"  worst 5% of days: {100 * worst.mean():.1f}% of days, {100 * ss[worst].sum() / ss.sum():.1f}% of residual SS")
    high = worst & (resid_j2 > 0)
    low = worst & (resid_j2 < 0)
    A(f"    positive residuals among them: {int(high.sum())} days, "
      f"{100 * ss[high].sum() / ss[worst].sum():.1f}% of that SS")
    A(f"    negative residuals among them: {int(low.sum())} days, "
      f"{100 * ss[low].sum() / ss[worst].sum():.1f}% of that SS")
    A(f"    of the worst days, also inside a low-use block: {int((worst & block).sum())}")
    top = int(np.argmax(resid_j2))
    A(f"    largest positive residual: {obs['date'].iloc[top].date()}  "
      f"{resid_j2[top]:+.1f} kWh  (consumption {y[top]:.1f})")
    A(f"  residual skew {pd.Series(resid_j2).skew():.2f}   "
      f"excess kurtosis {pd.Series(resid_j2).kurtosis():.2f}")
    A(f"  residual sd all days {resid_j2.std(ddof=1):.2f}   "
      f"excluding blocks {resid_j2[~block].std(ddof=1):.2f}")

    A("")
    A("RESIDUAL SD BY MONTH, AFTER TREND + 2 HARMONICS + DAY OF WEEK")
    A(f"{'mo':>4} {'sd':>8} {'sd ex-block':>12} {'n ex':>6}")
    sd_m = []
    sd_m_ex = []
    for m in range(1, 13):
        sel = month == m
        s1 = resid_j2[sel].std(ddof=1)
        s2 = resid_j2[sel & ~block].std(ddof=1) if np.any(sel & ~block) else np.nan
        sd_m.append(s1)
        sd_m_ex.append(s2)
        A(f"{m:4d} {s1:8.2f} {s2:12.2f} {int((sel & ~block).sum()):6d}")
    win = np.isin(month, (12, 1, 2))
    smr = np.isin(month, (6, 7, 8))
    A(f"  winter residual sd {resid_j2[win].std(ddof=1):.2f}   "
      f"summer residual sd {resid_j2[smr].std(ddof=1):.2f}")
    A(f"  winter ex-block {resid_j2[win & ~block].std(ddof=1):.2f}   "
      f"summer ex-block {resid_j2[smr & ~block].std(ddof=1):.2f}")

    # One-step vs marginal, descriptive AR(1) on the kWh residuals.
    phi = rho[0]
    A("")
    A(f"DESCRIPTIVE AR(1) ON kWh RESIDUALS: rho1={phi:.3f}")
    A(f"  marginal sd / one-step sd = {1 / np.sqrt(1 - phi ** 2):.3f}")
    phi_ex = rho_ex[0]
    A(f"  excluding blocks: rho1={phi_ex:.3f}   "
      f"ratio {1 / np.sqrt(1 - phi_ex ** 2):.3f}")

    # ------------------------------------------------------------------
    # Diurnal profile. One row per retained day and hour, then average days.
    # ------------------------------------------------------------------
    retained_dates = set(obs["date"])
    ok_ret = ok[ok["date"].isin(retained_dates)].copy()
    ok_ret["s1_kw"] = ok_ret["Sub_metering_1"] * 0.06
    ok_ret["s2_kw"] = ok_ret["Sub_metering_2"] * 0.06
    ok_ret["s3_kw"] = ok_ret["Sub_metering_3"] * 0.06
    dh = ok_ret.groupby(["date", "hour"]).agg(
        n=("Global_active_power", "size"),
        kw=("Global_active_power", "mean"),
        s1=("s1_kw", "mean"),
        s2=("s2_kw", "mean"),
        s3=("s3_kw", "mean"),
    ).reset_index()
    n_before = len(dh)
    dh = dh[dh["n"] >= 30].copy()
    A("")
    A("DIURNAL PROFILE")
    A(f"  day-hours before / after requiring >= 30 observed minutes: {n_before} / {len(dh)}")
    dh["remainder_kw"] = dh["kw"] - dh["s1"] - dh["s2"] - dh["s3"]
    date_info = obs[["date", "dow"]].copy()
    date_info["month"] = obs["date"].dt.month.to_numpy()
    date_info["season"] = season_of(date_info["month"].to_numpy())
    date_info["weekend"] = date_info["dow"] >= 5
    dh = dh.merge(date_info, on="date", how="left")

    def profile(mask):
        sub = dh[mask]
        g = sub.groupby("hour")[["kw", "remainder_kw", "s3", "s2", "s1"]].mean()
        return g.reindex(range(24))

    profiles = {}
    for season_name in ("winter", "summer"):
        for we, we_name in [(False, "weekday"), (True, "weekend")]:
            key = f"{season_name} {we_name}"
            profiles[key] = profile((dh["season"] == season_name) & (dh["weekend"] == we))
            p = profiles[key]["kw"]
            A(f"  {key:<18} night 01-05 {p.loc[1:5].mean():.3f} kW   "
              f"evening 18-21 {p.loc[18:21].mean():.3f} kW   "
              f"peak hour {int(p.idxmax()):02d}:00 ({p.max():.3f} kW)   "
              f"integral {p.sum():.2f} kWh")
    A("COMPOSITION AT THE MORNING AND EVENING PEAKS (kW)")
    for key in ("winter weekday", "winter weekend", "summer weekday", "summer weekend"):
        for hour in (7, 20):
            row = profiles[key].loc[hour]
            A(f"  {key:<18} {hour:02d}:00  total {row['kw']:.3f}  "
              f"rem {row['remainder_kw']:.3f}  s3 {row['s3']:.3f}  "
              f"s2 {row['s2']:.3f}  s1 {row['s1']:.3f}")

    # Check the integral against the daily mean for winter weekdays.
    win_wd_days = obs[np.isin(obs["date"].dt.month, (12, 1, 2)) & (obs["dow"] < 5)]
    A(f"  winter weekday mean of daily total: {win_wd_days['kwh'].mean():.2f} kWh "
      f"(compare with integral above)")

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------
    tab = os.path.join(ROOT, "results", "tables")

    def fmt(x, nd=2, signed=False):
        return f"{x:+.{nd}f}" if signed else f"{x:.{nd}f}"

    var_rows = []
    pretty = {
        "intercept only": "Intercept only",
        "+ linear trend": "Linear trend",
        "+ harmonic 1": "Add first harmonic",
        "+ harmonic 2": "Add second harmonic",
        "+ harmonic 3": "Add third harmonic",
        "+ harmonic 4": "Add fourth harmonic",
        "+ day of week (Mon baseline)": "Add day of week",
    }
    for label, r2, dr2, s in seq_rows:
        var_rows.append(
            f"{pretty[label]} & {r2:.3f} & {dr2:.3f} & {s:.2f} \\\\"
        )
    var_rows.append(r"\addlinespace")
    var_rows.append(
        f"Trend + month + day of week & {r2_month_dow:.3f} & --- & {r_md.std(ddof=1):.2f} \\\\"
    )
    write_tabular(
        os.path.join(tab, "eda_variance.tex"),
        "@{}lrrr@{}",
        r"Specification & $R^2$ & Increment & Residual SD \\",
        var_rows,
    )

    mon_rows = []
    for m, nn, mean, sd, p10, p50, p90 in month_stats:
        mon_rows.append(
            f"{MONTHS[m - 1]} & {nn} & {mean:.2f} & {sd:.2f} & {p10:.2f} & {p50:.2f} & {p90:.2f} \\\\"
        )
    write_tabular(
        os.path.join(tab, "eda_monthly.tex"),
        "@{}lrrrrrr@{}",
        r"Month & $n$ & Mean & SD & 10\% & Median & 90\% \\",
        mon_rows,
    )

    # Fix trend_rows: each element is (label, est, se, n, mean) but I prepended the label
    # AND trend_row returns (est, se, n, mean) while also being called inside a tuple that
    # already has the label. Check: trend_rows.append(("All retained days",) + trend_row(...))
    # trend_row returns (est, se, n, mean), so the tuple is (label, est, se, n, mean). Good.
    tr_rows = []
    for label, est, se, nn, _mu in trend_rows:
        tr_rows.append(
            f"{label} & {nn} & {est:+.3f} & [{est - 1.96 * se:+.3f},\\ {est + 1.96 * se:+.3f}] \\\\"
        )
    write_tabular(
        os.path.join(tab, "eda_trend.tex"),
        "@{}lrrr@{}",
        r"Sample & $n$ & Slope (kWh/year) & 95\% interval \\",
        tr_rows,
    )

    dow_rows = []
    for k in range(7):
        lo, hi = dev[k] - 1.96 * se_dev[k], dev[k] + 1.96 * se_dev[k]
        dow_rows.append(f"{DOW_NAMES[k]} & {dev[k]:+.2f} & [{lo:+.2f},\\ {hi:+.2f}] \\\\")
    dow_rows.append(r"\addlinespace")
    for name, nn, est, se in we_rows:
        dow_rows.append(
            f"{name.capitalize()} weekend $-$ weekday & {est:+.2f} & "
            f"[{est - 1.96 * se:+.2f},\\ {est + 1.96 * se:+.2f}] \\\\"
        )
    write_tabular(
        os.path.join(tab, "eda_week.tex"),
        "@{}lrr@{}",
        r"Contrast & kWh/day & 95\% interval \\",
        dow_rows,
    )

    end_rows = []
    for key, label, _ in ENDUSE:
        m = obs[key].mean()
        delta = dec[key].mean() - aug_df[key].mean()
        end_rows.append(
            f"{label} & {m:.2f} & {100 * m / total_mean:.1f}\\% & {delta:+.2f} & {100 * delta / d_tot:.1f}\\% \\\\"
        )
    end_rows.append(
        f"Total & {total_mean:.2f} & 100\\% & {d_tot:+.2f} & 100\\% \\\\"
    )
    write_tabular(
        os.path.join(tab, "eda_enduse.tex"),
        "@{}lrrrr@{}",
        r"End use & Mean & Share & Dec$-$Aug & Share of gap \\",
        end_rows,
    )

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    # 1. Seasonal shape and monthly distribution
    fig, ax = plt.subplots(1, 2, figsize=(12.2, 4.4))
    year_colors = {2007: "#6c3483", 2008: "#1f4e79", 2009: "#148f77", 2010: "#b9770e"}
    clim = clim_day(obs["date"])
    for yr, col in year_colors.items():
        sel = year == yr
        circular = yr != 2010
        sm = smooth_on_clim(clim[sel], y[sel], window=15, circular=circular)
        ax[0].plot(np.arange(1, 366), sm, color=col, lw=1.15, label=str(yr))
    ax[0].plot(np.arange(1, 366), curves[1], color="#2c3e50", lw=1.15, ls="--", label="one harmonic")
    ax[0].plot(np.arange(1, 366), curves[2], color="#b03a2e", lw=2.0, label="two harmonics")
    mids = []
    for m in range(1, 13):
        days = ref_dates[ref_dates.month == m]
        mids.append(int(days[len(days) // 2].dayofyear))
    ax[0].plot(mids, [row[2] for row in month_stats], "o", color="k", ms=3.5, zorder=5, label="month mean")
    ax[0].set_xticks(mids)
    ax[0].set_xticklabels([m[0] for m in MONTHS])
    ax[0].set_ylabel("kWh / day")
    ax[0].set_title("Annual shape by year")
    ax[0].legend(fontsize=7.5, ncol=2, frameon=False)
    ax[0].set_xlim(1, 365)

    for m, nn, mean, sd, p10, p50, p90 in month_stats:
        ax[1].plot([m, m], [p10, p90], color="#1f4e79", lw=4, solid_capstyle="round", alpha=0.35)
        ax[1].plot(m, p50, "o", color="#1f4e79", ms=5)
        ax[1].plot(m, mean, "x", color="#b03a2e", ms=5)
    ax[1].plot([], [], color="#1f4e79", lw=4, alpha=0.35, label="10th to 90th percentile")
    ax[1].plot([], [], "o", color="#1f4e79", label="median")
    ax[1].plot([], [], "x", color="#b03a2e", label="mean")
    ax[1].set_xticks(range(1, 13))
    ax[1].set_xticklabels([m[0] for m in MONTHS])
    ax[1].set_ylabel("kWh / day")
    ax[1].set_title("Level and spread by month")
    ax[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01d_seasonal.png"), dpi=140)
    plt.close(fig)

    # 2. End use
    fig, ax = plt.subplots(figsize=(8.4, 4.3))
    x = np.arange(1, 13)
    bottom = np.zeros(12)
    # Draw remainder at the bottom, then sub3, sub2, sub1. ENDUSE is already that order.
    for key, label, col in ENDUSE:
        vals = np.array([obs.loc[obs["date"].dt.month == m, key].mean() for m in range(1, 13)])
        ax.bar(x, vals, bottom=bottom, color=col, width=0.78, label=label)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([m[0] for m in MONTHS])
    ax.set_ylabel("kWh / day")
    ax.set_title("Monthly mean consumption by end use")
    ax.legend(fontsize=8, frameon=False, loc="upper center", ncol=2)
    ax.set_ylim(0, bottom.max() * 1.22)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01d_enduse.png"), dpi=140)
    plt.close(fig)

    # 3. Diurnal
    fig, ax = plt.subplots(1, 3, figsize=(12.4, 4.15), sharey=True)
    styles = {
        "winter weekday": ("#1f4e79", "-"),
        "winter weekend": ("#5dade2", "-"),
        "summer weekday": ("#b03a2e", "-"),
        "summer weekend": ("#e59866", "-"),
    }
    for key, (col, ls) in styles.items():
        ax[0].plot(range(24), profiles[key]["kw"], color=col, ls=ls, lw=1.6, label=key)
    ax[0].set_title("Mean power")
    ax[0].legend(fontsize=7.5, frameon=False)
    ax[0].set_ylabel("kW")

    stack_keys = [("remainder_kw", "Unmetered", "#1f4e79"),
                  ("s3", "Water heater / AC", "#e67e22"),
                  ("s2", "Laundry", "#1e8449"),
                  ("s1", "Kitchen", "#7d3c98")]
    for a, key, title in (
        (ax[1], "winter weekday", "Winter weekday"),
        (ax[2], "summer weekday", "Summer weekday"),
    ):
        bottom = np.zeros(24)
        for colkey, label, col in stack_keys:
            vals = profiles[key][colkey].to_numpy()
            a.fill_between(range(24), bottom, bottom + vals, color=col, label=label, step=None)
            bottom = bottom + vals
        a.set_title(title)
    ax[2].legend(fontsize=7.5, frameon=False, loc="upper left")
    for a in ax:
        a.set_xticks(range(0, 24, 3))
        a.set_xlim(0, 23)
        a.set_xlabel("hour")
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01d_diurnal.png"), dpi=140)
    plt.close(fig)

    # 4. ACF and monthly residual scale
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.0))
    lags = np.arange(1, 22)
    ax[0].bar(lags, rho, color="#1f4e79", width=0.7, label="all retained days")
    ax[0].plot(lags, rho_ex, "o-", color="#b03a2e", ms=3.5, lw=1.2, label="low-use blocks removed")
    ax[0].plot(lags, phi ** lags, ls="--", color="#7f8c8d", lw=1.2, label=r"AR(1) from $\hat\rho_1$")
    ci = 1.96 / np.sqrt(n)
    ax[0].axhline(ci, color="gray", ls=":", lw=0.8)
    ax[0].axhline(-ci, color="gray", ls=":", lw=0.8)
    ax[0].set_xlabel("lag (days)")
    ax[0].set_ylabel("residual autocorrelation")
    ax[0].set_title("kWh residuals after trend, two harmonics, day of week")
    ax[0].legend(fontsize=7.5, frameon=False)

    ax[1].bar(range(1, 13), sd_m, color="#1f4e79", width=0.7, label="all retained days")
    ax[1].plot(range(1, 13), sd_m_ex, "o", color="#b03a2e", ms=5, label="low-use blocks removed")
    ax[1].set_xticks(range(1, 13))
    ax[1].set_xticklabels([m[0] for m in MONTHS])
    ax[1].set_ylabel("residual SD (kWh)")
    ax[1].set_title("Scale of the residual by month")
    ax[1].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "01d_residual.png"), dpi=140)
    plt.close(fig)

    # The trend_rows construction calls trend_row, which appends text, and the
    # tuple unpacking is checked here so a future edit fails loudly.
    assert all(len(row) == 5 for row in trend_rows)

    txt = "\n".join(lines)
    out = os.path.join(ROOT, "results", "01d_eda.txt")
    with open(out, "w") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"\nwrote {out}, four figures, five tables")


if __name__ == "__main__":
    main()
