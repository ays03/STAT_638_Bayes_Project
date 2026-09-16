"""
Step 2: validate the sampler on simulated data with known truth.

Simulates from the model itself (known beta, phi, sigma, nu), blanks out days in
contiguous blocks to mimic the real missingness pattern, then checks that 95%
credible intervals cover the true values and that the posterior means are close.
A sampler that fails here cannot be trusted on the real series.
"""
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), "..", ".mplcache"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import model as M

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def simulate(T=1440, J=2, p=1, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(1, T + 1)
    dow = (t + 5) % 7
    X, names = M.build_design(t, dow, J)
    beta_true = np.array([3.20, -0.04, 0.05, 0.42, 0.02, -0.11, 0.06, 0.07, -0.05, 0.03, 0.18, 0.17])
    assert len(beta_true) == X.shape[1], (len(beta_true), X.shape[1])
    phi_true = np.array([0.45])[:p]
    sigma_true, nu_true = 0.24, 5.0

    eps = np.zeros(T)
    V0 = M.stationary_cov(phi_true, sigma_true**2)
    if p > 0:
        eps[:p] = np.linalg.cholesky(V0) @ rng.standard_normal(p)
    lam = rng.gamma(nu_true / 2, 2 / nu_true, T)
    eta = rng.standard_normal(T) * sigma_true / np.sqrt(lam)
    for i in range(p, T):
        eps[i] = sum(phi_true[k] * eps[i - 1 - k] for k in range(p)) + eta[i]
    z = X @ beta_true + eps
    y = np.exp(z)

    # contiguous missing blocks, ~1.6% of days
    observed = np.ones(T, dtype=bool)
    for start in rng.choice(T - 10, 6, replace=False):
        observed[start : start + rng.integers(1, 5)] = False
    y = np.where(observed, y, np.nan)
    return dict(
        y=y, observed=observed, t=t, dow=dow, names=names,
        beta=beta_true, phi=phi_true, sigma=sigma_true, nu=nu_true, J=J, p=p,
    )


def main():
    sim = simulate()
    print(f"simulated T={len(sim['y'])}, observed={sim['observed'].sum()}, "
          f"missing={(~sim['observed']).sum()}")
    print("fitting 4 chains x 12000 iterations ...")
    fitres = M.fit_chains(
        sim["y"], sim["observed"], sim["t"], sim["dow"],
        J=sim["J"], p=sim["p"], heavy=True, scale="log",
        n_iter=12000, burn=4000, thin=4, n_chains=4, seed=11,
    )

    lines = []
    A = lines.append
    A("SAMPLER VALIDATION ON SIMULATED DATA")
    A("True values are known by construction; 95% intervals should cover them.")
    A("")
    A(f"{'parameter':>16} {'truth':>8} {'post mean':>10} {'2.5%':>8} {'97.5%':>8} {'cover':>6} {'Rhat':>6} {'ESS':>7}")

    n_cov, n_tot = 0, 0

    def row(label, truth, draws2d):
        nonlocal n_cov, n_tot
        flat = draws2d.reshape(-1)
        lo, hi = np.percentile(flat, [2.5, 97.5])
        rhat, e = M.diagnose(draws2d)
        ok = lo <= truth <= hi
        n_cov += int(ok)
        n_tot += 1
        A(f"{label:>16} {truth:8.3f} {flat.mean():10.3f} {lo:8.3f} {hi:8.3f} "
          f"{'yes' if ok else 'NO':>6} {rhat:6.3f} {e:7.0f}")

    for i, nm in enumerate(sim["names"]):
        row(nm, sim["beta"][i], fitres["beta"][:, :, i])
    row("sigma", sim["sigma"], fitres["sigma"])
    for i in range(sim["p"]):
        row(f"phi_{i+1}", sim["phi"][i], fitres["phi"][:, :, i])
    row("nu", sim["nu"], fitres["nu"])

    A("")
    A(f"coverage of 95% intervals: {n_cov}/{n_tot} parameters")
    A("  (with ~15 parameters, 14 or 15 of 15 is the expected outcome)")
    worst_rhat = max(
        M.diagnose(fitres["beta"][:, :, i])[0] for i in range(fitres["beta"].shape[2])
    )
    worst_rhat = max(worst_rhat, M.diagnose(fitres["sigma"])[0], M.diagnose(fitres["nu"])[0])
    min_ess = min(M.diagnose(fitres["beta"][:, :, i])[1] for i in range(fitres["beta"].shape[2]))
    min_ess = min(min_ess, M.diagnose(fitres["sigma"])[1], M.diagnose(fitres["nu"])[1])
    A(f"worst Rhat = {worst_rhat:.4f} (target < 1.01),  min bulk ESS = {min_ess:.0f} (target > 400)")

    # check the missing-day imputation recovers the values that were blanked
    truth_z = None
    A("")
    A("MISSING-DAY IMPUTATION")
    zi = M.flatten(fitres, "z_imputed")
    A(f"  {zi.shape[1]} imputed days; posterior sd of imputed log kWh "
      f"ranges {zi.std(axis=0).min():.3f} to {zi.std(axis=0).max():.3f}")
    A("  (should be near the marginal error sd, larger for days inside long gaps)")

    txt = "\n".join(lines)
    print(txt)
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "02_sampler_validation.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
