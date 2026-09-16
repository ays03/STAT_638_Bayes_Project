"""
Bayesian harmonic regression with AR(p), heavy-tailed errors, and missing days.

MODEL
    Response on a chosen scale: z_t = g(y_t), with g = log or identity.
    Defined on a gap-free daily grid t = 1, ..., T, so the error process exists
    on every calendar day whether or not y_t was observed.

        z_t = x_t' beta + eps_t
        eps_t = sum_{k=1}^{p} phi_k eps_{t-k} + eta_t
        eta_t ~ N(0, sigma^2 / lambda_t),   lambda_t ~ Gamma(nu/2, nu/2)

    Marginalising lambda_t gives eta_t ~ t_nu(0, sigma^2). Setting heavy=False
    fixes lambda_t = 1 (Gaussian errors); setting p = 0 gives independent errors.

    x_t contains: intercept, centred linear trend in years, J harmonic pairs at
    period 365 days, and six day-of-week dummies (Monday is the baseline).

    (eps_1, ..., eps_p) is drawn from the exact stationary distribution of the
    AR(p), so the chain starts stationary and no burn-in of the error process
    is needed.

PRIORS
    beta        ~ N(0, diag(tau^2)), tau set by response scale (weakly informative)
    sigma       ~ half-t(3, 0, A), via the standard inverse-gamma auxiliary variable
    phi         ~ Uniform over the stationary region
    nu - 2      ~ Exponential(mean 10)
    z_t missing ~ implied by the AR process (imputed, not dropped)

SAMPLING
    Gibbs for beta (conjugate normal), sigma^2 and its auxiliary xi
    (inverse-gamma), lambda (gamma), and the missing z_t (normal).
    Random-walk Metropolis for phi (rejecting outside the stationary region)
    and for log(nu - 2).
"""
import warnings

import numpy as np
from scipy.special import gammaln
from scipy.stats import norm as _norm
from scipy.stats import t as _tdist

# macOS Accelerate's BLAS sets IEEE status flags inside its SIMD kernels, so
# numpy reports spurious divide/overflow/invalid warnings from matmul even for
# well-conditioned inputs with entirely finite results. Verified harmless.
warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

# ----------------------------------------------------------------------------
# design matrix
# ----------------------------------------------------------------------------


def build_design(t, dow, J, t_center=None):
    """t: 1-based integer day index. dow: 0=Monday. Returns (X, names)."""
    t = np.asarray(t, dtype=float)
    if t_center is None:
        t_center = t.mean()
    cols = [np.ones_like(t), (t - t_center) / 365.0]
    names = ["intercept", "trend_per_year"]
    for j in range(1, J + 1):
        cols += [np.sin(2 * np.pi * j * t / 365.0), np.cos(2 * np.pi * j * t / 365.0)]
        names += [f"sin_h{j}", f"cos_h{j}"]
    for k in range(1, 7):
        cols.append((np.asarray(dow) == k).astype(float))
        names.append(f"dow_{['Mon','Tue','Wed','Thu','Fri','Sat','Sun'][k]}")
    return np.column_stack(cols), names


# ----------------------------------------------------------------------------
# AR(p) helpers
# ----------------------------------------------------------------------------


def is_stationary(phi):
    """Roots of 1 - phi_1 B - ... - phi_p B^p outside the unit circle."""
    p = len(phi)
    if p == 0:
        return True
    companion = np.zeros((p, p))
    companion[0, :] = phi
    if p > 1:
        companion[1:, :-1] = np.eye(p - 1)
    return np.max(np.abs(np.linalg.eigvals(companion))) < 1.0 - 1e-8


def stationary_cov(phi, sigma2):
    """Covariance of (eps_1, ..., eps_p) under the stationary AR(p)."""
    p = len(phi)
    if p == 0:
        return np.zeros((0, 0))
    if p == 1:
        return np.array([[sigma2 / (1.0 - phi[0] ** 2)]])
    # solve the discrete Lyapunov equation V = A V A' + Q for the companion form
    A = np.zeros((p, p))
    A[0, :] = phi
    A[1:, :-1] = np.eye(p - 1)
    Q = np.zeros((p, p))
    Q[0, 0] = sigma2
    # vec(V) = (I - A kron A)^{-1} vec(Q)
    M = np.eye(p * p) - np.kron(A, A)
    V = np.linalg.solve(M, Q.reshape(-1)).reshape(p, p)
    return 0.5 * (V + V.T)


def ar_resid(eps, phi):
    """eta_t for t = p..T-1 (0-based), i.e. eps_t - sum phi_k eps_{t-k}."""
    p = len(phi)
    if p == 0:
        return eps.copy()
    out = eps[p:].copy()
    for k in range(1, p + 1):
        out -= phi[k - 1] * eps[p - k : len(eps) - k]
    return out


def quasi_difference(M, phi):
    """Apply (1 - phi_1 B - ... - phi_p B^p) down the rows of M."""
    p = len(phi)
    if p == 0:
        return M.copy()
    out = M[p:].copy()
    for k in range(1, p + 1):
        out -= phi[k - 1] * M[p - k : len(M) - k]
    return out


def marginal_sd_ratio(phi):
    """sd of the stationary AR(p) divided by the one-step innovation sd."""
    if len(phi) == 0:
        return 1.0
    return float(np.sqrt(stationary_cov(np.asarray(phi), 1.0)[0, 0]))


# ----------------------------------------------------------------------------
# sampler
# ----------------------------------------------------------------------------

PRIOR_TAU = {
    # (intercept sd, other-coefficient sd) by response scale
    "log": (5.0, 2.0),
    "identity": (100.0, 40.0),
}
PRIOR_SIGMA_A = {"log": 1.0, "identity": 20.0}


def fit(
    y,
    observed,
    t,
    dow,
    J=1,
    p=1,
    heavy=True,
    scale="log",
    n_iter=20000,
    burn=5000,
    thin=5,
    seed=1,
    ll_index=None,
    n_paths=400,
    verbose=False,
):
    """Run one chain. y and observed are on the full daily grid (length T);
    y may contain NaN wherever observed is False.

    ll_index fixes the days at which the pointwise log-likelihood is recorded.
    Pass the same index set to every model so WAIC compares like with like;
    if None, the model's own maximal usable set is used.

    Returns a dict of posterior draws (rows = saved iterations)."""
    rng = np.random.default_rng(seed)
    T = len(y)
    observed = np.asarray(observed, dtype=bool)
    X, names = build_design(t, dow, J)
    k = X.shape[1]

    z = np.where(observed, np.log(y) if scale == "log" else y, np.nan)
    # initialise missing days at the mean of the observed response
    z = np.where(np.isnan(z), np.nanmean(z), z)
    miss_idx = np.where(~observed)[0]

    tau_int, tau_oth = PRIOR_TAU[scale]
    prior_prec = np.diag(np.r_[1.0 / tau_int**2, np.full(k - 1, 1.0 / tau_oth**2)])
    sigma_A = PRIOR_SIGMA_A[scale]
    nu_sigma = 3.0
    nu_exp_mean = 10.0

    # --- initial values from OLS ---
    beta = np.linalg.lstsq(X[observed], z[observed], rcond=None)[0]
    eps = z - X @ beta
    sigma2 = float(np.var(eps[observed], ddof=k))
    xi = 1.0
    phi = np.full(p, 0.3 / max(p, 1)) if p > 0 else np.zeros(0)
    lam = np.ones(T)
    nu = 8.0

    # adaptive random-walk scales
    step_phi = 0.05 * np.ones(p) if p > 0 else np.zeros(0)
    step_nu = 0.3
    acc_phi, acc_nu, n_prop = 0, 0, 0

    keep = np.arange(burn, n_iter, thin)
    n_keep = len(keep)
    out = {
        "beta": np.zeros((n_keep, k)),
        "sigma": np.zeros(n_keep),
        "phi": np.zeros((n_keep, p)),
        "nu": np.zeros(n_keep),
        "eps_last": np.zeros((n_keep, max(p, 1))),
        "z_imputed": np.zeros((n_keep, len(miss_idx))),
        "loglik": None,
    }
    # pointwise log-likelihood for WAIC: days that are observed and whose p
    # preceding calendar days are also observed
    if ll_index is not None:
        ll_idx = np.asarray(ll_index)
    elif p > 0:
        ll_idx = np.array(
            [i for i in range(p, T) if observed[i] and all(observed[i - kk] for kk in range(1, p + 1))]
        )
    else:
        ll_idx = np.where(observed)[0]
    out["loglik"] = np.zeros((n_keep, len(ll_idx)))
    out["ll_idx"] = ll_idx
    out["names"] = names

    # thinned error paths, for one-step-ahead residuals and anomaly detection
    path_every = max(1, n_keep // max(n_paths, 1))
    path_slots = len(range(0, n_keep, path_every))
    out["eps_path"] = np.zeros((path_slots, T), dtype=np.float32)
    out["lam_path"] = np.zeros((path_slots, T), dtype=np.float32)
    out["path_pars"] = np.zeros((path_slots, 3 + p), dtype=np.float64)  # sigma, nu, spare, phi
    save, psave = 0, 0

    def log_prior_nu(v):
        return -(v - 2.0) / nu_exp_mean if v > 2.0 else -np.inf

    def eta_and_init(eps_, phi_):
        return ar_resid(eps_, phi_), eps_[:p]

    for it in range(n_iter):
        # ---------------- beta | phi, sigma2, lambda, z ----------------
        Xt = quasi_difference(X, phi)
        zt = quasi_difference(z.reshape(-1, 1), phi).ravel()
        w = lam[p:] / sigma2
        prec = prior_prec + (Xt * w[:, None]).T @ Xt
        rhs = (Xt * w[:, None]).T @ zt
        if p > 0:
            V0 = stationary_cov(phi, sigma2)
            P0 = np.linalg.inv(V0)
            X0, z0 = X[:p], z[:p]
            prec = prec + X0.T @ P0 @ X0
            rhs = rhs + X0.T @ P0 @ z0
        L = np.linalg.cholesky(prec)
        mean = np.linalg.solve(L.T, np.linalg.solve(L, rhs))
        beta = mean + np.linalg.solve(L.T, rng.standard_normal(k))
        eps = z - X @ beta

        # ---------------- sigma2 | rest ----------------
        eta = ar_resid(eps, phi)
        ss = float(np.sum(lam[p:] * eta**2))
        n_eff = T - p
        if p > 0:
            # the initial block contributes eps_{1:p}' (sigma2 V)^{-1} eps_{1:p},
            # and V is proportional to sigma2, so it scales out cleanly
            Vunit = stationary_cov(phi, 1.0)
            ss += float(eps[:p] @ np.linalg.solve(Vunit, eps[:p]))
            n_eff += p
        shape = 0.5 * (nu_sigma + n_eff)
        rate = nu_sigma / (2.0 * xi) + 0.5 * ss
        sigma2 = float(rate / rng.gamma(shape, 1.0))
        xi = float((1.0 / sigma_A**2 + nu_sigma / sigma2) / rng.gamma(0.5 * (nu_sigma + 1), 1.0))

        # ---------------- lambda | rest ----------------
        if heavy:
            eta = ar_resid(eps, phi)
            shape_l = 0.5 * (nu + 1.0)
            rate_l = 0.5 * (nu + eta**2 / sigma2)
            lam[p:] = rng.gamma(shape_l, 1.0 / rate_l)
            lam[:p] = 1.0
        else:
            lam[:] = 1.0

        # ---------------- nu | lambda  (Metropolis on log(nu-2)) ----------------
        if heavy:
            n_prop += 1
            lo = lam[p:]
            prop = 2.0 + np.exp(np.log(nu - 2.0) + step_nu * rng.standard_normal())

            sum_log_lo, sum_lo, m_lo = np.sum(np.log(lo)), np.sum(lo), len(lo)

            def ll_nu(v):
                a = 0.5 * v
                return m_lo * (a * np.log(a) - gammaln(a)) + (a - 1) * sum_log_lo - a * sum_lo

            log_r = (
                ll_nu(prop) + log_prior_nu(prop) + np.log(prop - 2.0)
                - ll_nu(nu) - log_prior_nu(nu) - np.log(nu - 2.0)
            )
            if np.log(rng.random()) < log_r:
                nu = prop
                acc_nu += 1

        # ---------------- phi | rest  (Metropolis) ----------------
        if p > 0:
            def log_post_phi(ph):
                if not is_stationary(ph):
                    return -np.inf
                e = ar_resid(eps, ph)
                Vunit = stationary_cov(ph, 1.0)
                sign, logdet = np.linalg.slogdet(Vunit)
                q0 = float(eps[:p] @ np.linalg.solve(Vunit, eps[:p]))
                return (
                    -0.5 * np.sum(lam[p:] * e**2) / sigma2
                    - 0.5 * logdet
                    - 0.5 * q0 / sigma2
                )

            cur = log_post_phi(phi)
            prop = phi + step_phi * rng.standard_normal(p)
            if np.log(rng.random()) < log_post_phi(prop) - cur:
                phi = prop
                acc_phi += 1

        # ---------------- missing z_t | rest ----------------
        if len(miss_idx) > 0:
            for i in miss_idx:
                # quadratic in eps_i: collect derivative and curvature
                d_terms, s_terms = [], []
                if i >= p:
                    d_terms.append(1.0)
                    s_terms.append(i)
                for kk in range(1, p + 1):
                    if i + kk < T and i + kk >= p:
                        d_terms.append(-phi[kk - 1])
                        s_terms.append(i + kk)
                if not d_terms:
                    continue
                e_all = ar_resid(eps, phi)
                grad, curv = 0.0, 0.0
                for d, s in zip(d_terms, s_terms):
                    grad += lam[s] * e_all[s - p] * d
                    curv += lam[s] * d * d
                if i < p:  # also inside the stationary initial block
                    Vunit = stationary_cov(phi, 1.0)
                    Pi = np.linalg.inv(Vunit)
                    grad += float(Pi[i] @ eps[:p])
                    curv += float(Pi[i, i])
                if curv <= 0:
                    continue
                cond_mean_eps = eps[i] - grad / curv
                eps[i] = cond_mean_eps + np.sqrt(sigma2 / curv) * rng.standard_normal()
                z[i] = X[i] @ beta + eps[i]

        # ---------------- adapt during burn-in ----------------
        if it < burn and (it + 1) % 200 == 0 and n_prop > 0:
            if p > 0:
                r = acc_phi / 200.0
                step_phi *= np.exp((r - 0.3) * 1.0)
                step_phi = np.clip(step_phi, 1e-4, 1.0)
                acc_phi = 0
            if heavy:
                r = acc_nu / 200.0
                step_nu *= np.exp((r - 0.3) * 1.0)
                step_nu = float(np.clip(step_nu, 1e-3, 3.0))
                acc_nu = 0
            n_prop = 0

        # ---------------- save ----------------
        if it >= burn and (it - burn) % thin == 0 and save < n_keep:
            out["beta"][save] = beta
            out["sigma"][save] = np.sqrt(sigma2)
            if p > 0:
                out["phi"][save] = phi
                out["eps_last"][save, :p] = eps[T - p : T][::-1]  # [eps_T, eps_{T-1}, ...]
            else:
                out["eps_last"][save, 0] = 0.0
            out["nu"][save] = nu if heavy else np.inf
            out["z_imputed"][save] = z[miss_idx]
            # one-step-ahead conditional log density, marginal over lambda
            e = ar_resid(eps, phi) if p > 0 else eps
            ll = _log_t(e[ll_idx - p] if p > 0 else e[ll_idx], np.sqrt(sigma2), nu, heavy)
            if scale == "log":
                ll = ll - z[ll_idx]  # Jacobian: density on the kWh scale
            out["loglik"][save] = ll
            if save % path_every == 0 and psave < path_slots:
                out["eps_path"][psave] = eps
                out["lam_path"][psave] = lam
                out["path_pars"][psave, 0] = np.sqrt(sigma2)
                out["path_pars"][psave, 1] = nu if heavy else np.inf
                if p > 0:
                    out["path_pars"][psave, 3:] = phi
                psave += 1
            save += 1

        if verbose and (it + 1) % 5000 == 0:
            print(f"    iter {it+1}/{n_iter}  sigma={np.sqrt(sigma2):.3f} "
                  f"phi={np.round(phi,3)} nu={nu:.1f}")

    out["miss_idx"] = miss_idx
    out["X"] = X
    out["scale"] = scale
    out["J"], out["p"], out["heavy"] = J, p, heavy
    return out


def _log_t(e, sigma, nu, heavy):
    """Log density of a scaled t (or normal if heavy is False). sigma and nu may
    be scalars or arrays broadcasting against e, so the same routine serves both
    the sampler (scalar parameters, vector residuals) and the predictive
    mixture (vector parameters, one residual per draw)."""
    if not heavy or not np.all(np.isfinite(nu)):
        return -0.5 * np.log(2 * np.pi) - np.log(sigma) - 0.5 * (e / sigma) ** 2
    return (
        gammaln(0.5 * (nu + 1))
        - gammaln(0.5 * nu)
        - 0.5 * np.log(np.pi * nu)
        - np.log(sigma)
        - 0.5 * (nu + 1) * np.log1p((e / sigma) ** 2 / nu)
    )


# ----------------------------------------------------------------------------
# multi-chain driver
# ----------------------------------------------------------------------------


def fit_chains(y, observed, t, dow, n_chains=4, seed=1, **kw):
    chains = [fit(y, observed, t, dow, seed=seed + 100 * c, **kw) for c in range(n_chains)]
    merged = {}
    for key in ("beta", "sigma", "phi", "nu", "eps_last", "z_imputed", "loglik"):
        merged[key] = np.stack([c[key] for c in chains])  # (chain, draw, ...)
    for key in ("eps_path", "lam_path", "path_pars"):
        merged[key] = np.concatenate([c[key] for c in chains], axis=0)
    for key in ("names", "ll_idx", "miss_idx", "X", "scale", "J", "p", "heavy"):
        merged[key] = chains[0][key]
    return merged


def common_ll_index(observed, max_lag=2):
    """Days that are observed and whose max_lag preceding days are also observed.
    Use one shared set across models so WAIC is computed on identical terms."""
    observed = np.asarray(observed, dtype=bool)
    T = len(observed)
    return np.array(
        [i for i in range(max_lag, T) if observed[i] and all(observed[i - k] for k in range(1, max_lag + 1))]
    )


def one_step_diagnostics(merged, z_obs, observed):
    """From the saved error paths, compute per-day one-step-ahead standardised
    residuals and posterior predictive p-values under the fitted model.

    Returns dict with arrays of length T (NaN on days not evaluable)."""
    eps = merged["eps_path"].astype(np.float64)  # (S, T)
    pars = merged["path_pars"]
    sigma, nu = pars[:, 0], pars[:, 1]
    p = merged["p"]
    phi = pars[:, 3:] if p > 0 else np.zeros((len(pars), 0))
    S, T = eps.shape

    pred = np.zeros((S, T))  # one-step-ahead prediction of eps_t
    for k in range(1, p + 1):
        pred[:, k:] += phi[:, [k - 1]] * eps[:, : T - k]
    e = eps - pred  # innovation
    std = e / sigma[:, None]

    heavy = merged["heavy"]
    if heavy:
        cdf = _tdist.cdf(std, df=nu[:, None])
    else:
        cdf = _norm.cdf(std)

    ok = np.asarray(observed, dtype=bool).copy()
    for k in range(1, p + 1):
        ok[k:] &= np.asarray(observed, dtype=bool)[: T - k]
    ok[:p] = False

    res = {
        "std_resid_mean": np.where(ok, std.mean(axis=0), np.nan),
        "std_resid_sd": np.where(ok, std.std(axis=0, ddof=1), np.nan),
        "ppp": np.where(ok, cdf.mean(axis=0), np.nan),
        "lam_mean": merged["lam_path"].astype(np.float64).mean(axis=0),
        "evaluable": ok,
    }
    return res


def flatten(merged, key):
    """Collapse the (chain, draw, ...) axes into one. Handles zero-width
    trailing dimensions, which occur when p = 0 leaves phi empty."""
    a = merged[key]
    return a.reshape((a.shape[0] * a.shape[1],) + a.shape[2:])


# ----------------------------------------------------------------------------
# convergence diagnostics (rank-normalised split-Rhat and bulk ESS)
# ----------------------------------------------------------------------------


def _rank_normalise(x):
    shape = x.shape
    flat = x.reshape(-1)
    ranks = flat.argsort().argsort() + 1
    return _norm.ppf((ranks - 0.375) / (len(flat) + 0.25)).reshape(shape)


def split_rhat(draws):
    """draws: (chain, iter)."""
    m, n = draws.shape
    half = n // 2
    s = np.concatenate([draws[:, :half], draws[:, half : 2 * half]], axis=0)
    m2, n2 = s.shape
    means = s.mean(axis=1)
    W = s.var(axis=1, ddof=1).mean()
    B = n2 * means.var(ddof=1)
    if W <= 0:
        return np.nan
    var_hat = ((n2 - 1) * W + B) / n2
    return float(np.sqrt(var_hat / W))


def ess(draws):
    """Bulk effective sample size via Geyer's initial positive sequence."""
    m, n = draws.shape
    if np.allclose(draws, draws.flat[0]):
        return float(m * n)
    acov = np.zeros((m, n))
    for c in range(m):
        x = draws[c] - draws[c].mean()
        f = np.fft.rfft(x, 2 * n)
        acov[c] = np.fft.irfft(f * np.conj(f), 2 * n)[:n] / n
    chain_mean = draws.mean(axis=1)
    mean_var = acov[:, 0].mean() * n / (n - 1)
    var_plus = mean_var * (n - 1) / n + (chain_mean.var(ddof=1) if m > 1 else 0)
    rho = np.zeros(n)
    rho[0] = 1.0
    rho_hat_even, rho_hat_odd = 1.0, 1.0 - (mean_var - acov[:, 1].mean()) / var_plus
    tau = rho_hat_even + 2 * rho_hat_odd
    kk = 1
    while kk + 2 < n - 2 and (rho_hat_even + rho_hat_odd) > 0:
        rho_hat_even = 1.0 - (mean_var - acov[:, kk + 1].mean()) / var_plus
        rho_hat_odd = 1.0 - (mean_var - acov[:, kk + 2].mean()) / var_plus
        if rho_hat_even + rho_hat_odd > 0:
            tau += 2 * (rho_hat_even + rho_hat_odd)
        kk += 2
    tau = max(tau, 1.0 / np.log10(max(m * n, 11)))
    return float(m * n / tau)


def diagnose(draws):
    """draws: (chain, iter). Returns (rhat on rank-normalised draws, bulk ESS)."""
    rn = _rank_normalise(draws)
    return split_rhat(rn), ess(rn)


# ----------------------------------------------------------------------------
# WAIC
# ----------------------------------------------------------------------------


def waic(loglik):
    """loglik: (chain, draw, n_obs) pointwise log densities."""
    ll = loglik.reshape(-1, loglik.shape[-1])
    lppd_i = np.log(np.mean(np.exp(ll - ll.max(axis=0)), axis=0)) + ll.max(axis=0)
    pwaic_i = ll.var(axis=0, ddof=1)
    elpd_i = lppd_i - pwaic_i
    n = ll.shape[1]
    return {
        "elpd_waic": float(elpd_i.sum()),
        "p_waic": float(pwaic_i.sum()),
        "waic": float(-2 * elpd_i.sum()),
        "se": float(np.sqrt(n * elpd_i.var(ddof=1))),
        "n": n,
        "elpd_i": elpd_i,
    }


# ----------------------------------------------------------------------------
# forecasting
# ----------------------------------------------------------------------------


def forecast(
    merged, t_future, dow_future, t_center, eps_init=None, n_draws=None,
    draw_idx=None, seed=7,
):
    """Simulate the posterior predictive distribution of y on future days.

    eps_init, if given, is (N, p) with columns [eps_origin, eps_origin-1, ...],
    overriding the error state saved at the end of the fitted series. This is
    what makes rolling-origin forecasts possible without refitting: given beta,
    the error at any day with an observed response is just z_t - x_t'beta.

    Returns a dict with
        y      (N, H) predictive draws on the kWh scale
        mu     (N, H) systematic component only, kWh scale
        loc    (N, H) conditional mean of the response on the MODELLING scale,
                      before the final innovation is added. Together with
                      (sigma, nu) this gives an exact mixture representation of
                      the predictive density, so log scores need no smoothing.
        sigma, nu, scale
    """
    rng = np.random.default_rng(seed)
    beta = flatten(merged, "beta")
    sigma = flatten(merged, "sigma")
    phi = flatten(merged, "phi")
    nu = flatten(merged, "nu")
    eps_last = flatten(merged, "eps_last")
    N = beta.shape[0]
    sel = np.arange(N)
    if draw_idx is not None:
        sel = np.asarray(draw_idx)
    elif n_draws is not None and n_draws < N:
        sel = rng.choice(N, n_draws, replace=False)
    if len(sel) != N:
        beta, sigma, phi, nu, eps_last = beta[sel], sigma[sel], phi[sel], nu[sel], eps_last[sel]
        N = len(sel)
    if eps_init is not None and len(eps_init) != N:
        raise ValueError(
            f"eps_init has {len(eps_init)} rows but {N} draws are in use; build it "
            "from the same draw_idx passed here."
        )

    p = merged["p"]
    Xf, _ = build_design(t_future, dow_future, merged["J"], t_center=t_center)
    H = len(t_future)
    mu = beta @ Xf.T  # (N, H)

    if p > 0:
        hist = eps_last[:, :p].copy() if eps_init is None else np.asarray(eps_init, float).copy()
    else:
        hist = np.zeros((N, 1))

    z = np.zeros((N, H))
    loc = np.zeros((N, H))
    heavy = merged["heavy"]
    for h in range(H):
        pred_eps = np.sum(phi * hist, axis=1) if p > 0 else np.zeros(N)
        loc[:, h] = mu[:, h] + pred_eps
        lam = rng.gamma(nu / 2.0, 2.0 / nu) if heavy else np.ones(N)
        eta = rng.standard_normal(N) * sigma / np.sqrt(lam)
        e_new = pred_eps + eta
        if p > 0:
            hist = np.column_stack([e_new, hist[:, :-1]])
        z[:, h] = mu[:, h] + e_new

    out = dict(sigma=sigma, nu=nu, scale=merged["scale"], loc=loc, heavy=heavy, sel=sel)
    if merged["scale"] == "log":
        out["y"], out["mu"] = np.exp(z), np.exp(mu)
    else:
        out["y"], out["mu"] = z, mu
    return out


def log_pred_density(fc, y_actual):
    """Exact-mixture log predictive density of the observed kWh values.

    The final innovation is integrated analytically (it is t_nu or normal given
    the draw), so this is a Rao-Blackwellised estimate rather than a kernel
    smooth of the predictive sample. Returns an array of length H."""
    y_actual = np.asarray(y_actual, float)
    H = fc["loc"].shape[1]
    out = np.full(H, np.nan)
    z_act = np.log(y_actual) if fc["scale"] == "log" else y_actual
    for h in range(H):
        if not np.isfinite(z_act[h]):
            continue
        ll = _log_t(z_act[h] - fc["loc"][:, h], fc["sigma"], fc["nu"], fc["heavy"])
        if fc["scale"] == "log":
            ll = ll - z_act[h]  # Jacobian to the kWh scale
        m = ll.max()
        out[h] = m + np.log(np.mean(np.exp(ll - m)))
    return out


def pit(fc, y_actual):
    """Probability integral transform: the predictive CDF evaluated at the
    actual value, computed exactly from the mixture rather than from sample
    ranks. Well-calibrated forecasts give a uniform histogram."""
    y_actual = np.asarray(y_actual, float)
    z_act = np.log(y_actual) if fc["scale"] == "log" else y_actual
    H = fc["loc"].shape[1]
    out = np.full(H, np.nan)
    for h in range(H):
        if not np.isfinite(z_act[h]):
            continue
        s = (z_act[h] - fc["loc"][:, h]) / fc["sigma"]
        out[h] = float(np.mean(_tdist.cdf(s, df=fc["nu"]) if fc["heavy"] else _norm.cdf(s)))
    return out


def crps(samples, y_actual):
    """Exact CRPS of an ensemble forecast, per column of samples (N, H)."""
    y_actual = np.asarray(y_actual, float)
    N, H = samples.shape
    out = np.full(H, np.nan)
    for h in range(H):
        if not np.isfinite(y_actual[h]):
            continue
        x = np.sort(samples[:, h])
        term1 = np.mean(np.abs(x - y_actual[h]))
        # E|X - X'| = (2/N^2) * sum_i (2i - N + 1) * x_(i)   with i from 0
        i = np.arange(N)
        term2 = 2.0 * np.sum((2 * i - N + 1) * x) / (N * N)
        out[h] = term1 - 0.5 * term2
    return out
