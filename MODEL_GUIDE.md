# Model guide: Bayesian forecasting of household electricity demand

This guide collects the model-related explanations in one place. It describes
the data entering the model, the equations and every symbol, the priors,
estimation, model comparison, forecasts, findings, limitations, and recommended
next steps. The full technical derivations and reproducible results remain in
`report.tex`, `src/model.py`, and `results/`.

## 1. What is being modelled?

- The data come from one household near Paris, from December 2006 through
  November 2010.
- The raw file records average active power every minute in kilowatts (kW).
- A minute lasts \(1/60\) hour, so daily energy consumption is

  \[
  y_d=\frac{1}{60}\sum_{m\in d}\mathrm{GAP}_m ,
  \]

  where:
  - \(d\) is a calendar day;
  - \(m\in d\) means all recorded minutes in day \(d\);
  - \(\mathrm{GAP}_m\) is global active power during minute \(m\), in kW; and
  - \(y_d\) is total energy used that day, in kilowatt-hours (kWh).

- The factor \(1/60\) converts one minute of kW into kWh. This is a daily sum,
  not a daily mean.
- The two partial boundary days are removed. A remaining day is treated as
  observed only when at least 95% of its 1,440 minutes are present.
- The final grid has \(T=1{,}440\) calendar days: 1,417 observed days and 23
  missing days. Missing days stay in the grid so the time-series relationship
  remains continuous.

The scientific goal is to describe long-term, yearly, and weekly patterns;
forecast future daily demand with uncertainty; and identify unusual days or
periods.

## 2. General model

The complete model family is

\[
z_t=g(y_t),
\]

\[
z_t=\mathbf{x}_t^\top\boldsymbol{\beta}+\varepsilon_t,
\]

\[
\mathbf{x}_t^\top\boldsymbol{\beta}
=\beta_0
+\beta_1\frac{t-\bar t}{365}
+\sum_{j=1}^{J}\left[
  a_j\sin\left(\frac{2\pi jt}{365}\right)
  +b_j\cos\left(\frac{2\pi jt}{365}\right)
\right]
+\sum_{k=2}^{7}\delta_k
  \mathbf{1}\{\operatorname{dow}_t=k\},
\]

\[
\varepsilon_t
=\sum_{\ell=1}^{p}\phi_\ell\varepsilon_{t-\ell}+\eta_t,
\qquad
\eta_t\sim t_\nu(0,\sigma^2).
\]

In plain language:

> daily consumption = calendar-based prediction + a serially correlated
> leftover, with occasional large surprises allowed.

### How to read the notation

- \(=\): “equals.”
- \(+\) and \(-\): add and subtract.
- A fraction such as \(a/b\): divide \(a\) by \(b\).
- A subscript, as in \(y_t\): identifies a day, lag, harmonic, or category.
- A superscript 2, as in \(\sigma^2\): square the quantity.
- Bold symbols such as \(\mathbf{x}\) or \(\boldsymbol{\beta}\): vectors
  containing several values.
- \(\sum\): add a sequence of terms over the stated index range.
- \(\prod\): multiply a sequence of terms over the stated index range.
- \(\in\): “belongs to,” as in a minute belonging to a day.
- \(\approx\): “approximately equals.”
- \(<,>,\le,\ge\): less than, greater than, less than or equal to, and greater
  than or equal to.
- \(P(A\mid B)\): probability of event \(A\), given information \(B\).
- \(E(X)\), \(\operatorname{Var}(X)\), and \(\operatorname{SD}(X)\): mean,
  variance, and standard deviation of \(X\).
- \([L,U]\) following an estimate: lower and upper endpoints of an interval.
  In the reported parameter results these are 95% posterior credible
  intervals.

### The response and time symbols

- \(t\): day number on the gap-free calendar, \(t=1,\ldots,T\).
- \(T\): number of calendar days, here \(1{,}440\).
- \(y_t\): observed or latent electricity consumption on day \(t\), in kWh.
- \(g(\cdot)\): response transformation.
- \(z_t\): response on the modelling scale.
  - Identity model: \(g(y_t)=y_t\), so \(z_t=y_t\).
  - Log model: \(g(y_t)=\log y_t\), so \(z_t=\log y_t\).
- \(\bar t\): average day index. Centering time at \(\bar t\) reduces
  correlation between the trend coefficient and intercept.

The headline model M7 uses the identity transformation, so all of its
coefficients are directly on the kWh scale.

### The regression symbols

- \(\mathbf{x}_t\): vector of known predictors for day \(t\): a constant,
  centred time, seasonal sine/cosine terms, and weekday indicators.
- \(\boldsymbol{\beta}\): vector containing all regression coefficients.
- \({}^\top\): transpose. Thus
  \(\mathbf{x}_t^\top\boldsymbol{\beta}\) is the dot product of predictors and
  coefficients, producing one systematic prediction.
- \(\beta_0\): intercept. In the fitted design it uses Monday as the weekday
  reference, the trend is zero at the record midpoint, and seasonal terms are
  represented separately. It is therefore not, by itself, the overall
  average-day consumption.
- \(\beta_1\): linear change per year, because
  \((t-\bar t)/365\) measures time in years.
- \(365\): annual period required by the project.

The reported **average-day underlying level at the record midpoint** is formed
by adding the average weekday offset to the intercept. For M7 it is 25.96 kWh,
whereas the raw Monday-reference intercept posterior mean is 23.82 kWh.

### Annual harmonics

- \(J\): number of harmonic pairs.
- \(j\): harmonic number, \(j=1,\ldots,J\).
- \(j=1\): one full wave every 365 days.
- \(j=2\): two waves every 365 days, or one every 182.5 days. Combined with the
  first harmonic, it allows an asymmetric annual shape.
- \(a_j\): coefficient multiplying harmonic \(j\)'s sine term.
- \(b_j\): coefficient multiplying harmonic \(j\)'s cosine term.
- \(\pi\): the constant pi.
- \(\sin(\cdot)\), \(\cos(\cdot)\): periodic functions used to represent
  seasonality.

The sine and cosine coefficients can be converted into an easier
amplitude-and-phase form:

\[
a_j\sin(\omega_jt)+b_j\cos(\omega_jt)
=A_j\cos(\omega_jt-\psi_j),
\]

\[
\omega_j=\frac{2\pi j}{365},
\qquad
A_j=\sqrt{a_j^2+b_j^2},
\qquad
\psi_j=\operatorname{atan2}(a_j,b_j).
\]

- \(\omega_j\): angular frequency of harmonic \(j\).
- \(A_j\): amplitude, the maximum distance of that individual wave from zero.
- \(2A_j\): peak-to-trough swing of that individual wave.
- \(\psi_j\): phase, which determines the date of the peak.
- \(\operatorname{atan2}\): an angle function that uses the signs of both
  coefficients to put the phase in the correct quadrant.
- \(\sqrt{\cdot}\): square root.

The first and second harmonics must be added before finding the range of the
full seasonal curve. Therefore the full seasonal range is not generally
\(2A_1+2A_2\).

### Weekday effects

- \(\operatorname{dow}_t\): day of week for day \(t\).
- \(k\): weekday category.
- \(\delta_k\): effect of weekday \(k\) relative to Monday.
- \(\mathbf{1}\{\operatorname{dow}_t=k\}\): indicator function. It equals 1
  when day \(t\) is weekday \(k\), and 0 otherwise.
- Monday has no separate coefficient and is the reference category. Tuesday
  through Sunday each have one coefficient.

For presentation, the seven weekday effects are recentered around their weekly
mean. These recentered values answer the more intuitive question, “How far is
this weekday above or below an average day?”

### Error process

- \(\varepsilon_t\): residual or leftover on day \(t\), after removing the
  systematic trend, season, and weekday prediction.
- \(p\): autoregressive order, or number of earlier residuals used.
- \(\ell\): lag index, \(\ell=1,\ldots,p\).
- \(\phi_\ell\): persistence coefficient at lag \(\ell\).
- \(\eta_t\): new innovation on day \(t\), after accounting for previous
  residuals.
- \(\sum\): “add all listed terms.”

For M7, \(p=1\), so the error equation simplifies to

\[
\varepsilon_t=\phi\varepsilon_{t-1}+\eta_t.
\]

This says that some of yesterday's unexplained high or low use carries into
today, followed by a new surprise. The fitted \(\phi\) is about 0.375, so that
memory fades quickly.

### Heavy-tailed innovations

\[
\eta_t\sim t_\nu(0,\sigma^2).
\]

- \(\sim\): “is distributed as.”
- \(t_\nu\): Student-\(t\) distribution with \(\nu\) degrees of freedom.
- \(0\): innovations are centred at zero.
- \(\sigma\): innovation scale, measuring the typical size of a new surprise.
- \(\sigma^2\): squared scale parameter.
- \(\nu\): degrees of freedom controlling tail thickness. Smaller values allow
  more extreme observations; as \(\nu\) becomes large, the distribution
  approaches a normal distribution.

Strictly, under this parameterization \(\sigma\) is a Student-\(t\) **scale**,
not its exact marginal standard deviation. When \(\nu>2\),

\[
\operatorname{SD}(\eta_t)
=\sigma\sqrt{\frac{\nu}{\nu-2}}.
\]

The project output informally labels \(\sigma\) as the one-step innovation
standard deviation; “innovation scale” is the mathematically precise name.

The Student-\(t\) is implemented as a normal scale mixture:

\[
\eta_t\mid\lambda_t
\sim N\left(0,\frac{\sigma^2}{\lambda_t}\right),
\qquad
\lambda_t\sim
\operatorname{Gamma}\left(\frac{\nu}{2},\frac{\nu}{2}\right).
\]

- \(\mid\): “conditional on” or “given.”
- \(N(\mu,v)\): normal distribution with mean \(\mu\) and variance \(v\).
- \(\lambda_t\): latent day-specific precision or outlier weight.
- \(\operatorname{Gamma}(\alpha,\beta)\): gamma distribution, using shape
  \(\alpha\) and rate \(\beta\) here.
- A normal day generally has \(\lambda_t\) near 1.
- A small \(\lambda_t\) makes \(\sigma^2/\lambda_t\) large, allowing an unusual
  day without inflating uncertainty for every ordinary day.
- Integrating out \(\lambda_t\) produces the Student-\(t\) distribution above.

### Missing days

The 23 missing \(y_t\) values are not set to zero and the dates are not
deleted. Their \(z_t\) and \(\varepsilon_t\) values are sampled from the same
model. This preserves the AR relationship across gaps and carries missing-data
uncertainty into parameter estimates and forecasts.

## 3. Headline model M7

The model used for the main scientific conclusions is

\[
y_t=\mathbf{x}_t^\top\boldsymbol{\beta}+\varepsilon_t,
\]

\[
\varepsilon_t=\phi\varepsilon_{t-1}+\eta_t,
\qquad
\eta_t\sim t_\nu(0,\sigma^2),
\]

with:

- response measured directly in kWh;
- a centred linear trend;
- \(J=2\) annual harmonic pairs;
- six weekday indicators, with Monday as the reference;
- an AR(1) residual process; and
- heavy-tailed Student-\(t\) innovations.

M7 was selected as the headline model because it was the best of the seven
Bayesian candidates in the held-out 12-month forecast, and because its
predictive intervals were substantially sharper and better calibrated than the
log-scale versions. It is the best tested model in this candidate set, not a
claim that no other model could improve it.

## 4. Priors

The model is Bayesian, so unknown parameters receive probability distributions
before seeing the data. For M7:

\[
\beta_0\sim N(0,100^2),
\]

\[
\beta_r\sim N(0,40^2)
\quad\text{for every other regression coefficient }r,
\]

\[
\sigma\sim\operatorname{half\text{-}t}_3(0,20),
\]

\[
\phi\sim\operatorname{Uniform}(-1,1),
\]

\[
\nu-2\sim\operatorname{Exponential}(\text{mean }10).
\]

Meanings:

- The normal priors for regression coefficients are deliberately wide on the
  kWh scale.
- A half-\(t\) prior keeps \(\sigma>0\) while allowing a broad range of
  plausible scales.
- For AR(1), \(-1<\phi<1\) is the stationary region. Stationarity prevents
  shocks from growing without bound.
- The prior on \(\nu-2\) forces \(\nu>2\), so the Student-\(t\) variance exists,
  while still permitting heavy tails.
- Missing responses have the distribution implied by the regression and AR
  process rather than a separate arbitrary prior.

The log-scale candidates use smaller prior scales because their coefficients
are naturally much smaller: 5 for the intercept, 2 for other coefficients,
and 1 for the half-\(t\) scale prior.

## 5. Candidate models and why they were compared

- **M1:** one harmonic, independent Gaussian errors, log response. Naive
  baseline.
- **M2:** one harmonic, AR(1) Student-\(t\) errors, log response. Tests serial
  dependence and heavy tails.
- **M3:** two harmonics, AR(1) Student-\(t\) errors, log response. Tests the
  required second harmonic.
- **M4:** three harmonics, AR(1) Student-\(t\) errors, log response. Tests a
  more flexible annual curve.
- **M5:** two harmonics, AR(2) Student-\(t\) errors, log response. Tests longer
  residual memory.
- **M6:** two harmonics, AR(1) Gaussian errors, log response. Isolates the value
  of heavy tails.
- **M7:** two harmonics, AR(1) Student-\(t\) errors, identity/kWh response.
  Tests whether the log transformation helps.

WAIC preferred M5 in sample, but the held-out forecasts rejected AR(2). The
extra lag was fitting long absence blocks that did not generalize. Held-out
testing also showed that direct kWh modelling beat the log transformation.
Forecast performance, rather than WAIC alone, therefore determined the
headline model.

Adding a third harmonic improved the log-scale model, but the project did not
fit the corresponding third-harmonic kWh model. That remains an important
extension.

## 6. How the posterior was estimated and checked

The posterior distribution was fitted with a custom Markov chain Monte Carlo
(MCMC) sampler:

- Gibbs updates for regression coefficients, innovation scale, latent
  \(\lambda_t\) weights, and missing responses;
- random-walk Metropolis updates for AR coefficients and
  \(\log(\nu-2)\);
- exact stationary initialization of the first AR residuals;
- four chains of 40,000 iterations for the full-data comparison;
- 10,000 burn-in iterations discarded per chain;
- every fifth remaining iteration retained, producing 24,000 posterior draws.

Important diagnostics:

- \(\widehat R\) compares within-chain and between-chain behavior. Values near
  1 indicate convergence.
- ESS is effective sample size after accounting for dependence between MCMC
  draws. Larger is better.
- For M7, the worst \(\widehat R\) was 1.0069 and minimum bulk ESS was 1,072,
  meeting the project targets.
- A simulated-data validation recovered 14 of 15 known parameters inside 95%
  credible intervals, which is consistent with expected coverage.

## 7. Forecast equations

For one posterior draw, future systematic demand is

\[
\mu_{T+h}=\mathbf{x}_{T+h}^\top\boldsymbol{\beta}.
\]

For M7, future residuals are generated recursively:

\[
\varepsilon_{T+h}
=\phi\varepsilon_{T+h-1}+\eta_{T+h},
\]

\[
y_{T+h}=\mu_{T+h}+\varepsilon_{T+h}.
\]

Repeating this over many posterior draws creates the posterior predictive
distribution. It includes:

- uncertainty in coefficients;
- uncertainty in \(\phi,\sigma,\nu\);
- future random innovations; and
- uncertainty propagated through the AR process.

If the innovation variance is denoted \(v_\eta\), the conditional AR(1)
forecast variance at horizon \(h\), with parameters fixed, is

\[
\operatorname{Var}(\varepsilon_{T+h}\mid\varepsilon_T)
=v_\eta\sum_{r=0}^{h-1}\phi^{2r}.
\]

For Student-\(t\) innovations,

\[
v_\eta=\sigma^2\frac{\nu}{\nu-2}.
\]

As \(h\) grows, the standard-deviation ratio relative to one step approaches

\[
\frac{1}{\sqrt{1-\phi^2}}.
\]

With \(\phi\approx0.375\), the ceiling is about 1.08. Therefore knowledge of
the latest residual makes near-term intervals about 8% narrower, but almost
all of that advantage disappears within four to seven days.

## 8. How forecasts were evaluated

The final twelve months were held out and never used to estimate the model
parameters.

- **Single-origin design:** forecast the entire holdout year from 30 November
  2009.
- **Rolling-origin design:** repeatedly forecast horizons 1–60 after updating
  the current residual state with newly observed consumption.

Scores:

- **Log predictive density (LPD):** rewards a predictive distribution that is
  both sharp and located correctly. Higher is better.
- **Continuous ranked probability score (CRPS):** measures predictive error
  while using the whole forecast distribution. Lower is better and its unit is
  kWh.
- **Coverage:** fraction of observations inside nominal 50% or 95% predictive
  intervals. Good calibration means observed coverage is close to the stated
  percentage.
- **PIT:** predictive cumulative probability at the observed value. A roughly
  uniform PIT distribution indicates calibration.

M7's rolling 95% intervals covered 94.5% at one day and 95.2% at seven days.
Their mean width rose from about 28.9 kWh at one day to about 31.4 kWh at seven
days, then remained almost flat.

## 9. Main findings from M7

All brackets below are 95% posterior credible intervals.

### Overall level and trend

- Average-day underlying level at the record midpoint:
  **25.96 [25.41, 26.50] kWh/day**.
- Trend: **-0.186 [-0.667, 0.292] kWh/year**.
- Percentage trend: **-0.72% [-2.57%, 1.13%] per year**.
- Posterior probability of a decline: **0.775**.
- Because the interval includes zero, there is no conclusive long-term rise or
  decline.

### Annual pattern

- First-harmonic amplitude: **8.08 [7.31, 8.85] kWh**.
- First-harmonic peak-to-trough swing: about **16.15 kWh**.
- First-harmonic peak: approximately **14 January**.
- Second-harmonic amplitude: **2.43 [1.68, 3.19] kWh**.
- Amplitude ratio \(A_2/A_1\): **0.302 [0.204, 0.405]**.
- Probability that \(A_2>0.15A_1\): **0.999**.
- Combined seasonal range: **18.25 [16.48, 20.06] kWh**.

The second harmonic clearly matters. It creates a sharper winter peak and a
longer, flatter summer low.

### Weekly pattern

Effects recentered around the weekly mean:

- Monday: -2.15 kWh.
- Tuesday: -0.41 kWh.
- Wednesday: -0.10 kWh.
- Thursday: -2.75 kWh.
- Friday: -1.12 kWh.
- Saturday: +3.50 kWh.
- Sunday: +3.02 kWh.
- Weekend minus weekday: **4.56 [3.71, 5.42] kWh/day**.
- Posterior probability that weekend demand exceeds weekday demand:
  **1.0000**.

### Error process

- Innovation scale \(\sigma\): **5.59 [5.25, 5.94] kWh**.
- AR coefficient \(\phi\): **0.375 [0.323, 0.426]**.
- Degrees of freedom \(\nu\): **5.81 [4.38, 7.81]**.

These values show modest short-term persistence and clearly heavier tails than
a Gaussian model.

### Unusual days and periods

- The one-step diagnostic is

  \[
  p_t=P(y_t^{\mathrm{rep}}\le y_t\mid\mathbf{y}),
  \]

  where \(y_t^{\mathrm{rep}}\) is a value replicated from the fitted model,
  \(\mathbf{y}\) is the complete observed dataset, and \(p_t\) is the
  posterior predictive probability that a replicate would be no larger than
  the observed day. Values near 0 indicate unusually low use; values near 1
  indicate unusually high use.
- One-step posterior predictive checks found 31 extreme days at the 1% tails:
  23 unusually high and 8 unusually low.
- A one-step check can miss a sustained absence because AR(1) adapts after its
  first few days.
- A separate 14-day rolling-mean posterior predictive check identified 13
  sustained anomalous episodes.
- The clearest low periods resemble summer vacancies; high winter episodes
  are consistent with cold spells.

These are interpretations supported by timing and load shape, not observed
occupancy or weather measurements.

## 10. What the model tells us operationally

- Winter versus summer is the largest predictable difference.
- Weekends require more electricity than weekdays for this household.
- There is not enough evidence to plan around a long-term upward or downward
  trend.
- Yesterday's consumption improves forecasts mainly for the next few days.
- Beyond about one week, the forecast is effectively a seasonal-plus-weekday
  distribution.
- Individual-day forecasts remain wide because one household's behavior is
  intrinsically variable.
- Forecast intervals should be used for planning; a point forecast alone hides
  most of the relevant uncertainty.

## 11. Important limitations

- This is one household, not a population or electricity grid. Results should
  not be generalized without new data.
- There is no temperature variable. Calendar harmonics represent average
  weather but cannot predict a specific cold snap.
- There is no occupancy variable. Vacations are treated as outliers rather
  than a separate away-from-home state.
- Weekday effects are additive and constant, although exploratory work shows a
  larger weekend premium in winter.
- The innovation scale is constant, although residual variation is larger in
  winter than summer.
- Trend and season are mildly confounded because the record is just short of
  four complete years.
- The required 365-day period ignores leap-year drift.
- Conditional WAIC is not fully reliable for selecting forecasting models with
  autocorrelated data; the holdout result is more trustworthy.
- A simple day-of-year climatology beat every Bayesian candidate on the
  single-origin one-year forecast. Its advantage partly reflects similar
  August absences in training and test years, but it also shows that the
  two-harmonic annual curve is not flexible enough for the best long-range
  forecast.

## 12. Recommended next steps

1. Add daily Paris-area temperature; this is the highest-value missing
   predictor.
2. Fit a kWh-scale model with more harmonics, including \(J=3\) through
   \(J=6\), because only the log-scale third-harmonic model was tested.
3. Compare versions with and without the linear trend for long-horizon
   forecasts.
4. Add a seasonal innovation scale so winter can be more variable than summer.
5. Allow the weekend effect to interact with season.
6. Consider a latent home/away state for multi-week absences.
7. Include the day-of-year climatology in the same rolling-origin evaluation
   used for Bayesian candidates.
8. Use daily rather than every-third-day rolling origins and score all horizons
   on a common target-date set.

## 13. Common distinctions

- **Coefficient vs effect:** \(\delta_k\) is relative to Monday; the displayed
  weekday profile is recentered relative to the weekly mean.
- **Residual vs innovation:** \(\varepsilon_t\) includes persistence from prior
  days; \(\eta_t\) is the new one-day shock.
- **Credible interval vs predictive interval:** a credible interval describes
  uncertainty about a parameter; a predictive interval describes uncertainty
  about a future observation and is much wider.
- **Innovation scale vs marginal variability:** \(\sigma\) controls one new
  shock; AR persistence increases longer-run residual variability.
- **In-sample fit vs forecast skill:** WAIC summarizes fitted predictive
  density; the held-out year directly measures forecasting performance.
- **Best candidate vs final possible model:** M7 is the headline model among
  M1–M7, while the climatology result and untested kWh models show clear room
  for improvement.

## 14. Where to verify each part

- Core equations, sampler, forecasting, diagnostics: `src/model.py`
- Candidate model definitions and WAIC: `src/03_fit_models.py`
- Held-out evaluation: `src/04_holdout.py`
- Derived effects and anomaly checks: `src/05_derived.py`
- Posterior summaries: `results/03_convergence.txt`
- WAIC comparison: `results/03_waic.txt`
- Forecast results: `results/04_holdout.txt`
- Interpreted posterior quantities: `results/05_derived.txt`
- Full technical narrative: `report.tex`
