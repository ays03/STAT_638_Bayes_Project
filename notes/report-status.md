# Report status (as of Nov 9)

Reviewed `Report.Rmd` at commit `1717f0d` (Aysegul's Nov 7 Section 3 update) against the rubric in `Bayes Project Instructions.docx`. Updated the same day after Section 7 was added.

**Short version.** Sections 1–5 are drafted, and Sections 4 and 5 are essentially finished. Section 7 is done. Sections 6, 8 and the AI Use Statement have not been started. The report currently does **not** knit, and the draft is already about 3,790 words against the 2,500-word cap before the two missing sections are added.

---

## Where each section stands

| Rubric section | Status | Words now / target |
|---|---|---|
| 1. Scientific question and data | Drafted, far too long | 889 / 200–250 |
| 2. Bayesian model (ours is "EDA and Bayesian Model") | Drafted, needs cleanup | ~1,010 / 250–300 |
| 3. Prior distributions | Half done; prior values are wrong | 753 / 250–300 |
| 4. Posterior computation and diagnostics | Done (one caveat) | 293 / 250–300 |
| 5. Posterior results and interpretation | Done; answers Q1–Q6 | 397 / 350–400 |
| 6. Sensitivity analysis and model checking | Not started | 0 / 250–300 |
| 7. Alternative analysis and comparison | Done (Table 7, Figure 10) | 448 / 400–450 |
| 8. Innovation or extension, plus "why Bayesian" | Not started | 0 / 200–250 |
| AI Use Statement | Missing | 0 / max 50 |
| References | Placeholder | (not counted) |

Word counts are prose only (code chunks and captions excluded, inline math counted as one word), so treat them as approximate.

**Word budget.** Sections 1–5 are allowed at most 1,550 words and are at about 3,340. Section 7 is at 448. Sections 6 and 8 need another 450–550. Roughly 1,800 words still have to come out of Sections 1–3.

---

## Problems in the drafted sections

### Section 1 (Scientific question and data)

- Most of the length is the data cleaning and transformation write-up (1.3.2 and 1.3.3). Supplementary material does not count toward the limit, so moving those to an appendix is the easiest cut.
- 1.2 "Why Is This a Bayesian Problem" overlaps with the "Why Bayesian?" paragraph that now ends Section 7, so it can be cut.

### Section 2 (EDA and Bayesian model)

- Two leftover notes are written as top-level headers, so they render as huge headings in the HTML:
  - line 483: `# **(IDK if we still need this part???)**However, the sustained absence periods ...`
  - line 593: `# HERE!`
- Section 2.3 (Figures 3 and 4) is still on the kWh and log scales from before we switched to √kWh.
- The model equation in 2.4 never states that the response is √kWh, or that the final model uses two harmonics (J = 2) and AR(1) (p = 1).
- 2.4 says the t errors are there to account for the absence periods, but Section 5 says the t errors downweight 20 isolated days and *not* the absence blocks. One of the two needs to change.

### Section 3 (Prior distributions)

- **Wrong prior values.** The table (line 631, "Kilowatt-hour prior") and the text give N(0, 100²), N(0, 40²), and half-t with scale 20. Those are the kWh-scale priors. The final √kWh model actually uses (`src/model.R` lines 37–38):
  - intercept: N(0, 20²)
  - trend, harmonics, day of week: N(0, 5²)
  - σ: half-t, 3 df, scale 2
- Two overlapping versions of the section are pasted together: the "First, … Fifth" walkthrough and the shorter version below the figure. Only one should stay.
- Line 669 "Fourth," is blank; that is where the AR coefficient prior paragraph should go.
- The "Fifth" paragraph has a broken `[(` fragment and repeats itself.
- Aysegul's "New Suggestion" paragraph (line 659) promises a prior sensitivity check on Bayes factors and the tail parameter. The report does not compute any Bayes factors, and no prior-sensitivity fits exist yet. Either run that check in Section 6 or reword the promise.

### Section 4 (Posterior computation)

- The text says the sampler uses a half-t prior on σ. Line 328 of `src/model.R` has `rate <- nu_sigma / (2 * xi) + 0.5 * ss`, which should be `nu_sigma / xi + 0.5 * ss`. With the extra 2, the σ² step and the ξ step do not match, so the implied prior is not a half-t.
- The data term dominates that update, so fixing it moves the σ posterior by under 0.1%, and nothing reported should change. Fixing line 328 and rerunning `src/06_fit_sqrt.R` (about 12 minutes) makes the Section 3 and 4 claims true.

### Section 5 (Results)

- Content is done. It refers ahead to Section 6, which does not exist yet.

---

## Why the report does not knit

Confirmed by knitting a copy of the repo. There are three blockers, in the order knitr hits them:

1. **Duplicate chunk labels.** `fig1` and `fig2` are each used twice. The old chunks marked "DON'T NEED THIS CHUNK ANYMORE" at lines 174 and 486 still have the same labels. Deleting those two chunks fixes it.
2. **Figure 2 crashes** with `'x' and 'y' lengths differ`. The Figure 2 chunk expects `y` to be the 1,417 observed days from the Figure 1 chunk, but the chunks in between reset `y <- d$kwh` (lines 367 and 412), which is all 1,440 calendar days. Resetting `y <- as.numeric(obs$kwh)` at the top of the Figure 2 chunk fixes it.
3. **`kableExtra` must be installed** for Table 3 (line 588). Run `install.packages("kableExtra")`, or drop the `kable_styling()` call.

With those three fixed, the report knits in a few seconds, because the MCMC results are cached in `results/06_sqrt/`.

The `Report.html` lying around locally is from Oct 31, so it is stale.

---

## How the rendered HTML looks

- `knitr::opts_chunk$set(echo = TRUE)` (line 8) shows code everywhere: 9 visible code blocks, about 105 lines of R.
- 6 blocks of raw R output appear: `## [1] 0.4` (twice), `## [1] 57`, `quartz_off_screen 2` (twice, from `dev.off()`), and the `here() starts at ...` message.
- The rubric says routine code should be hidden and raw output avoided. Setting `echo = FALSE` in the setup chunk (or `code_folding: hide` in the YAML), wrapping `dev.off()` in `invisible()`, and deleting the bare print lines would cover it.
- The Box–Cox profile likelihood plot (line 446) has no caption or figure number.
- Table numbering: there are two "Table 2"s (innovations in 2.2 and priors in 3) and no Table 4.

### References

- Currently just two entries under a plain `# REFERENCES:` heading.
- `ref.bib` exists but is not hooked into the YAML header (no `bibliography:` field).
- `ref.bib` has `chib1993`, `geweke1993`, and `vehtari2021` each twice.
- Works cited in the text are missing from it, for example Gelman (2006), Huang and Wand (2013), and the UCI dataset itself.

---

## Material we already have for the missing sections

### Section 6 (Sensitivity and model checking)

- Already run, in `results/06_sqrt/`: the one-harmonic fit (`J1`), the AR(2) fit (`ar2`), and the kWh-scale fit (`kwh`), with paired WAIC comparisons in `summary.rds`. These cover sensitivity to the seasonal terms, the AR order, and the response scale.
- Holdout coverage by horizon is already computed (used in Figure 8).
- Still needed:
  - the prior sensitivity run that Section 3 now promises;
  - the seasonal coverage check of the predictive intervals that Section 2.2 says will be "in Section 6".

### Section 7 (Alternative analysis): done

- `src/08_frequentist.R` fits the final √kWh model (trend, two harmonics, day of week, AR(1) with t errors) by maximum likelihood with no priors, in about 5 seconds. It writes `results/08_freq/` (Table 7 and `summary.rds`) and `R_figures/10_freq_vs_bayes.png` (Figure 10).
- Missing days: the likelihood conditions on the first day of each observed run, so the 23 missing days drop out (1,408 terms) instead of being imputed.
- kWh-scale confidence intervals come from simulating the estimator's asymptotic normal distribution (Krinsky and Robb, 1986). Forecasts are scored on the same 113 rolling origins as Section 5.
- Result: every typical-day estimate agrees with the posterior mean within 0.4 kWh, and interval widths are within 2%. The trend p-value is 0.55, against P(decline) = 0.76. The second harmonic is kept by both (LRT p < 0.001). Plug-in forecast intervals are about 2% narrower, with 97.2% versus 97.6% coverage.
- The section ends with the "Why is the Bayesian approach justified?" answer.
- **Rerun note:** the Bayesian column of Table 7 is read from `results/06_sqrt/summary.rds`. If the Bayesian fit is rerun (for example after the σ-prior fix), delete `results/08_freq/summary.rds` and knit again, or run `Rscript src/08_frequentist.R`.

### Section 8 (Innovation or extension)

The model already contains three items from the rubric's list of acceptable extensions. Pick one (or two) and write up motivation, what changed, and what we learned:

- t errors as an alternative likelihood for heavy tails;
- Bayesian imputation of the 23 missing days (data augmentation in the Gibbs sampler);
- the 14-day posterior predictive anomaly check (Figure 9) as a problem-specific predictive quantity.

Section 7 already answers "Why is the Bayesian approach justified?", so Section 8 can point to it rather than repeat it, and 1.2 can be cut.

### AI Use Statement

Not written. Maximum 50 words: which tools, how they were used, and how we verified the output.
