# messieR

**M**issingness, **E**rrors-based, **S**ample **S**ize **I**nference for
**E**SM **R**esearch.

`messieR` simulates person-specific VAR(1) processes under the
departures that real experience-sampling (ESM) data actually shows, and
provides three inferential procedures for the fitted model: the
parametric Wald interval, a model-based i.i.d. residual bootstrap (with
bias correction), and a nonparametric stationary block bootstrap. It
exists to answer a design question: *which departures from the textbook
power analysis actually cost you, and how much?* The accompanying thesis
answers that question for non-Gaussian innovations (very little), series
length (a lot), and missingness (about a third of your power at 30%
non-response).

This is research code accompanying an MSc thesis at KU Leuven
(Statistics & Data Science). It is a single-file core with a light dependency 
footprint: base R, stats, and parallel for the analysis scripts, plus MASS for 
drawing Gaussian innovations. It is not on CRAN.

## Installation

Today, source the core file directly:

``` r
source("messieR_core.R")
```

Once the repository is public and split into standard package layout:

``` r
# install.packages("devtools")
devtools::install_github("<username>/messieR")
```

## Quick start

Simulate a bivariate VAR(1) with skewed innovations, knock 30% of the
second variable out under a MAR mechanism, fit, and run all inference
procedures:

``` r
source("messieR_core.R")
set.seed(1)

phi <- matrix(c(0.5, 0.1,
                0.1, 0.5), 2, byrow = TRUE)

Y  <- sim_var1(120, phi, sigma = diag(2))    # complete series, T = 120
Ym <- inject_mar(Y, rate = 0.30)             # 30% MAR on variable 2,
                                             # driven by lagged variable 1

fit <- fit_var1(Ym)                          # gap-aware OLS
fit$phi_hat                                  # point estimates
res <- boot_var1_all(Ym, B = 999)            # Wald + residual + stationary
res[res$coef == "phi12", ]                   # cross-lagged inference
```

`boot_var1_all()` returns one row per procedure per coefficient with the
estimate, interval, p-value, and (for the bootstrap procedures) the
bias-corrected estimate `est_bc`. Supplying `true_phi` (simulation
studies) adds coverage indicators.

## Three things to know before trusting results

**1. Fitting is gap-aware, and that only works if your gaps are
visible.** `fit_var1()` lags the series first and drops incomplete pairs
second, so a transition enters the regression only when both endpoints
are observed and adjacent. This is the transition-level analogue of
listwise deletion, and the order of operations matters: deleting rows
*before* lagging splices non-adjacent observations into false
transitions and attenuates the dynamics. Some public ESM archives record
unanswered prompts as **absent rows rather than NA rows**. If yours
does, reindex against the recorded prompt schedule (and break
transitions at day boundaries) before fitting; in one openESM dataset,
skipping this step understated the autoregressive coefficient by about a
tenth and the cross-lagged coefficient by about an eighth.

**2. The residual bootstrap respects your missingness pattern.**
`boot_residual()` rebuilds each bootstrap series from a fully observed
starting row and then reimposes the observed NA mask, so every bootstrap
refit faces the same gap structure the data did. Comparisons across
missingness mechanisms therefore isolate the mechanism itself.

**3. Two inferential conventions coexist, on purpose.** Interval-type
columns (`ci_lo`, `ci_hi`, `covered`) come from the studentized
(percentile-t) bootstrap interval; the `pval` and `reject_05` columns
come from a percentile bootstrap p-value. These are related but not
identical constructions, and for a badly calibrated resampler they can
disagree. Report them separately, as the thesis does; do not treat
`1 - covered` and `reject_05` as interchangeable.

## Function overview

| Group | Functions |
|------------------------------------|------------------------------------|
| Simulation | `sim_var1()`, `is_stationary()`, `resample_pool()` |
| Missingness injection | `inject_mcar()`, `inject_mar()` |
| Estimation | `fit_var1()` |
| Bootstrap inference | `boot_var1_all()`, `boot_residual()`, `boot_stationary()`, `boot_moving_block()` |
| Tuning defaults | `default_block_length()` (5.03 T\^(1/4)), `default_stationary_p()` |
| Monte Carlo error | `mcse_proportion()`, `mcse_mean()`, `mcse_bias()`, `mcse_mse()` |

`sim_var1()` chooses its innovation source by precedence: a supplied
`residual_pool` (empirical innovations, resampled i.i.d. or in blocks)
beats a supplied `innov_fn` (any generator, e.g. skew-normal or t),
which beats `sigma` (Gaussian). Every series discards a `burn_in`
(default 1000) so the retained observations sit at the stationary
distribution.

A full function-by-function reference with arguments, values, and
details is in `messieR-manual.pdf` (source: `messieR-manual.tex`). Once
the file is split into package layout, `roxygen2::document()`
regenerates native help so `?fit_var1` and `??bootstrap` work in R;
until then the manual PDF is the reference.

## Reproducing the thesis

The analysis scripts call the core and nothing else:

-   `00_revol_validation.R` — external benchmark against the published
    minimum sample sizes of Revol, Lafit & Ceulemans (2024). Check A
    runs in minutes and is the primary validation.
-   `01_run_grid.R` — the 486-condition production grid (fully crossed;
    resumable; per-replicate seeding, so results are bit-identical
    across machines and across interrupted runs).
-   `02_aggregate_figures.R` — health checks, MCSE-aware aggregation,
    size-corrected power, thesis tables, base-R figures.
-   `03_openesm_context.R` — catalogue metadata used to anchor the T
    grid.
-   `04_openesm_demo.R` — the worked real-data demonstration on three
    openESM datasets (two-phase: inspect, then fit).
-   `05_reindex_check.R` — measures the absent-rows hazard by refitting
    one archive under both indexing conventions.

## Citing

Until a software paper exists, cite the thesis:

> Tumuluru, K. (2026). *Realistic simulation-based power analysis for
> person-specific VAR(1) models in ESM research.* MSc thesis, KU Leuven.

Related work this builds on: Revol, Lafit & Ceulemans (2024), *Behavior
Research Methods* 56:7152–7167 (the PAA framework and Shiny application
this machinery is designed to feed); Siepe & Kloft's openESM database
(the real-data demonstrations).

## Roadmap

Package split with `roxygen2` docs and `testthat` tests; an Rmd vignette
built from the openESM demo; empirical-residual condition libraries
drawn from openESM; integration hooks for the PAA Shiny application by
Jordan Revol. Contributions and bug reports are welcome once the
repository is public.

## License

MIT (see `LICENSE`).
