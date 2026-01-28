# Return-Distributions-R

Fit a broad set of univariate return distributions to asset return series, compare models by logLik/AIC/BIC, and optionally plot fitted densities. Includes a simulation harness to validate fitting and model selection.

## Requirements

- R (tested with R 4.4.x)
- Core packages: `sn`, `tseries`, `fBasics`, `nortest`
- Optional (only if plotting those overlays): `logcondens`, `logspline`

## Files

- `xreturns_dist.r`: main fitting script
- `xsim_dist.r`: simulation-based model selection test harness
- `return_utils.r`: return computation + summary stats helpers
- `dist_utils.r`: distribution log-pdf/moments/metrics utilities
- `plot_utils.r`: density plot helpers

## Input Format

CSV with a `Date` column and one or more asset price columns:

```csv
Date,SPY,VIX
1993-01-29,24.2414,12.42
1993-02-01,24.4138,12.33
```

## Usage

```sh
Rscript xreturns_dist.r prices.csv
```

The script expects the first argument to be the CSV path unless you set `infile` in the script.

### Simulation harness

```sh
Rscript xsim_dist.r
Rscript xsim_dist.r 50 5000   # 50 sims per model, 5000 observations
```

`xsim_dist.r` simulates returns from selected distributions, fits the chosen models directly, and reports selection counts and accuracy for logLik/AIC/BIC.

## Key Options (xreturns_dist.r)

Edit at the top of `xreturns_dist.r`:

- `return_type`: `"log"` or `"simple"`
- `ret_scale`: multiply returns (e.g., 100 for percent)
- `fixed_nu`: set degrees of freedom for t-family models (`NA` to estimate)
- `models_to_fit`: vector of model names to fit
- `max_models`: optional cap on number of models fitted
- `do_density_plots`: generate per-asset density plots
- `do_kde_plot`, `do_logcondens_plot`, `do_logspline_plot`: add nonparametric overlays

Environment overrides (useful for automation):

- `XRETURNS_MODELS`: comma-separated list of models
- `XRETURNS_MAX_MODELS`: integer cap or `NA`
- `XRETURNS_DENSITY_PLOTS`: `TRUE`/`FALSE`

Example:

```sh
set XRETURNS_MODELS=Normal,Logistic,SkewNormal,SkewT
set XRETURNS_DENSITY_PLOTS=FALSE
Rscript xreturns_dist.r spy_vix.csv
```

## Models Implemented

### Symmetric / location-scale families
- Normal
- Logistic
- Laplace
- Cauchy
- Hyperbolic secant (Sech)
- Generalized secant hyperbolic (GSH)
- Generalized error distribution (GED)
- Student t (T)

### Skewed or asymmetric families
- Asymmetric Laplace (ALaplace)
- Skewed GED (SGED; Fernandez-Steel skewing)
- Skewed GSH (SGSH)
- NEF-GHS (natural exponential family based on GSH)

### Azzalini / t-family skewing
- Azzalini skew-normal (SkewNormal)
- Azzalini skew-t (SkewT)
- Fernandez-Steel skew-normal (FSSkewNormal)
- Fernandez-Steel skew-t (FSSkewT)
- Jones-Faddy skew-t (JFSkewT)
- Non-central t (NCT)

### Other families
- EGB2 (McDonald)
- NIG (normal-inverse Gaussian)
- Variance-Gamma (VG)
- Hyperbolic (Hyperbolic)
- Generalized hyperbolic (GH)
- Champernowne (with `lambda` controlling sech/logistic limits)
- Normal-Laplace
- Generalized t (GT; McDonald–Newey)
- Skewed generalized t (SGT; Theodossiou)

## Output

- Return summary table (excess kurtosis)
- Fit tables per model with parameters, moments (where defined), logLik/AIC/BIC
- AIC/BIC summary tables across assets and models
- Per-asset model ranking tables (logLik/AIC/BIC)
- Average model ranking table across assets (deltas + average moments)
- Optional density plots saved as `<SYMBOL>.png`

## Notes

- All reported kurtosis values are **excess** kurtosis.
- For t-family models, moments require sufficiently large `nu`.
- Some moments are computed numerically and may emit warnings; likelihoods remain valid.
- Many fits use `optim(..., method = "BFGS")` with numerical gradients.

## Example

```sh
Rscript xreturns_dist.r spy_vix.csv
```

