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
- Symmetric Student t (T)

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

Sample output from
```sh
Rscript xreturns_dist.r spy_vix.csv
```
is
```
file: spy_vix.csv
rows: 8300  date range: 1993-01-29 to 2026-01-20
return_type: log
ret_scale: 100
fixed_nu: 5
return summary stats by asset (sample; kurtosis is excess)
 asset    n  median   mean     sd    skew    kurt      min     max   jb_p dago_p   ad_p
   SPY 8299  0.0682 0.0401 1.1738 -0.2488 11.4002 -11.5886 13.5578 0.0000 0.0000 0.0000
   VIX 8299 -0.4309 0.0058 6.8684  0.9537  6.6975 -44.2449 76.8245 0.0000 0.0000 0.0000

normal fit table
 asset   mean     sd      logLik        AIC        BIC
   SPY 0.0401 1.1738 -13105.2173 26214.4345 26228.4823
   VIX 0.0058 6.8684 -27766.8640 55537.7280 55551.7758

logistic fit table
 asset      mu      s    mean     sd   skew kurtosis      logLik        AIC        BIC
   SPY  0.0655 0.5724  0.0655 1.0382 0.0000   1.2000 -12239.0531 24482.1061 24496.1539
   VIX -0.2611 3.5329 -0.2611 6.4079 0.0000   1.2000 -27219.1925 54442.3850 54456.4328

EGB2 fit table
 asset      mu  sigma      p      q   mean     sd    skew kurtosis      logLik        AIC        BIC
   SPY  0.1178 0.0565 0.0697 0.0770 0.0401 1.0970 -0.2087   2.9823 -11960.1751 23928.3503 23956.4458
   VIX -1.4425 1.5187 0.4044 0.3015 0.0058 6.6662  0.5073   2.5707 -27089.1122 54186.2245 54214.3200

NIG fit table
 asset      mu  delta  alpha    beta   mean     sd    skew kurtosis      logLik        AIC        BIC
   SPY  0.1282 0.7486 0.5671 -0.0663 0.0401 1.1608 -0.5397   7.5035 -11910.6923 23829.3846 23857.4802
   VIX -1.3459 6.2424 0.1455  0.0308 0.0058 6.7785  0.6738   3.9843 -27069.9892 54147.9785 54176.0740

AIC summary table
 asset best_model AIC.Normal AIC.Logistic   AIC.EGB2    AIC.NIG
   SPY        NIG 26214.4345   24482.1061 23928.3503 23829.3846
   VIX        NIG 55537.7280   54442.3850 54186.2245 54147.9785

BIC summary table
 asset best_model BIC.Normal BIC.Logistic   BIC.EGB2    BIC.NIG
   SPY        NIG 26228.4823   24496.1539 23956.4458 23857.4802
   VIX        NIG 55551.7758   54456.4328 54214.3200 54176.0740

model ranking table: SPY
 logLik_model AIC_model BIC_model   dLogLik      dAIC      dBIC
          NIG       NIG       NIG    0.0000    0.0000    0.0000
         EGB2      EGB2      EGB2   49.4828   98.9657   98.9657
     Logistic  Logistic  Logistic  328.3608  652.7215  638.6737
       Normal    Normal    Normal 1194.5250 2385.0499 2371.0021

model ranking table: VIX
 logLik_model AIC_model BIC_model  dLogLik      dAIC      dBIC
          NIG       NIG       NIG   0.0000    0.0000    0.0000
         EGB2      EGB2      EGB2  19.1230   38.2460   38.2460
     Logistic  Logistic  Logistic 149.2033  294.4065  280.3587
       Normal    Normal    Normal 696.8748 1389.7495 1375.7018

average model ranking table
    model  dLogLik      dAIC      dBIC k    mean     sd   skew kurtosis time_sec
      NIG   0.0000    0.0000    0.0000 4  0.0229 3.9697 0.0670   5.7439   3.0600
     EGB2  34.3029   68.6058   68.6058 4  0.0230 3.8816 0.1493   2.7765   2.7200
 Logistic 238.7820  473.5640  459.5162 2 -0.0978 3.7231 0.0000   1.2000   0.0800
   Normal 945.6999 1887.3997 1873.3520 2  0.0230 4.0211 0.0000   0.0000   0.0200

Total time elapsed: 7.26 seconds
```

