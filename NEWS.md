# dda 0.2.0

---

### New features

-	`dda.bagging` performs bootstrap aggregation (bagging) of `dda.indep`, `dda.resdist`, and `dda.vardist` objects to evaluate the stability of direction dependence decisions, with accompanying `print` and `summary` methods. Test statistics are aggregated with `agg.stat`; p-values are pooled with the asymptotically exact harmonic mean p-value (Wilson, 2019).
-	`dda.indep`, `cdda.indep`, and `dda.resdist` now include a `robust` argument applying Siegel's (1982) repeated median estimation to the causally competing models.
-	In `dda.resdist`, skewness and kurtosis differences under `prob.trans = TRUE` are now computed as target minus alternative, so differences > 0 suggest the target model under both settings of `prob.trans` (previously, differences < 0 suggested the target model under `prob.trans = TRUE`).
-	Joint higher moment statistics now carry short labels in `print`, `summary`, and `plot`: `CoSk` (co-skewness difference), `CoKurt` (co-kurtosis difference), `HS.CoSk` (Hyvarinen-Smith co-skewness), `HS.CoKurt` (Hyvarinen-Smith co-kurtosis, `dda.resdist` only), `CC.CoKurt` (Chen-Chan co-kurtosis), and `Rtanh` (Hyvarinen-Smith hyperbolic tangent). The `stat` argument of `plot` for `cdda.vardist` objects takes these labels; the previous values (`"coskew"`, `"cokurt"`, `"rhs"`, `"rcc"`, `"rtanh"`) still work.
-	The `print` method of `dda.vardist` now reports bootstrap CIs in two tables: marginal higher moment differences (skewness and kurtosis) and joint higher moment differences (`CoSk`, `HS.CoSk`, `CoKurt`, `CC.CoKurt`, and `Rtanh`). The reported statistics are unchanged.
-	Argument, element, and class names now use dot case (e.g., `agg.stat`, `trim.prob`, `aggregated.stats`, `decision.proportions`, `breusch.pagan`, `call.info`, and class `dda.bagging`). The OLS summary of bagged objects is now the S3 generic `ols.summary` (previously, `summary_ols`) with its own `print` method.

### Bug fixes

-	In `dda.indep`, `boot.type = "bca"` now falls back to percentile intervals with a warning when the acceleration constant cannot be calculated (previously, the function stopped).
-	In `cdda.indep`, HSIC, dCor, and MI differences now compare the target model of the target panel with the alternative model of the alternative panel at each moderator value (Li et al., 2023), so differences > 0 suggest the target model (previously, each panel compared its own model with an in-panel reverse model, which favored the target of each panel). `summary` prints one table per statistic, and `plot` draws one panel titled "Target vs. Alternative". MI differences may be biased because the two models rest on differently purified variables and their joint entropies do not cancel.
-	In `cdda.vardist`, `conf.level` is now passed to the bootstrap CIs for continuous moderators (previously, 95% CIs were always returned).
-	In `cdda.indep` and `cdda.vardist`, results for factor moderators are now labeled in the order of the factor levels (previously, labels followed the order of appearance in the data and could be attached to the wrong results).
-	`cdda.indep` and `cdda.vardist` now return a clear error when a continuous moderator is used without `modval` (previously, "argument is of length zero").
-	`cdda.indep` and `cdda.vardist` now return "Specified moderator not found in the target model." for a moderator that is neither numeric nor a factor (previously, "object 'delete.mod' not found").
-	With `modval = "JN"`, `cdda.indep` and `cdda.vardist` now stop with an error when the Johnson-Neyman bounds leave no moderator values in the observed range.
-	For factor moderators, `cdda.indep` now passes on `parallelize` and `cores`.
-	The `summary` methods of `cdda.indep` and `cdda.vardist` now place the Target Model and Alternative Model labels over their columns and print all statistics with three decimals.
-	The `summary` method of `cdda.vardist` now prints the bootstrap CI header above the skewness differences (previously, above the Anscombe-Glynn tests, which are not bootstrapped).
-	The `summary` method of `cdda.vardist` now prints the D'Agostino and Anscombe-Glynn tests of both the predictor and the outcome in each panel (previously, only the predictor tests were printed).
-	Without `ylim`, the `plot` methods of `cdda.indep` and `cdda.vardist` now scale the y-axis to the plotted statistics and include zero (previously, a fixed padding of 0.1 flattened small statistics such as Hyvarinen-Smith tanh differences).
-	The Breusch-Pagan decision in `dda.bagging` now uses the standard Breusch-Pagan tests of the target and alternative models (previously, the robust tests).
-	Loading the package no longer resets the graphics layout or opens a graphics device (`Rplots.pdf`).
-	The help pages of `dda.vardist`, `dda.resdist`, and `cdda.vardist` now explain the statistic labels.
-	Various minor documentation changes and clarifications.


# dda 0.1.1

---

### Bug fixes

-	In `dda.indep`, the bootstrap HSIC method argument is now `hsic.method = "bootstrap"` (previously `hsic.method = "boot"` in 0.1.0).
-	`dda.vardist` and `dda.resdist` now consistently return error messages when the number of bootstrap replications (B) is too small for `boot.type = "bca"` (previously, behavior could be inconsistent).
-	Various minor documentation changes and clarifications.


# dda 0.1.0

---

### Initial release

- First CRAN release of `dda`.
- Includes five core `dda` functions and S3 generics where applicable (`print`, `summary`, `plot`).
- Documentation provided for all exported/user-facing functions.

# dda 0.0.0.9000

---

### Development

- Added NEWS.md using [newsmd](https://github.com/Dschaykib/newsmd) package.
