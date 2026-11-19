# dda 0.2.0

---

### New features
-	Under `prob.trans = TRUE`, `dda.resdist` now computes skewness and kurtosis differences as target minus alternative, so that differences > 0 suggest the target model under both `prob.trans = FALSE` and `prob.trans = TRUE` (previously, differences < 0 suggested the target model under `prob.trans = TRUE`).
-	The print method of `dda.vardist` now reports bootstrap CIs in two tables: marginal higher moment differences (skewness and kurtosis) and joint higher moment differences (co-skewness, Hyvarinen-Smith co-skewness, co-kurtosis, Chen-Chan co-kurtosis, and Hyvarinen-Smith tanh). The reported statistics are unchanged.
-	`dda.bagging` performs bootstrap aggregation (bagging) of `dda.indep`, `dda.resdist`, and `dda.vardist` objects to evaluate the stability of direction dependence decisions, with accompanying `print` and `summary` methods. Test statistics are aggregated with `agg.stat`; p-values are pooled with the asymptotically exact harmonic mean p-value (Wilson, 2019).
-	Argument, element, and class names now use dot case (e.g., `agg.stat`, `trim.prob`, `aggregated.stats`, `decision.proportions`, `breusch.pagan`, `call.info`, and class `dda.bagging`). The OLS summary of bagged objects is now `ols.summary` (previously `summary_ols`), an S3 generic whose method for `dda.bagging` objects returns an object of class `ols.summary` with its own `print` method.
-	`dda.indep`, `cdda.indep`, and `dda.resdist` now include a `robust` argument applying Siegel's (1982) repeated median estimation to the causally competing models.

### Bug fixes

-	In `dda.indep`, `boot.type = "bca"` now falls back to percentile intervals with a warning when the acceleration constant cannot be calculated (previously, the function stopped).
-	In `cdda.vardist`, `conf.level` is now passed to the bootstrap CIs for continuous moderators (previously, 95% CIs were always returned).
-	The `summary` methods of `cdda.indep` and `cdda.vardist` now align the Target Model and Alternative Model headers with their columns.
-	Without `ylim`, the `plot` methods of `cdda.indep` and `cdda.vardist` now scale the y-axis to the plotted statistics and include zero (previously, a fixed padding of 0.1 flattened small statistics such as Hyvarinen-Smith tanh differences).
-	Loading the package no longer resets the graphics layout or opens a graphics device (`Rplots.pdf`).
-	In `cdda.indep`, HSIC, dCor and MI difference statistics now compare the target model of the target panel with the alternative model of the alternative panel at each moderator value (Li et al., 2023), so differences > 0 suggest the target model. `summary` prints one table per statistic and `plot` draws one panel. Previously, each panel compared its own model with an in-panel reverse model, which favored the target of each panel. MI differences are computed the same way but may be biased, because the two models rest on differently purified variables and their joint entropies do not cancel; see Details in `?cdda.indep`.
-	In `cdda.indep` and `cdda.vardist`, results for factor moderators are now labelled in the order of the factor levels (previously, labels followed the order of appearance in the data and could be attached to the wrong results).
-	`cdda.indep` and `cdda.vardist` now return a clear error when a continuous moderator is used without `modval` (previously, "argument is of length zero").
-	`cdda.indep` and `cdda.vardist` now return "Specified moderator not found in the target model." for a moderator that is neither numeric nor a factor (previously, "object 'delete.mod' not found").
-	With `modval = "JN"`, `cdda.indep` and `cdda.vardist` now stop with an error when the Johnson-Neyman bounds leave no moderator values in the observed range.
-	For factor moderators, `cdda.indep` now passes `parallelize` and `cores` on.
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
