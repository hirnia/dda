# dda 0.2.0

---

### New features
-	Under `prob.trans = TRUE`, `dda.resdist` now computes skewness and kurtosis differences as target minus alternative, so that differences > 0 suggest the target model under both `prob.trans = FALSE` and `prob.trans = TRUE` (previously, differences < 0 suggested the target model under `prob.trans = TRUE`).
-	The print method of `dda.vardist` now reports bootstrap CIs in two tables: marginal higher moment differences (skewness and kurtosis) and joint higher moment differences (co-skewness, Hyvarinen-Smith co-skewness, co-kurtosis, Chen-Chan co-kurtosis, and Hyvarinen-Smith tanh). The reported statistics are unchanged.
-	`dda.bagging` performs bootstrap aggregation (bagging) of `dda.indep`, `dda.resdist`, and `dda.vardist` objects to evaluate the stability of direction dependence decisions, with accompanying `print` and `summary` methods.
-	`dda.indep`, `cdda.indep`, and `dda.resdist` now include a `robust` argument applying Siegel's (1982) repeated median estimation to the causally competing models.

### Bug fixes

-	In `dda.indep`, `boot.type = "bca"` now falls back to percentile intervals with a warning when the acceleration constant cannot be calculated (previously, the function stopped).
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
