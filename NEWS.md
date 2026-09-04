# dda 0.2.0

---

### New features

-	`dda.bagging` performs bootstrap aggregation (bagging) of `dda.indep`, `dda.resdist`, and `dda.vardist` objects to evaluate the stability of direction dependence decisions, with accompanying `print` and `summary` methods.
-	`dda.indep`, `cdda.indep`, and `dda.resdist` gain a `robust` argument applying Siegel's (1982) repeated median estimation to the causally competing models.

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
