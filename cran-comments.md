## Submission - dda 0.2.0

This is a feature update to the CRAN package 'dda' (current CRAN version 0.1.1).
All changes are backward-compatible.

## Summary of changes

* New dda.bagging(): bootstrap aggregation of dda.indep, dda.resdist, and
  dda.vardist objects to evaluate the stability of direction dependence
  decisions, with accompanying print() and summary() methods.
* New 'robust' argument for dda.indep(), cdda.indep(), and dda.resdist():
  Siegel (1982) repeated median estimation of the causally competing models.
* See NEWS.md for the full changelog.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.5.2
* win-builder: R-devel and R-release
* macOS: R-release

## Reverse dependencies

There are no reverse dependencies on CRAN. Confirmed with
revdepcheck::revdep_check() on 2026-09-02: 0 reverse dependencies
checked, 0 new problems.

## Maintainer

Wolfgang Wiedermann <wiedermannw@missouri.edu>
