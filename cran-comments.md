## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

The note reports two things:

* This is a new submission.
* Suggests 'dimsprepr', which is not in a mainstream repository yet. It is by
  the same maintainer and was submitted to CRAN on 2026-11-22. pondeR uses it
  only conditionally: examples, tests and the vignette skip it when it is absent.

## Test environments

* Local Windows 11, R 4.5.2, with `R CMD check --as-cran --run-donttest`
* Local Windows 11, R 4.5.2, with the suggested packages unavailable

## Reverse dependencies

There are no reverse dependencies.
