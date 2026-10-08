## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS (aarch64), R 4.5.3
* win-builder (devel and release): to be run before submission

## Notes for the reviewer

* All random numbers come from R's RNG (`R::norm_rand`, `R::unif_rand`, ...), so `set.seed()` and `seed =` make results reproducible. No OpenMP.
* Long-running examples are in `\donttest{}`.
