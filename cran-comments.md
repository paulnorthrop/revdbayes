## R CMD check results

0 errors | 0 warnings | 1 note

evdbayes is in Suggests: but it was archived on 2022-04-23.  revdbayes does not import functions from evdbayes.  It uses evdbayes only in 2 of the package vignettes and in tests.  requireNamespace() is used to make these tasks conditional on evdbayes being available. 

## Test environments

- Debian Linux, GCC (R-patched and R-devel) on R-hub
- Fedora Linux, GCC (R-devel) on R-hub
- Oracle Solaris 10, x86, 32 bit, R-release
- macOS 10.13.6 High Sierra, R-release on R-hub
- win-builder (R-devel and R-release)

## Downstream dependencies

The downstream dependencies of revdbayes (lax, lite, threshr, smovie, exdex, fitteR, mev, rust, distributions3) passed R CMD check.
