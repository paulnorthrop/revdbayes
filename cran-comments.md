A patch to fix the issues at https://cran.r-project.org/web/checks/check_results_revdbayes.html. Originally, there were ERRORs on r-release-macos-x86_64 and r-oldrel-macos-x86_64, stemming from the unit tests, but these seem to be false positives because they disappeared.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- macOS (R-release), ubuntu (R-oldrel, R-release, R-devel), windows (R-release) using the rcmdcheck package
- macOS builder (R-release)
- win-builder (R-devel, R-release and R-oldrelease)

## Downstream dependencies

Apart from warnings and notes unrelated to revdbayes, the downstream dependencies of revdbayes (lax, lite, threshr, smovie, exdex, fitteR, mev, rust, distributions3) passed R CMD check.
