A patch to fix the problems at https://cran.r-project.org/web/checks/check_results_revdbayes.html.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Debian Linux, clang (R-devel) on R-hub
- Fedora Linux, clang (R-devel) on R-hub
- macOS (R-release), ubuntu (R-oldrel, R-release, R-devel), windows (R-release) using the rcmdcheck package
- win-builder (R-devel, R-release and R-oldrelease)

## Downstream dependencies

Apart from warnings and notes unrelated to revdbayes, the downstream dependencies of revdbayes (lax, lite, threshr, smovie, exdex, fitteR, mev, rust, distributions3) passed R CMD check. (I will fix these problems for exdex, lax, rust, threshr, smovie shortly.)
