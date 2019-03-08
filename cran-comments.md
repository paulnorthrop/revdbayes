This is a resubmission of revdbayes 1.3.3, which was submitted on 28/2 but didn't make it to CRAN owing to a failing test in a strong dependency:

*** Changes to worse in reverse dependencies ***
Debian: <https://win-builder.r-project.org/incoming_pretest/revdbayes_1.3.3_20190228_173557/reverseDependencies/summary.txt>

I am also the maintainer for affected package, threshr.  I have fixed the problem (in threshr) and submitted a new version, threshr 1.0.1 to CRAN.  The check results for threshr 1.0.1 are now complete on all platforms.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Fedora Linux, clang, gfortran (on r-hub), R-devel 
- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- win-builder (R-devel and R-release)

## Downstream dependencies

The five downstream dependencies of revdbayes (mev, smovie, threshr, fitteR, rust) passed R CMD check.
