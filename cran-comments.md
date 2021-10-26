## Release summary

This if the first package release on CRAN.

## R CMD check results

0 errors | 0 warnings | 1 note

* There were no ERRORs or WARNINGs.

* There were 2 NOTEs:

  * This is a new release.

  * checking installed package size ... NOTE
    installed size is 16.2Mb
    sub-directories of 1Mb or more:
      libs  15.3Mb

It seems that on LINUX architectures, the CHECK returns one NOTE because the `libs` subdirectory is above the 1Mb threshold. 
However, it seems that this NOTE only appears under LINUX builds, but not under Windows or MacOS builds.
My understanding is that this inflation of the `libs` subdirectory is due to the use of `Rcpp`.
Indeed, some functions of the `psychtm` package have been written in C++ using `Rcpp`.
They are needed to perform model estimation for text mining.
Without the speed improvements gained from those C++ functions, this package would become impractical.

## Downstream dependencies

* There are currently no downstream dependencies for this package.
