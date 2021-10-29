## Release summary

* This is the first package release on CRAN.

## R CMD check results

* The package was tested with R CMD CHECK on:

  * Ubuntu Linux 20.04.1 LTS, R-release, GCC

  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit

  * Fedora Linux, R-devel, clang, gfortran

  * Debian Linux, R-devel, GCC ASAN/UBSAN

* In the tests from Windows, Fedora Linux, and Debian Linux above, there were no ERRORs or WARNINGs.
  
  * R CMD CHECK returned: 0 errors | 0 warnings | 1 note

* In the tests from Ubuntu Linux, there were no ERRORs or WARNINGS.

  * R CMD CHECK returned: 0 errors | 0 warnings | 2 notes

  * There were 2 NOTEs:
    
    * checking CRAN incoming feasibility ... NOTE
      Maintainer: 'Kenneth Wilcox <kwilcox3@nd.edu>'
      New submission

      * In response, this is a new release.

    * checking installed package size ... NOTE
      installed size is 16.2Mb
      sub-directories of 1Mb or more:
        libs  15.3Mb

      * In response, it seems that on Ubuntu Linux architectures, the CHECK
        returns one NOTE because the `libs` subdirectory is above the 1Mb
        threshold. However, it seems that this NOTE only appears under Linux
        builds, but not under Windows or MacOS builds. My understanding is that
        this inflation of the `libs` subdirectory is due to the use of `Rcpp`.
        Indeed, some functions of the `psychtm` package have been written in C++
        using `Rcpp`. They are needed to perform model estimation for text
        mining. Without the speed improvements gained from those C++ functions,
        this package would become impractical.

## Downstream dependencies

* There are currently no downstream dependencies for this package.
