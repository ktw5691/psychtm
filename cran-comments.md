## Release summary

* This is the first package release on CRAN.

## R CMD check

### Tests

* The package was tested with R CMD CHECK and passed with no ERRORs or WARNINGs
  on:

  * R-devel (2021-10-23 r81086; to become version 4.2.0)
    
    * Debian Linux, R-devel, GCC ASAN/UBSAN
    
    * Fedora Linux, R-devel, clang, gfortran
    
    * Windows Server 2008 R2 SP1, 32/64 bit
    
    * Ubuntu Linux 20.04.1 LTS, R-devel with rchk

  * R-release (2021-08-10; version 4.1.1)
  
    * Ubuntu Linux 20.04.1 LTS, GCC
    
    * Windows Server 2008 R2 SP1, 32/64 bit

  * R-oldrel (2021-03-31; version 4.0.5)
  
    * Windows Server 2008 R2 SP1, 32/64 bit

### Results

* In the tests from Windows, Fedora Linux, and Debian Linux above, there were no
  ERRORs or WARNINGs and R CMD CHECK returned:
  0 errors | 0 warnings | 1 note

* In the tests from Ubuntu Linux and Fedora (R-devel) Linux, there were no ERRORs or
  WARNINGS and R CMD CHECK returned: 0 errors | 0 warnings | 2 notes

  * There were 2 NOTEs:
    
    * checking CRAN incoming feasibility ... NOTE
      Maintainer: 'Kenneth Wilcox <kwilcox3@nd.edu>'
      New submission

      * In response, this is a new release.

    * checking installed package size ... NOTE
      installed size is 16.2Mb
      sub-directories of 1Mb or more:
        libs  15.3Mb

      * In response, it seems that on Ubuntu Linux (R 4.1.1) and Fedora
        Linux (R-devel) architectures, the
        CHECK returns one NOTE because the `libs` subdirectory is above the 1Mb
        threshold (installed size is 16.2Mb on Ubuntu and 8.0Mb on Fedora).
        However, it seems that this NOTE only appears under Linux builds, but
        not under Windows or MacOS builds. My understanding is that this
        inflation of the `libs` subdirectory is due to the use of `Rcpp`.
        Indeed, some functions of the `psychtm` package have been written in C++
        using `Rcpp`. They are needed to perform model estimation for text
        mining. Without the speed improvements gained from those C++ functions,
        this package would become impractical.

## Downstream dependencies

* There are currently no downstream dependencies for this package.
