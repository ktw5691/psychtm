
<!-- README.md is generated from README.Rmd. Please edit that file -->

psychtm Package for Text Mining in Psychological Research
=========================================================

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis Build
Status](https://travis-ci.com/ktw5691/psychtm.svg?branch=master)](https://travis-ci.com/ktw5691/psychtm)
[![covr
codecov](https://codecov.io/gh/ktw5691/psychtm/branch/master/graph/badge.svg)](https://codecov.io/gh/ktw5691/psychtm)

Package Features
----------------

-   Bayesian estimation of the SLDAX model and popular models subsumed
    by the SLDAX model, including the SLDA model, LDA model, and
    regression models

-   Estimate SLDAX, SLDA, and regression models with support for both
    continuous and dichotomous outcomes

-   Model comparison of supervised models using WAIC

-   Visualization of model results

How to Cite the Package
-----------------------

Wilcox, K. T. (2020). psychtm: Text mining for psychological research. R
package version 2020.5. Retrieved from
<https://github.com/ktw5691/psychtm>

Installation
------------

To install this package from source:

1.  Ensure that appropriate `C++` compilers are installed on your
    computer:

-   Windows users may need to install
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/). For easier
    command line use, be sure to select the option to install Rtools to
    their path.

-   Mac users will have to download
    [Xcode](https://itunes.apple.com/ca/app/xcode/id497799835?mt=12) and
    its related Command Line Tools (found within Xcode’s Preference Pane
    under Downloads/Components).

-   Most Linux distributions should already have up-to-date compilers
    (or if not they can be updated easily).

1.  Install the `psychtm` package:

-   `psychtm` is not currently available on CRAN, but can be installed
    from this repository: 3.1) If necessary, install the `devtools` R
    package

<!-- -->

    install.packages("devtools")

3.2) Install the latest stable development version from the Github
source code:

    devtools::install_github("ktw5691/psychtm")

3.3) Alternatively, if you are interested in the most recent (untested)
development snapshot:

    devtools::install_github("ktw5691/psychtm@devel")

Limitations
-----------

-   This package should be considered beta software and is still under
    active development: use at your own risk

-   This package uses a Gibbs sampling algorithm that currently can be
    memory-intensive for a large corpus; future developments are planned
    to better scale to larger corpora

Getting Help
------------

-   If you think you have found a bug, please [open
    issues](https://github.com/ktw5691/psychtm/issues) and provide a
    [minimal complete verifiable
    example](https://stackoverflow.com/help/mcve).
