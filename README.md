
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psychtm Package for Text Mining in Psychological Research

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis Build
Status](https://travis-ci.com/ktw5691/psychtm.svg?branch=master)](https://travis-ci.com/ktw5691/psychtm)
[![covr
codecov](https://codecov.io/gh/ktw5691/psychtm/branch/master/graph/badge.svg)](https://codecov.io/gh/ktw5691/psychtm)

## Package Features

-   Bayesian estimation of the SLDAX model and popular models subsumed
    by the SLDAX model, including the SLDA model, LDA model, and
    regression models

-   Estimate SLDAX, SLDA, and regression models with support for both
    continuous and dichotomous outcomes

-   Model fit assessment using coherence and exclusivity metrics

-   Model comparison of supervised models using WAIC

## How to Cite the Package

Wilcox, K. T., Jacobucci, R., Zhang, Z., Ammerman, B. A. (2021).
Supervised latent Dirichlet allocation with covariates: A Bayesian
structural and measurement model of text and covariates. *PsyArXiv*.
<https://doi.org/10.31234/osf.io/62tc3>

## Installation

### Ensure that appropriate `C++` compilers are installed on your computer:

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

### Install the `psychtm` package:

-   `psychtm` is not currently available on CRAN, but can be installed
    from this repository
-   If necessary, install the `devtools` R package

``` r
install.packages("devtools")
```

#### Option 1: Install the latest stable development version from the Github source code:

``` r
devtools::install_github("ktw5691/psychtm")
```

#### Option 2: If you are interested in the most recent (untested) development snapshot:

``` r
devtools::install_github("ktw5691/psychtm@devel")
```

## Limitations

-   This package should be considered beta software and is still under
    active development: use at your own risk
-   Documentation is also expected to be expanded and modified in future
    releases. See the **Getting Help** section below.
-   This package uses a Gibbs sampling algorithm that currently can be
    memory-intensive for a large corpus; future developments are planned
    to better scale to larger corpora

## Getting Help

-   If you think you have found a bug, please [open
    issues](https://github.com/ktw5691/psychtm/issues) and provide a
    [minimal complete verifiable
    example](https://stackoverflow.com/help/mcve).
