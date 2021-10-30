
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psychtm: A package for text mining in psychological research

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ktw5691/psychtm/workflows/R-CMD-check/badge.svg)](https://github.com/ktw5691/psychtm/actions)
![CRAN status](https://www.r-pkg.org/badges/version/psychtm) [![Codecov
test
coverage](https://codecov.io/gh/ktw5691/psychtm/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ktw5691/psychtm?branch=main)
<!-- badges: end -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/psychtm)](https://CRAN.R-project.org/package=psychtm) -->

The goal of `psychtm` is to make text mining models and methods
accessible for social science researchers, particularly within
psychology. This package allows users to

-   Estimate the SLDAX topic model and popular models subsumed by SLDAX,
    including SLDA, LDA, and regression models;

-   Obtain posterior inferences;

-   Assess model fit using coherence and exclusivity metrics.

## Installation

Once on CRAN, install the package as usual:

``` r
install.packages("psychtm")
```

Alternatively, you can install the most current development version:

-   If necessary, first install the `devtools` R package,

``` r
install.packages("devtools")
```

### Option 1: Install the latest stable version from Github

``` r
devtools::install_github("ktw5691/psychtm")
```

### Option 2: Install the latest development snapshot

``` r
devtools::install_github("ktw5691/psychtm@devel")
```

## Example

This is a basic example which shows you how to (1) prepare text
documents stored in a data frame; (2) fit a supervised topic model with
covariates (SLDAX); and (3) summarize the regression relationships from
the estimated SLDAX model.

``` r
library(psychtm)
library(lda) # Required if using `prep_docs()`

data(teacher_rate)  # Synthetic student ratings of instructors
docs_vocab <- prep_docs(teacher_rate, "doc")
vocab_len <- length(docs_vocab$vocab)
fit_sldax <- gibbs_sldax(rating ~ I(grade - 1),
                         data = teacher_rate,
                         docs = docs_vocab$documents,
                         V = vocab_len,
                         K = 2,
                         model = "sldax")
eta_post <- post_regression(fit_sldax)
```

``` r
summary(eta_post)
#> 
#> Iterations = 1:100
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 100 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD  Naive SE Time-series SE
#> I(grade - 1) -0.2656 0.007307 0.0007307      0.0007307
#> topic1        4.6165 0.122216 0.0122216      0.0804883
#> topic2        4.8189 0.034301 0.0034301      0.0034301
#> effect_t1    -0.2024 0.134106 0.0134106      0.0884898
#> effect_t2     0.2024 0.134106 0.0134106      0.0884898
#> sigma2        1.1422 0.028296 0.0028296      0.0028296
#> 
#> 2. Quantiles for each variable:
#> 
#>                  2.5%     25%     50%     75%    97.5%
#> I(grade - 1) -0.27849 -0.2711 -0.2659 -0.2601 -0.25175
#> topic1        4.34365  4.5709  4.6584  4.6945  4.76228
#> topic2        4.75032  4.7994  4.8181  4.8420  4.87593
#> effect_t1    -0.51412 -0.2639 -0.1828 -0.1086 -0.01216
#> effect_t2     0.01216  0.1086  0.1828  0.2639  0.51412
#> sigma2        1.08793  1.1245  1.1445  1.1599  1.20649
```

For a more detailed example of the key functionality of this package,
explore the vignette(s) for a good starting point:

``` r
browseVignettes("psychtm")
```

## How to Cite the Package

Wilcox, K. T., Jacobucci, R., Zhang, Z., Ammerman, B. A. (2021).
Supervised latent Dirichlet allocation with covariates: A Bayesian
structural and measurement model of text and covariates. *PsyArXiv*.
<https://doi.org/10.31234/osf.io/62tc3>

## Common Troubleshooting

Ensure that appropriate `C++` compilers are installed on your computer:

-   Mac users will have to download
    [Xcode](https://apps.apple.com/ca/app/xcode/id497799835?mt=12) and
    its related Command Line Tools (found within Xcode’s Preference Pane
    under Downloads/Components).

-   Windows users may need to install
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/). For easier
    command line use, be sure to select the option to install Rtools to
    their path.

-   Most Linux distributions should already have up-to-date compilers.

## Limitations

-   This package uses a Gibbs sampling algorithm that can be
    memory-intensive for a large corpus.

## Getting Help

If you think you have found a bug, please [open an
issue](https://github.com/ktw5691/psychtm/issues) and provide a [minimal
complete verifiable example](https://stackoverflow.com/help/mcve).
