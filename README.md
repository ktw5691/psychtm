
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psychtm: A package for text mining in psychological research

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ktw5691/psychtm/workflows/R-CMD-check/badge.svg)](https://github.com/ktw5691/psychtm/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/psychtm)](https://CRAN.R-project.org/package=psychtm)
[![Codecov test
coverage](https://codecov.io/gh/ktw5691/psychtm/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ktw5691/psychtm?branch=main)
[![](https://cranlogs.r-pkg.org/badges/grand-total/psychtm?color=blue)](https://CRAN.R-project.org/package=psychtm)
[![](https://cranlogs.r-pkg.org/badges/last-month/psychtm?color=blue)](https://CRAN.R-project.org/package=psychtm)
[![](https://cranlogs.r-pkg.org/badges/last-week/psychtm?color=blue)](https://CRAN.R-project.org/package=psychtm)
<!-- badges: end -->

The goal of `psychtm` is to make text mining models and methods
accessible for social science researchers, particularly within
psychology. This package allows users to

- Estimate the SLDAX topic model and popular models subsumed by SLDAX,
  including SLDA, LDA, and regression models;

- Obtain posterior inferences;

- Assess model fit using coherence and exclusivity metrics.

## Installation

You can install the package from CRAN as usual:

``` r
install.packages("psychtm")
```

Alternatively, you can install the most current development version:

- If necessary, first install the `devtools` R package,

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

set.seed(42) # For reproducibility

data(teacher_rate)  # Synthetic student ratings of instructors
docs_vocab <- prep_docs(teacher_rate, "doc")
vocab_len <- length(docs_vocab$vocab)
fit_sldax <- gibbs_sldax(rating ~ I(grade - 1),
                         data = teacher_rate,
                         docs = docs_vocab$documents,
                         V = vocab_len,
                         K = 2,
                         burn = 100L,
                         m = 600L,
                         model = "sldax")
eta_post <- post_regression(fit_sldax)
```

``` r
summary(eta_post)
#> 
#> Iterations = 101:600
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 500 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean       SD  Naive SE Time-series SE
#> I(grade - 1) -0.2640 0.007956 0.0003558      0.0003558
#> topic1        4.8323 0.029389 0.0013143      0.0015705
#> topic2        4.1385 0.167145 0.0074749      0.0256058
#> effect_t1     0.6938 0.180692 0.0080808      0.0258165
#> effect_t2    -0.6938 0.180692 0.0080808      0.0258165
#> sigma2        1.1347 0.027863 0.0012461      0.0011197
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%     25%     50%     75%   97.5%
#> I(grade - 1) -0.2790 -0.2690 -0.2642 -0.2584 -0.2489
#> topic1        4.7740  4.8134  4.8328  4.8525  4.8876
#> topic2        3.8256  4.0248  4.1407  4.2423  4.4778
#> effect_t1     0.3242  0.5785  0.6910  0.8139  1.0387
#> effect_t2    -1.0387 -0.8139 -0.6910 -0.5785 -0.3242
#> sigma2        1.0798  1.1179  1.1337  1.1508  1.1899
```

For a more detailed example of the key functionality of this package,
explore the vignette(s) for a good starting point:

``` r
browseVignettes("psychtm")
```

## How to Cite the Package

Wilcox, K. T., Jacobucci, R., Zhang, Z., Ammerman, B. A. (2023).
Supervised latent Dirichlet allocation with covariates: A Bayesian
structural and measurement model of text and covariates. *Psychological
Methods*. Advance online publication.
<https://doi.org/10.1037/met0000541>

## Common Troubleshooting

Ensure that appropriate `C++` compilers are installed on your computer:

- Mac users will have to download
  [Xcode](https://apps.apple.com/ca/app/xcode/id497799835?mt=12) and its
  related Command Line Tools (found within Xcode’s Preference Pane under
  Downloads/Components).

- Windows users may need to install
  [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/). For easier
  command line use, be sure to select the option to install Rtools to
  their path.

- Most Linux distributions should already have up-to-date compilers.

## Limitations

- This package uses a Gibbs sampling algorithm that can be
  memory-intensive for a large corpus.

## Getting Help

If you think you have found a bug, please [open an
issue](https://github.com/ktw5691/psychtm/issues) and provide a [minimal
complete verifiable example](https://stackoverflow.com/help/mcve).
