
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psychtm: A package for text mining in psychological research

<!-- badges: start -->

[![Project Status:
WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![covr
codecov](https://app.codecov.io/gh/ktw5691/psychtm/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ktw5691/psychtm)
[![R-CMD-check](https://github.com/ktw5691/psychtm/workflows/R-CMD-check/badge.svg)](https://github.com/ktw5691/psychtm/actions)
<!-- badges: end -->

The goal of psychtm is to make text mining models and methods accessible
for social science researchers, particularly those within psychology.
This package allows users to

-   Perform Bayesian estimation of the SLDAX model and popular models
    subsumed by the SLDAX model, including the SLDA model, LDA model,
    and regression models;

-   Estimate SLDAX, SLDA, and regression models with support for both
    continuous and dichotomous outcomes;

-   Obtain Bayesian posterior inferences;

-   Assess model fit using coherence and exclusivity metrics.

## Installation

You can install the development version of psychtm like so:

-   If necessary, first install the `devtools` R package

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
                         K = 3,
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
#>                  Mean       SD  Naive SE Time-series SE
#> I(grade - 1) -0.26447 0.008505 0.0008505      0.0005087
#> topic1        4.81556 0.068277 0.0068277      0.0143247
#> topic2        4.73025 0.071178 0.0071178      0.0131204
#> topic3        4.59128 0.230035 0.0230035      0.1336056
#> effect_t1     0.15480 0.142195 0.0142195      0.0628737
#> effect_t2     0.02683 0.164625 0.0164625      0.0784209
#> effect_t3    -0.18163 0.260353 0.0260353      0.1515511
#> sigma2        1.14072 0.024834 0.0024834      0.0024834
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%      25%       50%       75%   97.5%
#> I(grade - 1) -0.2826 -0.27008 -0.264742 -0.258232 -0.2487
#> topic1        4.7019  4.76530  4.816858  4.853650  4.9594
#> topic2        4.5897  4.68463  4.721874  4.786010  4.8427
#> topic3        4.0576  4.43925  4.659748  4.759582  4.8644
#> effect_t1    -0.1079  0.05395  0.135724  0.254635  0.4519
#> effect_t2    -0.2498 -0.08630 -0.009826  0.133714  0.3969
#> effect_t3    -0.7557 -0.36032 -0.105162 -0.005648  0.1440
#> sigma2        1.0924  1.12802  1.141213  1.155407  1.1908
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

### Common Troubleshooting

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
