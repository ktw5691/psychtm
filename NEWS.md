# psychtm 2021.1.0

* Initial release for [CRAN](https://CRAN.R-project.org/).

## Major changes

* Requires R packages:
  * Requires `coda` v0.4 or higher.
  * Requires `Rcpp` v0.11.0 or higher.
  * Requires `tibble` v2.1.3 or higher.
* Added a vignette to illustrate package use: closes (#47).
* Added example data set `teacher_rate` for demonstration in vignette.
* New function `prep_docs()` to prepare documents in a data frame for modeling with `gibbs_sldax()`.
* New function `post_regression()` to summarize regression relationships for objects of class `Mlr`, `Logistic`, or `Sldax`: closes (#7).

## Minor changes

* Added a `NEWS.md` file to track changes to the package.
* Suggests R packages:
  * Suggests `spelling`.
  * Suggests `knitr` v1.22 or higher.
  * Suggests `lda`
  * Suggests `testthat` v3.0.2 or higher.
  * Suggests `rmarkdown`.
* Uses `roxygen2` v7.1.2 with Markdown support for documentation.
* Expanded and improved documentation.
* Remove user-facing documentation in man/ for internal functions.
* Language for package now listed as `en-US`.
* S4 class definitions moved to a single file.
* S4 generic functions moved to a single file.
* S4 methods for `Sldax` objects moved to a separate file.
* Deprecated `gg_coef()`; this function will be removed in a future release.
* Use of `stop()` and `warning()` in the event of errors or warnings instead of `print()` and `cat()`.
* Use of `Rcpp_cout` and `Rcpp_cerr` for printing messages and errors from compiled functions.
* Improved README.
* Added CITATION information.

## Bug fixes

* Fixed namespace issues for importing and exporting packages, methods, and functions.
* Corrected bug when fitting an SLDAX model with only a single predictor where variable names were accidentally omitted.
* Fixed bug where `get_zbar.Sldax()` was missing a return statement.
* Makevars.win no longer asks for `OpenMP` as `OpenMP` is not currently used.

# psychtm 2020.9-alpha

## Non-breaking Changes

* Fixed bug when correcting label switching for an LDA model.

# psychtm 2020.8-alpha

## Non-breaking Changes

* Fixed bug when correcting label switching when manifest predictors for an SLDAX model were included where the manifest predictors' regression coefficients were not stored in the model object.

# psychtm 2020.6-alpha

## Non-breaking Changes

* Fixes to Travis CI build configuration file.
* Updated README.

# psychtm 2020.5-alpha

## Non-breaking Changes

* Minor documentation fix.

# psychtm 2020.4-alpha

## New Features

* *Coherence* [(Mimno, Wallach, Talley, Leenders, & McCallum, 2011)](https://dl.acm.org/doi/10.5555/2145432.2145462) and *exclusivity* [(Roberts, Stewart, & Airoldi, 2013)](https://scholar.princeton.edu/sites/default/files/bstewart/files/stmnips2013.pdf) metrics of topic model fit can now be computed:
    * `get_coherence()` computes the coherence score for each topic,
    * `get_exclusivity()` computes the exclusivity score for each topic.
* Label switching correction using the Stephens (2000) algorithm is now supported for SLDAX and other topic models: closes (#42).
    * `gibb_sldax()` gains the argument `correct_ls` (`TRUE` / `FALSE`),
    * A logical flag indicating whether a label switching correction was applied is now stored in `@extra$corrected_label_switching` for objects of class `Model`.
* Added options to sample the topic proportions and topic-word probabilities or just obtain a point estimate to `gibbs_sldax()`:
    * Argument `sample_theta` (default = `TRUE`) for the topic proportions,
    * Argument `sample_beta` (default = `TRUE`) for the topic-word probabilities.
* Added option to store the chain of sampled topic assignments to `gibbs_sldax()`.
    * Argument `return_assignments` (default = `FALSE`).
    * Keeping this `FALSE` is dramatically more memory efficient: closes (#37).
* Added a check to ensure that the number of observations supplied to `docs` and `data` arguments of `gibbs_sldax()` are equal since missing data in either argument is currently not supported: closes (#36).
* Extractor functions are now available to safely obtain the contents of object slots: closes (#1).
    * For example, to obtain the posterior samples of regression coefficients from a fitted SLDAX model called `my_fit`, use `eta(my_fit)`.

## Breaking Changes

* Fixed WAIC calculations: closes (#11).
* Created hierarchical S4 classes that inherit slots as needed from parents (may not be a breaking change for many users):
    * `Model`: a generic model class,
    * `Mlr`: a linear regression class which inherits from `Model`,
    * `Logistic`: a logistic regression class which inherits from `Model`,
    * `Sldax`: a topic model class which inherits from `Mlr` and `Logistic`.

## Non-breaking Changes

* Imports `label.switching` package.
* Changed default prior hyperparameters in `gibb_sldax()` for the Dirichlet priors on the topic proportions and on the topic-word probabilities to 1.0 (flat over the simplex).
* Minor fix in calculation for when to message user that burn-in period is finished if `verbose = TRUE` specified to `gibbs_sldax()`.

# psychtm v2020.3-alpha

## New Features

* Fixed error in `gibbs_*()` functions where interrupting the process (e.g., hitting the `Esc` key while running) and then trying to run these functions again would fail with the error `ERROR: there is already an InterruptableProgressMonitor instance defined`.

## Breaking Changes

* Requires `RcppProgress` v0.4.2.

## Non-breaking Changes

* None.

# psychtm v2020.2-alpha

## New Features

* Fix initialization error in `gibbs_logistic()` where regression coefficients could start with infinite values and the algorithm could not recover: closes (#39).

## Breaking Changes

* None.

## Non-breaking Changes

* See above in **New Features**.

# psychtm 2020.1-alpha

## New Features

* Better error handling for `gibbs_*()` functions.
* Increased maximum number of iterations attempted from 1,000 to 10,000 when trying to sample regression coefficients in `gibbs_sldax()` with `constrain_eta = TRUE` and now warns user if exceeded.

## Breaking Changes

* None.

## Non-breaking Changes

* Import `Rcpp`, `RcppArmadillo`, and `RcppProgress` packages.
* Split up C++ codebase into multiple files for easier maintenance going forward and faster compilation if using `ccache`.

# psychtm 2019.7-alpha

## New Features

* Added unit tests (still a work in progress to cover entire package).
* Documentation fixes for a number of functions: closes (#33).
* Faster computation for `get_zbar()`.
* Better argument checks for some functions.

## Breaking Changes

* Intermediate C++ functions are no longer exported: `.est_betak()`, `.est_thetad()`, `.count_topic_word()`, `.gibbs_mlr_cpp()`, `.gibbs_logistic_cpp()`, `.gibbs_sldax_cpp()`, `rmvnorm_cpp()`, `pwaic_d()`, `waic_d()`.
* `get_toptopics()` no longer requires a `Sldax` model object and instead needs a topic-proportions matrix for all documents ('theta'): closes (#32).
* Removed obsolete classes (`LDA`, `Slda`, `Sldalogit`) and reorganized `Sldax` class.

## Non-breaking Changes

* Changes to `est_beta()`/`est_theta()`: remove unwanted dim name in returned matrix; handle error in case where chain was shorter than 10 iterations.
* `gg_coef()` no longer automatically prints a plot if the result of `gg_coef()` is assigned to an object.
* Vectorized computation in `term_score()` for 2-3x speed improvement, especially for large matrices.
* Import `is()` from `methods` package.
