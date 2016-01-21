[![CRAN version](http://www.r-pkg.org/badges/version/BayesComm)](http://www.r-pkg.org/pkg/BayesComm)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/BayesComm)](http://www.r-pkg.org/pkg/BayesComm)
[![Research software impact](http://depsy.org/api/package/cran/BayesComm/badge.svg)](http://depsy.org/package/r/BayesComm)
[![Build Status](https://travis-ci.org/goldingn/BayesComm.svg?branch=master)](https://travis-ci.org/goldingn/BayesComm)

# BayesComm

BayesComm is an R package for fitting species interaction distribution models using Bayesian multivariate binary (probit) regression.
This is a work in progress and will be the subject of an upcoming paper. There are a number of updates planned - potentially including a change of name - so expect some fairly major changes.

There is a version of BayesComm [on CRAN](http://cran.r-project.org/web/packages/BayesComm/index.html), but this is the development version so will get the updates first.

### installing BayesComm

To install this development version of BayesComm you can use the ```install_github``` function in the [```devtools```](http://cran.r-project.org/web/packages/devtools/index.html) package, like this:

```{r}
# install devtools if you haven't already
# install.packages('devtools')

# load the package
library(devtools)

# install BayesComm from github (the version from goldingn's repo at least)
install_github('BayesComm', 'goldingn')

# and load it
library(BayesComm)
```

### reporting bugs
If you find a bug in the code or have suggestions for improvements, please let me know via the issues reporting system (button on the right of this page).
