
# SLAM

<!-- badges: start -->
<!-- badges: end -->

The `SLAM` package is designed to simulate fishery dynamics and data, and to apply an age-structured assessment model to simulated or real fishery data, for a short-lived semelparous species such as octopus.


## Installation
The `SLAM` package can be installed directly from GitHub. A couple of prerequisites are required to install `SLAM` from GitHub. RTools is required because the package uses [TMB](https://cran.r-project.org/package=TMB). The `devtools` package can check if RTools is installed on your system:

```
# Load `devtools` if it's installed
chk <- require('devtools')

# install `devtools` if necessary
if (!chk) {
  install.packages('devtools')
  library('devtools')
}
```

```
# Check if RTools is installed
if (!devtools::find_rtools()) {
  stop('RTools needs to be installed')
}
```
If you get an error message after running the code above, you need to install RTools. Instructions for installing RTools for Windows machines are available [here](https://cran.r-project.org/bin/windows/Rtools/rtools40.html). A quick web search will find similar instructions for  installing RTools on OS machines.

### Installing `SLAM`

Once RTools is installed, `SLAM` can be installed with:

``` 
devtools::install_github('blue-matter/SLAM')
```

## Help Documentation and User Guide

Help Documentation and a User Guide for the `SLAM` package are available on the [SLAM homepage](https://blue-matter.github.io/SLAM/).

