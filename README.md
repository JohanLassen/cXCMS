
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cXCMS

<!-- badges: start -->
<!-- badges: end -->

The goal of cXCMS is to …

### Installation

You can install the development version of cXCMS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("johanLassen/cXMS")
```

### Overview

The cluster xcms package is dependent on the regular xcms package. We
aim to make as few changes as possible to ease the use and maintenance
of the package.

``` r
library(cXCMS)
library(xcms)
```

Here’s the idea behind the package: The workflow makes the peak picking
and the peak integration parallel to (1) reduce memory consumption and
(2) improve support for multi-node high performance cluster processing.
The first achievement benefits all workflows by lowering memory
requirements, while the second achievement allows parallelization. All
of the processes are orchestrated by a simple metadata dataframe that we
will describe below.

##### Workflow

![workflow](./images/workflow.png)
