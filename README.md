
<!-- README.md is generated from README.Rmd. Please edit that file -->

# demdrainer

[![Travis build
status](https://travis-ci.org/bocinsky/demdrainer.svg?branch=master)](https://travis-ci.org/bocinsky/demdrainer)
[![Coverage
status](https://codecov.io/gh/bocinsky/demdrainer/branch/master/graph/badge.svg)](https://codecov.io/github/bocinsky/demdrainer?branch=master)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bocinsky/demdrainer")
```

## Drain a DEM

``` r
  library(magrittr)
  library(demdrainer)
  library(mapview)
  library(leafsync)
  
  dem <- FedData::get_ned(demdrainer::mcphee, label = "McPhee")
  
  drained_dem <- drain_dem(dem, label = "McPhee")
  
  sync(mapview(dem), mapview(drained_dem))
  
```

## Contributing

Please note that this project is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project
you agree to abide by its terms.
