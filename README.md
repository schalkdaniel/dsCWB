
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Actions
Status](https://github.com/schalkdaniel/dsCWB/workflows/R-CMD-check/badge.svg)](https://github.com/schalkdaniel/dsCWB/actions)
[![License: LGPL
v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![codecov](https://codecov.io/gh/schalkdaniel/dsCWB/branch/main/graph/badge.svg?token=0K9P2WBKNH)](https://codecov.io/gh/schalkdaniel/dsCWB)

# Component-wise boosting for DataSHIELD

The package provides functionality to conduct and visualize
component-wise boosting on decentralized data. The basis is the
DataSHIELD\](<https://www.datashield.org/>) infrastructure for
distributed computing. This package provides the calculation of the
[**component-wise
boosting**](https://www.tandfonline.com/doi/abs/10.1198/016214503000125).
Note that DataSHIELD uses an option `datashield.privacyLevel` to
indicate the minimal amount of numbers required to be allowed to share
an aggregated value of these numbers. Instead of setting the option, we
directly retrieve the privacy level from the
[`DESCRIPTION`](https://github.com/schalkdaniel/dsCWB/blob/master/DESCRIPTION)
file each time a function calls for it. This options is set to 5 by
default.

## Installation

At the moment, there is no CRAN version available. Install the
development version from GitHub:

``` r
remotes::install_github("schalkdaniel/dsCWB")
```

#### Register methods

It is necessary to register the assign and aggregate methods in the OPAL
administration. These methods are registered automatically when
publishing the package on OPAL (see
[`DESCRIPTION`](https://github.com/schalkdaniel/dsCWB/blob/main/DESCRIPTION)).

Note that the package needs to be installed at both locations, the
server and the analysts machine.
