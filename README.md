
# visPed R package

An R package that extends the kinship2 R package to plot pedigrees
affected by multiple cancer syndromes. Designed to be compatible with
the BayesMendel and associated packages.

## Installation

Use the latest github version by running

``` r
devtools::install_github("gavin-k-lee/visPed")
library(visPed)
```

## Usage

Pass a compatible pedigree (data frame) into the function.

``` r
visPed(pedigree)
```
