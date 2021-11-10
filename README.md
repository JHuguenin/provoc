
<!-- README.md is generated from README.Rmd. Please edit that file -->

# provoc <a href='https://github.com/JHuguenin/provoc'><img src='img/imgfile.png' align="right" height="138" /></a>

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/JHuguenin/provoc/workflows/R-CMD-check/badge.svg)](https://github.com/JHuguenin/provoc/actions) -->
<!-- [![R build status](https://github.com/JHuguenin/provoc/workflows/R-CMD-check/badge.svg)](https://github.com/mitchelloharawild/icons/actions?workflow=R-CMD-check) -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/provoc)](https://cran.r-project.org/package=provoc)
<!-- [![Coverage status](https://codecov.io/gh/mitchelloharawild/icons/branch/master/graph/badge.svg)](https://codecov.io/gh/mitchelloharawild/icon?branch=master) -->

<!-- badges: end -->

**P**erform a **R**apid **O**verview for **V**olatile **O**rganic
**C**ompounds.

The `provoc` package has been developed to support PTR-ToF-MS users in
their analyses. It has been designed for a quick import of data into R
and visualization of the first results in a few minutes. It
automatically detects peaks and provides a matrix for further analysis.
Some chemometrics functions are proposed.

It is still a young and wild package that will appreciate feedback and
new ideas for its development. Do not hesitate to contact the author to
get or provide help.

# Installation

The **development** version can be installed from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("jhuguenin/provoc")
```

The package requires the update of many dependencies:

-   `dygraphs` (&gt;= 1.1.0)  
-   `magrittr` (&gt;= 2.0.0)  
-   `MALDIquant` (&gt;= 1.19.0)  
-   `rhdf5` (&gt;= 2.34.0)  
-   `rmarkdown` (&gt;= 2.11.0)  
-   `scales` (&gt;= 1.1.0)  
-   `stringr` (&gt;= 1.4.0)  
-   `usethis` (&gt;= 2.0.0)  
-   `viridis` (&gt;= 0.6.0)  
-   `xts` (&gt;= 0.12.0)

# Usage

``` r
library(provoc)
```

Icons can be inserted inline using inline Icons can also be inserted
using usual R chunks.

``` r
# working directory
wd <- "~/R/data_test/miscalenous" # without final "/"
setwd(wd)

# import
sp <- import.h5(wd)
```

``` r
sp <- import.meta("meta_1")
```

``` r
saveRDS(sp$workflow, "workflow.rds")
wf <- readRDS("workflow.rds")
```

``` r
sp <- re.calc.T.para(sp)
sp <- re.init.T.para(sp)
```

``` r
# a dynamic plot :
dy.spectra(sel_sp = sp$mt$meta[sp$acq,"end"], new_color = FALSE)
# a standart plot :
fx.spectra(sel_sp = sp$mt$meta[sp$acq,"end"], pkm = 137, pkM = 137, leg = "l")
fx.spectra(sel_sp = 1, pkm = 59, pkM = 150)
```

``` r
kinetic.plot(M_num = M.Z(c(59, 137)), each_mass = TRUE,
                         group = "grp1", graph_type = "dy",
                         Y_exp = FALSE, time_format = "date")
# Liste des variables :
# M_num         les masses analysees. M.Z(c(69, 205, 157)) ou c(69.055, 205.158, 157.021)
# each_mass     un graphe pour chaque masse ou non. Logical TRUE or FALSE
# group         groupe avec le nom de la colonne en argument. Ex : "grp1". Or FALSE
# graph_type    choisi soit des graphe fixe "fx" au format tiff soit des graphe dynamique "dy" au format html
# Y_exp         L'ordonnée exponentielle. Logical. TRUE or FALSE
# time_format   L'abscisse est une duree ("time") ou une date ("date").
```
