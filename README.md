# cfda: Categorical Functional Data Analysis

[![R-CMD-check](https://github.com/modal-inria/cfda/workflows/R-CMD-check/badge.svg)](https://github.com/modal-inria/cfda/actions) [![codecov](https://codecov.io/gh/modal-inria/cfda/branch/master/graphs/badge.svg)](https://codecov.io/gh/modal-inria/cfda) 

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/cfda)](https://cran.r-project.org/package=cfda) [![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/cfda?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/cfda) [![Downloads](https://cranlogs.r-pkg.org/badges/cfda)](https://cran.rstudio.com/web/packages/cfda/index.html)

**cfda** provides functions for the analysis of categorical functional data. 

The main contribution is the computation of an optimal encoding (real functional variable) of each state of the categorical functional data.


## Installation

From CRAN:

``` r
install.packages("cfda")
```

From github:

``` r
library(remotes)
install_github("modal-inria/cfda", build_vignettes = TRUE)
```

## Vignette

Once the package is installed, a vignette describing the mathematical background and showing an example is available using the R command:

``` r
RShowDoc("cfda", package = "cfda")
```

## Credits

**cfda** is developed by Cristian Preda (Inria Lille, Université de Lille), Quentin Grimonprez (Inria Lille) and Vincent Vandewalle (Inria Lille, Université de Lille).

Copyright Inria - Université de Lille


## Citation

Preda C, Grimonprez Q, Vandewalle V. Categorical Functional Data Analysis. The cfda R Package. Mathematics. 2021; 9(23):3074. https://doi.org/10.3390/math9233074


``` bibtex
@Article{,
  title = {Categorical Functional Data Analysis. The cfda R Package},
  author = {{Preda} and {Cristian} and {Grimonprez} and {Quentin} and {Vandewalle} and {Vincent}},
  journal = {Mathematics},
  year = {2021},
  volume = {9},
  number = {23},
  url = {https://www.mdpi.com/2227-7390/9/23/3074},
  doi = {10.3390/math9233074},
}
```

## Licence

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.en.html) for more details.
