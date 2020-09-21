# cfda: Categorical Functional Data Analysis

[![Build Status](https://travis-ci.com/modal-inria/cfda.svg?branch=master)](https://travis-ci.com/modal-inria/cfda) [![Build status](https://ci.appveyor.com/api/projects/status/902s96okh97clt5q/branch/master?svg=true)](https://ci.appveyor.com/project/Quentin62/cfda/branch/master) [![codecov](https://codecov.io/gh/modal-inria/cfda/branch/master/graphs/badge.svg)](https://codecov.io/gh/modal-inria/cfda) 

**cfda** provides functions for the analysis of categorical functional data. 

The main contribution is the computation of an optimal encoding (real functional variable) of each state of the categorical functional data.


## Installation

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

## Licence

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU Affero General Public License](https://www.gnu.org/licenses/agpl-3.0.en.html) for more details.
