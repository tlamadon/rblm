# rblm

[![Testing multiarch installation](https://github.com/tlamadon/rblm/actions/workflows/test-multiarch-install.yml/badge.svg)](https://github.com/tlamadon/rblm/actions/workflows/test-multiarch-install.yml)
[![doc](https://img.shields.io/badge/doc-latest-blue)](https://tlamadon.github.io/rblm/index.html)

R package for estimation of complementarities in models with two-sided heterogeneity. You can find documentation for the different functions here [rblm documentation](https://tlamadon.github.io/rblm/index.html).

There is a simple example on simulated data here: [tutorial](https://tlamadon.github.io/rblm/articles/example.html).

Link to github repo: [github](https://github.com/tlamadon/rblm).

## Install from github using pak

To install directly the package, run the following:

    install.packages("pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s", .Platform$pkgType, R.Version()$os, R.Version()$arch))
    library(pak);
    
    # install our package for the main branch
    pak::pkg_install("tlamadon/rblm")
    
## Install from within package
    
First download the repository using git, then open the projet in Rstudio. Finally run the following:    
    
    install.packages("devtools") # this is in case you do not have devtools already
    require(devtools)
    document(".")
    install(".")
    
## Example

You should find the following working example `vignettes/example.Rmd`.
    
