# rblm

R package for estimation of complementarities in models with two-sided heterogeneity. You can find documentation for the different functions here [rblm documentation](https://tlamadon.github.io/rblm/index.html).

There is a simple example on simulated data here: [turorial](https://tlamadon.github.io/rblm/articles/example.html).

Link to github repo: [github](https://github.com/tlamadon/rblm).

## Install from github

To install directly the package, run the following:

    install.packages("devtools") # this is in case you do not have devtools already
    require(devtools)
    install_github("tlamadon/rblm")
    
## Install from within package
    
First download the repository using git, then open the projet in Rstudio. Finally run the following:    
    
    install.packages("devtools") # this is in case you do not have devtools already
    require(devtools)
    document(".")
    install(".")
    
## Example

You should find the following working example `vignettes/example.Rmd`.
    
